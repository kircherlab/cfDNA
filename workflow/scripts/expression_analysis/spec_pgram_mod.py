import pandas as pd
import numpy as np
from scipy import stats, signal, fft


def spec_pgram(
    x,
    xfreq=1,
    spans=None,
    kernel=None,
    taper=0.1,
    pad=0,
    fast=True,
    demean=False,
    detrend=True,
    minimal=True,
    option_summary=False,
    **kwargs
):
    """Computes the spectral density estimate using a periodogram.

    Args:
        x (numpy array): Univariate or multivariate time series.
        xfreq (optional): Number of samples per unit time. Defaults to 1.
        spans (optional): Sequence of spans for convoluted Daniell smoothers. Defaults to None.
        kernel (optional): Defines Kernel for smoothing. Defaults to None.
        taper (optional): Defines proportion for tapering start and end of series to avoud end-of-signal effects. Defaults to 0.1.
        pad (optional): Pads the provided series before computation, adding pad*(length of series) zeros at the end. Defaults to 0.
        fast (optional): [description]. Defaults to True.
        demean (optional): Demeans series. Defaults to False.
        detrend (optional): Detrends series. Defaults to True.
        minimal (optional): Returns only frequency and spectrum. Overrides option_summary. Defaults to True.
        option_summary (optional): Returns specified options alongside results. Defaults to False.

    Adapted from R's stats::spec.pgram.

    Based on versions at https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/spectrum.R and
    https://github.com/telmo-correa/time-series-analysis/blob/master/Python/spectrum.py
    """

    def spec_taper(x, p=0.1):
        """
        Apply a cosine-bell taper to a time series.
        
        Computes a tapered version of x, with tapering proportion p at each end of x.
    
        Adapted from R's stats::spec.taper.
        """
        p = np.r_[p]
        assert np.all((p >= 0) & (p < 0.5)), "'p' must be between 0 and 0.5"

        x = np.r_[x].astype("float64")
        original_shape = x.shape

        assert len(original_shape) <= 2, "'x' must have at most 2 dimensions"
        while len(x.shape) < 2:
            x = np.expand_dims(x, axis=1)

        nrow, ncol = x.shape
        if len(p) == 1:
            p = p * np.ones(ncol)
        else:
            assert (
                len(p) == nc
            ), "length of 'p' must be 1 or equal the number of columns of 'x'"

        for i in range(ncol):
            m = int(np.floor(nrow * p[i]))
            if m == 0:
                continue
            w = 0.5 * (1 - np.cos(np.pi * np.arange(1, 2 * m, step=2) / (2 * m)))
            x[:, i] = np.r_[w, np.ones(nrow - 2 * m), w[::-1]] * x[:, i]

        x = np.reshape(x, original_shape)
        return x

    def daniell_window_modified(m):
        """ Single-pass modified Daniell kernel window.
        
        Weight is normalized to add up to 1, and all values are the same, other than the first and the
        last, which are divided by 2.
        """

        def w(k):
            return np.where(
                np.abs(k) < m, 1 / (2 * m), np.where(np.abs(k) == m, 1 / (4 * m), 0)
            )

        return w(np.arange(-m, m + 1))

    def daniell_window_convolve(v):
        """ Convolved version of multiple modified Daniell kernel windows.
        
        Parameter v should be an iterable of m values.
        """

        if len(v) == 0:
            return np.r_[1]

        if len(v) == 1:
            return daniell_window_modified(v[0])

        return signal.convolve(
            daniell_window_modified(v[0]), daniell_window_convolve(v[1:])
        )

    # Ensure we can store non-integers in x, and that it is a numpy object
    x = np.r_[x].astype("float64")
    original_shape = x.shape

    # Ensure correct dimensions
    assert len(original_shape) <= 2, "'x' must have at most 2 dimensions"
    while len(x.shape) < 2:
        x = np.expand_dims(x, axis=1)

    # N/N0 = number of rows
    N, nser = x.shape
    N0 = N

    # Ensure only one of spans, kernel is provided, and build the kernel window if needed
    assert (spans is None) or (
        kernel is None
    ), "must specify only one of 'spans' or 'kernel'"
    if spans is not None:
        kernel = daniell_window_convolve(np.floor_divide(np.r_[spans], 2))

    # Detrend and/or demean the series
    if detrend:
        t = np.arange(N) - (N - 1) / 2
        sumt2 = N * (N ** 2 - 1) / 12
        x -= (
            np.repeat(np.expand_dims(np.mean(x, axis=0), 0), N, axis=0)
            + np.outer(np.sum(x.T * t, axis=1), t / sumt2).T
        )
    elif demean:
        x -= np.mean(x, axis=0)

    # Compute taper and taper adjustment variables
    x = spec_taper(x, taper)
    u2 = 1 - (5 / 8) * taper * 2
    u4 = 1 - (93 / 128) * taper * 2

    # Pad the series with copies of the same shape, but filled with zeroes
    if pad > 0:
        x = np.r_[x, np.zeros((np.int(pad * x.shape[0]), x.shape[1]))]
        N = x.shape[0]

    # Further pad the series to accelerate FFT computation
    if fast:
        newN = fft.next_fast_len(N, True)
        x = np.r_[x, np.zeros((newN - N, x.shape[1]))]
        N = newN

    # Compute the Fourier frequencies (R's spec.pgram convention style)
    Nspec = int(np.floor(N / 2))
    freq = (np.arange(Nspec) + 1) * xfreq / N

    # Translations to keep same row / column convention as stats::mvfft
    xfft = fft.fft(x.T).T

    # Compute the periodogram for each i, j
    pgram = np.empty((N, nser, nser), dtype="complex")
    for i in range(nser):
        for j in range(nser):
            pgram[:, i, j] = xfft[:, i] * np.conj(xfft[:, j]) / (N0 * xfreq)
            pgram[0, i, j] = 0.5 * (pgram[1, i, j] + pgram[-1, i, j])

    if kernel is None:
        # Values pre-adjustment
        df = 2
        bandwidth = np.sqrt(1 / 12)
    else:

        def conv_circular(signal, kernel):
            """
            Performs 1D circular convolution, in the same style as R::kernapply,
            assuming the kernel window is centered at 0.
            """
            pad = len(signal) - len(kernel)
            half_window = int((len(kernel) + 1) / 2)
            indexes = range(-half_window, len(signal) - half_window)
            orig_conv = np.real(
                fft.ifft(fft.fft(signal) * fft.fft(np.r_[np.zeros(pad), kernel]))
            )
            return orig_conv.take(indexes, mode="wrap")

        # Convolve pgram with kernel with circular conv
        for i in range(nser):
            for j in range(nser):
                pgram[:, i, j] = conv_circular(pgram[:, i, j], kernel)

        df = 2 / np.sum(kernel ** 2)
        m = (len(kernel) - 1) / 2
        k = np.arange(-m, m + 1)
        bandwidth = np.sqrt(np.sum((1 / 12 + k ** 2) * kernel))

    df = df / (u4 / u2 ** 2) * (N0 / N)
    bandwidth = bandwidth * xfreq / N

    # Remove padded results
    pgram = pgram[1 : (Nspec + 1), :, :]

    spec = np.empty((Nspec, nser))
    for i in range(nser):
        spec[:, i] = np.real(pgram[:, i, i])

    if minimal == True:
        spec = spec / u2
        spec = spec.squeeze()

        results = {
            "freq": freq,
            "spec": spec,
        }

    else:
        if nser == 1:
            coh = None
            phase = None
        else:
            coh = np.empty((Nspec, int(nser * (nser - 1) / 2)))
            phase = np.empty((Nspec, int(nser * (nser - 1) / 2)))
            for i in range(nser):
                for j in range(i + 1, nser):
                    index = int(i + j * (j - 1) / 2)
                    coh[:, index] = np.abs(pgram[:, i, j]) ** 2 / (
                        spec[:, i] * spec[:, j]
                    )
                    phase[:, index] = np.angle(pgram[:, i, j])

        spec = spec / u2
        spec = spec.squeeze()

        if option_summary == True:
            results = {
                "freq": freq,
                "spec": spec,
                "coherency": coh,
                "phase": phase,
                "kernel": kernel,
                "degrees of freedom": df,
                "bandwidth": bandwidth,
                "n.used": N,
                "orig.n": N0,
                "taper": taper,
                "pad": pad,
                "detrend": detrend,
                "demean": demean,
                "method": "Raw Periodogram"
                if kernel is None
                else "Smoothed Periodogram",
            }
        else:
            results = {
                "freq": freq,
                "spec": spec,
                "coherency": coh,
                "phase": phase,
                "kernel": kernel,
                "degrees of freedom": df,
                "bandwidth": bandwidth,
                "n.used": N,
                "orig.n": N0,
            }

    return results
