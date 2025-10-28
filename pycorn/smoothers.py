import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def butter_lowpass(cutoff: float, fs: float, order: int = 5)-> tuple[np.ndarray, np.ndarray]:
    """
    Design a low-pass Butterworth filter.

    Parameters
    ----------
    cutoff : float
        Cutoff frequency of the filter [Hz].
    fs : float
        Sampling frequency of the signal [Hz].
    order : int, optional
        Order of the Butterworth filter. Higher values yield a sharper cutoff.
        Default is 5.

    Returns
    -------
    b : ndarray
        Numerator (feedforward) coefficients of the filter.
    a : ndarray
        Denominator (feedback) coefficients of the filter.

    Notes
    -----
    This function uses `scipy.signal.butter` to compute the filter coefficients.
    It is typically used to design filters for chromatographic signal smoothing.
    """
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data: np.ndarray, cutoff: float, fs: float, order: int = 5) -> np.ndarray:
    """
    Apply a low-pass Butterworth filter to smooth chromatographic data.

    Parameters
    ----------
    data : array_like
        Input signal to be filtered (e.g., raw chromatogram intensity values).
    cutoff : float
        Cutoff frequency of the filter [Hz].
    fs : float
        Sampling frequency of the signal [Hz].
    order : int, optional
        Order of the Butterworth filter. Higher values yield a sharper cutoff.
        Default is 5.

    Returns
    -------
    y : ndarray
        Filtered signal of the same shape as `data`.

    Notes
    -----
    This function first designs a Butterworth low-pass filter using
    `butter_lowpass`, then applies it to the data via `scipy.signal.lfilter`.
    It helps remove high-frequency noise while preserving chromatographic peak
    shapes.
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y  


# --- Derivative thresholding ---
def spike_remover_derivative(x, y, k=10):
    """
    Detect and replace sharp single-sample spikes using derivative thresholding.

    Parameters
    ----------
    x : ndarray
        Independent variable (e.g., time or distance).
    y : ndarray
        Dependent variable (signal).
    k : float
        Threshold multiplier for spike detection; higher means less sensitive spike detection.

    Returns
    -------
    y_clean : ndarray
        Signal with spikes replaced by linear interpolation.
    """
    dy = np.diff(y)
    med_abs_dy = np.median(np.abs(dy))
    threshold = k * med_abs_dy
    spike_idx = np.where(np.abs(dy) > threshold)[0]

    y_clean = y.copy()
    mask = np.ones_like(y, dtype=bool)
    mask[spike_idx] = False
    mask[spike_idx + 1] = False  # handle both sides of spike

    interp = interp1d(x[mask], y[mask], kind='linear', fill_value="extrapolate")
    y_clean[~mask] = interp(x[~mask])
    return y_clean


# --- Hampel/MAD-based filtering ---
def spike_remover_hampel(x, y, window_size=5, n_sigma=6):
    """
    Replace outliers based on the Hampel identifier using local Median Absolute Deviation (MAD).

    Parameters
    ----------
    x : ndarray
        Independent variable.
    y : ndarray
        Signal.
    window_size : int
        Half window size for local neighborhood; higher means less sensitive spike detection (larger window, thus the spike must be of greater intensity)
    n_sigma : float
        Threshold multiplier of MAD to flag outliers; higher means less sensitive spike detection (higher threshold for detection)

    Returns
    -------
    y_clean : ndarray
        Signal with outliers replaced by linear interpolation.
    """
    y_clean = y.copy()
    n = len(y)
    for i in range(window_size, n - window_size):
        window = y[i - window_size:i + window_size + 1]
        median = np.median(window)
        mad = np.median(np.abs(window - median))
        if np.abs(y[i] - median) > n_sigma * mad:
            y_clean[i] = np.nan

    nan_mask = np.isnan(y_clean)
    if np.any(nan_mask):
        interp = interp1d(x[~nan_mask], y_clean[~nan_mask], kind='linear', fill_value="extrapolate")
        y_clean[nan_mask] = interp(x[nan_mask])
    return y_clean