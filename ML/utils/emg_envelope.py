"""
emg_envelope.py
---------------
EMG envelope extraction pipeline for sEMG-to-torque modelling.

Implements the standard 4-step processing pipeline:
    1. Bandpass filter  — remove motion artefacts and high-frequency noise
    2. Full-wave rectify — |x(t)|
    3. Low-pass filter  — extract the linear envelope
    4. Normalize ÷ max  — scale to [0, 1] using the global max across trials

"""

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt


# Default parameters 

# Known EMG column names across different FLB files:
#   Default naming:  gm, gl, sol, ta
#   File-header naming: ta_emg, mg_emg, sol_emg, lg_emg
# The user can always override by passing emg_columns explicitly.
DEFAULT_EMG_COLUMNS = ['gm', 'gl', 'sol', 'ta']
DEFAULT_FS = 1000.0        # Hz — sampling frequency used in REKLAB experiments
DEFAULT_BP_LOW = 30.0      # Hz — bandpass lower cutoff
DEFAULT_BP_HIGH = 300.0    # Hz — bandpass upper cutoff
DEFAULT_LP_CUTOFF = 2.0    # Hz — lowpass cutoff for envelope (Can also use 1.5 Hz)
DEFAULT_ORDER = 4          # Butterworth filter order

_EMG_KEYWORDS = {'emg', 'gm', 'gl', 'sol', 'ta', 'mg', 'lg'}
_NON_EMG_NAMES = {'time', 'pos', 'tq', 'position', 'torque'}


def detect_emg_columns(df):
    """
    Auto-detect EMG column names in a trial DataFrame.

    Identifies columns by checking for known EMG-related keywords
    (e.g. 'emg', 'gm', 'sol') while excluding non-EMG columns
    (e.g. 'time', 'pos', 'torque').

    Parameters
    ----------
    df : pd.DataFrame
        A trial DataFrame as returned by ``read_flb``.

    Returns
    -------
    list of str
        Column names identified as EMG channels.
    """
    emg_cols = []
    for col in df.columns:
        col_lower = col.lower()
        if col_lower in _NON_EMG_NAMES:
            continue
        if any(kw in col_lower for kw in _EMG_KEYWORDS):
            emg_cols.append(col)
    return emg_cols



# Single-signal envelope extraction

def extract_envelope(signal, fs=DEFAULT_FS,
                     bp_low=DEFAULT_BP_LOW, bp_high=DEFAULT_BP_HIGH,
                     lp_cutoff=DEFAULT_LP_CUTOFF, order=DEFAULT_ORDER,
                     demean=True,
                     return_intermediates=False):
    """
    Extract the linear envelope from a single raw EMG signal.

    Steps:
        0. Demean (remove DC offset) — matches MATLAB ``demean.m``
        1. Bandpass filter (zero-phase Butterworth)
        2. Full-wave rectification  |x(t)|
        3. Lowpass filter (zero-phase Butterworth)

    Parameters
    ----------
    signal : array-like
        1-D array of raw EMG samples.
    fs : float
        Sampling frequency in Hz (default 1000).
    bp_low : float
        Bandpass lower cutoff in Hz (default 30).
    bp_high : float
        Bandpass upper cutoff in Hz (default 300).
    lp_cutoff : float
        Lowpass cutoff in Hz for envelope extraction (default 2).
    order : int
        Butterworth filter order (default 4).
    demean : bool
        If True, subtract mean before filtering (default True).
        Required for wireless EMG sensors with DC offset.
    return_intermediates : bool
        If True, also return the bandpass-filtered+rectified signal
        (``rectified``) alongside the envelope.

    Returns
    -------
    envelope : np.ndarray
        Un-normalized linear envelope (same length as input).
    rectified : np.ndarray
        Bandpass-filtered and rectified signal — only returned when
        ``return_intermediates=True``.
    """
    signal = np.asarray(signal, dtype=np.float64)

    # Step 0: Remove DC offset (matches MATLAB demean.m: y = x - mean(x))
    if demean:
        signal = signal - np.mean(signal)

    nyq = 0.5 * fs

    # Step 1: Bandpass filter
    sos_bp = butter(order, [bp_low / nyq, bp_high / nyq], btype='band', output='sos')
    filtered = sosfiltfilt(sos_bp, signal)

    # Step 2: Full-wave rectification
    rectified = np.abs(filtered)

    # Step 3: Lowpass filter → linear envelope
    sos_lp = butter(order, lp_cutoff / nyq, btype='low', output='sos')
    envelope = sosfiltfilt(sos_lp, rectified)

    # Clip: sosfiltfilt can produce tiny negative values at edges after rectification
    envelope = np.clip(envelope, 0, None)

    if return_intermediates:
        return envelope, rectified
    return envelope


# Batch processing + normalization

def compute_mvc_max(mvc_trials, emg_columns=None, fs=DEFAULT_FS,
                    bp_low=DEFAULT_BP_LOW, bp_high=DEFAULT_BP_HIGH,
                    lp_cutoff=DEFAULT_LP_CUTOFF, order=DEFAULT_ORDER,
                    demean=True):
    """
    Compute per-channel envelope max values from MVC trials.

    Extracts the linear envelope from each MVC trial and returns the peak
    value per channel across all MVC trials.  Pass the result as
    ``max_values`` to ``process_trials`` or ``normalize_with_max`` so that
    all normalization is relative to MVC rather than to the trial set itself.

    Parameters
    ----------
    mvc_trials : list of pd.DataFrame
        MVC trial DataFrames as returned by ``read_flb``.
    emg_columns : list of str, optional
        EMG column names to process (default: auto-detected).
    fs : float
        Sampling frequency in Hz.
    bp_low, bp_high : float
        Bandpass cutoffs in Hz.
    lp_cutoff : float
        Lowpass cutoff in Hz for envelope extraction.
    order : int
        Butterworth filter order.

    Returns
    -------
    max_values : dict
        ``{column_name: mvc_peak}`` — peak envelope value per channel across
        all MVC trials.  Pass directly to ``process_trials`` or
        ``normalize_with_max``.
    """
    if emg_columns is None:
        emg_columns = detect_emg_columns(mvc_trials[0])

    max_values = {}
    for col in emg_columns:
        col_max = 0.0
        for df in mvc_trials:
            env = extract_envelope(df[col].values, fs=fs,
                                   bp_low=bp_low, bp_high=bp_high,
                                   lp_cutoff=lp_cutoff, order=order,
                                   demean=demean)
            col_max = max(col_max, env.max())
        max_values[col] = col_max if col_max > 0 else 1.0

    print(f"MVC max computed from {len(mvc_trials)} MVC trial(s):")
    for col in emg_columns:
        print(f"  {col} MVC max: {max_values[col]:.6f}")

    return max_values


def process_trials(trials, emg_columns=None, max_values=None, fs=DEFAULT_FS,
                   bp_low=DEFAULT_BP_LOW, bp_high=DEFAULT_BP_HIGH,
                   lp_cutoff=DEFAULT_LP_CUTOFF, order=DEFAULT_ORDER,
                   demean=True):
    """
    Process all trials for one subject: extract envelopes and normalize.

    For each EMG column in each trial DataFrame, this function:
        1. Extracts the linear envelope via ``extract_envelope``
        2. Determines the normalization max per channel — either from
           pre-computed MVC max values (``max_values`` argument) or from the
           global max across the provided trials
        3. Normalizes each envelope to [0, 1]

    New columns are added to each DataFrame:
        - ``{col}_env``      — un-normalized envelope
        - ``{col}_env_norm`` — normalized envelope [0, 1]

    Parameters
    ----------
    trials : list of pd.DataFrame
        Trial DataFrames as returned by ``read_flb``.  Each must contain the
        columns listed in *emg_columns*.
    emg_columns : list of str, optional
        EMG column names to process (default: auto-detected).
    max_values : dict, optional
        Pre-computed per-channel max values (e.g. from ``compute_mvc_max``).
        If provided, these are used for normalization instead of the max
        computed across *trials*.
    fs : float
        Sampling frequency in Hz.
    bp_low, bp_high : float
        Bandpass cutoffs in Hz.
    lp_cutoff : float
        Lowpass cutoff in Hz for envelope extraction.
    order : int
        Butterworth filter order.

    Returns
    -------
    trials : list of pd.DataFrame
        Same DataFrames with ``*_env`` and ``*_env_norm`` columns added.
    max_values : dict
        ``{column_name: max_used}`` — the max values actually used for
        normalization (either MVC-derived or computed from trials).
    """
    if emg_columns is None:
        emg_columns = detect_emg_columns(trials[0])

    # Step 1: Extract envelopes
    for df in trials:
        for col in emg_columns:
            env, rect = extract_envelope(df[col].values, fs=fs,
                                         bp_low=bp_low, bp_high=bp_high,
                                         lp_cutoff=lp_cutoff, order=order,
                                         demean=demean,
                                         return_intermediates=True)
            df[f'{col}_rect'] = rect
            df[f'{col}_env'] = env

    # Step 2: Determine max per channel
    if max_values is None:
        # Fall back to global max across the provided trials
        max_values = {}
        for col in emg_columns:
            global_max = max(df[f'{col}_env'].max() for df in trials)
            max_values[col] = global_max if global_max > 0 else 1.0

    # Step 3: Normalize
    for df in trials:
        for col in emg_columns:
            df[f'{col}_env_norm'] = df[f'{col}_env'] / max_values[col]

    print(f"Processed {len(trials)} trials × {len(emg_columns)} EMG channels")
    print(f"  Bandpass: {bp_low}–{bp_high} Hz | Lowpass: {lp_cutoff} Hz | "
          f"Order: {order}")
    for col in emg_columns:
        print(f"  {col} max used: {max_values[col]:.6f}")

    return trials, max_values


def normalize_with_max(trials, max_values, emg_columns=None, fs=DEFAULT_FS,
                       bp_low=DEFAULT_BP_LOW, bp_high=DEFAULT_BP_HIGH,
                       lp_cutoff=DEFAULT_LP_CUTOFF, order=DEFAULT_ORDER,
                       demean=True):
    """
    Extract envelopes and normalize trials using pre-computed max values.

    Use this to normalize a **test set** with max values obtained from
    ``compute_mvc_max`` (MVC-based normalization) or from ``process_trials``.

    Parameters
    ----------
    trials : list of pd.DataFrame
        Trial DataFrames (e.g. test-set trials).
    max_values : dict
        ``{column_name: max}`` as returned by ``compute_mvc_max`` or
        ``process_trials``.
    emg_columns : list of str, optional
        EMG column names to process (default: ['gm', 'gl', 'sol', 'ta']).
    fs, bp_low, bp_high, lp_cutoff, order :
        Filter parameters (should match those used for training).

    Returns
    -------
    trials : list of pd.DataFrame
        Same DataFrames with ``*_env`` and ``*_env_norm`` columns added.
    """
    if emg_columns is None:
        emg_columns = detect_emg_columns(trials[0])

    for df in trials:
        for col in emg_columns:
            env, rect = extract_envelope(df[col].values, fs=fs,
                                         bp_low=bp_low, bp_high=bp_high,
                                         lp_cutoff=lp_cutoff, order=order,
                                         demean=demean,
                                         return_intermediates=True)
            df[f'{col}_rect'] = rect
            df[f'{col}_env'] = env
            df[f'{col}_env_norm'] = env / max_values[col]

    print(f"Normalized {len(trials)} trials using pre-computed max values")

    return trials
