#!/usr/bin/env python3
"""
Estimate time lag between two EDF files using cross-correlation.

This script compares signals from two EDF recordings (e.g., from different devices
recording the same session) and estimates the constant time offset between them.
"""

import numpy as np
import scipy.signal as signal
import argparse
import sys
from datetime import datetime, timedelta
import json

# Import from existing modules
from edf_viewer import load_edf


def resample_signal(signal_data, original_fs, target_fs):
    """Resample signal to target sample rate using anti-aliased resampling."""
    if original_fs == target_fs:
        return signal_data

    num_samples = int(len(signal_data) * target_fs / original_fs)
    return signal.resample(signal_data, num_samples)


def preprocess_signal(signal_data):
    """Apply z-score normalization to signal."""
    # Z-score normalize
    signal_data = (signal_data - np.mean(signal_data)) / np.std(signal_data)
    return signal_data


def compute_crosscorr(sig1, sig2, fs, max_lag):
    """
    Compute cross-correlation and return lag that maximizes correlation.

    Returns
    -------
    lag : float
        Lag in seconds (negative = sig2 ahead, positive = sig2 behind)
    correlation : float
        Peak correlation coefficient
    corr_function : ndarray
        Full cross-correlation array
    lag_times : ndarray
        Lag time array corresponding to corr_function
    """
    # Cross-correlation (FFT-based)
    corr = signal.correlate(sig2, sig1, mode='full', method='auto')
    corr = corr / len(sig1)  # Normalize

    # Lag array
    lags = np.arange(-len(sig1)+1, len(sig2))
    lag_times = lags / fs

    # Restrict to search range
    mask = np.abs(lag_times) <= max_lag
    corr_search = corr[mask]
    lags_search = lags[mask]
    lag_times_search = lag_times[mask]

    # Find peak
    peak_idx_search = np.argmax(corr_search)
    peak_correlation = corr_search[peak_idx_search]

    # Parabolic interpolation for sub-sample refinement
    if 0 < peak_idx_search < len(corr_search) - 1:
        y1 = corr_search[peak_idx_search - 1]
        y2 = corr_search[peak_idx_search]
        y3 = corr_search[peak_idx_search + 1]

        delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3) if (y1 - 2*y2 + y3) != 0 else 0
        refined_lag = (lags_search[peak_idx_search] + delta) / fs
    else:
        refined_lag = lag_times_search[peak_idx_search]

    return refined_lag, peak_correlation, corr, lag_times


def bootstrap_lag_ci(sig1, sig2, fs, max_lag, window_size=30, n_boot=300):
    """
    Estimate confidence interval for lag using bootstrap resampling.

    Parameters
    ----------
    sig1, sig2 : ndarray
        Preprocessed signals (same length, same sample rate)
    fs : float
        Sample rate
    max_lag : float
        Maximum lag to search (seconds)
    window_size : float
        Window size for segmentation (seconds)
    n_boot : int
        Number of bootstrap iterations

    Returns
    -------
    ci_lower, ci_upper : float
        95% confidence interval bounds
    boot_lags : ndarray
        Array of bootstrap lag estimates
    """
    # Segment signals into windows
    window_samples = int(window_size * fs)
    overlap_samples = window_samples // 2

    # Create sliding windows
    windows1 = []
    windows2 = []
    start = 0
    while start + window_samples <= len(sig1):
        windows1.append(sig1[start:start + window_samples])
        windows2.append(sig2[start:start + window_samples])
        start += overlap_samples

    n_windows = len(windows1)
    if n_windows < 3:
        # Not enough windows for bootstrap, return None
        return None, None, None

    rng = np.random.default_rng(seed=42)
    boot_lags = np.zeros(n_boot)

    for b in range(n_boot):
        # Resample windows with replacement
        boot_indices = rng.choice(n_windows, size=n_windows, replace=True)

        # Concatenate windows
        boot_sig1 = np.concatenate([windows1[i] for i in boot_indices])
        boot_sig2 = np.concatenate([windows2[i] for i in boot_indices])

        # Compute lag for this bootstrap sample
        lag, _, _, _ = compute_crosscorr(boot_sig1, boot_sig2, fs, max_lag)
        boot_lags[b] = lag

    # Compute 95% CI
    ci_lower = np.percentile(boot_lags, 2.5)
    ci_upper = np.percentile(boot_lags, 97.5)

    return ci_lower, ci_upper, boot_lags


def analyze_stability(sig1, sig2, fs, global_lag, window_size=60, max_lag=0.2):
    """
    Analyze stability of correlation over time using sliding windows.

    Parameters
    ----------
    sig1, sig2 : ndarray
        Preprocessed signals (same length, same sample rate)
    fs : float
        Sample rate
    global_lag : float
        The global lag estimate to use as reference (seconds)
    window_size : float
        Window size for analysis (seconds)
    max_lag : float
        Maximum lag deviation to search around global_lag (seconds)

    Returns
    -------
    window_times : ndarray
        Center time of each window (seconds)
    correlations : ndarray
        Peak correlation coefficient for each window
    lag_residuals : ndarray
        Optimal lag - global_lag for each window (seconds)
    """
    window_samples = int(window_size * fs)
    step_samples = window_samples // 2  # 50% overlap

    # We'll extract windows from the unaligned signals and compute the full lag for each
    # Then compute residuals as (window_lag - global_lag)

    # Sliding window analysis
    window_times = []
    correlations = []
    lag_residuals = []

    # Maximum search range for cross-correlation
    # We search in ±(|global_lag| + max_lag) to capture the expected lag ± allowed deviation
    search_limit = abs(global_lag) + max_lag + 1.0

    start = 0
    while start + window_samples <= min(len(sig1), len(sig2)):
        # Extract window from SAME indices in both signals
        win1 = sig1[start:start + window_samples]
        win2 = sig2[start:start + window_samples]

        # Find optimal lag for this window (search in wide range)
        lag_window, corr_window, _, _ = compute_crosscorr(win1, win2, fs, search_limit)

        # Only accept this window if the lag is within max_lag of global_lag
        residual = lag_window - global_lag
        if abs(residual) <= max_lag:
            # Record results
            window_center = (start + window_samples / 2) / fs
            window_times.append(window_center)
            correlations.append(corr_window)
            lag_residuals.append(residual)

        start += step_samples

    return np.array(window_times), np.array(correlations), np.array(lag_residuals)


def plot_stability(window_times, correlations, lag_residuals, results):
    """Generate stability analysis plots."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Plot 1: Correlation coefficient over time
    ax1 = axes[0]
    ax1.plot(window_times / 60, correlations, 'o-', linewidth=1, markersize=3)
    ax1.axhline(results['peak_correlation'], color='red', linestyle='--',
                alpha=0.5, label=f'Global peak: {results["peak_correlation"]:.3f}')
    ax1.set_ylabel('Correlation Coefficient')
    ax1.set_title('Correlation Quality Over Time')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim([max(0, min(correlations) - 0.1), min(1, max(correlations) + 0.1)])

    # Plot 2: Lag residual over time
    ax2 = axes[1]
    ax2.plot(window_times / 60, lag_residuals * 1000, 'o-', linewidth=1, markersize=3, color='orange')
    ax2.axhline(0, color='red', linestyle='--', alpha=0.5, label='Global lag estimate')
    ax2.set_xlabel('Time (minutes)')
    ax2.set_ylabel('Lag Residual (ms)')
    ax2.set_title('Lag Drift Over Time (deviation from global estimate)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Add interpretation text
    lag_std = np.std(lag_residuals) * 1000  # ms
    lag_trend = np.polyfit(window_times, lag_residuals * 1000, 1)[0]  # ms/s

    fig.text(0.5, 0.02,
             f'Lag variability: {lag_std:.2f} ms (std) | '
             f'Drift rate: {lag_trend:.3f} ms/s ({lag_trend * 60:.2f} ms/min)',
             ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.show()


def print_lag_report(results, quiet=False):
    """Print formatted lag estimation report."""
    if quiet:
        print(f"{results['lag_estimate']:.3f}")
        return

    print("=== EDF Lag Estimation ===")
    print()
    print("Files:")
    print(f"  File 1 (reference): {results['file1']}")
    print(f"  File 2 (to align):  {results['file2']}")
    print()
    print(f"Channel: {results['channel1']}", end="")
    if results['channel2'] != results['channel1']:
        print(f" (file 1), {results['channel2']} (file 2)")
    else:
        print(" (both files)")

    print(f"Sample rates: {results['sample_rate1']:.0f} Hz (file 1), ", end="")
    print(f"{results['sample_rate2']:.0f} Hz (file 2)", end="")
    if results['resampled']:
        print(f" → resampled to {results['resampled_rate']:.0f} Hz")
    else:
        print()

    # Show preprocessing info
    print("Preprocessing: z-score normalized")

    print()
    print(f"Time overlap: {results['overlap_duration']:.1f} s ({results['overlap_duration']/60:.1f} min)")
    print(f"Analysis duration: {results['analysis_duration']:.1f} s")
    print(f"Max lag search range: ±{results['max_lag']:.1f} s")
    print()
    print("=== Results ===")
    print()
    print(f"Estimated lag: {results['lag_estimate']:+.3f} s")

    if results['bootstrap_enabled'] and results['ci_lower'] is not None:
        print(f"  95% CI: [{results['ci_lower']:+.3f}, {results['ci_upper']:+.3f}] s")
        print(f"  Bootstrap iterations: {results['n_boot']}")

    print()
    print("Interpretation:")
    if results['lag_estimate'] < 0:
        print(f"  • File 2 is {abs(results['lag_estimate']):.3f} seconds AHEAD of File 1")
        print(f"  • Events in File 2 occur {abs(results['lag_estimate']):.3f} seconds EARLIER")
    elif results['lag_estimate'] > 0:
        print(f"  • File 2 is {results['lag_estimate']:.3f} seconds BEHIND File 1")
        print(f"  • Events in File 2 occur {results['lag_estimate']:.3f} seconds LATER")
    else:
        print("  • Files are already aligned (lag ≈ 0)")

    print(f"  • To align File 2 with File 1, shift by {results['lag_estimate']:+.3f} seconds")

    print()
    print("Cross-correlation quality:")
    print(f"  • Peak correlation: {results['peak_correlation']:.3f} (max possible: 1.0)")

    if results['peak_correlation'] < 0.3:
        print("  ⚠ WARNING: Low correlation - signals may not be related")
    elif results['peak_correlation'] < 0.6:
        print("  ⚠ NOTE: Moderate correlation - consider checking channel/filter settings")

    print()
    print("=== Verification ===")
    print()
    print("To visualize the alignment:")
    print(f"  python edf_viewer.py {results['file1']} {results['file2']} --shift {results['lag_estimate']:.3f} --duration 60")
    print()


def plot_diagnostics(results):
    """Generate diagnostic plots."""
    import matplotlib.pyplot as plt

    # Plot 1: Cross-correlation function
    plt.figure(figsize=(12, 4))
    plt.plot(results['lag_times'], results['corr_function'], linewidth=0.5)
    plt.axvline(results['lag_estimate'], color='red', linestyle='--',
                label=f'Peak lag: {results["lag_estimate"]:.3f}s')
    plt.axhline(0, color='k', linewidth=0.5)
    plt.xlabel('Lag (s)')
    plt.ylabel('Cross-correlation')
    plt.title('Cross-Correlation Function')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    # Plot 2: Bootstrap distribution (if available)
    if results['bootstrap_enabled'] and results['boot_lags'] is not None:
        plt.figure(figsize=(8, 4))
        plt.hist(results['boot_lags'], bins=50, alpha=0.7, edgecolor='black')
        plt.axvline(results['lag_estimate'], color='red', linestyle='--',
                   label=f'Estimate: {results["lag_estimate"]:.3f}s')
        plt.axvline(results['ci_lower'], color='orange', linestyle=':',
                   label=f'95% CI: [{results["ci_lower"]:.3f}, {results["ci_upper"]:.3f}]s')
        plt.axvline(results['ci_upper'], color='orange', linestyle=':')
        plt.xlabel('Lag (s)')
        plt.ylabel('Frequency')
        plt.title(f'Bootstrap Distribution (n={results["n_boot"]})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

    plt.show()


def create_parser():
    """Create argument parser."""
    parser = argparse.ArgumentParser(
        description='Estimate time lag between two EDF files using cross-correlation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (same channel in both files)
  python estimate_lag.py file1.edf file2.edf --channel F7

  # Use different channels from each file (if channel names differ)
  python estimate_lag.py file1.edf file2.edf --channel1 F7 --channel2 F8

  # Disable bootstrap for speed
  python estimate_lag.py file1.edf file2.edf --channel F7 --no-bootstrap

  # Full analysis with plots
  python estimate_lag.py file1.edf file2.edf --channel F7 --plot

  # Analyze correlation stability over time (detect clock drift)
  python estimate_lag.py file1.edf file2.edf --channel F7 --stability-plot

  # Custom stability window size (30s windows)
  python estimate_lag.py file1.edf file2.edf --channel F7 --stability-plot --stability-window 30

  # Tighter stability search range (±100ms) for well-aligned signals
  python estimate_lag.py file1.edf file2.edf --channel F7 --stability-plot --stability-max-lag 0.1

  # Save results to JSON
  python estimate_lag.py file1.edf file2.edf --channel F7 --output results.json

  # Minimal output for scripting
  python estimate_lag.py file1.edf file2.edf --channel F7 --quiet
        """
    )

    parser.add_argument('file1', help='First EDF file (reference)')
    parser.add_argument('file2', help='Second EDF file (to be aligned)')

    # Channel selection (either --channel OR both --channel1 and --channel2 required)
    parser.add_argument('--channel', help='Channel name to use from both files')
    parser.add_argument('--channel1', help='Channel from file 1 (use with --channel2)')
    parser.add_argument('--channel2', help='Channel from file 2 (use with --channel1)')

    # Processing options
    parser.add_argument('--duration', type=float,
                       help='Use only first N seconds (default: all overlap)')
    parser.add_argument('--max-lag', type=float, default=120,
                       help='Maximum lag to search in seconds (default: 120)')

    # Bootstrap options (default: enabled)
    parser.add_argument('--no-bootstrap', action='store_true',
                       help='Disable bootstrap confidence interval (faster)')
    parser.add_argument('--n-boot', type=int, default=300,
                       help='Number of bootstrap iterations (default: 300)')
    parser.add_argument('--window-size', type=float, default=30,
                       help='Bootstrap window size in seconds (default: 30)')

    # Output options
    parser.add_argument('--plot', action='store_true',
                       help='Show diagnostic plots')
    parser.add_argument('--stability-plot', action='store_true',
                       help='Show stability analysis (correlation quality and lag drift over time)')
    parser.add_argument('--stability-window', type=float, default=60,
                       help='Window size for stability analysis in seconds (default: 60)')
    parser.add_argument('--stability-max-lag', type=float, default=0.2,
                       help='Max lag deviation to search in stability analysis in seconds (default: 0.2)')
    parser.add_argument('--quiet', action='store_true',
                       help='Minimal output (just print lag value)')
    parser.add_argument('--output', help='Save results to JSON file')

    return parser


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()

    # Validate channel arguments
    if args.channel:
        # Using same channel for both files
        if args.channel1 or args.channel2:
            print("Error: Cannot use --channel with --channel1/--channel2", file=sys.stderr)
            sys.exit(1)
        channel1 = args.channel
        channel2 = args.channel
    elif args.channel1 and args.channel2:
        # Using different channels for each file
        channel1 = args.channel1
        channel2 = args.channel2
    elif args.channel1 or args.channel2:
        # Only one of channel1/channel2 provided
        print("Error: --channel1 and --channel2 must be used together", file=sys.stderr)
        sys.exit(1)
    else:
        # No channel arguments provided
        print("Error: Must specify either --channel or both --channel1 and --channel2", file=sys.stderr)
        sys.exit(1)

    try:
        # Load EDF files
        if not args.quiet:
            print(f"Loading {args.file1}...")
        data1 = load_edf(args.file1)

        if not args.quiet:
            print(f"Loading {args.file2}...")
        data2 = load_edf(args.file2)

        # Validate channels exist
        if channel1 not in data1['signals']:
            print(f"Error: Channel '{channel1}' not found in {args.file1}", file=sys.stderr)
            print(f"Available channels: {', '.join(data1['signals'].keys())}", file=sys.stderr)
            sys.exit(1)

        if channel2 not in data2['signals']:
            print(f"Error: Channel '{channel2}' not found in {args.file2}", file=sys.stderr)
            print(f"Available channels: {', '.join(data2['signals'].keys())}", file=sys.stderr)
            sys.exit(1)

        # Extract signals and metadata
        sig1 = data1['signals'][channel1]
        sig2 = data2['signals'][channel2]
        fs1 = data1['sample_rates'][channel1]
        fs2 = data2['sample_rates'][channel2]

        # Find time overlap
        start1 = data1['start_time']
        start2 = data2['start_time']
        dur1 = len(sig1) / fs1
        dur2 = len(sig2) / fs2
        end1 = start1 + timedelta(seconds=dur1)
        end2 = start2 + timedelta(seconds=dur2)

        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)

        if overlap_start >= overlap_end:
            print("Error: No time overlap between recordings!", file=sys.stderr)
            print(f"  {args.file1}: {start1} to {end1}", file=sys.stderr)
            print(f"  {args.file2}: {start2} to {end2}", file=sys.stderr)
            sys.exit(1)

        overlap_duration = (overlap_end - overlap_start).total_seconds()

        # Extract overlapping segments
        offset1 = (overlap_start - start1).total_seconds()
        offset2 = (overlap_start - start2).total_seconds()

        idx1_start = int(offset1 * fs1)
        idx2_start = int(offset2 * fs2)
        idx1_end = idx1_start + int(overlap_duration * fs1)
        idx2_end = idx2_start + int(overlap_duration * fs2)

        sig1_overlap = sig1[idx1_start:idx1_end]
        sig2_overlap = sig2[idx2_start:idx2_end]

        # Resample to common rate (lower of the two)
        target_fs = min(fs1, fs2)
        resampled = (fs1 != fs2)

        if fs1 != target_fs:
            sig1_overlap = resample_signal(sig1_overlap, fs1, target_fs)
        if fs2 != target_fs:
            sig2_overlap = resample_signal(sig2_overlap, fs2, target_fs)

        # Truncate to requested duration if specified
        analysis_duration = overlap_duration
        if args.duration and args.duration < overlap_duration:
            n_samples = int(args.duration * target_fs)
            sig1_overlap = sig1_overlap[:n_samples]
            sig2_overlap = sig2_overlap[:n_samples]
            analysis_duration = args.duration

        # Preprocess signals (z-score normalization only)
        sig1_proc = preprocess_signal(sig1_overlap)
        sig2_proc = preprocess_signal(sig2_overlap)

        # Compute cross-correlation
        lag, peak_corr, corr_func, lag_times = compute_crosscorr(
            sig1_proc, sig2_proc, target_fs, args.max_lag
        )

        # Bootstrap confidence interval (if enabled)
        bootstrap_enabled = not args.no_bootstrap
        ci_lower, ci_upper, boot_lags = None, None, None

        if bootstrap_enabled:
            if not args.quiet:
                print("Computing bootstrap confidence interval...")
            ci_lower, ci_upper, boot_lags = bootstrap_lag_ci(
                sig1_proc, sig2_proc, target_fs, args.max_lag,
                args.window_size, args.n_boot
            )

            if ci_lower is None and not args.quiet:
                print("Warning: Not enough data for bootstrap CI (need >3 windows)", file=sys.stderr)

        # Prepare results dictionary
        results = {
            'file1': args.file1,
            'file2': args.file2,
            'channel1': channel1,
            'channel2': channel2,
            'sample_rate1': fs1,
            'sample_rate2': fs2,
            'resampled': resampled,
            'resampled_rate': target_fs,
            'overlap_duration': overlap_duration,
            'analysis_duration': analysis_duration,
            'max_lag': args.max_lag,
            'lag_estimate': lag,
            'peak_correlation': peak_corr,
            'bootstrap_enabled': bootstrap_enabled,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'n_boot': args.n_boot if bootstrap_enabled else None,
            'boot_lags': boot_lags,
            'corr_function': corr_func,
            'lag_times': lag_times,
            'timestamp': datetime.now().isoformat()
        }

        # Print report
        print_lag_report(results, quiet=args.quiet)

        # Save to JSON if requested
        if args.output:
            # Remove non-serializable arrays
            json_results = {k: v for k, v in results.items()
                          if k not in ['boot_lags', 'corr_function', 'lag_times']}
            with open(args.output, 'w') as f:
                json.dump(json_results, f, indent=2)
            if not args.quiet:
                print(f"Results saved to {args.output}")

        # Plot diagnostics if requested
        if args.plot:
            plot_diagnostics(results)

        # Stability analysis if requested
        if args.stability_plot:
            if not args.quiet:
                print("Computing stability analysis...")

            window_times, correlations, lag_residuals = analyze_stability(
                sig1_proc, sig2_proc, target_fs, lag,
                window_size=args.stability_window,
                max_lag=args.stability_max_lag
            )

            if len(window_times) == 0:
                print("Error: No valid windows found in stability analysis!", file=sys.stderr)
                print(f"  The lag in all windows deviated by more than {args.stability_max_lag}s from the global lag.", file=sys.stderr)
                print(f"  Try increasing --stability-max-lag", file=sys.stderr)
            else:
                if not args.quiet:
                    print(f"\nStability Analysis ({args.stability_window}s windows, ±{args.stability_max_lag}s search range):")
                    print(f"  Windows analyzed: {len(window_times)}")
                    print(f"  Correlation range: {np.min(correlations):.3f} to {np.max(correlations):.3f}")
                    print(f"  Lag variability (std): {np.std(lag_residuals) * 1000:.2f} ms")
                    lag_trend = np.polyfit(window_times, lag_residuals * 1000, 1)[0]
                    print(f"  Drift rate: {lag_trend:.3f} ms/s ({lag_trend * 60:.2f} ms/min)")

                    # Warn if residuals are large relative to search range
                    max_residual = np.max(np.abs(lag_residuals))
                    if max_residual > args.stability_max_lag * 0.8:
                        print(f"  ⚠ WARNING: Max lag residual ({max_residual:.3f}s) is close to search range limit")
                        print(f"            Consider increasing --stability-max-lag")

                plot_stability(window_times, correlations, lag_residuals, results)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
