#!/usr/bin/env python3
"""
Apply time lag correction to an EDF file.

This script adjusts the start_time metadata and trims signal data to align
with whole-second boundaries (EDF limitation). Used after estimate_lag.py
to create time-aligned EDF files.

Example:
    # Find lag between two files
    python estimate_lag.py file1.edf file2.edf --channel F7 --quiet
    # Output: -15.234

    # Apply correction to file2
    python apply_lag_correction.py file2.edf -15.234 -o file2_corrected.edf
"""

import argparse
import sys
from datetime import timedelta
import numpy as np
import pyedflib


def apply_lag_correction(input_path, lag_seconds, output_path, verbose=False):
    """
    Apply lag correction to an EDF file.

    Parameters
    ----------
    input_path : str
        Path to input EDF file
    lag_seconds : float
        Lag correction in seconds (from estimate_lag.py)
        Negative: file is ahead, shift backward
        Positive: file is behind, shift forward
    output_path : str
        Path to output corrected EDF file
    verbose : bool
        Print detailed information
    """
    if verbose:
        print(f"Loading {input_path}...")

    # Load input EDF
    f_in = pyedflib.EdfReader(input_path)

    try:
        n_channels = f_in.signals_in_file

        # Read metadata
        old_start_time = f_in.getStartdatetime()

        # Read all signals and their metadata
        signals = []
        headers = []

        for i in range(n_channels):
            signal = f_in.readSignal(i)
            signals.append(signal)

            # Get signal header info
            header = {
                'label': f_in.getLabel(i),
                'dimension': f_in.getPhysicalDimension(i),
                'sample_frequency': f_in.getSampleFrequency(i),
                'physical_max': f_in.getPhysicalMaximum(i),
                'physical_min': f_in.getPhysicalMinimum(i),
                'digital_max': f_in.getDigitalMaximum(i),
                'digital_min': f_in.getDigitalMinimum(i),
                'prefilter': f_in.getPrefilter(i),
                'transducer': f_in.getTransducer(i)
            }
            headers.append(header)

        # Read annotations
        # readAnnotations() returns (onsets, durations, labels) as three separate lists
        annotation_onsets, annotation_durations, annotation_labels = f_in.readAnnotations()
        annotations = list(zip(annotation_onsets, annotation_durations, annotation_labels))

    finally:
        f_in.close()

    # Calculate new start time
    new_start_time_exact = old_start_time + timedelta(seconds=lag_seconds)

    # Extract fractional seconds (EDF only supports 1-second resolution)
    fractional_seconds = new_start_time_exact.microsecond / 1e6

    # Truncate to whole seconds for EDF header
    new_start_time_edf = new_start_time_exact.replace(microsecond=0)

    if verbose:
        print(f"\nTime adjustment:")
        print(f"  Original start time: {old_start_time}")
        print(f"  Lag correction: {lag_seconds:+.3f} s")
        print(f"  New start time (exact): {new_start_time_exact}")
        print(f"  New start time (EDF): {new_start_time_edf}")
        print(f"  Fractional seconds: {fractional_seconds:.6f} s")

    # Trim signals to align with whole-second boundary
    # EDF can only represent integer seconds, so we skip forward by the fractional part
    trimmed_signals = []
    max_trim = 0

    for i, signal in enumerate(signals):
        fs = headers[i]['sample_frequency']
        samples_to_trim = round(fractional_seconds * fs)

        if samples_to_trim > len(signal):
            print(f"Error: Cannot trim {samples_to_trim} samples from channel {headers[i]['label']} (only {len(signal)} samples)",
                  file=sys.stderr)
            sys.exit(1)

        trimmed_signal = signal[samples_to_trim:]
        trimmed_signals.append(trimmed_signal)

        max_trim = max(max_trim, samples_to_trim)

        if verbose:
            print(f"  Channel {headers[i]['label']}: trimmed {samples_to_trim} samples ({samples_to_trim/fs:.6f} s)")

    # Adjust annotations
    # Annotations are referenced to start_time, so we recalculate based on absolute times
    adjusted_annotations = []

    for onset, duration, label in annotations:
        # Original absolute time
        abs_time_original = old_start_time + timedelta(seconds=onset)

        # Corrected absolute time
        abs_time_corrected = abs_time_original + timedelta(seconds=lag_seconds)

        # New onset relative to new start_time
        new_onset = (abs_time_corrected - new_start_time_edf).total_seconds()

        # Only keep annotations that fall within the corrected recording
        if new_onset >= 0:
            adjusted_annotations.append((new_onset, duration, label))

            if verbose and len(adjusted_annotations) <= 5:
                print(f"  Annotation '{label}': {onset:.3f}s → {new_onset:.3f}s")

    if verbose and len(annotations) > 0:
        dropped = len(annotations) - len(adjusted_annotations)
        if dropped > 0:
            print(f"  Dropped {dropped} annotations that fell before new start time")
        print(f"  Total annotations: {len(adjusted_annotations)}")

    # Write corrected EDF file
    if verbose:
        print(f"\nWriting corrected EDF to {output_path}...")

    f_out = pyedflib.EdfWriter(output_path, n_channels, file_type=pyedflib.FILETYPE_EDFPLUS)

    try:
        # Set start time
        f_out.setStartdatetime(new_start_time_edf)

        # Set signal headers
        for i in range(n_channels):
            f_out.setSignalHeader(i, headers[i])

        # Write signal data
        f_out.writeSamples(trimmed_signals)

        # Write annotations
        for onset, duration, label in adjusted_annotations:
            f_out.writeAnnotation(onset, duration, label)

    finally:
        f_out.close()

    if verbose:
        print(f"\nSuccess! Corrected file saved to {output_path}")
        print(f"Total samples trimmed: {max_trim}")
        print(f"Correction error: {abs(samples_to_trim/fs - fractional_seconds):.6f} s (due to integer sample trimming)")


def main():
    parser = argparse.ArgumentParser(
        description='Apply time lag correction to an EDF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Find lag using estimate_lag.py, then apply correction
  LAG=$(python estimate_lag.py file1.edf file2.edf --channel F7 --quiet)
  python apply_lag_correction.py file2.edf $LAG -o file2_corrected.edf

  # Apply negative lag (file is ahead, shift backward)
  python apply_lag_correction.py recording.edf -15.234 -o recording_corrected.edf

  # Apply positive lag (file is behind, shift forward)
  python apply_lag_correction.py recording.edf 8.567 -o recording_corrected.edf

  # Verbose output
  python apply_lag_correction.py input.edf -15.234 -o output.edf --verbose

Notes:
  - EDF format only supports 1-second resolution for start_time
  - Fractional seconds are handled by trimming signal samples
  - Annotations are adjusted to maintain absolute timing
  - Lag value typically comes from estimate_lag.py output
        """
    )

    parser.add_argument('input', help='Input EDF file')
    parser.add_argument('lag', type=float,
                       help='Lag correction in seconds (negative=shift backward, positive=shift forward)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output corrected EDF file')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed processing information')

    args = parser.parse_args()

    try:
        apply_lag_correction(args.input, args.lag, args.output, args.verbose)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
