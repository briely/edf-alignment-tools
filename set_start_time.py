#!/usr/bin/env python3
"""
Set or adjust the start time of an EDF file.

This script can either:
1. Set the start time to an absolute datetime
2. Adjust the start time by a lag offset (e.g., from estimate_lag.py)

In both cases, the script adjusts start_time metadata and trims signal data to align
with whole-second boundaries (EDF limitation).

Examples:
    # Set absolute start time
    python set_start_time.py input.edf --set-time "2024-03-15 14:30:45" -o output.edf

    # Adjust by lag offset (from estimate_lag.py)
    python set_start_time.py input.edf --adjust-by -15.234 -o output.edf
"""

import argparse
import sys
from datetime import datetime, timedelta
import numpy as np
import pyedflib


def set_start_time(input_path, new_start_time, output_path, verbose=False):
    """
    Set the start time of an EDF file to an absolute datetime.

    Parameters
    ----------
    input_path : str
        Path to input EDF file
    new_start_time : datetime
        New start time to set
    output_path : str
        Path to output EDF file
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
        annotation_onsets, annotation_durations, annotation_labels = f_in.readAnnotations()
        annotations = list(zip(annotation_onsets, annotation_durations, annotation_labels))

    finally:
        f_in.close()

    # Calculate the effective lag
    lag_seconds = (new_start_time - old_start_time).total_seconds()

    # Extract fractional seconds (EDF only supports 1-second resolution)
    fractional_seconds = new_start_time.microsecond / 1e6

    # Truncate to whole seconds for EDF header
    new_start_time_edf = new_start_time.replace(microsecond=0)

    if verbose:
        print(f"\nTime adjustment:")
        print(f"  Original start time: {old_start_time}")
        print(f"  New start time (exact): {new_start_time}")
        print(f"  New start time (EDF): {new_start_time_edf}")
        print(f"  Effective adjustment: {lag_seconds:+.3f} s")
        print(f"  Fractional seconds: {fractional_seconds:.6f} s")

    # Trim signals to align with whole-second boundary
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
        if max_trim > 0:
            print(f"Trimming error: {abs(samples_to_trim/fs - fractional_seconds):.6f} s (due to integer sample trimming)")


def adjust_start_time(input_path, lag_seconds, output_path, verbose=False):
    """
    Adjust the start time of an EDF file by a lag offset.

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

    # Load input EDF to get current start time
    f_in = pyedflib.EdfReader(input_path)
    try:
        old_start_time = f_in.getStartdatetime()
    finally:
        f_in.close()

    # Calculate new start time
    new_start_time = old_start_time + timedelta(seconds=lag_seconds)

    # Use the set_start_time function
    set_start_time(input_path, new_start_time, output_path, verbose)


def main():
    parser = argparse.ArgumentParser(
        description='Set or adjust the start time of an EDF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Set absolute start time
  python set_start_time.py input.edf --set-time "2024-03-15 14:30:45" -o output.edf

  # Set start time with microsecond precision
  python set_start_time.py input.edf --set-time "2024-03-15 14:30:45.123456" -o output.edf

  # Adjust by lag offset (from estimate_lag.py)
  LAG=$(python estimate_lag.py file1.edf file2.edf --channel F7 --quiet)
  python set_start_time.py file2.edf --adjust-by $LAG -o file2_corrected.edf

  # Apply negative lag (file is ahead, shift backward)
  python set_start_time.py recording.edf --adjust-by -15.234 -o recording_corrected.edf

  # Apply positive lag (file is behind, shift forward)
  python set_start_time.py recording.edf --adjust-by 8.567 -o recording_corrected.edf

  # Verbose output
  python set_start_time.py input.edf --adjust-by -15.234 -o output.edf --verbose

Notes:
  - EDF format only supports 1-second resolution for start_time
  - Fractional seconds are handled by trimming signal samples
  - Annotations are adjusted to maintain absolute timing
  - Must use either --set-time OR --adjust-by (mutually exclusive)
        """
    )

    parser.add_argument('input', help='Input EDF file')

    # Mutually exclusive group for setting vs adjusting time
    time_group = parser.add_mutually_exclusive_group(required=True)
    time_group.add_argument('--set-time', type=str, metavar='DATETIME',
                           help='Set start time to absolute datetime (format: "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD HH:MM:SS.ffffff")')
    time_group.add_argument('--adjust-by', type=float, metavar='SECONDS',
                           help='Adjust start time by lag in seconds (negative=shift backward, positive=shift forward)')

    parser.add_argument('-o', '--output', required=True,
                       help='Output EDF file')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed processing information')

    args = parser.parse_args()

    try:
        if args.set_time:
            # Parse the datetime string
            # Try with microseconds first, then without
            try:
                new_start_time = datetime.strptime(args.set_time, "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                try:
                    new_start_time = datetime.strptime(args.set_time, "%Y-%m-%d %H:%M:%S")
                except ValueError:
                    print(f"Error: Invalid datetime format. Use 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DD HH:MM:SS.ffffff'",
                          file=sys.stderr)
                    sys.exit(1)

            set_start_time(args.input, new_start_time, args.output, args.verbose)

        elif args.adjust_by is not None:
            adjust_start_time(args.input, args.adjust_by, args.output, args.verbose)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
