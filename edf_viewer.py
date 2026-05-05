#!/usr/bin/env python3
"""CLI tool to load and visualize EDF files."""

import argparse
import sys
import numpy as np
import pyedflib
import matplotlib.pyplot as plt
from typing import Dict, List, Optional


def load_edf(edf_path: str) -> Dict:
    """
    Load all channels from an EDF file.

    Parameters
    ----------
    edf_path : str
        Path to EDF file

    Returns
    -------
    data : dict
        Dictionary containing:
        - 'signals': dict mapping channel labels to signal arrays
        - 'sample_rates': dict mapping channel labels to sample rates
        - 'start_time': datetime of recording start
        - 'annotations': list of (onset, duration, label) tuples
    """
    f = pyedflib.EdfReader(edf_path)

    try:
        n_channels = f.signals_in_file

        # Load all channels
        signals = {}
        sample_rates = {}

        for i in range(n_channels):
            label = f.getLabel(i)
            signal = f.readSignal(i)
            sample_rate = f.getSampleFrequency(i)

            signals[label] = signal
            sample_rates[label] = sample_rate

        # Get recording start time
        start_time = f.getStartdatetime()

        # Get annotations
        annotations = f.readAnnotations()

        # Get header metadata
        # These return bytes in pyedflib, so decode them
        def decode_field(field):
            if isinstance(field, bytes):
                return field.decode('latin-1').strip()
            return str(field).strip() if field else ''

        patient_id = decode_field(f.getPatientName())
        patient_code = decode_field(f.getPatientCode())
        patient_additional = decode_field(f.getPatientAdditional())
        recording_id = decode_field(f.getRecordingAdditional())
        admin_code = decode_field(f.getAdmincode())
        technician = decode_field(f.getTechnician())
        equipment = decode_field(f.getEquipment())

        # Get file type
        file_type_code = f.filetype
        file_type_map = {
            0: 'EDF',
            1: 'EDF+C (Continuous)',
            2: 'EDF+D (Discontinuous)',
            3: 'BDF',
            4: 'BDF+C (Continuous)',
            5: 'BDF+D (Discontinuous)'
        }
        file_type = file_type_map.get(file_type_code, f'Unknown ({file_type_code})')

        return {
            'signals': signals,
            'sample_rates': sample_rates,
            'start_time': start_time,
            'annotations': annotations,
            'metadata': {
                'file_type': file_type,
                'file_type_code': file_type_code,
                'patient_id': patient_id,
                'patient_code': patient_code,
                'patient_additional': patient_additional,
                'recording_id': recording_id,
                'admin_code': admin_code,
                'technician': technician,
                'equipment': equipment
            }
        }
    finally:
        f.close()


def plot_edf_channels(
    data: Dict,
    duration_s: Optional[float] = None,
    offset_s: float = 0.0,
    channels: Optional[List[str]] = None,
    show_annotations: bool = True,
    annotations_all_channels: bool = False,
    custom_title: Optional[str] = None
) -> None:
    """
    Plot EDF channels as time series.

    Parameters
    ----------
    data : dict
        Data dictionary from load_edf()
    duration_s : float, optional
        Duration to plot in seconds (default: plot all)
    offset_s : float
        Start offset in seconds (default: 0)
    channels : list of str, optional
        Specific channels to plot (default: all channels)
    show_annotations : bool
        Whether to show annotation markers (default: True)
    annotations_all_channels : bool
        Whether to show annotations on all subplots (default: False, only first)
    custom_title : str, optional
        Custom title for the plot (default: auto-generated)
    """
    signals = data['signals']
    sample_rates = data['sample_rates']
    annotations = data['annotations']

    # Determine which channels to plot
    if channels is None:
        channel_list = list(signals.keys())
    else:
        # Validate requested channels exist
        channel_list = []
        for ch in channels:
            if ch in signals:
                channel_list.append(ch)
            else:
                print(f"Warning: Channel '{ch}' not found, skipping")

        if len(channel_list) == 0:
            print("Error: No valid channels to plot")
            return

    # Create subplots
    n_plots = len(channel_list)
    fig, axes = plt.subplots(n_plots, 1, figsize=(14, 2 * n_plots), sharex=True)

    if n_plots == 1:
        axes = [axes]

    # Plot each channel
    for i, ch_label in enumerate(channel_list):
        ax = axes[i]
        signal = signals[ch_label]
        fs = sample_rates[ch_label]

        # Create time axis
        time = np.arange(len(signal)) / fs

        # Apply offset and duration constraints
        start_idx = int(offset_s * fs)
        if duration_s is not None:
            end_idx = int((offset_s + duration_s) * fs)
        else:
            end_idx = len(signal)

        # Ensure indices are within bounds
        start_idx = max(0, min(start_idx, len(signal) - 1))
        end_idx = max(start_idx + 1, min(end_idx, len(signal)))

        time_plot = time[start_idx:end_idx]
        signal_plot = signal[start_idx:end_idx]

        # Plot signal
        ax.plot(time_plot, signal_plot, linewidth=0.5, color='black')
        ax.set_ylabel(f"{ch_label}\n({fs:.0f} Hz)")
        ax.grid(True, alpha=0.3)

        # Add annotations if requested
        # Show on first subplot only (default), or all subplots if annotations_all_channels=True
        if show_annotations and (annotations_all_channels or i == 0) and len(annotations[0]) > 0:
            for onset, duration, label in zip(*annotations):
                if offset_s <= onset <= (offset_s + (duration_s or time[-1])):
                    ax.axvline(onset, color='red', alpha=0.5, linestyle='--', linewidth=1)
                    ax.text(onset, ax.get_ylim()[1], label,
                           rotation=90, verticalalignment='top', fontsize=8)

    # Set x-axis label
    axes[-1].set_xlabel("Time (s)")

    # Set title
    if custom_title:
        axes[0].set_title(custom_title)
    else:
        start_time_str = data['start_time'].strftime('%Y-%m-%d %H:%M:%S')
        duration_str = f"{duration_s}s" if duration_s else "all"
        axes[0].set_title(f"EDF Recording (Start: {start_time_str}, Duration: {duration_str})")

    plt.tight_layout()
    plt.show()


def plot_multiple_edf_channels(
    data_list: List[Dict],
    file_names: List[str],
    duration_s: Optional[float] = None,
    offset_s: float = 0.0,
    channels: Optional[List[str]] = None,
    show_annotations: bool = True,
    time_shifts: Optional[List[float]] = None,
    annotations_all_channels: bool = False,
    custom_title: Optional[str] = None
) -> None:
    """
    Plot multiple EDF files overlaid on the same axes, aligned by recording timestamp.

    Parameters
    ----------
    data_list : list of dict
        List of data dictionaries from load_edf()
    file_names : list of str
        List of file names for legend labels
    duration_s : float, optional
        Duration to plot in seconds (default: plot overlapping region)
    offset_s : float
        Start offset in seconds from the earliest recording start (default: 0)
    channels : list of str, optional
        Specific channels to plot (default: all unique channels across all files)
    show_annotations : bool
        Whether to show annotation markers (default: True)
    time_shifts : list of float, optional
        Manual time shifts in seconds for each file (starting from second file).
        Negative values shift backward (events occur earlier).
    annotations_all_channels : bool
        Whether to show annotations on all subplots with labels (default: False, only first)
    custom_title : str, optional
        Custom title for the plot (default: auto-generated)
    """
    if len(data_list) == 0:
        print("Error: No data to plot")
        return

    # Prepare time shifts (first file always has 0 shift)
    if time_shifts is None:
        time_shifts = [0.0] * len(data_list)
    else:
        # Prepend 0 for first file if not included
        if len(time_shifts) == len(data_list) - 1:
            time_shifts = [0.0] + list(time_shifts)
        elif len(time_shifts) != len(data_list):
            raise ValueError(f"Number of time shifts ({len(time_shifts)}) must match number of files ({len(data_list)}) or be one less")

    # Calculate absolute time ranges for each recording (with manual shifts applied)
    import datetime
    recording_info = []
    for idx, (data, file_name, shift) in enumerate(zip(data_list, file_names, time_shifts)):
        start_time = data['start_time'] + datetime.timedelta(seconds=shift)
        # Get duration from first channel
        first_channel = list(data['signals'].keys())[0]
        duration = len(data['signals'][first_channel]) / data['sample_rates'][first_channel]
        end_time = start_time + datetime.timedelta(seconds=duration)

        recording_info.append({
            'file_name': file_name,
            'start_time': start_time,
            'end_time': end_time,
            'duration': duration,
            'shift': shift
        })

    # Find the overlapping time interval
    earliest_start = min(info['start_time'] for info in recording_info)
    latest_start = max(info['start_time'] for info in recording_info)
    earliest_end = min(info['end_time'] for info in recording_info)
    latest_end = max(info['end_time'] for info in recording_info)

    # Check if there's any overlap
    overlap_start = latest_start
    overlap_end = earliest_end

    if overlap_start >= overlap_end:
        print("Error: No overlapping time interval between recordings!")
        print("\nRecording time ranges:")
        for info in recording_info:
            print(f"  {info['file_name']}: {info['start_time']} to {info['end_time']} ({info['duration']:.1f}s)")
        return

    # Calculate time range to plot
    # Use earliest start as reference (t=0)
    overlap_duration = (overlap_end - overlap_start).total_seconds()
    total_span = (latest_end - earliest_start).total_seconds()

    print(f"\nTime alignment:")
    print(f"  Reference time (t=0): {earliest_start}")
    print(f"  Total time span: {total_span:.1f}s")
    print(f"  Overlapping interval: {overlap_duration:.1f}s ({overlap_start} to {overlap_end})")

    for info in recording_info:
        offset = (info['start_time'] - earliest_start).total_seconds()
        shift_str = f" (shift: {info['shift']:+.1f}s)" if info['shift'] != 0 else ""
        print(f"  {info['file_name']}: starts at t={offset:.1f}s, duration={info['duration']:.1f}s{shift_str}")

    # Collect all unique channel names across all files
    all_channels = set()
    for data in data_list:
        all_channels.update(data['signals'].keys())

    # Determine which channels to plot
    if channels is None:
        channel_list = sorted(list(all_channels))
    else:
        channel_list = []
        for ch in channels:
            if ch in all_channels:
                channel_list.append(ch)
            else:
                print(f"Warning: Channel '{ch}' not found in any file, skipping")

        if len(channel_list) == 0:
            print("Error: No valid channels to plot")
            return

    # Create subplots
    n_plots = len(channel_list)
    fig, axes = plt.subplots(n_plots, 1, figsize=(14, 2 * n_plots), sharex=True)

    if n_plots == 1:
        axes = [axes]

    # Color cycle for different files
    colors = plt.cm.tab10(np.linspace(0, 1, len(data_list)))
    line_styles = ['-', '-', '-.', ':']

    # Plot each channel
    for i, ch_label in enumerate(channel_list):
        ax = axes[i]

        # Plot data from each file
        for file_idx, (data, file_name, info) in enumerate(zip(data_list, file_names, recording_info)):
            signals = data['signals']
            sample_rates = data['sample_rates']

            # Skip if channel not in this file
            if ch_label not in signals:
                continue

            signal = signals[ch_label]
            fs = sample_rates[ch_label]

            # Calculate time offset for this recording relative to earliest_start
            time_offset = (info['start_time'] - earliest_start).total_seconds()

            # Create time axis (relative to earliest_start)
            time = time_offset + np.arange(len(signal)) / fs

            # Apply offset and duration constraints
            plot_start = offset_s
            if duration_s is not None:
                plot_end = offset_s + duration_s
            else:
                plot_end = total_span

            # Find indices that fall within plot range
            mask = (time >= plot_start) & (time <= plot_end)
            time_plot = time[mask]
            signal_plot = signal[mask]

            if len(time_plot) == 0:
                continue

            # Plot signal with color and style
            label = f"{file_name} ({fs:.0f} Hz)" if i == 0 else None
            linestyle = line_styles[file_idx % len(line_styles)]
            ax.plot(time_plot, signal_plot, linewidth=0.5,
                   color=colors[file_idx], linestyle=linestyle,
                   label=label, alpha=0.8)

        ax.set_ylabel(f"{ch_label}")
        ax.grid(True, alpha=0.3)

        # Add legend to first subplot
        if i == 0:
            ax.legend(loc='upper right', fontsize=8)

        # Add annotations from all files if requested
        # Show on first subplot only (default), or all subplots if annotations_all_channels=True
        if show_annotations and (annotations_all_channels or i == 0):
            for file_idx, (data, info) in enumerate(zip(data_list, recording_info)):
                annotations = data['annotations']
                if len(annotations[0]) > 0:
                    time_offset = (info['start_time'] - earliest_start).total_seconds()
                    for onset, dur, ann_label in zip(*annotations):
                        abs_onset = time_offset + onset
                        if offset_s <= abs_onset <= (offset_s + (duration_s or total_span)):
                            # Draw vertical line
                            ax.axvline(abs_onset, color=colors[file_idx], alpha=0.3,
                                     linestyle='--', linewidth=1)

                            # Add text label if showing on all channels
                            if annotations_all_channels:
                                # Offset labels vertically for different files to avoid overlap
                                y_pos = ax.get_ylim()[1] - (file_idx * 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]))
                                ax.text(abs_onset, y_pos, ann_label,
                                       rotation=90, verticalalignment='top',
                                       fontsize=7, color=colors[file_idx], alpha=0.8)

    # Set x-axis label
    axes[-1].set_xlabel(f"Time (s) relative to {earliest_start.strftime('%Y-%m-%d %H:%M:%S')}")

    # Set title
    if custom_title:
        axes[0].set_title(custom_title)
    else:
        duration_str = f"{duration_s}s" if duration_s else f"all ({total_span:.1f}s)"
        axes[0].set_title(f"Multi-File EDF Comparison ({len(data_list)} files, Duration: {duration_str})")

    plt.tight_layout()
    plt.show()


def print_edf_info(data: Dict) -> None:
    """
    Print summary information about the EDF file.

    Parameters
    ----------
    data : dict
        Data dictionary from load_edf()
    """
    signals = data['signals']
    sample_rates = data['sample_rates']
    annotations = data['annotations']
    start_time = data['start_time']
    metadata = data.get('metadata', {})

    print(f"Recording start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Print header metadata if available
    if metadata:
        print(f"\nEDF Header Metadata:")
        print(f"  File Type: {metadata.get('file_type', 'Unknown')}")
        print(f"\n  Local Patient Identification:")
        print(f"    Patient Name: {metadata.get('patient_id') or '(empty)'}")
        print(f"    Patient Code: {metadata.get('patient_code') or '(empty)'}")
        print(f"    Patient Additional: {metadata.get('patient_additional') or '(empty)'}")
        print(f"\n  Local Recording Identification:")
        print(f"    Recording Additional: {metadata.get('recording_id') or '(empty)'}")
        print(f"    Admin Code: {metadata.get('admin_code') or '(empty)'}")
        print(f"    Technician: {metadata.get('technician') or '(empty)'}")
        print(f"    Equipment: {metadata.get('equipment') or '(empty)'}")

    print(f"\nChannels ({len(signals)}):")

    for label in signals.keys():
        signal = signals[label]
        fs = sample_rates[label]
        duration = len(signal) / fs

        print(f"  {label}:")
        print(f"    Sample rate: {fs:.2f} Hz")
        print(f"    Samples: {len(signal)}")
        print(f"    Duration: {duration:.2f} s")
        print(f"    Range: [{signal.min():.2f}, {signal.max():.2f}]")

    # Print annotations
    if len(annotations[0]) > 0:
        print(f"\nAnnotations ({len(annotations[0])}):")
        for onset, duration, label in zip(*annotations):
            print(f"  {onset:.2f}s: {label} (duration: {duration:.2f}s)")
    else:
        print("\nNo annotations found")


def create_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="EDF File Viewer - Load and visualize EDF recordings",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # View all channels
  python edf_viewer.py recording.edf

  # Show file info only
  python edf_viewer.py recording.edf --info

  # Plot first 60 seconds
  python edf_viewer.py recording.edf --duration 60

  # Plot specific time window
  python edf_viewer.py recording.edf --offset 30 --duration 60

  # Plot specific channels
  python edf_viewer.py recording.edf --channels F7,F8

  # Show annotations on all channels with labels
  python edf_viewer.py recording.edf --annotations-all-channels

  # Set custom title
  python edf_viewer.py recording.edf --title "Baseline Recording - Subject 01"

  # Compare two EDF files
  python edf_viewer.py file1.edf file2.edf --duration 30

  # Compare with time shift (shift file2 backward by 20.5 seconds)
  python edf_viewer.py file1.edf file2.edf --shift -20.5 --duration 30

  # Compare files with custom title and annotations
  python edf_viewer.py file1.edf file2.edf --annotations-all-channels \
      --title "Device Comparison - Session 2026-04-21" --duration 60
        """
    )

    parser.add_argument(
        "edf_file",
        nargs='+',
        help="Path to EDF file(s) - multiple files will be overlaid"
    )

    parser.add_argument(
        "--info", "-i",
        action="store_true",
        help="Print file info and exit (don't plot)"
    )

    parser.add_argument(
        "--duration", "-d",
        type=float,
        metavar="SECONDS",
        help="Duration to plot in seconds (default: plot all)"
    )

    parser.add_argument(
        "--offset", "-o",
        type=float,
        default=0.0,
        metavar="SECONDS",
        help="Start offset in seconds (default: 0)"
    )

    parser.add_argument(
        "--channels", "-c",
        metavar="CH1,CH2,...",
        help="Comma-separated list of channels to plot (default: all)"
    )

    parser.add_argument(
        "--no-annotations",
        action="store_true",
        help="Hide annotation markers"
    )

    parser.add_argument(
        "--annotations-all-channels",
        action="store_true",
        help="Show annotations on all channel subplots (with labels). "
             "By default, annotations are only shown on the first subplot."
    )

    parser.add_argument(
        "--title",
        metavar="TEXT",
        help="Custom title for the plot (default: auto-generated from file info)"
    )

    parser.add_argument(
        "--shift",
        metavar="SECONDS",
        help="Time shift in seconds for aligning files. Single value shifts all files after the first, "
             "or comma-separated values (one per file, starting with second file). "
             "Negative values make events occur earlier. Example: --shift -20.5"
    )

    return parser


def main():
    """CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()

    try:
        # Convert to list if single file
        edf_files = args.edf_file if isinstance(args.edf_file, list) else [args.edf_file]

        # Load all EDF files
        data_list = []
        file_names = []

        for edf_file in edf_files:
            print(f"Loading EDF file: {edf_file}")
            data = load_edf(edf_file)
            print(f"Loaded successfully\n")

            # Print info for each file
            print(f"=== {edf_file} ===")
            print_edf_info(data)
            print()

            data_list.append(data)
            # Extract just the filename for legend
            import os
            file_names.append(os.path.basename(edf_file))

        # Print summary of all available channels if multiple files
        if len(data_list) > 1:
            print("=== Available Channels Across All Files ===")
            all_channels = set()
            for data in data_list:
                all_channels.update(data['signals'].keys())

            for ch_name in sorted(all_channels):
                # Show which files contain this channel
                sources = []
                for i, data in enumerate(data_list):
                    if ch_name in data['signals']:
                        fs = data['sample_rates'][ch_name]
                        sources.append(f"{file_names[i]} ({fs:.0f} Hz)")
                print(f"  {ch_name}: {', '.join(sources)}")
            print()

        # Plot if not info-only mode
        if not args.info:
            print("Generating plot...")

            # Parse channel list if provided
            channels = None
            if args.channels:
                channels = [ch.strip() for ch in args.channels.split(',')]

            # Parse time shifts if provided
            time_shifts = None
            if args.shift:
                try:
                    if ',' in args.shift:
                        # Multiple shifts provided
                        time_shifts = [float(s.strip()) for s in args.shift.split(',')]
                    else:
                        # Single shift - apply to all files after the first
                        shift_val = float(args.shift)
                        time_shifts = [shift_val] * (len(data_list) - 1)
                except ValueError:
                    print(f"Error: Invalid shift value: {args.shift}", file=sys.stderr)
                    sys.exit(1)

            # Use multi-file plot if multiple files, otherwise single-file plot
            if len(data_list) > 1:
                plot_multiple_edf_channels(
                    data_list,
                    file_names,
                    duration_s=args.duration,
                    offset_s=args.offset,
                    channels=channels,
                    show_annotations=not args.no_annotations,
                    time_shifts=time_shifts,
                    annotations_all_channels=args.annotations_all_channels,
                    custom_title=args.title
                )
            else:
                plot_edf_channels(
                    data_list[0],
                    duration_s=args.duration,
                    offset_s=args.offset,
                    channels=channels,
                    show_annotations=not args.no_annotations,
                    annotations_all_channels=args.annotations_all_channels,
                    custom_title=args.title
                )

    except FileNotFoundError as e:
        print(f"Error: File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
