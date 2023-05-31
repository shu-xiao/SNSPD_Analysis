#!/usr/bin/env python3

from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import errno
from enum import Enum
import json

# User defined functions
from Timing_Analyzer import *
from tdmsUtils import *
from plotUtils import *
import config

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-d',default="./plots",type=str,help='output directory')
parser.add_argument('--report','-r',default=100,type=int,help='report every x events')
parser.add_argument('--debug_report','-b',default=10000,type=int,help='report every x events')
parser.add_argument('--display_report','-p',default=10000,type=int,help='report every x events')
parser.add_argument('--subset','-s',default=-1,type=int,help='Process a subset of data. -1 = all')
args = parser.parse_args()

def debugPrint(string):
    if (config.DEBUG): print(string)

def SingleTDMS_analysis(in_filename):
    with TdmsFile.open(in_filename) as tdms_file:
        # Read Meta Data (Basic information)
        metadata = tdms_file.properties
        metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
        print(metadata_df)
        totalEvents = metadata_df.loc[metadata_df['metaKey'] == 'Total Events', 'metaValue'].iloc[0]
        # Read Groups and Channels
        Read_Groups_and_Channels(tdms_file)

        channel_sum = 0.0
        channel_length = 0
        nPass=0
        ranges=[]
        # Sideband region numpy
        prePulse_mean, postPulse_mean, prePulse_stdev, postPulse_stdev, prePulse_range, postPulse_range, prePulse_integral, postPulse_integral = [], [], [], [], [], [], [], []
        # Signal region numpy
        ch1_pulse_spline_ranges, ch1_pulse_diff_ranges, ch1_pulse_amplitudes, ch1_pulse_arrivalTs, ch1_pulse_riseTs, ch1_pulse_spline_integrals = [], [], [], [], [], []
        # Start Looping through events
        for event, chunk in enumerate(tdms_file.data_chunks()):
            # Choose a subset of the whole data to do the analysis. -1 = run All
            if (event == args.subset ): break
            # Loop progress
            if ((event+1)%args.report==0): print (f"==========Processing {event}/{totalEvents} event==========")
            config.DEBUG = True if (event+1)%args.debug_report==0 else False
            config.DISPLAY = True if (event+1)%args.display_report==0 else False
            # Skip chunk larger than totalEvents
            if (event > int(totalEvents)-1): break
            # Skip first event
            if (event == 0 ): continue
            # Read ch1 into np array
            ch1 = chunk['ADC Readout Channels']['ch1']._data()
            ch2 = chunk['ADC Readout Channels']['ch2']._data()
            if (config.DISPLAY): event_display_2ch(ch1,ch2,f'Waveform')

            ch1_diff = np.diff(ch1)
            ch2_diff = np.diff(ch2)
            # Create a spline interpolation function for the data
            x_index = np.arange(len(ch1))
            ch1_spline = CubicSpline(x_index, ch1)
            ch2_spline = CubicSpline(x_index, ch2)
            # Get ch2 (trigger) arrival times
            ch2_turning_pedestals, ch2_turning_peaks = Get_turning_times(ch2_spline, 0.4, 0, len(ch2), 'Fall', config.DEBUG)
            for ipulse, (ch2_turning_pedestal, ch2_turning_peak) in enumerate(zip(ch2_turning_pedestals, ch2_turning_peaks)):
                # Skip last pulse due to distortion of the oscilloscop at the boundary
                if (ipulse >= config.NpulsePerTrigger-1): continue
                # Skip unreasonable turning points
                if ( ch2_turning_peak.x < 0 or ch2_turning_pedestal.x < 0 ): continue
                # Define time of arrival at the 50% level of the falling slope
                ch2_arrivalT = Get_Function_Arrival(ch2_spline, (ch2_turning_pedestal.y+ch2_turning_peak.y)/2, ch2_turning_pedestal.x, ch2_turning_peak.x)
                # Define signal pulse region
                Pulse_startT =  int(ch2_arrivalT) + 210
                Pulse_endT =  int(ch2_arrivalT) + 250
                # Define pre-pulse (sideband) region
                prePulse_startT =  int(ch2_arrivalT) + 10
                prePulse_endT =  int(ch2_arrivalT) + 180
                # Define post-pulse (sideband) region
                postPulse_startT =  int(ch2_arrivalT) + 300
                postPulse_endT =  int(ch2_arrivalT) + 800
                event_display_2ch(ch1[prePulse_startT:postPulse_endT], ch1_diff[prePulse_startT:postPulse_endT], f'Waveform#{event}_pulse{ipulse}')
                # Sideband characteristic
                prePulse_mean.append(np.mean(ch1[prePulse_startT:prePulse_endT])) # mean
                postPulse_mean.append(np.mean(ch1[postPulse_startT:postPulse_endT]))
                prePulse_stdev.append(np.std(ch1[prePulse_startT:prePulse_endT])) # stdev
                postPulse_stdev.append(np.std(ch1[postPulse_startT:postPulse_endT]))
                prePulse_range.append(np.ptp(ch1[prePulse_startT:prePulse_endT])) # max - min
                postPulse_range.append(np.ptp(ch1[postPulse_startT:postPulse_endT]))
                prePulse_integral.append(np.sum(ch1[prePulse_startT:prePulse_endT])) # max - min
                postPulse_integral.append(np.sum(ch1[postPulse_startT:postPulse_endT]))

                # Pulse pre-selection using sideband region
                if (prePulse_range[-1] < 0.057 or prePulse_stdev[-1] < 0.013 or postPulse_range[-1] < 0.075 or postPulse_stdev[-1] < 0.014):
                    debugPrint(f'Event{event}_Pulse{ipulse} pass preselection')
                    # Pulse region
                    ch1_pulse = ch1[Pulse_startT:Pulse_endT]
                    ch1_pulse_xIndex = np.arange(len(ch1_pulse))
                    # Cubic Spline Fit
                    ch1_pulse_spline = CubicSpline(ch1_pulse_xIndex, ch1_pulse)
                    # Pulse spline range
                    ch1_pulse_spline_range = Get_FunctionMax(ch1_pulse_spline, 7, 25).y - Get_FunctionMin(ch1_pulse_spline, 7, 25).y
                    ch1_pulse_spline_ranges.append(ch1_pulse_spline_range)
                    # Pulse spline integral
                    ch1_pulse_spline_integral, error = Get_function_integral(ch1_pulse_spline, 7, 25)
                    ch1_pulse_spline_integrals.append(ch1_pulse_spline_integral)
                    # Derivative of pulse region
                    ch1_pulse_diff = ch1_diff[Pulse_startT:Pulse_endT]
                    ch1_pulse_diff_xIndex = np.arange(len(ch1_pulse_diff))
                    ch1_pulse_diff_spline = ch1_pulse_spline.derivative()
                    # Diff spline range
                    ch1_pulse_diff_range = Get_FunctionMax(ch1_pulse_diff_spline, 8, 19).y - Get_FunctionMin(ch1_pulse_diff_spline, 8, 19).y
                    ch1_pulse_diff_ranges.append(ch1_pulse_diff_range)

                    debugPrint(f'Pulse range = {ch1_pulse_spline_range}, Diff range = {ch1_pulse_diff_range}')
                    # Event Selection
                    if (ch1_pulse_spline_range > 0.025 and ch1_pulse_diff_range > 0.01):
                        debugPrint(f'Event{event}_Pulse{ipulse} Pass event selection')
                        # Find turning point
                        ch1_pulse_turning_pedestals, ch1_pulse_turning_peaks = Get_turning_times(ch1_pulse_spline, 0.04, 6, 25, 'Rise', config.DEBUG)
                        if (len(ch1_pulse_turning_peaks)>0):
                            # Get pulse amplitude --> Defined as range between pulse rising turning points
                            ch1_pulse_amplitude = ch1_pulse_turning_peaks[0].y - ch1_pulse_turning_pedestals[0].y
                            ch1_pulse_amplitudes.append(ch1_pulse_amplitude)
                            # Get 50% pulse amplitude level
                            ch1_pulse_10 = ch1_pulse_turning_peaks[0].y*0.1 + ch1_pulse_turning_pedestals[0].y*0.9
                            ch1_pulse_50 = ch1_pulse_turning_peaks[0].y*0.5 + ch1_pulse_turning_pedestals[0].y*0.5
                            ch1_pulse_90 = ch1_pulse_turning_peaks[0].y*0.9 + ch1_pulse_turning_pedestals[0].y*0.1
                            # Get Arrival time
                            ch1_pulse_arrivalT = Get_Function_Arrival(ch1_pulse_spline, ch1_pulse_50, ch1_pulse_turning_pedestals[0].x, ch1_pulse_turning_peaks[0].x) + Pulse_startT - ch2_arrivalT
                            ch1_pulse_arrivalTs.append(ch1_pulse_arrivalT)
                            # Get Rise time
                            ch1_pulse_riseT = Get_Function_RiseFall_Range(ch1_pulse_spline, ch1_pulse_10, ch1_pulse_90, ch1_pulse_turning_pedestals[0].x, ch1_pulse_turning_peaks[0].x)
                            ch1_pulse_riseTs.append(ch1_pulse_riseT)

                            debugPrint(f'Pulse amplitude = {ch1_pulse_amplitude}, arrival Time = {ch1_pulse_arrivalT}, rise Time = {ch1_pulse_riseT}')
                        else:
                            print(f'Abnormal Event{event}_Pulse{ipulse}. Pass Event Selection, but can\'t find turning points')
                            # event_display_2ch(ch1_pulse_diff, ch1_pulse, f'Wavform#{event}_pulse{ipulse}')
                        # Create a check point for amplitude
                        # if (ch1_pulse_amplitude < 0):
                        #     print('Abnormal Event{event}_Pulse{ipulse}. Pulse amplitude is negative')
                        #     exit()
                    else:
                        debugPrint(f'Event{event}_Pulse{ipulse} Fail event selection')
                        # event_display_2ch(ch1_pulse_diff, ch1_pulse, f'Wavform#{event}_pulse{ipulse}')

                    # ch1_pulse_diff_turning_pedestals, ch1_pulse_diff_turning_peaks = Get_turning_times(ch1_pulse_diff_spline, 0.02, 0, 'Rise', config.DEBUG)
                    # display_spline_fit(ch1_pulse_spline, ch1_pulse_xIndex)
                    if (config.DISPLAY): event_display_2ch(ch1_pulse, ch1_pulse_diff, f'Wavform#{event}_pulse{ipulse}')
                else:
                    debugPrint (f'Event{event}_Pulse{ipulse} fail preselection ')

        # Plots

        # Define plotDir
        if(in_filename.find('.txt')!=0):
            basename = in_filename.rsplit('/',1)[1].split('.txt')[0]
            plotDir = args.outputDir + '/' + basename
        else:
            in_filename.rsplit('/',1)[1].split('.tdms')[0]
            plotDir = args.outputDir + '/' + basename
        # make plotdir
        try:
            os.makedirs(plotDir)
        except OSError as e:
            if e.errno == errno.EEXIST:
                print('Plot directory exists.')
            else:
                raise

        # Create root filen
        hfile = ROOT.TFile(f'{plotDir}/{basename}.root', 'RECREATE', 'analysis histograms of {basename} measurements' )

        # Sideband Region
        Sideband_results={}
        # Plot two histograms side-by-side
        Sideband_results['mean']     = plot_2histo(prePulse_mean, postPulse_mean, 50, -0.01, 0.01, 'prePulse', 'postPulse', 'Sideband mean', f'{plotDir}/sideband_mean.png')
        Sideband_results['stdev']    = plot_2histo(prePulse_stdev, postPulse_stdev, 50, 0, 0.05, 'prePulse', 'postPulse', 'Sideband stdev', f'{plotDir}/sideband_stdev.png')
        Sideband_results['range']    = plot_2histo(prePulse_range, postPulse_range, 50, 0, 0.2, 'prePulse', 'postPulse', 'Sideband range', f'{plotDir}/sideband_range.png')
        Sideband_results['integral'] = plot_2histo(prePulse_integral, postPulse_integral, 50, -1, 1, 'prePulse', 'postPulse', 'Sideband integral', f'{plotDir}/sideband_integral.png')
        # Signal Region
        Signal_results={}
        # plot_2histo(ch1_pulse_amplitudes, ch1_pulse_ranges, 50, 0, 0.2, 'Spline amplitude', 'Range', 'Pulse amplitude')
        Signal_results['amplitude'] = plot_histo(ch1_pulse_amplitudes, 256, -0.25, 0.25, 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude.png')
        Signal_results['range']     = plot_histo(ch1_pulse_diff_ranges, 50, 0, 0.1, 'Voltage [V]', 'Signal Region differentiate range',f'{plotDir}/signal_diff.png')
        Signal_results['arrivalT']  = plot_histo(ch1_pulse_arrivalTs, 50, 220, 230, 'Time [index]', 'Pulse arrival time',f'{plotDir}/signal_arrivalT.png')
        Signal_results['riseT']     = plot_histo(ch1_pulse_riseTs, 50, 0, 7, 'Time [index]', 'Pulse rise time',f'{plotDir}/signal_riseT.png')
        Signal_results['integral']  = plot_histo(ch1_pulse_spline_integrals, 50, -1, 1, 'Voltage [V]', 'Signal Region Integral',f'{plotDir}/signal_integral.png')
        plot_histo(ch1_pulse_amplitudes, 50, 0, 0.5, 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude_zoomin_largerbin.png')

        # root histograms
        hist_amplitudes  = plot_histo_root(ch1_pulse_amplitudes       , 50, 0, 0.5,   'signal_amplitude', 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude_root.png')
        hist_diff_ranges = plot_histo_root(ch1_pulse_diff_ranges      , 50, 0, 0.1,   'signal_diff', 'Voltage [V]', 'Signal Region differentiate range',f'{plotDir}/signal_diff_root.png')
        hist_arrivalTs   = plot_histo_root(ch1_pulse_arrivalTs        , 50, 220, 230, 'signal_arrivalT', 'Time [index]', 'Pulse arrival time',f'{plotDir}/signal_arrivalT_root.png')
        hist_riseTs      = plot_histo_root(ch1_pulse_riseTs           , 50, 0, 7,     'signal_riseT', 'Time [index]', 'Pulse rise time',f'{plotDir}/signal_riseT_root.png')
        hist_integrals   = plot_histo_root(ch1_pulse_spline_integrals , 50, -1, 1,    'signal_integral', 'Voltage [V]', 'Signal Region Integral',f'{plotDir}/signal_integral_root.png')

        SDE = len(ch1_pulse_amplitudes) / len(ch1_pulse_spline_integrals)
        print(f'SDE for {in_filename} : {SDE}')

        hfile.Write()
        hfile.Close()

        return SDE, Sideband_results, Signal_results

if __name__ == "__main__":

    for in_filename in args.in_filenames:
        SDE, Sideband_results, Signal_results = SingleTDMS_analysis(in_filename)
        with open('test.txt','a') as f:
            f.write(f'{in_filename}\n')
            f.write(f'SDE={SDE}\n')
            # json_str = json.dumps(Signal_results, indent=None)
            # json_str = json_str.replace('],', ',\n')
            # f.write(f'{json_str}\n')
            # json_str = json.dumps(Sideband_results, indent=None)
            # json_str = json_str.replace('],', ',\n')
            # f.write(f'{json_str}\n')
