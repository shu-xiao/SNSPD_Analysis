#!/usr/bin/env python

from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
np.seterr(all='ignore')
from enum import Enum
import json
import math

# User defined functions
from Timing_Analyzer import *
from tdmsUtils import *
from plotUtils import *
import config

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-d',default="./plots/test",type=str,help='output directory')
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
        vertical_range = float(metadata_df.loc[metadata_df['metaKey'] == 'vertical range Sig', 'metaValue'].iloc[0])
        # Read Groups and Channels
        Read_Groups_and_Channels(tdms_file)

        channel_sum = 0.0
        channel_length = 0
        nPass=0
        ranges=[]
        # Sideband region numpy
        prePulse_mean, postPulse_mean, prePulse_stdev, postPulse_stdev, prePulse_range, postPulse_range, prePulse_integral, postPulse_integral = [], [], [], [], [], [], [], []
        # Signal region simple array
        pulseRanges = []
        # Signal region numpy
        chSig_pulse_spline_ranges, chSig_pulse_diff_ranges, chSig_pulse_amplitudes, chSig_pulse_arrivalTs, chSig_pulse_arrivalTs_220, chSig_pulse_riseTs, chSig_pulse_spline_integrals = [], [], [], [], [], [], []
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
            # initialize parameter
            chSig_pulse_arrivalT_0 = -1
            # Read chSig into np array
            chSig = chunk['ADC Readout Channels']['chSig']._data()
            chTrig = chunk['ADC Readout Channels']['chTrig']._data()
            if (config.DISPLAY): event_display_2ch(chSig,chTrig,f'Waveform', 0.02)

            chSig_diff = np.diff(chSig)
            chTrig_diff = np.diff(chTrig)
            # Create a spline interpolation function for the data
            x_index = np.arange(len(chSig))
            chSig_spline = CubicSpline(x_index, chSig)
            chTrig_spline = CubicSpline(x_index, chTrig)
            # Get chTrig (trigger) arrival times
            chTrig_turning_pedestals, chTrig_turning_peaks = Get_turning_times(chTrig_spline, 0.1, 0, len(chTrig), 'Fall', config.DEBUG)
            # chTrig_turning_pedestals, chTrig_turning_peaks = Get_turning_times(chTrig_spline, 0.4, 0, len(chTrig), 'Fall', config.DEBUG)
            for ipulse, (chTrig_turning_pedestal, chTrig_turning_peak) in enumerate(zip(chTrig_turning_pedestals, chTrig_turning_peaks)):
                # Skip last pulse due to distortion of the oscilloscop at the boundary
                if (ipulse >= config.NpulsePerTrigger-1): continue
                # Skip unreasonable turning points
                if ( chTrig_turning_peak.x < 0 or chTrig_turning_pedestal.x < 0 ): continue
                # Define time of arrival at the 50% level of the falling slope
                # chTrig_arrivalT = Get_Function_Arrival(chTrig_spline, (chTrig_turning_pedestal.y+chTrig_turning_peak.y)/2, chTrig_turning_pedestal.x, chTrig_turning_peak.x)
                chTrig_arrivalT = Get_Function_Arrival(chTrig_spline, chTrig_turning_pedestal.y-0.1, chTrig_turning_pedestal.x, chTrig_turning_peak.x)
                if (chTrig_arrivalT<0) : continue
                # print(chTrig_arrivalT, int(chTrig_arrivalT))
                # Define signal pulse region
                Pulse_startT =  int(chTrig_arrivalT) + 205
                Pulse_endT =  int(chTrig_arrivalT) + 250
                # Define pre-pulse (sideband) region
                prePulse_startT =  int(chTrig_arrivalT) + 10
                prePulse_endT =  int(chTrig_arrivalT) + 180
                # Define post-pulse (sideband) region
                postPulse_startT =  int(chTrig_arrivalT) + 500
                postPulse_endT =  int(chTrig_arrivalT) + 800
                event_display(chSig[prePulse_startT:postPulse_endT], f'Waveform#{event}_pulse{ipulse}')

                # Sideband characteristic
                prePulse_mean.append(np.mean(chSig[prePulse_startT:prePulse_endT])) # mean
                postPulse_mean.append(np.mean(chSig[postPulse_startT:postPulse_endT]))
                prePulse_stdev.append(np.std(chSig[prePulse_startT:prePulse_endT])) # stdev
                postPulse_stdev.append(np.std(chSig[postPulse_startT:postPulse_endT]))
                prePulse_range.append(np.ptp(chSig[prePulse_startT:prePulse_endT])) # max - min
                postPulse_range.append(np.ptp(chSig[postPulse_startT:postPulse_endT]))
                prePulse_integral.append(np.sum(chSig[prePulse_startT:prePulse_endT])) # max - min
                postPulse_integral.append(np.sum(chSig[postPulse_startT:postPulse_endT]))
                pulseRanges.append(np.ptp(chSig[Pulse_startT:Pulse_endT]))
                # Pulse pre-selection using sideband region
                # if (prePulse_range[-1] < 0.057 or prePulse_stdev[-1] < 0.013 or postPulse_range[-1] < 0.075 or postPulse_stdev[-1] < 0.014):
                if (prePulse_range[-1] < 0.1 and prePulse_stdev[-1] < 0.02 and postPulse_range[-1] < 0.1 and postPulse_stdev[-1] < 0.02):
                    debugPrint(f'Event{event}_Pulse{ipulse} pass preselection')
                    # Pulse region
                    chSig_pulse = chSig[Pulse_startT:Pulse_endT]
                    chSig_pulse_xIndex = np.arange(len(chSig_pulse))
                    # Cubic Spline Fit
                    chSig_pulse_spline = CubicSpline(chSig_pulse_xIndex, chSig_pulse)
                    # Pulse spline range
                    # chSig_pulse_spline_range = Get_FunctionMax(chSig_pulse_spline, 7, 25).y - Get_FunctionMin(chSig_pulse_spline, 7, 25).y
                    chSig_pulse_spline_range = Get_FunctionMax(chSig_pulse_spline, 6, 32).y - Get_FunctionMin(chSig_pulse_spline, 6, 32).y
                    chSig_pulse_spline_ranges.append(chSig_pulse_spline_range)
                    # Pulse spline integral
                    # chSig_pulse_spline_integral, error = Get_function_integral(chSig_pulse_spline, 7, 25)
                    # chSig_pulse_spline_integrals.append(chSig_pulse_spline_integral)
                    # Derivative of pulse region
                    chSig_pulse_diff = chSig_diff[Pulse_startT:Pulse_endT]
                    chSig_pulse_diff_xIndex = np.arange(len(chSig_pulse_diff))
                    chSig_pulse_diff_spline = chSig_pulse_spline.derivative()
                    # Diff spline range
                    chSig_pulse_diff_range = Get_FunctionMax(chSig_pulse_diff_spline, 8, 19).y - Get_FunctionMin(chSig_pulse_diff_spline, 8, 19).y
                    chSig_pulse_diff_ranges.append(chSig_pulse_diff_range)

                    debugPrint(f'Pulse range = {chSig_pulse_spline_range}, Diff range = {chSig_pulse_diff_range}')
                    # Event Selection
                    if (chSig_pulse_spline_range > 0.1):
                    # if (chSig_pulse_spline_range > 0.002 and chSig_pulse_diff_range > 0.001):
                    # if (chSig_pulse_spline_range > 0.02 and chSig_pulse_diff_range > 0.01):
                        debugPrint(f'Event{event}_Pulse{ipulse} Pass event selection')
                        # simple range
                        # chSig_pulse_range = np.ptp(chSig_pulse[6:26])
                        chSig_pulse_range = np.ptp(chSig_pulse[6:32])
                        pulseRanges.append(chSig_pulse_range)
                        # Get pulse amplitude --> Defined as range between pulse rising turning points
                        # chSig_pulse_amplitude = chSig_pulse_turning_peaks[0].y - chSig_pulse_turning_pedestals[0].y
                        chSig_pulse_amplitude = chSig_pulse_spline_range
                        chSig_pulse_amplitudes.append(chSig_pulse_amplitude)
                        # Find turning point
                        chSig_pulse_turning_pedestals, chSig_pulse_turning_peaks = Get_turning_times(chSig_pulse_spline, 0.002, 6, 25, 'Rise', config.DEBUG)
                        # chSig_pulse_turning_pedestals, chSig_pulse_turning_peaks = Get_turning_times(chSig_pulse_spline, 0.03, 6, 25, 'Rise', config.DEBUG)
                        # chSig_pulse_turning_pedestals, chSig_pulse_turning_peaks = Get_turning_times(chSig_pulse_spline, 0.02, 6, 25, 'Rise', config.DEBUG)
                        if (len(chSig_pulse_turning_peaks)>0):
                            if ( abs(chSig_pulse_amplitude-chSig_pulse_range)/chSig_pulse_range > 0.5 ):
                                print(chSig_pulse_amplitude, chSig_pulse_range)
                                event_display(chSig_pulse,f'Waveform#{event}_pulse{ipulse}')
                            # Get 50% pulse amplitude level
                            chSig_pulse_10 = chSig_pulse_turning_peaks[0].y*0.1 + chSig_pulse_turning_pedestals[0].y*0.9
                            chSig_pulse_50 = chSig_pulse_turning_peaks[0].y*0.5 + chSig_pulse_turning_pedestals[0].y*0.5
                            chSig_pulse_90 = chSig_pulse_turning_peaks[0].y*0.9 + chSig_pulse_turning_pedestals[0].y*0.1
                            # Get Arrival time
                            chSig_pulse_arrivalT = Get_Function_Arrival(chSig_pulse_spline, chSig_pulse_50, chSig_pulse_turning_pedestals[0].x, chSig_pulse_turning_peaks[0].x) + Pulse_startT - chTrig_arrivalT  #int(chTrig_arrivalT) + 205
                            if (ipulse==0): chSig_pulse_arrivalT_0 = chSig_pulse_arrivalT
                            if (chSig_pulse_arrivalT > 0): chSig_pulse_arrivalTs.append(chSig_pulse_arrivalT)
                            if (ipulse>0 and chSig_pulse_arrivalT > 0 and chSig_pulse_arrivalT_0 > 219.4 and chSig_pulse_arrivalT_0 < 220.6): chSig_pulse_arrivalTs_220.append(chSig_pulse_arrivalT)
                            # Get Rise time
                            chSig_pulse_riseT = Get_Function_RiseFall_Range(chSig_pulse_spline, chSig_pulse_10, chSig_pulse_90, chSig_pulse_turning_pedestals[0].x, chSig_pulse_turning_peaks[0].x)
                            if (chSig_pulse_riseT > 0): chSig_pulse_riseTs.append(chSig_pulse_riseT)

                            debugPrint(f'Pulse amplitude = {chSig_pulse_amplitude}, arrival Time = {chSig_pulse_arrivalT}, rise Time = {chSig_pulse_riseT}')
                        else:
                            print(f'Abnormal Event{event}_Pulse{ipulse}. Pass Event Selection, but can\'t find turning points')
                            # event_display_2ch(chSig_pulse_diff, chSig_pulse, f'Wavform#{event}_pulse{ipulse}')
                        # Create a check point for amplitude
                        # if (chSig_pulse_amplitude < 0):
                        #     print('Abnormal Event{event}_Pulse{ipulse}. Pulse amplitude is negative')
                        #     exit()
                    else:
                        debugPrint(f'Event{event}_Pulse{ipulse} Fail event selection')
                        # event_display_2ch(chSig_pulse_diff, chSig_pulse, f'Wavform#{event}_pulse{ipulse}')

                    # chSig_pulse_diff_turning_pedestals, chSig_pulse_diff_turning_peaks = Get_turning_times(chSig_pulse_diff_spline, 0.02, 0, 'Rise', config.DEBUG)
                    # display_spline_fit(chSig_pulse_spline, chSig_pulse_xIndex)
                    if (config.DISPLAY): event_display_2ch(chSig_pulse, chSig_pulse_diff, f'Wavform#{event}_pulse{ipulse}',0.02)
                else:
                    debugPrint (f'Event{event}_Pulse{ipulse} fail preselection ')

        # Directories
        if(in_filename.find('.txt')!=-1):
            basename = in_filename.rsplit('/',1)[1].split('.txt')[0]
        else:
            basename = in_filename.rsplit('/',1)[1].split('.tdms')[0]

        baseDir = in_filename.split('/')[-2]
        plotDir = args.outputDir + '/' + baseDir + '/' + basename
        createDir(baseDir)
        createDir(plotDir)

        # Create root filen
        hfile = ROOT.TFile(f'{plotDir}/{basename}.root', 'RECREATE', 'analysis histograms of {basename} measurements' )

        # Plots
        # Sideband Region
        Sideband_results={}
        # Plot two histograms side-by-side
        Sideband_results['mean']     = plot_2histo(prePulse_mean, postPulse_mean, 50, -0.0025, 0.0025, 'prePulse', 'postPulse', 'Sideband mean', f'{plotDir}/sideband_mean.png')
        Sideband_results['stdev']    = plot_2histo(prePulse_stdev, postPulse_stdev, 50, 0, 0.005, 'prePulse', 'postPulse', 'Sideband stdev', f'{plotDir}/sideband_stdev.png')
        Sideband_results['range']    = plot_2histo(prePulse_range, postPulse_range, 100, 0, 0.05, 'prePulse', 'postPulse', 'Sideband range', f'{plotDir}/sideband_range.png')
        Sideband_results['integral'] = plot_2histo(prePulse_integral, postPulse_integral, 50, -1, 1, 'prePulse', 'postPulse', 'Sideband integral', f'{plotDir}/sideband_integral.png')
        # # # Signal Region
        # # Signal_results={}
        # # plot_2histo(chSig_pulse_amplitudes, chSig_pulse_ranges, 50, 0, 0.2, 'Spline amplitude', 'Range', 'Pulse amplitude')
        # # Signal_results['amplitude'] = plot_histo(chSig_pulse_amplitudes, 256, -0.25, 0.25, 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude.png')
        # # Signal_results['range']     = plot_histo(chSig_pulse_diff_ranges, 50, 0, 0.1, 'Voltage [V]', 'Signal Region differentiate range',f'{plotDir}/signal_diff.png')
        # # Signal_results['arrivalT']  = plot_histo(chSig_pulse_arrivalTs, 100, 217, 221, 'Time [index]', 'Pulse arrival time',f'{plotDir}/signal_arrivalT.png')
        # # Signal_results['riseT']     = plot_histo(chSig_pulse_riseTs, 50, 0, 7, 'Time [index]', 'Pulse rise time',f'{plotDir}/signal_riseT.png')
        # # Signal_results['integral']  = plot_histo(chSig_pulse_spline_integrals, 50, -1, 1, 'Voltage [V]', 'Signal Region Integral',f'{plotDir}/signal_integral.png')
        # # plot_histo(chSig_pulse_amplitudes, 50, 0, 0.5, 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude_zoomin_largerbin.png')

        # # root histograms
        hist_prePulse_mean = plot_histo_root(prePulse_mean, 50, -0.0025, 0.0025, 'prePulse_mean', 'Voltage [V]', 'PrePulse mean', f'{plotDir}/sideband_prepulse_mean.png')
        hist_postPulse_mean = plot_histo_root(postPulse_mean, 50, -0.01, 0.01, 'postPulse_mean', 'Voltage [V]', 'PostPulse mean', f'{plotDir}/sideband_postpulse_mean.png')
        hist_prePulse_stdev = plot_histo_root(prePulse_stdev, 50, 0, 0.005, 'prePulse_stdev', 'Voltage [V]', 'PrePulse stdev', f'{plotDir}/sideband_prepulse_stdev.png')
        hist_postPulse_stdev = plot_histo_root(postPulse_stdev, 50, 0, 0.005, 'postPulse_stdev', 'Voltage [V]', 'PostPulse stdev', f'{plotDir}/sideband_postpulse_stdev.png')
        hist_prePulse_range = plot_histo_root(prePulse_range, 100, 0, 0.05, 'prePulse_range', 'Voltage [V]', 'PrePulse range', f'{plotDir}/sideband_prepulse_range.png')
        hist_postPulse_range = plot_histo_root(postPulse_range, 100, 0, 0.05, 'postPulse_range', 'Voltage [V]', 'PostPulse range', f'{plotDir}/sideband_postpulse_range.png')


        vertical_range = 1
        beg = 218
        end = 222

        # # beg = 220
        # # end = 224
        hist_amplitudes        = plot_histo_root(chSig_pulse_amplitudes       , 256*4, 0,   0.05*4, 'signal_amplitude', 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude_root.png')
        # hist_amplitudes        = plot_histo_root(chSig_pulse_amplitudes       , 256*2, 0,   vertical_range*2, 'signal_amplitude', 'Voltage [V]', 'Pulse amplitude',f'{plotDir}/signal_amplitude_root.png')
        hist_range             = plot_histo_root(pulseRanges                , 256*2, 0,   vertical_range*2, 'signal_range', 'Voltage [V]', 'Pulse range',f'{plotDir}/signal_range_root.png')
        hist_range_0p1mVbinned = plot_histo_root(pulseRanges                , 1000,   0,   0.1,             'signal_range_0p1mVbinned', 'Voltage [V]', 'Pulse range',f'{plotDir}/signal_range_root_0p1mVbinned.png')
        hist_diff_ranges       = plot_histo_root(chSig_pulse_diff_ranges      , 50,    0,   0.1,              'signal_diff', 'Voltage [V]', 'Signal Region differentiate range',f'{plotDir}/signal_diff_root.png')
        hist_arrivalTs         = plot_histo_root(chSig_pulse_arrivalTs        , 40,    beg, end,              'signal_arrivalT', 'Time [index]', 'Pulse arrival time',f'{plotDir}/signal_arrivalT_root.png')
        hist_arrivalTs_220     = plot_histo_root(chSig_pulse_arrivalTs_220    , 40,    beg, end,              'signal_arrivalT_220', 'Time [index]', 'Pulse arrival time',f'{plotDir}/signal_arrivalT_220_root.png')
        hist_riseTs            = plot_histo_root(chSig_pulse_riseTs           , 50,    0,   7,                'signal_riseT', 'Time [index]', 'Pulse rise time',f'{plotDir}/signal_riseT_root.png')
        hist_integrals         = plot_histo_root(chSig_pulse_spline_integrals , 50,    -1,  1,                'signal_integral', 'Voltage [V]', 'Signal Region Integral',f'{plotDir}/signal_integral_root.png')

        range_mean, range_meanE, range_std, range_stdE = fit_histo_gaus(hist_range, 0, vertical_range*2, 'fit_signal_range', 'Voltage [V]', 'Pulse range fit', f'{plotDir}/signal_range_fit.png')
        arrivalT_mean, arrivalT_meanE, arrivalT_std, arrivalT_stdE = fit_histo_gaus(hist_arrivalTs, 219.5, 221, 'fit_signal_arrivalT', 'Time [index]', 'Pulse arrival time fit', f'{plotDir}/signal_arrivalT_fit.png')
        # arrivalT_mean, arrivalT_meanE, arrivalT_std, arrivalT_stdE = fit_histo_gaus(hist_arrivalTs, 220, 224, 'fit_signal_arrivalT', 'Time [index]', 'Pulse arrival time fit', f'{plotDir}/signal_arrivalT_fit.png')
        SDE = len(pulseRanges) / len(chSig_pulse_spline_integrals)
        SDE_error = math.sqrt((SDE*(1-SDE))/len(chSig_pulse_spline_integrals))

        with open(f'{args.outputDir}/{baseDir}/test.txt','a') as f:
            f.write(f'{basename} ; Total Events : {len(chSig_pulse_spline_integrals)} ; SDE={SDE:.3f} pm {SDE_error:.3f}\n')
            f.write(f'Range Mean fit={range_mean:.5f} pm {range_meanE:.5f} ; Range std fit={range_std:.5f} pm {range_stdE} ; fit Time jitter={arrivalT_std:.3f} pm {arrivalT_stdE:.3f}\n')
            f.write(f'Range Mean={np.mean(pulseRanges):.5f} pm {np.std(pulseRanges)/math.sqrt(len(pulseRanges)):.5f} ; Range std={np.std(pulseRanges):.5f}\n')

        print(f'{basename} ; Total Events : {len(chSig_pulse_spline_integrals)} ; SDE={SDE:.3f} pm {SDE_error:.3f}\n')
        print(f'Range Mean fit={range_mean:.5f} pm {range_meanE:.5f} ; Range std fit={range_std:.5f} pm {range_stdE} ; fit Time jitter={arrivalT_std:.3f} pm {arrivalT_stdE:.3f}\n')
        print(f'Range Mean={np.mean(pulseRanges):.5f} pm {np.std(pulseRanges)/math.sqrt(len(pulseRanges)):.5f} ; Range std={np.std(pulseRanges):.5f}\n')

        # root_numpy.array2root(np.array(chSig_pulse_amplitudes, dtype=float), f'{plotDir}/dataset.root' , 'fTree')

        hfile.Write()
        hfile.Close()
        ROOT.gROOT.GetListOfFiles().Remove(hfile)

if __name__ == "__main__":

    createDir(args.outputDir)

    for in_filename in args.in_filenames:
        # baseDir = in_filename.split('/')[-2]
        # Temperature = in_filename.split('uW_')[1].split('K_')[0]
        # Voltage = in_filename.split('/')[-1].split('mV')[0]
        # print(f'Temperature : {Temperature}, Bias Voltage : {Voltage}')
        # SDE, range_average, range_std = SingleTDMS_analysis(in_filename)
        SingleTDMS_analysis(in_filename)
