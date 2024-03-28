#!/usr/bin/env python3

from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import errno
from enum import Enum
import json
import math
from scipy.signal import find_peaks

# User defined functions
from utils.Timing_Analyzer import *
from utils.tdmsUtils import *
from utils.plotUtils import *
# import config.config

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-d',default="./plots/test/",type=str,help='output directory')
parser.add_argument('--report','-r',default=100,type=int,help='report every x events')
parser.add_argument('--debug_report','-b',default=10000,type=int,help='report every x events')
parser.add_argument('--display_report','-p',default=10000,type=int,help='report every x events')
parser.add_argument('--subset','-s',default=-1,type=int,help='Process a subset of data. -1 = all')
args = parser.parse_args()

def debugPrint(string):
    if (args.debug_report): print(string)

def SingleTDMS_CW_analysis(in_filename):
    with TdmsFile.open(in_filename) as tdms_file:
        # Read Meta Data (Basic information)
        metadata = tdms_file.properties
        metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
        print(metadata_df)
        totalEvents = metadata_df.loc[metadata_df['metaKey'] == 'Total Events', 'metaValue'].iloc[0]
        timeWindow = float(metadata_df.loc[metadata_df['metaKey'] == 'Time window', 'metaValue'].iloc[0]) * 1E-9
        # Read Groups and Channels
        Read_Groups_and_Channels(tdms_file)
        # Initialize arrays
        counts, pulseRanges = np.array([]), np.array([])
        # Set threshold
        threshold = 0.3
        # Start Looping through events
        for event, chunk in enumerate(tdms_file.data_chunks()):
            # Choose a subset of the whole data to do the analysis. -1 = run All
            if (event == args.subset ): break
            # Loop progress
            if ((event+1)%args.report==0): print (f"==========Processing {event}/{totalEvents} event==========")
            args.debug_report = True if (event+1)%args.debug_report==0 else False
            args.display_report = True if (event+1)%args.display_report==0 else False
            # Skip chunk larger than totalEvents
            if (event > int(totalEvents)-1): break
            # Read chSig into np array
            chSig = chunk['ADC Readout Channels']['chSig']._data()
            # Initialize counts
            count=0
            # Find peaks
            peaks, _ = find_peaks(chSig, height=threshold, distance=300)

            # Count peaks
            num_peaks = len(peaks)
            print(f"Event {event} Number of peaks above threshold:", num_peaks)
            counts = np.append(counts, num_peaks)

            # peak ranges
            for peak in peaks:
                pulseRange = np.ptp(chSig[peak:peak+250])
                pulseRanges = np.append(pulseRanges, pulseRange)

            # Plot waveform with peaks
            if (args.display_report):
                info = r'$T=4.6K,\quad V_{Bias}=1.7V,\quad 100 \mu W,\quad 532nm\, CW\, laser$'
                plt.title(info, fontsize = 13, loc='right')
                plt.plot(chSig, label='data')
                plt.plot(peaks, chSig[peaks], "x", label='Found peaks')
                plt.xlabel('Index [0.4ns], Gate width 2ms',fontsize=15)
                plt.ylabel('ADC [V]',fontsize=15)
                plt.ylim(-0.06,0.06)
                plt.legend()
                plt.tight_layout()
                # plt.title('Waveform with Peaks')
                plt.show()

        # Results
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
        hist_counts  = plot_histo_root(counts      , 20, 0, np.max(counts), 'counts', f'Photon Counts / {timeWindow}s', '',f'{plotDir}/counts.png')
        hist_range   = plot_histo_root(pulseRanges , 20, 0.02, np.max(pulseRanges), 'signal_range', 'Voltage [V]', 'Pulse range',f'{plotDir}/signal_range.png')

        PC_mean = np.mean(counts)
        PC_meanE = np.std(counts)/math.sqrt(len(counts))
        PCR_mean = PC_mean/timeWindow
        PCR_meanE = PC_meanE/timeWindow

        with open(f'{args.outputDir}/{baseDir}/test.txt','a') as f:
            f.write(f'{basename} ; Window : {timeWindow}s ; Total Events : {len(counts)} ;')
            f.write(f'PCR mean={PCR_mean:.3f} pm {PCR_meanE:.3f}\n')
            f.write(f'Range Mean={np.mean(pulseRanges):.5f} pm {np.std(pulseRanges)/math.sqrt(len(pulseRanges)):.5f} ; Range std={np.std(pulseRanges):.5f}\n')

        print(f'{basename} ; Window : {timeWindow}s ; Total Events : {len(counts)} ;')
        print(f'PC mean={PC_mean:.3f} pm {PC_meanE:.3f} ; PCR mean={PCR_mean:.3f} pm {PCR_meanE:.3f}')
        print(f'Range Mean={np.mean(pulseRanges):.5f} pm {np.std(pulseRanges)/math.sqrt(len(pulseRanges)):.5f} ; Range std={np.std(pulseRanges):.5f}\n')

if __name__ == "__main__":

    try:
        os.makedirs(args.outputDir)
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('Plot directory exists.')
        else:
            raise

    for in_filename in args.in_filenames:
        SingleTDMS_CW_analysis(in_filename)
