#!/usr/bin/env python

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

# User defined functions
from Timing_Analyzer import *
from tdmsUtils import *
from plotUtils import *
import config

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
    if (config.DEBUG): print(string)

def SingleTDMS_CW_analysis(in_filename):
    with TdmsFile.open(in_filename) as tdms_file:
        # Read Meta Data (Basic information)
        metadata = tdms_file.properties
        metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
        print(metadata_df)
        totalEvents = metadata_df.loc[metadata_df['metaKey'] == 'Total Events', 'metaValue'].iloc[0]
        vertical_range = float(metadata_df.loc[metadata_df['metaKey'] == 'Ch1 vertical range', 'metaValue'].iloc[0])
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
        ch1_pulse_spline_ranges, ch1_pulse_diff_ranges, ch1_pulse_amplitudes, ch1_pulse_arrivalTs, ch1_pulse_arrivalTs_220, ch1_pulse_riseTs, ch1_pulse_spline_integrals = [], [], [], [], [], [], []
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
            # Read ch1 into np array
            ch1 = chunk['ADC Readout Channels']['ch1']._data()
            ch1_diff = np.diff(ch1)
            if (config.DISPLAY): event_display_2ch(ch1,ch1_diff,f'Waveform',0.1)
            # Initialize counts
            count=0


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
