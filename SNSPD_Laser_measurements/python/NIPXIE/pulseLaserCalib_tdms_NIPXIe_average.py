#!/usr/bin/env python3

from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
np.seterr(all='ignore')
from enum import Enum
import json
import math

# User defined functions
from ..utils.Timing_Analyzer import *
from ..utils.tdmsUtils import *
from ..utils.plotUtils import *
from ..config import config

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-d',default="./plots/test",type=str,help='output directory')
parser.add_argument('--report','-r',default=100,type=int,help='report every x events')
parser.add_argument('--debug_report','-b',default=10000,type=int,help='report every x events')
parser.add_argument('--display_report','-p',default=10000,type=int,help='report every x events')
parser.add_argument('--subset','-s',default=-1,type=int,help='Process a subset of data. -1 = all')
parser.add_argument('--txtfilename','-t',default="test",type=str,help='results txt file')
args = parser.parse_args()

def debugPrint(string):
    if (config.DEBUG): print(string)

def Pulse_analysis(chSig, pulseCount, event):
    event_display(chSig/pulseCount, f'Waveform#{event}')
    prePulse["std"].Fill(np.std(chSig[config.prePulse_startT:config.prePulse_endT]/pulseCount))
    prePulse["mean"].Fill(np.mean(chSig[config.prePulse_startT:config.prePulse_endT]/pulseCount))
    prePulse["range"].Fill(np.ptp(chSig[config.prePulse_startT:config.prePulse_endT]/pulseCount))
    posPulse["std"].Fill(np.std(chSig[config.postPulse_startT:config.postPulse_endT]/pulseCount))
    posPulse["mean"].Fill(np.mean(chSig[config.postPulse_startT:config.postPulse_endT]/pulseCount))
    posPulse["range"].Fill(np.ptp(chSig[config.postPulse_startT:config.postPulse_endT]/pulseCount))

    # Pulse region
    chSig_pulse = chSig[config.Pulse_startT:config.Pulse_endT]/pulseCount
    chSig_pulse_xIndex = np.arange(len(chSig_pulse))
    event_display(chSig_pulse, f'Waveform#{event}')
    print(event, np.ptp(chSig_pulse))
    Pulse["range"].Fill(np.ptp(chSig_pulse))


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

        # initialize parameter
        chSig_average = np.zeros(1000)
        pulseCount = 0

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

            # Read chSig into np array
            chSig = chunk['ADC Readout Channels']['chSig']._data()
            chTrig = chunk['ADC Readout Channels']['chTrig']._data()
            event_display_2ch(chSig,chTrig,f'Waveform', 0.02)

            x_index = np.arange(len(chSig))
            chTrig_spline = CubicSpline(x_index, chTrig)
            # Get chTrig (trigger) arrival times
            chTrig_turning_pedestals, chTrig_turning_peaks, chTrig_turning_ranges = Get_turning_times(chTrig_spline, 0.1, 0, len(chTrig), 'Fall', config.DEBUG)
            for ipulse, (chTrig_turning_pedestal, chTrig_turning_peak) in enumerate(zip(chTrig_turning_pedestals, chTrig_turning_peaks)):
                debugPrint(f'==========Event{event}_Pulse{ipulse}==========')
                # Skip last pulse due to distortion of the oscilloscop at the boundary
                if (ipulse >= config.NpulsePerTrigger-1): continue
                # Skip unreasonable turning points
                if ( chTrig_turning_peak.x < 0 or chTrig_turning_pedestal.x < 0 ): continue
                # Define time of arrival at the 50% level of the falling slope
                chTrig_arrivalT = Get_Function_Arrival(chTrig_spline, chTrig_turning_pedestal.y-0.1, chTrig_turning_pedestal.x, chTrig_turning_peak.x)
                if (chTrig_arrivalT<0) : continue
                chSig_average = np.add(chSig_average,chSig[int(chTrig_arrivalT):int(chTrig_arrivalT) + 1000])
                pulseCount = pulseCount + 1

            if (pulseCount>30):
                Pulse_analysis(chSig_average, pulseCount, event)
                chSig_average = np.zeros(1000)
                pulseCount = 0
            else:
                continue

if __name__ == "__main__":

    createDir(args.outputDir)

    # Sideband region statistics collection
    prePulse, posPulse, Pulse = {}, {}, {}
    prePulse["std"] = ROOT.TH1F("prePulse_std", "prePulse_std", 50, 0, 0.01)
    posPulse["std"] = ROOT.TH1F("posPulse_std", "posPulse_std", 50, 0, 0.01)
    prePulse["mean"] = ROOT.TH1F("prePulse_mean", "prePulse_mean", 50, -0.005, 0.005)
    posPulse["mean"] = ROOT.TH1F("posPulse_mean", "posPulse_mean", 50, -0.005, 0.005)
    prePulse["range"] = ROOT.TH1F("prePulse_range", "prePulse_range", 100, 0, 0.03)
    posPulse["range"] = ROOT.TH1F("posPulse_range", "posPulse_range", 100, 0, 0.03)
    Pulse["range"] = ROOT.TH1F("Pulse_range", "Pulse_range", 1024, 0, 0.5)

    for in_filename in args.in_filenames:
        SingleTDMS_analysis(in_filename)
        c1 = ROOT.TCanvas()
        Pulse["range"].Draw("HIST")
        c1.SaveAs("pulseRange.png")
        for key in prePulse:
            prePulse[key].Draw("HIST")
            c1.SaveAs(f"{prePulse[key].GetName()}.png")
        for key in posPulse:
            posPulse[key].Draw("HIST")
            c1.SaveAs(f"{posPulse[key].GetName()}.png")
