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
parser.add_argument('--debug_report','-b',action="store_true",help='report every x events')
parser.add_argument('--display_report','-p',action="store_true",help='report every x events')
parser.add_argument('--subset','-s',default=-1,type=int,help='Process a subset of data. -1 = all')
parser.add_argument('--txtfilename','-t',default="test",type=str,help='results txt file')
args = parser.parse_args()


def debugPrint(string):
    if (config.DEBUG): print(string)

def Advanced_pulse_analysis(data, trigT, event):
    data_xIndex = np.arange(len(data))
    # Cubic Spline Fit
    data_spline = CubicSpline(data_xIndex, data)
    # Find turning point
    data_turning_pedestals, data_turning_peaks, data_turning_ranges = Get_turning_times(data_spline, 0.02, 0, len(data), 'Rise', config.DEBUG)
    if (len(data_turning_peaks)<1):
        print(f'Abnormal Event{event}. Pass Event Selection, but can\'t find turning points')
        return
    # Get pulse amplitude --> Defined as range between pulse rising turning points
    data_amplitude = max(data_turning_ranges)
    Pulse["amplitude_avg"].Fill(data_amplitude)
    imax = data_turning_ranges.index(data_amplitude)

    # Get 50% pulse amplitude level
    data_10 = data_turning_peaks[imax].y*0.1 + data_turning_pedestals[imax].y*0.9
    data_50 = data_turning_peaks[imax].y*0.5 + data_turning_pedestals[imax].y*0.5
    data_90 = data_turning_peaks[imax].y*0.9 + data_turning_pedestals[imax].y*0.1
    # Get Arrival time
    data_arrivalT = Get_Function_Arrival(data_spline, data_50, data_turning_pedestals[imax].x, data_turning_peaks[imax].x) + config.Pulse_startT + trigT  #int(chTrig_arrivalT) + 205
    Pulse["arrivalT"].Fill(data_arrivalT)
    # Get Rise time
    data_riseT = Get_Function_RiseFall_Range(data_spline, data_10, data_90, data_turning_pedestals[imax].x, data_turning_peaks[imax].x)
    Pulse["riseT"].Fill(data_riseT)

    debugPrint(f'Pulse amplitude = {data_amplitude:.4f}, arrival Time = {data_arrivalT:.4f}, rise Time = {data_riseT:.4f}')
    display_spline_fit(data_spline, data_xIndex)

def Simple_pulse_analysis(data, event, prePulse, posPulse, Pulse):
    event_display(data, f'Waveform#{event}')
    if (event < 10):
        for i in range(1000):
            Pulse_avg_display.SetPoint(i,i,data[i])

    prePulse["std"].Fill(np.std(data[config.prePulse_startT:config.prePulse_endT]))
    prePulse["mean"].Fill(np.mean(data[config.prePulse_startT:config.prePulse_endT]))
    prePulse["range"].Fill(np.ptp(data[config.prePulse_startT:config.prePulse_endT]))
    posPulse["std"].Fill(np.std(data[config.postPulse_startT:config.postPulse_endT]))
    posPulse["mean"].Fill(np.mean(data[config.postPulse_startT:config.postPulse_endT]))
    posPulse["range"].Fill(np.ptp(data[config.postPulse_startT:config.postPulse_endT]))

    # Pulse region
    data_pulse = data[config.Pulse_startT:config.Pulse_endT]
    data_range = np.ptp(data_pulse)
    Pulse["range"].Fill(data_range)
    debugPrint(f'Range = {data_range:.5f}')
    event_display(data_pulse, f'Waveform#{event}')

def SingleTDMS_analysis(in_filename):
    with TdmsFile.open(in_filename) as tdms_file:
        # Read Meta Data (Basic information)
        metadata = tdms_file.properties
        metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
        print(metadata_df)
        totalEvents = int(metadata_df.loc[metadata_df['metaKey'] == 'Total Events', 'metaValue'].iloc[0])
        recordlength = int(metadata_df.loc[metadata_df['metaKey'] == 'record length', 'metaValue'].iloc[0])
        vertical_range = float(metadata_df.loc[metadata_df['metaKey'] == 'vertical range Sig', 'metaValue'].iloc[0])
        # Read Groups and Channels
        Read_Groups_and_Channels(tdms_file)
        chSig_total = tdms_file['ADC Readout Channels']['chSig']
        chTrig_total = tdms_file['ADC Readout Channels']['chTrig']

        # initialize parameter
        chSig_average = np.zeros(1000)
        chTrig_arrivalT_average = 0
        pulseCount = 0

        for event in range(totalEvents):
            # Choose a subset of the whole data to do the analysis. -1 = run All
            if (event == args.subset ): break
            # Loop progress
            if ((event+1)%args.report==0): print (f"==========Processing {event}/{totalEvents} event==========")
            # Read chSig into np array
            chSig = chSig_total[event * recordlength:(event+1) * recordlength]
            chTrig = chTrig_total[event * recordlength:(event+1) * recordlength]
            event_display_2ch(chSig,chTrig,f'Waveform', 0.02)
            # Get chTrig (trigger) arrival times
            x_index = np.arange(len(chTrig))
            chTrig_spline = CubicSpline(x_index, chTrig)
            chTrig_turning_pedestals, chTrig_turning_peaks, chTrig_turning_ranges = Get_turning_times(chTrig_spline, 0.1, 0, len(chTrig), 'Fall', config.DEBUG)
            # Loop over laser pulse
            for ipulse, (chTrig_turning_pedestal, chTrig_turning_peak) in enumerate(zip(chTrig_turning_pedestals, chTrig_turning_peaks)):
                debugPrint(f'==========Event{event}_Pulse{ipulse}==========')
                # Skip last pulse due to distortion of the oscilloscop at the boundary
                if (ipulse >= config.NpulsePerTrigger-1): continue
                # Skip unreasonable turning points
                if ( chTrig_turning_peak.x < 0 or chTrig_turning_pedestal.x < 0 ): continue
                # Define time of arrival at the 50% level of the falling slope
                chTrig_arrivalT = Get_Function_Arrival(chTrig_spline, chTrig_turning_pedestal.y-0.1, chTrig_turning_pedestal.x, chTrig_turning_peak.x)
                if (chTrig_arrivalT<0) : continue
                if (event==1 and ipulse == 3):
                    for i in range(1000):
                        Pulse_display.SetPoint(i,i,chSig[int(chTrig_arrivalT) + i])
                chSig_average = np.add(chSig_average,chSig[int(chTrig_arrivalT):int(chTrig_arrivalT) + 1000])
                chTrig_arrivalT_average = chTrig_arrivalT_average + ( chTrig_arrivalT - int(chTrig_arrivalT) )
                pulseCount = pulseCount + 1
            # Analysis after averaging pulses
            if (pulseCount>20):
                Simple_pulse_analysis(chSig_average/pulseCount, event)
                Advanced_pulse_analysis(chSig_average/pulseCount, chTrig_arrivalT_average/pulseCount, event)
                chSig_average = np.zeros(1000)
                chTrig_arrivalT_average = 0
                pulseCount = 0
            else:
                continue

if __name__ == "__main__":

    in_filename = args.in_filenames[0]
    if (args.debug_report==True): config.DEBUG = True
    if (args.display_report==True): config.DISPLAY = True

    # Make Directories
    if(in_filename.find('.txt')!=-1):
        basename = in_filename.rsplit('/',1)[1].split('.txt')[0]
    else:
        basename = in_filename.rsplit('/',1)[1].split('.tdms')[0]
    baseDir = in_filename.split('/')[-2]
    plotDir = args.outputDir + '/' + baseDir + '/' + basename
    createDir(args.outputDir)
    createDir(baseDir)
    createDir(plotDir)

    # Create root filen
    hfile = ROOT.TFile(f'{plotDir}/{basename}.root', 'RECREATE', 'analysis histograms of {basename} measurements' )

    # Histogram collection
    prePulse, posPulse, Pulse = {}, {}, {}
    prePulse["std"] = ROOT.TH1F("prePulse_std", "prePulse_std", 50, 0, 0.01)
    posPulse["std"] = ROOT.TH1F("posPulse_std", "posPulse_std", 50, 0, 0.01)
    prePulse["mean"] = ROOT.TH1F("prePulse_mean", "prePulse_mean", 50, -0.005, 0.005)
    posPulse["mean"] = ROOT.TH1F("posPulse_mean", "posPulse_mean", 50, -0.005, 0.005)
    prePulse["range"] = ROOT.TH1F("prePulse_range", "prePulse_range", 100, 0, 0.05)
    posPulse["range"] = ROOT.TH1F("posPulse_range", "posPulse_range", 100, 0, 0.05)

    prePulse_avg["std"] = ROOT.TH1F("prePulse_std", "prePulse_std", 50, 0, 0.01)
    posPulse_avg["std"] = ROOT.TH1F("posPulse_std", "posPulse_std", 50, 0, 0.01)
    prePulse_avg["mean"] = ROOT.TH1F("prePulse_mean", "prePulse_mean", 50, -0.005, 0.005)
    posPulse_avg["mean"] = ROOT.TH1F("posPulse_mean", "posPulse_mean", 50, -0.005, 0.005)
    prePulse_avg["range"] = ROOT.TH1F("prePulse_range", "prePulse_range", 100, 0, 0.05)
    posPulse_avg["range"] = ROOT.TH1F("posPulse_range", "posPulse_range", 100, 0, 0.05)

    Pulse["range"] = ROOT.TH1F("Pulse_range", "Pulse Range", 1024, 0, 0.5)
    Pulse["range"] = ROOT.TH1F("Pulse_range", "Pulse Range Average Over 20 pulses", 1024, 0, 0.5)
    Pulse["amplitude"] = ROOT.TH1F("Pulse_amplitude", "Pulse Amplitude Average Over 20 pulses", 1024, 0, 0.5)
    Pulse["arrivalT"] = ROOT.TH1F("Pulse_arrivalT", "Pulse arrival time", 100, 441, 445)
    Pulse["riseT"] = ROOT.TH1F("Pulse_riseT", "Pulse rising time", 100, 0, 6)
    Pulse_display = ROOT.TGraph()
    Pulse_avg_display = ROOT.TGraph()
    Pulse_display.SetName("Pulse_display")
    Pulse_avg_display.SetName("Pulse_avg_display")

    #################### Start Analysis ####################
    SingleTDMS_analysis(in_filename)
    #################### End Analysis ####################

    ###############################################
    #################### Plots ####################
    ###############################################
    c1 = ROOT.TCanvas()
    for key in Pulse:
        Pulse[key].Draw("HIST")
        c1.SaveAs(f"{plotDir}/{Pulse[key].GetName()}.png")
    for key in prePulse:
        prePulse[key].Draw("HIST")
        c1.SaveAs(f"{plotDir}/{prePulse[key].GetName()}.png")
    for key in posPulse:
        posPulse[key].Draw("HIST")
        c1.SaveAs(f"{plotDir}/{posPulse[key].GetName()}.png")

    Pulse_display.Write()
    Pulse_avg_display.Write()
    hfile.Write()
    hfile.Close()
    ROOT.gROOT.GetListOfFiles().Remove(hfile)
