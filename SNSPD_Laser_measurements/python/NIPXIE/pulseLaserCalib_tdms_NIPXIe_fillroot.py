#!/usr/bin/env python3

 #!/usr/bin/env python3

# printf '%s\n' /Volumes/T7/SNSPD/20230821/4.7K/*uW/*-1250mV.tdms | parallel --progress --jobs 3 python3 -m python.NIPXIE.pulseLaserCalib_tdms_NIPXIe_average {} --doSingle -d plots/20230821/20231011_1/.
# python3 -m python.NIPXIE.pulseLaserCalib_tdms_NIPXIe_fillroot -d plots/20240221/20231011_1/.

from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from array import array
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
parser.add_argument('--avgCount','-a',default=20,type=int,help='average pulse counts')
parser.add_argument('--report','-r',default=100,type=int,help='report every x events')
parser.add_argument('--debug_report','-b',action="store_true",help='report every x events')
parser.add_argument('--display_report','-p',action="store_true",help='report every x events')
parser.add_argument('--doAdvanced',action="store_true",help='do single pulse analysis')
parser.add_argument('--doAverage',action="store_true",help='do average pulse analysis')
parser.add_argument('--dryRun',action="store_true",help='dry run')
parser.add_argument('--subset','-s',default=-1,type=int,help='Process a subset of data. -1 = all')
parser.add_argument('--txtfilename','-t',default="test",type=str,help='results txt file')
args = parser.parse_args()


def debugPrint(string):
    if (config.DEBUG): print(string)

def Sideband_selection():
    if pre_range[0] < config.cut_preRange and pos_range[0] < config.cut_posRange and pre_std[0] < config.cut_preStd and pos_std[0] < config.cut_posStd:
        return True
    else:
        return True

def Advanced_pulse_analysis(data, trigT, event):
    data_xIndex = np.arange(len(data))
    # Cubic Spline Fit
    data_spline = CubicSpline(data_xIndex, data)
    # Find turning point
    data_turning_pedestals, data_turning_peaks, data_turning_ranges = Get_turning_times(data_spline, threshold, 0, len(data), 'Rise', config.DEBUG)
    if (len(data_turning_peaks)<1):
        # print(f'Abnormal Event{event}. Pass Event Selection, but can\'t find turning points')
        return
    # Get pulse amplitude --> Defined as range between pulse rising turning points
    pulse_amplitude[0] = max(data_turning_ranges)
    imax = data_turning_ranges.index(pulse_amplitude)

    # Get 50% pulse amplitude level
    data_10 = data_turning_peaks[imax].y*0.1 + data_turning_pedestals[imax].y*0.9
    data_50 = data_turning_peaks[imax].y*0.5 + data_turning_pedestals[imax].y*0.5
    data_90 = data_turning_peaks[imax].y*0.9 + data_turning_pedestals[imax].y*0.1
    # Get Arrival time
    pulse_arrivalT[0] = Get_Function_Arrival(data_spline, data_50, data_turning_pedestals[imax].x, data_turning_peaks[imax].x) + config.Pulse_startT + trigT  #int(chTrig_arrivalT) + 205
    # Get Rise time
    pulse_riseT[0] = Get_Function_RiseFall_Range(data_spline, data_10, data_90, data_turning_pedestals[imax].x, data_turning_peaks[imax].x)
    # display_spline_fit(data_spline, data_xIndex)
    debugPrint(f'Pulse amplitude = {pulse_amplitude:.4f}, arrival Time = {pulse_arrivalT:.4f}, rise Time = {pulse_riseT:.4f}')

def Simple_pulse_analysis(data, event, ipulse):
    # event_display(data, f'Waveform#{event}_{ipulse}')
    pre_std[0] = (np.std(data[config.prePulse_startT:config.prePulse_endT]))
    # pos_mean[0] = (np.mean(data[config.prePulse_startT:config.prePulse_endT]))
    pre_range[0] = (np.ptp(data[config.prePulse_startT:config.prePulse_endT]))
    pos_std[0] = (np.std(data[config.postPulse_startT:config.postPulse_endT]))
    # pre_mean[0] = (np.mean(data[config.postPulse_startT:config.postPulse_endT]))
    pos_range[0] = (np.ptp(data[config.postPulse_startT:config.postPulse_endT]))
    # Pulse region
    pulse_rise_range[0] = data[config.Pulse_rise_endT] - data[config.Pulse_startT]
    pulse_fall_range[0] = data[config.Pulse_rise_endT] - data[config.Pulse_endT]
    pulse_pre_range = (pulse_fall_range[0]-pre_range[0])/pre_range[0]
    debugPrint(f'Rise Range = {pulse_rise_range[0]:.5f}, Fall Range = {pulse_fall_range[0]:.5f}, Pre-range = {pre_range[0]:.5f}, (pulse-pre)/pre = {pulse_pre_range:.5f}')

def Common_mode_analysis(chSig_average, data):
    chSig_average = np.add(chSig_average,data)
    return chSig_average

def Find_Trigger_time_predefined():
    # chTrig_arrivalTs = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
    # chTrig_arrivalTs = [98]
    chTrig_arrivalTs = [58]
    return chTrig_arrivalTs

def Find_Trigger_time_splineFit(chTrig):
    chTrig_arrivalTs = []
    # Get chTrig (trigger) arrival times
    x_index = np.arange(len(chTrig))
    chTrig_spline = CubicSpline(x_index, chTrig)
    chTrig_turning_pedestals, chTrig_turning_peaks, chTrig_turning_ranges = Get_turning_times(chTrig_spline, 0.2, 0, len(chTrig), 'Fall', config.DEBUG)
    # Loop over laser pulse
    for ipulse, (chTrig_turning_pedestal, chTrig_turning_peak) in enumerate(zip(chTrig_turning_pedestals, chTrig_turning_peaks)):
        # Skip last pulse due to distortion of the oscilloscop at the boundary
        if (ipulse >= config.NpulsePerTrigger-1): continue
        # Skip unreasonable turning points
        if ( chTrig_turning_peak.x < 0 or chTrig_turning_pedestal.x < 0 ): continue
        # Define time of arrival at the 50% level of the falling slope
        chTrig_arrivalT = Get_Function_Arrival(chTrig_spline, chTrig_turning_pedestal.y-0.1, chTrig_turning_pedestal.x, chTrig_turning_peak.x)
        if (chTrig_arrivalT<0) : continue
        chTrig_arrivalTs.append(chTrig_arrivalT)
    return chTrig_arrivalTs

def average_plots(chSig_average, pulseCount):
    c1 = ROOT.TCanvas()
    Pulse_avg_display = ROOT.TGraph()
    Pulse_avg_display.SetName(f"Pulse_avg_{args.avgCount}_display")
    for i in range(len(chSig_average)): Pulse_avg_display.SetPoint(i,i,chSig_average[i]/outtree.GetEntries())
    Pulse_avg_display.SetMarkerStyle(4)
    Pulse_avg_display.SetMarkerSize(0.5)
    Pulse_avg_display.Draw("ALP")
    # c1.SaveAs(f"{outDir}/{Pulse_avg_display.GetName()}.png")
    Pulse_avg_display.Write()

###############################################
#################### Main  ####################
###############################################
def SingleTDMS_analysis():
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
        # chSig_average = np.zeros(avg_buffer)
        chSig_average = np.zeros(recordlength)
        pulseCount = 0
        # Start Loop
        print (f"==========Start Looping==========")
        for event in range(totalEvents):
            # if (args.dryRun): break
            # Choose a subset of the whole data to do the analysis. -1 = run All
            if (event == args.subset): break
            if (outtree.GetEntries() == config.totalTreeEvents): break
            # Loop progress
            if ((event)%args.report==0): print (f"==========Processing {event}/{totalEvents} event==========")
            # Read chSig into np array
            chSig = chSig_total[event * recordlength:(event+1) * recordlength]
            chTrig = chTrig_total[event * recordlength:(event+1) * recordlength]
            # event_display_2ch(chSig,chTrig,f'Waveform', 0.02)
            # chSig_diff = np.diff(chSig)
            # event_display_2ch(chSig,chSig_diff,f'Waveform', 0.04)
            #
            # Find Trigger timeing
            chTrig_arrivalTs = Find_Trigger_time_splineFit(chTrig) if args.doAdvanced else Find_Trigger_time_predefined()
            for ipulse, chTrig_arrivalT in enumerate(chTrig_arrivalTs):
                debugPrint(f'==========Event{event}_Pulse{ipulse}==========')
                pulseCount = pulseCount + 1
                # Record simple signal and sideband ranges into histogram
                # Simple_pulse_analysis(chSig[int(chTrig_arrivalT):int(chTrig_arrivalT) + avg_buffer], event, ipulse)
                Simple_pulse_analysis(chSig, event, ipulse)
                # chSig_average = Common_mode_analysis(chSig_average, chSig[int(chTrig_arrivalT):int(chTrig_arrivalT) + avg_buffer])
                # Record an example pulse waveform
                # if (event==2 and ipulse == 4):
                #     for i in range(avg_buffer): Pulse_display.SetPoint(i,i,chSig[int(chTrig_arrivalT) + i])
                # Do advanced analysis (Rising time, timing jitter, sophisticated amplitude)
                if (args.doAdvanced):
                    Advanced_pulse_analysis(chSig[int(chTrig_arrivalT):int(chTrig_arrivalT) + avg_buffer], chTrig_arrivalT - int(chTrig_arrivalT), event)
                if Sideband_selection():
                    debugPrint("pass sideband selection")
                    outtree.Fill()
                    chSig_average = Common_mode_analysis(chSig_average, chSig)
                    event_display_2ch(chSig,chTrig,f'Waveform', 0.02)
                else:
                    debugPrint("fail sideband selection")
                    debugPrint(f"{pre_std[0]}, {pos_std[0]}, {pre_range[0]}, {pos_range[0]}")
        ########## End Reading TDMS file ##########
        print (f"==========End Looping==========")
    # Output some numbers
    print(f"TotalEvents:{totalEvents}, TriggerPulse_Count:{pulseCount}, PassSideband_Count: {outtree.GetEntries()}")
    # Plots
    average_plots(chSig_average, pulseCount)


if __name__ == "__main__":

    if (args.debug_report==True): config.DEBUG = True
    if (args.display_report==True): config.DISPLAY = True

    ########## Init ##########
    in_filename = args.in_filenames[0]
    basename = in_filename.rsplit('/',1)[1].split('.tdms')[0]
    baseDir = in_filename.split('/')[-2]
    outDir = args.outputDir + '/' + baseDir + '/' + basename
    createDir(outDir)
    # Create root filen
    outfile = ROOT.TFile(f'{outDir}/{basename}.root', 'RECREATE', f'analysis histograms of {basename} measurements' )
    print (f'{outDir}/{basename}.root')
    outtree = ROOT.TTree("Result_tree","Pulse laser analysis results")
    # Define variables for branch
    pre_std, pos_std, pre_mean, pos_mean, pre_range, pos_range = array('f',[0]),array('f',[0]),array('f',[0]),array('f',[0]),array('f',[0]),array('f',[0])
    pulse_rise_range, pulse_fall_range, pulse_amplitude, pulse_arrivalT, pulse_riseT, pulse_pre_range = array('f',[0]),array('f',[0]),array('f',[0]),array('f',[0]),array('f',[0]),array('f',[0])
    # init branches
    # outtree.Branch('pre_std', pre_std, 'pre_std/F')
    # outtree.Branch('pos_std', pos_std, 'pos_std/F')
    # outtree.Branch('pre_mean', pre_mean, 'pre_mean/F')
    # outtree.Branch('pos_mean', pos_mean, 'pos_mean/F')
    outtree.Branch('pre_range', pre_range, 'pre_range/F')
    # outtree.Branch('pos_range', pos_range, 'pos_range/F')
    outtree.Branch('pulse_rise_range', pulse_rise_range, 'pulse_rise_range/F')
    outtree.Branch('pulse_fall_range', pulse_fall_range, 'pulse_fall_range/F')
    # outtree.Branch('pulse_pre_range', pulse_pre_range, 'pulse_pre_range/F')
    # outtree.Branch('pulse_amplitude', pulse_amplitude, 'pulse_amplitude/F')
    # outtree.Branch('pulse_arrivalT', pulse_arrivalT, 'pulse_arrivalT/F')
    # outtree.Branch('pulse_riseT', pulse_riseT, 'pulse_riseT/F')
    ########## End Init ##########

    # Start Analysis
    SingleTDMS_analysis()
    # plots()
    # End Analysis

    outtree.Write()
    outfile.Close()
