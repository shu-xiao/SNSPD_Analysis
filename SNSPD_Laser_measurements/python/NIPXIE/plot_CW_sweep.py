#!/usr/bin/env python3

import ROOT
from array import array
import argparse
import matplotlib.pyplot as plt
import numpy as np
import math

from ..utils.plotUtils import *
from ..utils.tdmsUtils import *
from ..utils.osUtils import *

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-o',default="./output",type=str,help='output directory')
parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
parser.add_argument('--debug','-d',action="store_true",help='debug mode')
args = parser.parse_args()

def project(tree, h, var, cut, title="", xtit="", ytit="", outDir="plots/test/", saveTitle="", save=False):
    # print (f'projecting var: {var}, cut: {cut} from tree: {tree.GetName()} into hist: {h.GetName()}')
    tree.Project(h.GetName(),var,cut)
    if (save):
        c_hist = ROOT.TCanvas()
        h.GetXaxis().SetTitle(xtit)
        h.GetYaxis().SetTitle(ytit)
        h.SetTitle(title)
        h.Draw()
        c_hist.SaveAs(f"{outDir}/{h.GetName}.png")

def get_info(in_filename):
    laser_power = in_filename.split('uW/')[0].split('/')[-1]
    bias_voltage = in_filename.split('mV')[0].split('/')[-1]
    bias_current = in_filename.split('nA')[0].split('_')[-1]
    return laser_power, bias_voltage, bias_current

def calculate_tree(in_filename,laser_power, bias_voltage, bias_current):
    plotDir= in_filename.rsplit("/",1)[0]
    infile = ROOT.TFile.Open(in_filename)
    intree = infile.Get('SNSPD_data')

    # initialize histo
    h_count = ROOT.TH1F("h_count","h_count",20,0,20)
    h_pulseRange = ROOT.TH1F("h_pulseRange","h_pulseRange",100,0,3)

    # Project variables to histos
    project(intree,h_count,"pulseCount","",f"{laser_power}uW_{bias_voltage}mV_{bias_current}uA","Pulse Counts","Event",plotDir,"h_count",True)
    project(intree,h_pulseRange,"pulseRange","",f"{laser_power}uW_{bias_voltage}mV_{bias_current}uA","Pulse Range (V)","Event",plotDir,"h_range",True)

    # Calculate
    count = h_count.GetMean()
    pulse_range = h_pulseRange.GetMean()
    return count, pulse_range

def sort_bias(bias, var):
    bias_array = np.array(bias)
    var_array = np.array(var)
    sorted_bias_array = np.sort(bias_array)
    sorted_var_array= var[bias_array.argsort()]
    return sorted_bias_array, sorted_var_array

def Compare_bias_var(bias, var, title="graph", xtit="Bias Current (#muA)",ytit=""):
    c1 = ROOT.TCanvas()
    outfile.cd()
    graph = ROOT.TGraph()
    st_bias, st_var = sort_bias(bias,var)
    for i, (b,v) in enumerate(zip(st_bias,st_var)): graph.SetPoint(i,b,v)
    graph.Draw("AP")
    graph.SetName(title)
    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(xtit)
    graph.GetYaxis().SetTitle(ytit)
    graph.Write()

def plots():
    biases, counts, pulse_ranges=[],[],[]
    for in_filename in args.in_filenames:
        # bias = float(in_filename.split("mV_")[1].split("nA")[0])/1000
        laser_power, bias_voltage, bias_current = get_info(in_filename)
        count, pulse_range = calculate_tree(in_filename,laser_power, bias_voltage, bias_current)
        biases.append(bias_current)
        counts.append(count)
        pulse_ranges.append(pulse_range)
        print(f"{bias_current}uA: Counts: {count}, {pulse_range*1000:.1f}mV+-{pulse_range_err*1000:.2f}mV")

    # Plots
    Compare_bias_var(biases,counts,title="g_count",ytit="Pulse Count Efficiency (%)")
    Compare_bias_var(biases,pulse_ranges,title="g_pulse_range",ytit="Pulse range mean (V)")
    Compare_bias_var(biases,pre_ranges,title="g_pre_range",ytit="Pre range mean (V)")

if __name__ == "__main__":
    laser_power = args.in_filenames[0].split('uW/')[0].split('/')[-1]
    baseDir = args.in_filenames[0].split('uW/')[0]
    outDir = baseDir + "uW/"
    createDir(outDir)
    outfile = ROOT.TFile(f'{outDir}/plot_{laser_power}uW.root', 'RECREATE', f'plots for laser_power {laser_power}uW' )
    # Compare plots
    plots()
    # plot_multiHistograms()
    # plot_DE_polar()
    print(f'Outfile: {outDir}/plot_{laser_power}uW.root')
    outfile.Write()
    outfile.Close()
