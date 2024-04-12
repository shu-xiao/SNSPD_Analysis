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
        c_hist.SaveAs(f"{outDir}/{saveTitle}.png")

def color(i):
    colorwheel = [416, 600, 800, 632, 880, 432, 616, 860, 820, 900, 420, 620, 820, 652, 1000, 452, 636, 842, 863, 823]
    # colorindex = int(i/11) + int(i%11)
    return colorwheel[i]

def plot_multiHistograms():
    c1 = ROOT.TCanvas("c1","c1",900,600)
    leg = ROOT.TLegend(0.6,0.45,0.85,0.85)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0);
    outputDir = args.in_filenames[0].rsplit('/',2)[0]
    # Sort file
    params, files = array('d'), []
    for in_filename in args.in_filenames:
        # param = float(in_filename.split('/')[-1].split('uW')[0])
        param = float(in_filename.split('uW')[0].split('/')[-1])
        # if (param>0.9):continue
        params.append(param)
        files.append(in_filename)
    sorted_data = sorted(zip(params, files))
    sorted_params, sorted_files = zip(*sorted_data)

    h_amplitude = {}
    means, stds, resos, photons, intensity = array('d'),array('d'),array('d'),array('d'),array('d')
    for ifile, (param, in_filename) in enumerate(zip(sorted_params, sorted_files)):
        infile = ROOT.TFile.Open(in_filename)
        h_amplitude[ifile] = infile.Get('Pulse_range')
        mean = h_amplitude[ifile].GetMean()
        std = h_amplitude[ifile].GetStdDev()
        photon = (mean/std) * (mean/std)
        means.append(mean)
        stds.append(std)
        resos.append(std / mean)
        photons.append(photon)
        intensity.append(param)
        h_amplitude[ifile].SetLineColor(color(ifile))
        h_amplitude[ifile].SetDirectory(0)
        h_amplitude[ifile].Scale(1/h_amplitude[ifile].GetEntries());
        leg.AddEntry(h_amplitude[ifile], f"{param}uW", "l")
        if (ifile==0): h_amplitude[ifile].Draw("HIST")
        else: h_amplitude[ifile].Draw("HISTsame")
        infile.Close()

    # Plot
    leg.Draw()
    h_amplitude[0].GetYaxis().SetTitleOffset(0.7)
    h_amplitude[0].GetYaxis().SetTitle("Normalized Entries / 0.195mV")
    h_amplitude[0].GetXaxis().SetTitle("Signal Amplitude [V]")
    # h_amplitude[param].GetXaxis().SetRangeUser(0,0.06)
    h_amplitude[0].SetTitle("")
    h_amplitude[0].SetDirectory(0)
    c1.SaveAs(f"{outputDir}/multiAmplitude.png")

    gmean = ROOT.TGraph(len(intensity),intensity,means)
    gmean.Draw("AP")
    # gmean.GetXaxis().SetRangeUser(0,8)
    gmean.SetTitle("mean")
    c1.SaveAs(f"{outputDir}/multiMean.png")

    print(intensity,stds)
    gstd = ROOT.TGraph(len(intensity),intensity,stds)
    gstd.Draw("AP")
    gstd.GetXaxis().SetRangeUser(0,8)
    gstd.SetTitle("StdDev")
    c1.SaveAs(f"{outputDir}/multiStd.png")

    greso = ROOT.TGraph(len(intensity),intensity,resos)
    greso.Draw("AP")
    greso.GetXaxis().SetRangeUser(0,8)
    greso.SetTitle("relative resolution")
    c1.SaveAs(f"{outputDir}/multiReso.png")

    gphoton = ROOT.TGraph(len(intensity),intensity,photons)
    gphoton.Draw("AP")
    gphoton.Fit("pol1")
    gphoton.GetXaxis().SetRangeUser(0,8)
    gphoton.SetMinimum(0)
    gphoton.SetTitle("photon number")
    c1.SaveAs(f"{outputDir}/multiPhotons.png")

def plot_DE_polar():
    radians, DEs, Pulse_ranges = [], [], []
    for in_filename in args.in_filenames:
        degree = float(in_filename.split('degrees')[0].split('/')[-1])*2
        radians.append(np.radians(degree))
        infile = ROOT.TFile.Open(in_filename)
        h_pulse_range = infile.Get('Pulse_range')
        h_sb_range = infile.Get('prePulse_range')
        DE = (h_pulse_range.Integral() / h_sb_range.Integral()) * 100
        Pulse_range = h_pulse_range.GetMean()
        print(f"{degree:.0f}: {DE:.1f}%, {Pulse_range:.2f}V")
        DEs.append(DE)
        Pulse_ranges.append(Pulse_range)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(radians, DEs, marker='o', markersize=6, fillstyle='none', alpha=0.75)
    ax.grid(True)
    ax.set_title("Detection efficiency (%)", va='bottom')
    plt.savefig('DE.png', format='png')
    plt.show()

    fig1, ax1 = plt.subplots(subplot_kw={'projection': 'polar'})
    ax1.plot(radians, Pulse_ranges, marker='o', markersize=6, fillstyle='none', alpha=0.75)
    ax1.grid(True)
    ax1.set_title("Pulse range (V)", va='bottom')
    plt.savefig('Pulse_range.png', format='png')
    plt.show()

def get_info(in_filename):
    laser_power = in_filename.split('uW/')[0].split('/')[-1]
    bias_voltage = int(in_filename.split('mV')[0].split('/')[-1])
    bias_current = int(in_filename.split('nA')[0].split('_')[-1])
    # infile = ROOT.TFile.Open(in_filename)
    # inputDataName = infile.Get('inputDataName').GetTitle()
    # with TdmsFile.open(inputDataName) as tdms_file:
    #     # Read Meta Data (Basic information)
    #     metadata = tdms_file.properties
    #     metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
    #     laser_power = float(metadata_df.loc[metadata_df['metaKey'] == 'Laser Power (uW)', 'metaValue'].iloc[0])
    #     bias_voltage = float(metadata_df.loc[metadata_df['metaKey'] == 'Bias Voltage (mV)', 'metaValue'].iloc[0])
    #     bias_current = float(metadata_df.loc[metadata_df['metaKey'] == 'Bias Current (uA)', 'metaValue'].iloc[0])
    return laser_power, bias_voltage, bias_current

def calculate_tree(in_filename):
    plotDir= in_filename.rsplit("/",1)[0]
    basename = in_filename.rsplit('/',1)[1].split('.root')[0]
    infile = ROOT.TFile.Open(in_filename)
    intree = infile.Get('Result_tree')

    # initialize histo
    nbin, range_min, range_max= 50, -0.5, 1.5
    h_pulse_fall_range = ROOT.TH1F("h_pulse_fall_range","h_pulse_fall_range",nbin,range_min,range_max)
    h_pre_range = ROOT.TH1F("h_pre_range","h_pre_range",100,0.,0.3)
    h_eff = ROOT.TH1F("h_eff","h_eff",2,0,2)
    h_diff = ROOT.TH1F("h_diff","h_diff",100,0,0.3)

    # Project variables to histos
    project(intree,h_pulse_fall_range,"pulse_fall_range","",basename,"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range",True)
    project(intree,h_pre_range,"pre_range","",basename,"pre_range (V)","Event",plotDir,"h_pre_range",True)
    project(intree,h_eff,"1","pulse_fall_range>0.1",basename,"Pulse detected","Event",plotDir,"h_eff",True)

    # Calculate
    eff = h_eff.Integral()/intree.GetEntries()
    pre_range = h_pre_range.GetMean()
    pulse_range = h_pulse_fall_range.GetMean()
    try:
        pulse_range_error = h_pulse_fall_range.GetRMS()/math.sqrt(h_pulse_fall_range.Integral())
    except ZeroDivisionError:
        pulse_range_error = 0
    return eff, pulse_range, pulse_range_error, pre_range

def Compare_bias_var(bias, var, title="graph", xtit="Bias Current (#muA)",ytit=""):
    c1 = ROOT.TCanvas()
    outfile.cd()
    graph = ROOT.TGraph()
    for i, (b,v) in enumerate(zip(bias,var)): graph.SetPoint(i,b,v)
    graph.Draw("AP")
    graph.SetName(title)
    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(xtit)
    graph.GetYaxis().SetTitle(ytit)
    graph.Write()

def plots():
    biases, effs, pulse_ranges, pulse_range_errs, pre_ranges=[],[],[],[],[]
    for in_filename in args.in_filenames:
        # bias = float(in_filename.split("mV_")[1].split("nA")[0])/1000
        laser_power, bias_voltage, bias_current = get_info(in_filename)
        print(bias_current)
        eff, pulse_range, pulse_range_err, pre_range = calculate_tree(in_filename)
        biases.append(bias_current)
        effs.append(eff)
        pulse_ranges.append(pulse_range)
        pulse_range_errs.append(pulse_range_err)
        pre_ranges.append(pre_range)
        print(f"{bias_current}uA: {eff*100:.1f}%, {pulse_range*1000:.1f}mV+-{pulse_range_err*1000:.2f}mV")

    # Plots
    Compare_bias_var(biases,effs,title="g_eff",ytit="Pulse Count Efficiency (%)")
    Compare_bias_var(biases,pulse_ranges,title="g_pulse_range",ytit="Pulse range mean (V)")
    Compare_bias_var(biases,pre_ranges,title="g_pre_range",ytit="Pre range mean (V)")

if __name__ == "__main__":
    laser_power, bias_voltage, bias_current = get_info(args.in_filenames[0])
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
