#!/usr/bin/env python3

import ROOT
from array import array
import argparse
import matplotlib.pyplot as plt
import numpy as np
import math

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-o',default="./output",type=str,help='output directory')
parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
parser.add_argument('--debug','-d',action="store_true",help='debug mode')
args = parser.parse_args()

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

def gettree(in_filename):
    plotDir= in_filename.rsplit("/",1)[0]
    infile = ROOT.TFile.Open(in_filename)
    intree = infile.Get('Result_tree')

    # initialize histo
    h_pulse_fall_range = ROOT.TH1F("h_pulse_fall_range","h_pulse_fall_range",100,-0.1,0.3)
    h_pre_range = ROOT.TH1F("h_pre_range","h_pre_range",100,0.,0.3)
    h_diff = ROOT.TH1F("h_diff","h_diff",100,0,0.3)

    # Draw
    c1 = ROOT.TCanvas()
    intree.Draw("pulse_rise_range")
    c1.SaveAs(f"{plotDir}/pulse_range.pdf")
    intree.Draw("pulse_fall_range>>h_pulse_fall_range")
    c1.SaveAs(f"{plotDir}/pulse_fall_range.pdf")
    intree.Draw("pre_range>>h_pre_range")
    c1.SaveAs(f"{plotDir}/pre_range.pdf")
    # intree.Draw("pulse_fall_range>>h_diff","(pulse_fall_range-pre_range)/pre_range > 0.0")
    intree.Draw("pulse_fall_range>>h_diff","pulse_fall_range > 0.02")
    c1.SaveAs(f"{plotDir}/pulse_range_cut.pdf")
    eff = h_diff.Integral()/intree.GetEntries()
    pulse_amplitude = h_pulse_fall_range.GetMean()
    try:
        pulse_amplitude_error = h_pulse_fall_range.GetRMS()/math.sqrt(h_pulse_fall_range.Integral())
    except ValueError as ve:
        pulse_amplitude_error = 0
    pre_range = h_pre_range.GetMean()
    return eff, pulse_amplitude, pulse_amplitude_error, pre_range

def plots():
    biases, effs, amps, amp_errs, pre_ranges=[],[],[],[],[]
    for in_filename in args.in_filenames:
        bias = float(in_filename.split("mV_")[1].split("nA")[0])/1000
        print(bias)
        eff, amp, amp_err, pre_range= gettree(in_filename)
        biases.append(bias)
        effs.append(eff)
        amps.append(amp)
        amp_errs.append(amp_err)
        pre_ranges.append(pre_range)
        print(f"{bias}uA: {eff*100:.1f}%, {amp*1000:.1f}mV+-{amp_err*1000:.2f}mV")
    # maxEff = max(effs)
    maxEff = 1
    g_eff = ROOT.TGraph()
    g_amp = ROOT.TGraphErrors()
    g_pre = ROOT.TGraph()
    for i, (bias,eff,amp,amp_err,pre_range) in enumerate(zip(biases,effs,amps,amp_errs,pre_ranges)):
        g_eff.SetPoint(i,bias,eff/maxEff)
        g_amp.SetPoint(i,bias,amp)
        g_amp.SetPointError(i, 0, amp_err)
        g_pre.SetPoint(i,bias,pre_range)


    plot_Dir = args.in_filenames[0].rsplit("/",2)[0]
    outfile = ROOT.TFile(f'{plot_Dir}/plots.root', 'RECREATE', f'plots' )
    c2 = ROOT.TCanvas()
    g_eff.Draw("AP")
    g_eff.SetName(f"g_eff")
    g_eff.GetXaxis().SetTitle("Bias Current (#muA)")
    g_eff.GetYaxis().SetTitle("Pulse Detection Efficiency")
    c2.SaveAs(f"{plot_Dir}/eff.png")
    g_eff.Write()
    g_amp.Draw("AP")
    g_amp.SetName(f"g_amp")
    g_amp.GetXaxis().SetTitle("Bias Current (#muA)")
    g_amp.GetYaxis().SetTitle("Pulse amplitude")
    c2.SaveAs(f"{plot_Dir}/amp.png")
    g_pre.Write()
    g_pre.Draw("AP")
    g_pre.SetName(f"g_pre_range")
    g_pre.GetXaxis().SetTitle("Bias Current (#muA)")
    g_pre.GetYaxis().SetTitle("Pre pulse range")
    c2.SaveAs(f"{plot_Dir}/pre_range.png")
    g_pre.Write()
    outfile.Write()
    outfile.Close()

if __name__ == "__main__":
    # plot_multiHistograms()
    # plot_DE_polar()
    plots()
