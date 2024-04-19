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
from ..utils.fitfunctions import alt_expo

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-o',default="./output",type=str,help='output directory')
parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
parser.add_argument('--debug','-d',action="store_true",help='debug mode')
parser.add_argument('--saveHist','-s',action="store_true",help='save histograms')
parser.add_argument('--sweepBias',action="store_true",help='sweep bias mode')
parser.add_argument('--sweepAll',action="store_true",help='Run all plots for each sweep variable')
parser.add_argument('--fit','-f',action="store_true",help='do fit')
args = parser.parse_args()

def project(tree, hist, var, cut, title="", xtit="", ytit="", outDir="plots/test/", saveTitle="", xmin=-1, xmax=-1):
    # print (f'projecting var: {var}, cut: {cut} from tree: {tree.GetName()} into hist: {h.GetName()}')
    tree.Project(hist.GetName(),var,cut)
    # move overflow to last bin
    nbins = hist.GetNbinsX();
    lastBin = hist.GetBinContent(nbins);
    overflow = hist.GetBinContent(nbins+1);
    firstBin = hist.GetBinContent(1);
    underflow = hist.GetBinContent(0);
    hist.SetBinContent(nbins,lastBin+overflow);
    hist.SetBinContent(nbins+1, 0.0);
    hist.SetBinContent(1,firstBin+underflow);
    hist.SetBinContent(0, 0.0);
    if (args.saveHist):
        c_hist = ROOT.TCanvas()
        hist.GetXaxis().SetTitle(xtit)
        hist.GetYaxis().SetTitle(ytit)
        hist.SetTitle(title)
        hist.Draw()
        if (xmax!=-1): hist.GetXaxis().SetRangeUser(xmin,xmax)
        c_hist.SaveAs(f"{outDir}/{saveTitle}.png")

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

def rebin(h, title="", xtit="", ytit="", outDir="plots/test/", saveTitle="", save=False):
    h_new = h.Clone()
    for i in range(h.GetNbinsX()):
        previousbin = h.GetBinContent(i)/2
        thisbin = h.GetBinContent(i+1)/2
        nextbin = h.GetBinContent(i+2)/2
        newbin = previousbin + nextbin
        h_new.SetBinContent(i,newbin)
    if (save):
        c_hist = ROOT.TCanvas()
        h_new.GetXaxis().SetTitle(xtit)
        h_new.GetYaxis().SetTitle(ytit)
        h_new.SetTitle(title)
        h_new.Draw()
        c_hist.SaveAs(f"{outDir}/{saveTitle}.png")
    return h_new

def get_info(in_filename):
    laser_power = int(in_filename.split('nW/')[0].split('/')[-1])
    bias_voltage = int(in_filename.split('mV')[0].split('_')[-1])
    bias_current = int(in_filename.split('uA')[0].split('_')[-1])
    # infile = ROOT.TFile.Open(in_filename)
    # inputDataName = infile.Get('inputDataName').GetTitle()
    # with TdmsFile.open(inputDataName) as tdms_file:
    #     # Read Meta Data (Basic information)
    #     metadata = tdms_file.properties
    #     metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
    #     laser_power = float(metadata_df.loc[metadata_df['metaKey'] == 'Laser Power (uW)', 'metaValue'].iloc[0])
    #     bias_voltage = float(metadata_df.loc[metadata_df['metaKey'] == 'Bias Voltage (mV)', 'metaValue'].iloc[0])
    #     bias_current = float(metadata_df.loc[metadata_df['metaKey'] == 'Bias Current (nA)', 'metaValue'].iloc[0])
    return laser_power, bias_voltage, bias_current

def sort_bias(bias, var):
    bias_array = np.array(bias)
    var_array = np.array(var)
    sorted_bias_array = np.sort(bias_array)
    sorted_var_array= var_array[bias_array.argsort()]
    return sorted_bias_array, sorted_var_array

def Graph_sweep_bias(power, bvs, bcs, stat, title="graph", ytit="", ymin=0, ymax=1, plotDir="./"):
    # outfile.cd()
    c1 = ROOT.TCanvas()
    leg = ROOT.TLegend(0.6,0.8,0.88,0.88)
    leg.SetNColumns(4)
    graph = ROOT.TGraph()
    index=0
    for ibv, bv in enumerate(bvs):
        key = str(power) + 'uW_' + str(bv) + 'mV'
        try:
            value = stat[key]
            graph.SetPoint(index,bcs[ibv],value)
            index+=1
        except KeyError:
            pass
    graph.Draw("ALP PMC PLC")
    graph.GetXaxis().SetTitle("Bias Current (#muA)")
    graph.GetYaxis().SetTitle(ytit)
    graph.GetYaxis().SetRangeUser(ymin,ymax)
    graph.SetMarkerStyle(20)
    leg.AddEntry(graph,f'{float(power/1000):.1f}#muW','lp')
    leg.Draw()
    c1.SaveAs(f"{plotDir}/{title}_sweep_bias_current.png")

def Graph_sweep_bias_err(power, bvs, bcs, stat, stat_err, title="graph", ytit="", ymin=0, ymax=1, plotDir="./"):
    # outfile.cd()
    c1 = ROOT.TCanvas()
    leg = ROOT.TLegend(0.6,0.8,0.88,0.88)
    leg.SetNColumns(4)
    graph = ROOT.TGraphErrors()
    index=0
    for ibv, bv in enumerate(bvs):
        key = str(power) + 'uW_' + str(bv) + 'mV'
        try:
            value = stat[key]
            graph.SetPoint(index,bcs[ibv],value)
            graph.SetPointError(index,0,stat_err[key])
            index+=1
        except KeyError:
            pass
    graph.Draw("ALP PMC PLC")
    graph.GetXaxis().SetTitle("Bias Current (#muA)")
    graph.GetYaxis().SetTitle(ytit)
    graph.GetYaxis().SetRangeUser(ymin,ymax)
    graph.SetMarkerStyle(20)
    leg.AddEntry(graph,f'{float(power/1000):.1f}#muW','lp')
    leg.Draw()
    c1.SaveAs(f"{plotDir}/{title}_sweep_bias_current.png")

def Graph_sweep(powers, bvs, bcs, stat, title="graph", ytit="", ymin=0, ymax=1):
    # outfile.cd()
    graphs_sweep_power, graphs_sweep_bv = {}, {}

    c1_power = ROOT.TCanvas()
    leg_power = ROOT.TLegend(0.15,0.8,0.9,0.98)
    leg_power.SetNColumns(4)
    for ibv, bv in enumerate(bvs):
        graphs_sweep_power[bv] = ROOT.TGraph()
        graphs_sweep_power[bv].GetXaxis().SetTitle("Laser Power (#muW)")
        graphs_sweep_power[bv].GetYaxis().SetTitle(ytit)
        index=0
        for ipow, power in enumerate(powers):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                value = stat[key]
                graphs_sweep_power[bv].SetPoint(index,float(power/1000),value)
                index+=1
            except KeyError:
                continue
        if (ibv==0):
            graphs_sweep_power[bv].Draw("ALP PMC PLC")
            graphs_sweep_power[bv].GetYaxis().SetRangeUser(ymin,ymax)
        else: graphs_sweep_power[bv].Draw("LPSame PMC PLC")
        leg_power.AddEntry(graphs_sweep_power[bv],f'{bv}mV','lp')
        graphs_sweep_power[bv].SetMarkerStyle(ibv+20)
    leg_power.Draw()
    c1_power.SaveAs(f"{title}_sweep_power.png")

    c1_bv = ROOT.TCanvas()
    leg_bv = ROOT.TLegend(0.15,0.8,0.9,0.98)
    leg_bv.SetNColumns(4)
    for ipow, power in enumerate(powers):
        graphs_sweep_bv[power] = ROOT.TGraph()
        graphs_sweep_bv[power].GetXaxis().SetTitle("Bias Current (#muA)")
        graphs_sweep_bv[power].GetYaxis().SetTitle(ytit)
        index=0
        for ibv, bv in enumerate(bvs):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                value = stat[key]
                graphs_sweep_bv[power].SetPoint(index,bcs[ibv],value)
                index+=1
            except KeyError:
                pass
        if (ipow==0):
            graphs_sweep_bv[power].Draw("ALP PMC PLC")
            graphs_sweep_bv[power].GetYaxis().SetRangeUser(ymin,ymax)
        else: graphs_sweep_bv[power].Draw("LPSame PMC PLC")
        leg_bv.AddEntry(graphs_sweep_bv[power],f'{float(power/1000):.1f}#muW','lp')
        graphs_sweep_bv[power].SetMarkerStyle(ipow+20)
    leg_bv.Draw()
    c1_bv.SaveAs(f"{title}_sweep_bias_current.png")
    return graphs_sweep_power, graphs_sweep_bv

def Graph_sweep_var_err(bias, var, var_err, title="graph", xtit="Bias Current (#mnA)",ytit=""):
    c1 = ROOT.TCanvas()
    outfile.cd()
    graph = ROOT.TGraphErrors()
    for i, (b,v,e) in enumerate(zip(bias,var,var_err)):
        graph.SetPoint(i,b,v)
        graph.SetPointError(i,0,e)
    graph.Draw("AP")
    graph.SetName(title)
    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(xtit)
    graph.GetYaxis().SetTitle(ytit)
    graph.Write()

def fit_graph_efficiency(graph,rangemin,rangemax,title,xTitle,saveTitle):
    fit = ROOT.TF1("fit","expo",rangemin, rangemax)
    fit.SetLineWidth(2)
    graph.Fit("fit","R0")
    # Draw hist
    c1 = ROOT.TCanvas()
    ROOT.gStyle.SetStatY(0.85);
    ROOT.gStyle.SetStatX(0.85);
    ROOT.gStyle.SetStatW(0.2);
    ROOT.gStyle.SetStatH(0.2);
    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(xTitle)
    graph.Draw("AP")
    ROOT.gPad.Update()
    alt_fit = ROOT.TF1("alt_fit",alt_expo,0,2.5,2)
    alt_fit.SetParameter(0,fit.GetParameter(0))
    alt_fit.SetParameter(1,fit.GetParameter(1))
    alt_fit.Draw("same")
    c1.SaveAs(saveTitle)

def multi_histo_canvas(powers,bvs,bcs,histos,plotDir):
    c_multi_bv = {}
    for ibv, bv in enumerate(bvs):
        c_multi_bv[bv] = ROOT.TCanvas(f"c_multi_{bv}",f"c_multi_{bv}",1800,900)
        c_multi_bv[bv].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(powers)/5) if int(len(powers)%5==0) else int(len(powers)/5)+1
        c_multi_bv[bv].Divide(cx,cy,0,0)
        index=0
        for ipow, power in enumerate(powers):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                pad = c_multi_bv[bv].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]/1000:.1f}uA")
                histos[key].Draw()
                ROOT.gPad.Update()
                stat = histos[key].FindObject("stats")
                ROOT.gStyle.SetOptStat(1101)
                stat.SetY1NDC(0.5)
                stat.SetY2NDC(0.99)
                stat.SetX1NDC(0.55)
                stat.SetX2NDC(0.99)
                stat.SetStatFormat("6.2g")
                stat.SetFillStyle(0)
                ROOT.gPad.Modified()
                ROOT.gPad.Draw()
                index+=1
            except KeyError:
                print(key)
                continue
        c_multi_bv[bv].SaveAs(f"{plotDir}/histos_sweep_power_{bv}mV.png")

    c_multi_power = {}
    for ipower, power in enumerate(powers):
        c_multi_power[power] = ROOT.TCanvas(f"c_multi_{power}",f"c_multi_{power}",1800,900)
        c_multi_power[power].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
        c_multi_power[power].Divide(cx,cy,0,0)
        index=0
        for ibv, bv in enumerate(bvs):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                pad = c_multi_power[power].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]/1000:.1f}uA")
                histos[key].Draw()
                ROOT.gPad.Update()
                ROOT.gStyle.SetOptStat(1101)
                stat = histos[key].FindObject("stats")
                stat.SetY1NDC(0.55)
                stat.SetY2NDC(0.99)
                stat.SetX1NDC(0.55)
                stat.SetX2NDC(0.99)
                stat.SetStatFormat("6.2g")
                stat.SetFillStyle(0)
                ROOT.gPad.Modified()
                ROOT.gPad.Draw()
                index+=1
            except KeyError:
                print(key)
                continue
        c_multi_power[power].SaveAs(f"{plotDir}/histos_sweep_bv_{power}nW.png")

def multi_histo_canvas_sweep_bias(power,bvs,bcs,histos,title,plotDir):
    c_multi = ROOT.TCanvas("c_multi","c_multi",1800,900)
    c_multi_log = ROOT.TCanvas("c_multi_log","c_multi_log",1800,900)
    c_multi.SetFixedAspectRatio(True)
    c_multi_log.SetFixedAspectRatio(True)
    ROOT.gStyle.SetPadBorderMode(0)
    cx = 5
    cy = int(len(bvs)/5) if int(len(bvs)%5==0) else int(len(bvs)/5)+1
    c_multi.Divide(cx,cy,0,0)
    c_multi_log.Divide(cx,cy,0,0)
    for ibv, bv in enumerate(bvs):
        key = str(power) + 'uW_' + str(bv) + 'mV'
        try:
            histos[key].GetXaxis().SetTitle("")
            histos[key].GetYaxis().SetTitle("")
            histos[key].GetXaxis().SetLabelSize(0.1)
            histos[key].GetYaxis().SetLabelSize(0.1)
            histos[key].SetTitle("")
            histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]}uA")
            pad = c_multi.cd(ibv+1)
            histos[key].Draw()
            ROOT.gPad.Update()
            ROOT.gStyle.SetOptStat(1101)
            stat = histos[key].FindObject("stats")
            stat.SetY1NDC(0.55)
            stat.SetY2NDC(0.99)
            stat.SetX1NDC(0.55)
            stat.SetX2NDC(0.99)
            stat.SetStatFormat("6.2g")
            stat.SetFillStyle(0)
            ROOT.gPad.Modified()
            ROOT.gPad.Draw()
            pad = c_multi_log.cd(ibv+1)
            pad.SetLogy()
            histos[key].Draw()
            ROOT.gPad.Update()
            ROOT.gStyle.SetOptStat(1101)
            stat = histos[key].FindObject("stats")
            stat.SetY1NDC(0.55)
            stat.SetY2NDC(0.99)
            stat.SetX1NDC(0.55)
            stat.SetX2NDC(0.99)
            stat.SetStatFormat("6.2g")
            stat.SetFillStyle(0)
            ROOT.gPad.Modified()
            ROOT.gPad.Draw()
        except KeyError:
            print(key)
            continue
    c_multi.SaveAs(f"{plotDir}/histos_{title}_sweep_bv_{power}nW.png")
    c_multi_log.SaveAs(f"{plotDir}/histos_{title}_sweep_bv_{power}nW_log.png")

def multi_histo_canvas_fit(powers,bvs,bcs,histos,plotDir):
    c_multi_bv = {}
    for ibv, bv in enumerate(bvs):
        c_multi_bv[bv] = ROOT.TCanvas(f"c_multi_{bv}",f"c_multi_{bv}",1800,900)
        c_multi_bv[bv].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(powers)/5) if int(len(powers)%5==0) else int(len(powers)/5)+1
        c_multi_bv[bv].Divide(cx,cy,0,0)
        index=0
        for ipow, power in enumerate(powers):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                pad = c_multi_bv[bv].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]}uA")
                histos[key].Draw()
                ROOT.gPad.Update()
                stat = histos[key].FindObject("stats")
                ROOT.gStyle.SetOptStat("n")
                ROOT.gStyle.SetOptFit(0o0011)
                stat.SetY1NDC(0.5)
                stat.SetY2NDC(0.99)
                stat.SetX1NDC(0.55)
                stat.SetX2NDC(0.99)
                stat.SetStatFormat("6.2g")
                stat.DrawClone()
                pad.Modified();
                index+=1
            except KeyError:
                print(key)
                continue
        # c_multi_bv[bv].SetLogy()
        c_multi_bv[bv].SaveAs(f"{plotDir}/histos_sweep_power_{bv}mV_fit.png")

    c_multi_power = {}
    for ipower, power in enumerate(powers):
        c_multi_power[power] = ROOT.TCanvas(f"c_multi_{power}",f"c_multi_{power}",1800,900)
        c_multi_power[power].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
        c_multi_power[power].Divide(cx,cy,0,0)
        index=0
        for ibv, bv in enumerate(bvs):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                pad = c_multi_power[power].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]}uA")
                histos[key].Draw()
                ROOT.gPad.Update()
                ROOT.gStyle.SetOptStat("n")
                ROOT.gStyle.SetOptFit(0o0011)
                stat = histos[key].FindObject("stats")
                stat.SetY1NDC(0.55)
                stat.SetY2NDC(0.99)
                stat.SetX1NDC(0.55)
                stat.SetX2NDC(0.99)
                stat.SetStatFormat("6.2g")
                stat.DrawClone()
                pad.Modified();
                index+=1
            except KeyError:
                print(key)
                continue
        c_multi_power[power].SaveAs(f"{plotDir}/histos_sweep_bv_{power}nW_fit.png")

def multi_histo_canvas_fitlimit(powers,bvs,bcs,histos,integrals,plotDir):
    c_multi_bv = {}
    for ibv, bv in enumerate(bvs):
        c_multi_bv[bv] = ROOT.TCanvas(f"c_multi_{bv}",f"c_multi_{bv}",1800,900)
        c_multi_bv[bv].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(powers)/5) if int(len(powers)%5==0) else int(len(powers)/5)+1
        c_multi_bv[bv].Divide(cx,cy,0,0)
        index=0
        for ipow, power in enumerate(powers):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                pad = c_multi_bv[bv].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]}uA")
                histos[key].Draw()
                ROOT.gPad.Update()
                stat = histos[key].FindObject("stats")
                ROOT.gStyle.SetOptFit(0)
                ROOT.gStyle.SetOptStat("ne")
                stat.SetY1NDC(0.5)
                stat.SetY2NDC(0.99)
                stat.SetX1NDC(0.55)
                stat.SetX2NDC(0.99)
                stat.SetOptFit(0000)
                stat.AddText(f"NoSig={noSignals[key]*100:.0f}%")
                stat.SetStatFormat("6.2g")
                stat.DrawClone()
                pad.Modified();
                index+=1
            except KeyError:
                print(key)
                continue
            # c_multi_bv[bv].SetLogy()
        c_multi_bv[bv].SaveAs(f"{plotDir}/histos_sweep_power_{bv}mV_fitlimit.png")

    c_multi_power = {}
    for ipower, power in enumerate(powers):
        c_multi_power[power] = ROOT.TCanvas(f"c_multi_{power}",f"c_multi_{power}",1800,900)
        c_multi_power[power].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
        c_multi_power[power].Divide(cx,cy,0,0)
        index=0
        for ibv, bv in enumerate(bvs):
            key = str(power) + 'uW_' + str(bv) + 'mV'
            try:
                pad = c_multi_power[power].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{float(power/1000):.1f}uW_{bcs[ibv]}uA")
                histos[key].Draw()
                ROOT.gPad.Update()
                ROOT.gStyle.SetOptFit(0)
                ROOT.gStyle.SetOptStat("ne")
                stat = histos[key].FindObject("stats")
                stat.SetY1NDC(0.55)
                stat.SetY2NDC(0.99)
                stat.SetX1NDC(0.55)
                stat.SetX2NDC(0.99)
                stat.SetStatFormat("6.2g")
                stat.SetOptFit(0000)
                stat.AddText(f"NoSig={noSignals[key]*100:.0f}%")
                stat.DrawClone()
                pad.Modified();
                index+=1
            except KeyError:
                print(key)
                continue

def spectrum_sweep_bias(power, bvs, bcs, graphs, plotDir="./"):
    c_spectrum = ROOT.TCanvas("c_spectrum","c_spectrum",1800,900)
    c_spectrum.SetLeftMargin(0.08)
    c_spectrum.SetRightMargin(0.032)
    c_spectrum.SetTopMargin(0.032)
    c_spectrum.SetBottomMargin(0.12)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    leg = ROOT.TLegend(0.5,0.3,0.88,0.88)
    leg.SetNColumns(3)
    for ibv, bv in enumerate(bvs):
        key = str(power) + 'uW_' + str(bv) + 'mV'
        try:
            graphs[key].GetXaxis().SetTitle("Index (0.4ns)")
            graphs[key].GetYaxis().SetTitle("Voltage (V)")
            graphs[key].GetYaxis().SetTitleOffset(0.5)
            graphs[key].SetTitle("")
            graphs[key].GetYaxis().SetRangeUser(-0.05,0.35)
            graphs[key].GetXaxis().SetRangeUser(307,330)
            graphs[key].SetMarkerStyle(ibv+20)
            graphs[key].SetMarkerSize(1)
            if(ibv==0): graphs[key].Draw("APL PLC PMC")
            else: graphs[key].Draw("PLSame PLC PMC")
            leg.AddEntry(graphs[key],f'{bcs[ibv]}uA','pl')
        except KeyError:
            print(key)
            continue
    leg.SetHeader(f"Sweep Bias @ {float(power/1000):.1f}uW","C")
    leg.Draw()
    c_spectrum.SaveAs(f"{plotDir}/spectrums_sweep_bv_{power}nW.png")

def multi_spectrum_canvas_sweep_bias(power, bvs, bcs, graphs, plotDir="./"):
    c_spectrum = ROOT.TCanvas("c_spectrum","c_spectrum",1800,900)
    # c_spectrum.SetFixedAspectRatio(True)
    # ROOT.gStyle.SetPadBorderMode(0)
    cx = 6
    cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
    c_spectrum.Divide(cx,cy,0,0)
    for ibv, bv in enumerate(bvs):
        key = str(power) + 'uW_' + str(bv) + 'mV'
        try:
            graphs[key].GetXaxis().SetTitle("")
            graphs[key].GetYaxis().SetTitle("")
            graphs[key].GetXaxis().SetLabelSize(0.1)
            graphs[key].GetYaxis().SetLabelSize(0.1)
            graphs[key].GetYaxis().SetRangeUser(-0.05,0.35)
            graphs[key].GetXaxis().SetRangeUser(305,330)
            graphs[key].SetTitle(f"{float(power/1000):.1f}uW_{bcs[ibv]}uA")
            pad = c_spectrum.cd(ibv+1)
            graphs[key].Draw()
        except KeyError:
            print(key)
            continue
    c_spectrum.SaveAs(f"{plotDir}/multi_spectrums_sweep_bv_{power}nW.png")

def gettree():
    for in_filename in args.in_filenames:
        laser_power, bias_voltage, bias_current = get_info(in_filename)
        basename = str(laser_power) + 'uW_' + str(bias_voltage) + 'mV'
        plotDir= in_filename.rsplit("/",1)[0]
        infile = ROOT.TFile.Open(in_filename)
        # Get graphs
        spectrums[basename]=infile.Get('Pulse_avg_display')
        # Get Tree
        intree = infile.Get('Result_tree')
        # initialize histo
        nbin, range_min, range_max= 64, -0.5, 2
        h_pulse_fall_range = ROOT.TH1F(f"h_pulse_fall_range_{bias_current}",f"h_pulse_fall_range_{bias_current}",nbin,range_min,range_max)
        h_pulse_fall_range_limitfit = ROOT.TH1F(f"h_pulse_fall_range_{bias_current}_limitfit",f"h_pulse_fall_range_{bias_current}_limitfit",nbin,range_min,range_max)
        h_pulse_fall_tau = ROOT.TH1F("h_pulse_fall_tau","h_pulse_fall_tau",20,0,5)
        h_pre_range = ROOT.TH1F("h_pre_range","h_pre_range",100,0.,0.3)
        h_eff = ROOT.TH1F("h_eff","h_eff",2,0,2)
        h_diff = ROOT.TH1F("h_diff","h_diff",100,0,0.3)
        # Project variables to histos
        try:
            project(intree,h_pulse_fall_range,"pulse_fall_range","",basename,"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range")
        except AttributeError:
            continue
        project(intree,h_pulse_fall_range_limitfit,"pulse_fall_range","",basename,"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range")
        project(intree,h_pulse_fall_tau,"pulse_fall_tau","",basename,"pulse fall time constant (0.4ns)",f"Event",plotDir,"h_pulse_fall_tau")
        project(intree,h_pre_range,"pre_range","",basename,"pre_range (V)","Event",plotDir,"h_pre_range")
        project(intree,h_eff,"1","pulse_fall_range>0.1",basename,"Pulse detected","Event",plotDir,"h_eff")
        # Rebin
        # h_pulse_fall_range_rebin1 = rebin(h_pulse_fall_range,f'{basename}_rebin',"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range_rebin1",True)
        # h_pulse_fall_range_rebin2 = rebin(h_pulse_fall_range_rebin1,f'{basename}_rebin',"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range_rebin2",True)
        # h_pulse_fall_range_rebin3 = rebin(h_pulse_fall_range_rebin2,f'{basename}_rebin',"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range_rebin3",True)
        # Fit
        if (args.fit):
            mean, mean_error, std, std_error, integral = fit_histo_gaus(h_pulse_fall_range, -0.5, 2, f"fit_pulse_fall_range_{basename}", 'pulse_range (V)', "", f"{plotDir}/fit_pulse_fall_range.png")
            mean, mean_error, std, std_error, const, const_error, integral = fit_histo_gaus_limit(h_pulse_fall_range_limitfit, -0.5, 2, 0.003854, 0.003854, 0.03126, 0.03126, 0, 10000,  f"fit_pulse_fall_range_{basename}", 'pulse_range (V)', "", f"{plotDir}/fit_pulse_fall_range_limitfit.png")
            noSignals[basename] = (integral/0.0390625) / h_pulse_fall_range_limitfit.GetEntries()
        # Calculate
        eff = h_eff.Integral()/intree.GetEntries()
        pre_range = h_pre_range.GetMean()
        pulse_range = h_pulse_fall_range.GetMean()
        pulse_fall_tau = h_pulse_fall_tau.GetMean()
        pre_range_err = h_pulse_fall_range.GetRMS()
        pulse_range_error = h_pulse_fall_range.GetRMS()
        try:
            pulse_range_stderror = h_pulse_fall_range.GetRMS()/math.sqrt(h_pulse_fall_range.Integral())
        except ZeroDivisionError:
            pulse_range_stderror = 0
        # Append sweep variables
        if laser_power not in Pows:
            Pows.append(laser_power)
        if bias_voltage not in BVs:
            BVs.append(bias_voltage)
            BCs.append(bias_current)
        # Fill stats dict
        effs[basename] = eff
        pulse_ranges[basename] = pulse_range
        pulse_range_errs[basename] = pulse_range_error
        pre_ranges[basename] = pre_range
        pre_range_errs[basename] = pre_range_err
        pulse_fall_taus[basename] = pulse_fall_tau
        # Histograms
        h_pulse_fall_taus[basename] = h_pulse_fall_tau.Clone()
        h_pulse_fall_ranges[basename] = h_pulse_fall_range.Clone()
        h_pulse_fall_ranges_limitfit[basename] = h_pulse_fall_range_limitfit.Clone()
        h_pulse_fall_taus[basename].SetDirectory(0)
        h_pulse_fall_ranges[basename].SetDirectory(0)
        h_pulse_fall_ranges_limitfit[basename].SetDirectory(0)
        print(f"{basename}: {eff*100:.1f}%, {pulse_range*1000:.1f}mV+-{pulse_range_error*1000:.2f}mV")

def plots():
    # Sort sweep variables
    Pows.sort()
    BVs.sort()
    BCs.sort()
    # Graphs of stats vs sweep. variables
    if (args.sweepBias):
        sweepBias_plotDir = args.in_filenames[0].split('degrees/')[0] + 'degrees'
        Graph_sweep_bias(Pows[0],BVs,BCs,effs,title="g_eff",ytit="Pulse Detection Efficiency (%)",ymin=0,ymax=1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias_err(Pows[0],BVs,BCs,pulse_ranges,pulse_range_errs, title="g_pulse_range",ytit="Pulse range mean (V)",ymin=0,ymax=max(pulse_ranges.values())*1.6,plotDir=sweepBias_plotDir)
        Graph_sweep_bias(Pows[0],BVs,BCs,pre_ranges,title="g_pre_range",ytit="Pre range mean (V)",ymin=min(pre_ranges.values())*0.8,ymax=max(pre_ranges.values())*1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias(Pows[0],BVs,BCs,pulse_fall_taus,title="g_fall_tau",ytit="Fall time constant (index=0.4ns)",ymin=min(pulse_fall_taus.values())*0.8,ymax=max(pulse_fall_taus.values())*1.2,plotDir=sweepBias_plotDir)
        spectrum_sweep_bias(Pows[0],BVs,BCs,spectrums,sweepBias_plotDir)
        multi_spectrum_canvas_sweep_bias(Pows[0],BVs,BCs,spectrums,sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Pows[0],BVs,BCs,h_pulse_fall_ranges,"range",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Pows[0],BVs,BCs,h_pulse_fall_taus,"tau",sweepBias_plotDir)
    if (args.sweepAll):
        sweepAll_plotDir = args.in_filenames[0].split('nW/')[0].rsplit('/',1)[0]
        Graph_sweep(Pows,BVs,BCs,effs,title="g_eff",ytit="Pulse Detection Efficiency (%)",ymin=0,ymax=1.2,plotDir=sweepAll_plotDir)
        Graph_sweep(Pows,BVs,BCs,pulse_ranges,title="g_pulse_range",ytit="Pulse range mean (V)",ymin=0,ymax=max(pulse_ranges.values())*1.2,plotDir=sweepAll_plotDir)
        Graph_sweep(Pows,BVs,BCs,pre_ranges,title="g_pre_range",ytit="Pre range mean (V)",ymin=min(pre_ranges.values())*0.8,ymax=max(pre_ranges.values())*1.2,plotDir=sweepAll_plotDir)
        # Multi histograms
        if (args.fit):
            graphs_sweep_power_noSignal, graphs_sweep_bv_noSignal = Graph_sweep(Pows,BVs,BCs,noSignals,title="g_no_signal",ytit="No signal (%)",ymin=0,ymax=1.2,plotDir=sweepAll_plotDir)
            # fit_graph_efficiency(graphs_sweep_power_noSignal[800],0.7,2.5,"Expofit","Laser Power (#muW)",f"{plotDir}/expoFit.png")
            multi_histo_canvas_fit(Pows,BVs,BCs,h_pulse_fall_ranges,sweepAll_plotDir)
            multi_histo_canvas_fitlimit(Pows,BVs,BCs,h_pulse_fall_ranges_limitfit,noSignals,sweepAll_plotDir)
        else:
            multi_histo_canvas(Pows,BVs,BCs,h_pulse_fall_ranges,sweepAll_plotDir)

if __name__ == "__main__":
    ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
    Pows,BVs,BCs = [],[],[] # List for sweep variables
    effs, pulse_ranges, pulse_range_errs, pulse_fall_taus, pre_ranges, pre_range_errs, noSignals={},{},{},{},{},{},{} # List for stats
    h_pulse_fall_ranges,h_pulse_fall_ranges_limitfit,h_pulse_fall_taus={},{},{} # List of histos
    spectrums={} #List of graphs
    gettree() # loop over the input files
    plots() # Plot them together
