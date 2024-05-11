#!/usr/bin/env python3
# NI-5162 specs: https://www.ni.com/docs/en-US/bundle/pxie-5162-specs/page/specs.html

import ROOT
from array import array
import argparse
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

from ..utils.plotUtils import *
from ..utils.tdmsUtils import *
from ..utils.osUtils import *
from ..utils.fitfunctions import alt_expo

np.seterr(all='ignore')

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-o',default="./output",type=str,help='output directory')
parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
parser.add_argument('--debug','-d',action="store_true",help='debug mode')
parser.add_argument('--initPlot','-s',action="store_true",help='save histograms')
parser.add_argument('--sweepBias',action="store_true",help='sweep bias mode')
parser.add_argument('--sweepAll',action="store_true",help='Run all plots for each sweep variable')
parser.add_argument('--sweepPolarization',action="store_true",help='Run plots for polarization')
parser.add_argument('--fit','-f',action="store_true",help='do fit')
parser.add_argument('--res',action="store_true",help='do resistivie mode analysis')
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
    if (args.initPlot):
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

def calculate_pulse_photon(laser_power,repetition_rate,attenuation,wavelength):
    c = 299792458
    h = 6.62607015e-34
    SNSPD_power = laser_power*attenuation
    frequency = c / wavelength
    single_photon_energy = h*frequency
    photon_number = SNSPD_power/(repetition_rate*single_photon_energy)
    return int(photon_number)

def get_info(in_filename):
    # Read the JSON file into a DataFrame
    jsonfile = in_filename.replace('.root', '.json')
    df = pd.read_json(jsonfile, lines=True)
    laser_power        = float(df[df['metaKey'] == 'Laser Power (uW)']['metaValue'].iloc[0])
    bias_voltage       = int(float(df[df['metaKey'] == 'Bias Voltage (mV)']['metaValue'].iloc[0]))
    bias_current       = int(float(df[df['metaKey'] == 'Bias Current (uA)']['metaValue'].iloc[0]))
    vertical_range     = float(df[df['metaKey'] == 'vertical range Sig']['metaValue'].iloc[0])
    vertical_offset    = float(df[df['metaKey'] == 'vertical offset Sig']['metaValue'].iloc[0])
    repetition_rate    = float(df[df['metaKey'] == 'Repetition rate (kHz)']['metaValue'].iloc[0]) * 1e+3
    wavelength         = float(df[df['metaKey'] == 'Laser Wavelength (nm)']['metaValue'].iloc[0]) * 1e-9
    polarization       = int(df[df['metaKey'] == 'Polarization']['metaValue'].iloc[0])
    sample_temperature = df[df['metaKey'] == 'Sample Temperature (K)']['metaValue'].iloc[0]
    attenuation = 0.03
    photon_number = calculate_pulse_photon(laser_power*1e-6,repetition_rate,attenuation,wavelength)
    # infile = ROOT.TFile.Open(in_filename)
    # inputDataName = infile.Get('inputDataName').GetTitle()
    # with TdmsFile.open(inputDataName) as tdms_file:
    #     # Read Meta Data (Basic information)
    #     metadata = tdms_file.properties
    #     metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
    #     laser_power = float(metadata_df.loc[metadata_df['metaKey'] == 'Laser Power (uW)', 'metaValue'].iloc[0])
    #     bias_voltage = float(metadata_df.loc[metadata_df['metaKey'] == 'Bias Voltage (mV)', 'metaValue'].iloc[0])
    #     bias_current = float(metadata_df.loc[metadata_df['metaKey'] == 'Bias Current (nA)', 'metaValue'].iloc[0])
    return laser_power, bias_voltage, bias_current, photon_number, polarization, sample_temperature, vertical_range, vertical_offset

def getkey(photon_number,bias_voltage,polarization,sample_temperature):
    key = str(photon_number) + '_' + str(bias_voltage) + 'mV_' + str(polarization) + 'degrees_' + str(sample_temperature) + 'K'
    return key

def sort_bias(bias, var):
    bias_array = np.array(bias)
    var_array = np.array(var)
    sorted_bias_array = np.sort(bias_array)
    sorted_var_array= var_array[bias_array.argsort()]
    return sorted_bias_array, sorted_var_array

def graph2list(graph):
    list=[]
    for i in range(graph.GetN()):
        list.append(graph.GetY()[i])
    return list

def Graph_sweep_polarization(photon, temp, bvs, bcs, Polars, stat, title="graph", ytit="", ymin=0, ymax=1, plotDir="./"):
    # outfile.cd()
    # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    fig, ax = plt.subplots()
    vals,pols=[],[]
    for ibv, bv in enumerate(bvs):
        for ipolar, polar in enumerate(Polars):
            key = getkey(photon,bv,polar,temp)
            try:
                value = stat[key]
                vals.append(value)
                pols.append(float(polar))
            except KeyError:
                pass
        print(pols,vals)
        if (len(pols)>2):
            # pols_rad = np.deg2rad(pols)  # Convert degrees to radians
            ax.plot(pols, vals, marker='o', markersize=6, fillstyle='none', alpha=0.75, label=f'{bv}mV')
        vals.clear()
        pols.clear()
    ax.grid(True)
    ax.set_title("Detection efficiency (%)", va='bottom')
    ax.legend()
    # ax.set_thetagrids(np.arange(0, 360, 30)) # Set theta gridlines with 30 degree separation
    plt.savefig('Polar.png', format='png')
    plt.close()

def Graph_sweep_bias(photon, polar, temp, bvs, bcs, stat, title="graph", ytit="", ymin=0, ymax=1, plotDir="./"):
    # outfile.cd()
    c1 = ROOT.TCanvas()
    leg = ROOT.TLegend(0.6,0.8,0.88,0.88)
    leg.SetNColumns(4)
    graph = ROOT.TGraph()
    index=0
    for ibv, bv in enumerate(bvs):
        key = getkey(photon,bv,polar,temp)
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
    leg.AddEntry(graph,f'{photon} photons','lp')
    leg.Draw()
    c1.SaveAs(f"{plotDir}/{title}_sweep_bias_current.png")

def Graph_sweep_bias_err(photon, polar, temp, bvs, bcs, stat, stat_err, title="graph", ytit="", ymin=0, ymax=1, plotDir="./"):
    # outfile.cd()
    c1 = ROOT.TCanvas()
    leg = ROOT.TLegend(0.6,0.8,0.88,0.88)
    leg.SetNColumns(4)
    graph = ROOT.TGraphErrors()
    index=0
    for ibv, bv in enumerate(bvs):
        key = getkey(photon,bv,polar,temp)
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
    leg.AddEntry(graph,f'{photon} photons','lp')
    leg.Draw()
    c1.SaveAs(f"{plotDir}/{title}_sweep_bias_current.png")

def Graph_sweep(polar, temp, photons, bvs, bcs, stat, title="graph", ytit="", plotDir="./"):
    fig, ax = plt.subplots()
    yvals,xvals=[],[]
    for i, bv in enumerate(bvs):
        for photon in photons:
            key = getkey(photon,bv,polar,temp)
            try:
                value = stat[key]
                yvals.append(value)
                xvals.append(photon)
            except KeyError:
                continue
        ax.plot(xvals, yvals, marker=markers(i), markersize=6, fillstyle='none', alpha=0.75, label=f'{bv}mV')
        xvals.clear()
        yvals.clear()
    ax.grid(True)
    ax.set_xlabel('Photon number',fontsize=15)
    ax.set_ylabel(ytit,fontsize=15)
    ax.set_title(f'{title}_sweep_photon',fontsize=15)
    ax.legend(title='Bias Voltage')
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig(f"{plotDir}/{title}_sweep_photon.png", format='png')
    print(f"{plotDir}/{title}_sweep_photon.png")
    plt.close()

    fig1, ax1 = plt.subplots()
    for i, photon in enumerate(photons):
        for ibv, bv in enumerate(bvs):
            key = getkey(photon,bv,polar,temp)
            try:
                value = stat[key]
                yvals.append(value)
                xvals.append(bv)
            except KeyError:
                continue
        ax1.plot(xvals, yvals, marker=markers(i), markersize=6, fillstyle='none', alpha=0.75, label=f'{photon}')
        xvals.clear()
        yvals.clear()
    ax1.grid(True)
    ax1.set_xlabel(r'Bias Current ($\mu$A)',fontsize=15)
    ax1.set_ylabel(ytit,fontsize=15)
    ax1.set_title(f'{title}_sweep_bias',fontsize=15)
    ax1.set_yscale('log')
    ax1.legend(title='Photon number')
    plt.tight_layout()
    plt.savefig(f"{plotDir}/{title}_sweep_bias.png", format='png')
    print(f"{plotDir}/{title}_sweep_bias.png")
    plt.close()

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

def fit_histo(hist, rangemin, rangemax, name, xTitle, title, saveTitle):
    c1 = ROOT.TCanvas()
    fit = ROOT.TF1("fit","gausn",rangemin, rangemax)
    fit.SetLineWidth(1)
    if (args.res):
        fit.SetParameter(0,100)
        fit.SetParameter(1,0.1)
        fit.SetParameter(2,0.02)
    else:
        fit.SetParameter(0,100)
        fit.SetParameter(1,0.2)
        fit.SetParameter(2,0.01)
    # fit.SetParameter(1,(mean_min+mean_max)/2)
    fitResults = hist.Fit("fit",'IQRS')
    mean = fit.GetParameter(1)
    mean_error = fit.GetParError(1)
    std = fit.GetParameter(2)
    std_error = fit.GetParError(2)
    const = fit.GetParameter(0)
    const_error = fit.GetParError(0)
    integral = fit.Integral(rangemin, rangemax)
    integral /= hist.GetBinWidth(1)
    # Draw hist
    hist.SetTitle(title)
    hist.GetXaxis().SetTitle(xTitle)
    hist.Draw()
    ROOT.gPad.Update()
    st = hist.GetListOfFunctions().FindObject("stats")
    st.AddText(f"Integral={integral}")
    fit.Draw("same")
    if (args.initPlot): c1.SaveAs(saveTitle)
    return mean, mean_error, std, std_error, const, const_error, integral, fitResults

def multi_histo_canvas(polar, temp, photons,bvs,bcs,histos,xtit,plotDir):
    c_multi_bv = {}
    for ibv, bv in enumerate(bvs):
        c_multi_bv[bv] = ROOT.TCanvas(f"c_multi_{bv}",f"c_multi_{bv}",1800,900)
        c_multi_bv[bv].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(photons)/5) if int(len(photons)%5==0) else int(len(photons)/5)+1
        c_multi_bv[bv].Divide(cx,cy,0,0)
        index=0
        for ipow, photon in enumerate(photons):
            key = getkey(photon,bv,polar,temp)
            try:
                pad = c_multi_bv[bv].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle(xtit)
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{photon}photons_{bcs[ibv]}uA")
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
        c_multi_bv[bv].SaveAs(f"{plotDir}/histos_sweep_photon_{bv}mV.png")

    c_multi_photon = {}
    for iphoton, photon in enumerate(photons):
        c_multi_photon[photon] = ROOT.TCanvas(f"c_multi_{photon}",f"c_multi_{photon}",1800,900)
        c_multi_photon[photon].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
        c_multi_photon[photon].Divide(cx,cy,0,0)
        index=0
        for ibv, bv in enumerate(bvs):
            key = getkey(photon,bv,polar,temp)
            try:
                pad = c_multi_photon[photon].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle(xtit)
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{int(photon)}photons_{bcs[ibv]}uA")
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
        c_multi_photon[photon].SaveAs(f"{plotDir}/histos_sweep_bv_{photon}.png")

def multi_histo_canvas_sweep_bias(photon, polar, temp, bvs,bcs,histos,title,plotDir):
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
        key = getkey(photon,bv,polar,temp)
        try:
            histos[key].GetXaxis().SetTitle("")
            histos[key].GetYaxis().SetTitle("")
            histos[key].GetXaxis().SetLabelSize(0.1)
            histos[key].GetYaxis().SetLabelSize(0.1)
            histos[key].SetTitle("")
            histos[key].SetName(f"{int(photon)}photons_{bcs[ibv]}uA")
            pad = c_multi.cd(ibv+1)
            histos[key].Draw()
            ROOT.gPad.Update()
            ROOT.gStyle.SetOptStat("nemr")
            stat = histos[key].FindObject("stats")
            stat.SetY1NDC(0.55)
            stat.SetY2NDC(0.99)
            stat.SetX1NDC(0.55)
            stat.SetX2NDC(0.99)
            stat.AddText(f"Integral={pulse_integrals[key]*100:.0f}%")
            stat.SetStatFormat("6.2g")
            stat.SetFillStyle(0)
            ROOT.gPad.Modified()
            ROOT.gPad.Draw()
            pad = c_multi_log.cd(ibv+1)
            pad.SetLogy()
            histos[key].Draw()
            ROOT.gPad.Update()
            ROOT.gStyle.SetOptStat("nemr")
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
    c_multi.SaveAs(f"{plotDir}/histos_{title}_sweep_bv_{photon}.png")
    c_multi_log.SaveAs(f"{plotDir}/histos_{title}_sweep_bv_{photon}_log.png")

def multi_histo_canvas_fit(polar, temp, photons,bvs,bcs,histos,plotDir):
    c_multi_bv = {}
    for ibv, bv in enumerate(bvs):
        c_multi_bv[bv] = ROOT.TCanvas(f"c_multi_{bv}",f"c_multi_{bv}",1800,900)
        c_multi_bv[bv].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(photons)/5) if int(len(photons)%5==0) else int(len(photons)/5)+1
        c_multi_bv[bv].Divide(cx,cy,0,0)
        index=0
        for ipow, photon in enumerate(photons):
            key = getkey(photon,bv,polar,temp)
            try:
                pad = c_multi_bv[bv].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{photon}_{bcs[ibv]}uA")
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
        c_multi_bv[bv].SaveAs(f"{plotDir}/histos_sweep_photon_{bv}mV_fit.png")

    c_multi_photon = {}
    for iphoton, photon in enumerate(photons):
        c_multi_photon[photon] = ROOT.TCanvas(f"c_multi_{photon}",f"c_multi_{photon}",1800,900)
        c_multi_photon[photon].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
        c_multi_photon[photon].Divide(cx,cy,0,0)
        index=0
        for ibv, bv in enumerate(bvs):
            key = getkey(photon,bv,polar,temp)
            try:
                pad = c_multi_photon[photon].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{photon}photons_{bcs[ibv]}uA")
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
        c_multi_photon[photon].SaveAs(f"{plotDir}/histos_sweep_bv_{photon}_fit.png")

def multi_histo_canvas_fitlimit(polar, temp, photons,bvs,bcs,histos,integrals,plotDir):
    c_multi_bv = {}
    for ibv, bv in enumerate(bvs):
        c_multi_bv[bv] = ROOT.TCanvas(f"c_multi_{bv}",f"c_multi_{bv}",1800,900)
        c_multi_bv[bv].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(photons)/5) if int(len(photons)%5==0) else int(len(photons)/5)+1
        c_multi_bv[bv].Divide(cx,cy,0,0)
        index=0
        for ipow, photon in enumerate(photons):
            key = getkey(photon,bv,polar,temp)
            try:
                pad = c_multi_bv[bv].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{int(photon)}photons_{bcs[ibv]}uA")
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
        c_multi_bv[bv].SaveAs(f"{plotDir}/histos_sweep_photon_{bv}mV_fitlimit.png")

    c_multi_photon = {}
    for iphoton, photon in enumerate(photons):
        c_multi_photon[photon] = ROOT.TCanvas(f"c_multi_{photon}",f"c_multi_{photon}",1800,900)
        c_multi_photon[photon].SetFixedAspectRatio(True)
        ROOT.gStyle.SetPadBorderMode(0)
        cx = 5
        cy = int(len(bvs)/6) if int(len(bvs)%6==0) else int(len(bvs)/6)+1
        c_multi_photon[photon].Divide(cx,cy,0,0)
        index=0
        for ibv, bv in enumerate(bvs):
            key = getkey(photon,bv,polar,temp)
            try:
                pad = c_multi_photon[photon].cd(index+1)
                # pad.SetLogy()
                histos[key].GetXaxis().SetTitle("")
                histos[key].GetYaxis().SetTitle("")
                histos[key].GetXaxis().SetLabelSize(0.1)
                histos[key].GetYaxis().SetLabelSize(0.1)
                histos[key].SetTitle("")
                histos[key].SetName(f"{int(photon)}photons_{bcs[ibv]}uA")
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

def spectrum_sweep_bias(photon, polar, temp, bvs, bcs, graphs, plotDir="./"):
    xvals=[]
    for i in range(1000):
        xvals.append(i)
    fig, ax = plt.subplots()
    pad_xnum = 5
    fig_multi, ax_multi = plt.subplots(int(len(bvs)/pad_xnum+1), pad_xnum, figsize=(36, 6*(int(len(bvs)/pad_xnum)+1)), sharex=True, sharey=True)
    index=0
    for ibv, bv in enumerate(bvs):
        key = getkey(photon,bv,polar,temp)
        try:
            ax.plot(xvals, graphs[key], marker='o', markersize=6, fillstyle='none', alpha=0.75, label=f'{bv}mV')
            # ax multi
            ax1 = ax_multi.flatten()[index]
            ax1.plot(xvals, graphs[key], marker='o', markersize=6, fillstyle='none', alpha=0.75, label=f'{bv}mV')
            ax1.tick_params(axis='both', which='major', labelsize=20) # Set larger ticks and labels
            ax1.legend(fontsize=30)
            ax1.grid()
            if index // pad_xnum == int(len(bvs)/pad_xnum):  ax1.set_xlabel('Index (0.4ns)', fontsize=30) # Add x and y titles only to the leftmost and bottom-most subplot
            if index % pad_xnum == 0:  ax1.set_ylabel('Voltage (V)', fontsize=30)
            ax1.set_xlim(300,340)
            ax1.legend()
            index+=1
        except KeyError:
            print(key)
            continue
    fig_multi.tight_layout()
    fig_multi.savefig(f"{plotDir}/spectrums_multi.png", format='png')
    print(f"{plotDir}/spectrums_multi.png")

    ax.grid(True)
    ax.set_xlabel('Index (0.4ns)',fontsize=15)
    ax.set_ylabel('Voltage (V)',fontsize=15)
    ax.set_xlim(300,340)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{plotDir}/spectrums.png", format='png')
    print(f"{plotDir}/spectrums.png")

def gettree():
    for i, in_filename in enumerate(args.in_filenames):
        print (f"{i}/{len(args.in_filenames)}: {in_filename}")
        laser_power, bias_voltage, bias_current, photon_number, polarization, sample_temperature, vertical_range, vertical_offset = get_info(in_filename)
        # if (bias_voltage > 530): continue
        basename = getkey(photon_number,bias_voltage,polarization,sample_temperature)
        plotDir= in_filename.rsplit("/",1)[0]
        infile = ROOT.TFile.Open(in_filename)
        # Get graphs
        Pulse_avg_display=infile.Get('Pulse_avg_display')
        # Get Tree
        intree = infile.Get('Result_tree')
        # Initialize variables
        if (vertical_range==0.05): nbin=pow(2,6)
        elif (vertical_range==0.1): nbin=pow(2,6.6)
        elif (vertical_range>0.2 and vertical_range<5): nbin=pow(2,7)
        else:
            print("Invalid vertical range!!!")
            break
        range_min, range_max= 0, vertical_range
        binsize = float((range_max-range_min)/nbin)
        if (args.res):
            fit_range_min = range_min
            fit_range_max = range_max
        else:
            fit_range_min = 0.1
            fit_range_max = range_max
        # initialize histo
        h_pulse_fall_range = ROOT.TH1F(f"h_pulse_fall_range_{bias_current}",f"h_pulse_fall_range_{bias_current}",nbin,range_min,range_max)
        h_pulse_fall_range_limitfit = ROOT.TH1F(f"h_pulse_fall_range_{bias_current}_limitfit",f"h_pulse_fall_range_{bias_current}_limitfit",nbin,range_min,range_max)
        h_pulse_fall_tau = ROOT.TH1F("h_pulse_fall_tau","h_pulse_fall_tau",20,0,5)
        h_pulse_FWHM = ROOT.TH1F("h_pulse_FWHM","h_pulse_FWHM",20,0,5)
        h_pre_range = ROOT.TH1F("h_pre_range","h_pre_range",15,0.02,0.17)
        h_eff = ROOT.TH1F("h_eff","h_eff",2,0,2)
        h_diff = ROOT.TH1F("h_diff","h_diff",100,0,0.3)
        # initialize graph
        g_pulse_fall_range = ROOT.TGraph()
        # Project variables to histos
        try:
            project(intree,h_pulse_fall_range,"pulse_fall_range_ptp","",basename,"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range")
        except AttributeError:
            continue
        project(intree,h_pulse_fall_range_limitfit,"pulse_fall_range","",basename,"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range")
        project(intree,h_pulse_fall_tau,"pulse_fall_tau","pulse_fall_tau>0",basename,"pulse fall time constant (0.4ns)",f"Event",plotDir,"h_pulse_fall_tau")
        project(intree,h_pulse_FWHM,"pulse_FWHM","pulse_FWHM>0 && pulse_FWHM<100",basename,"pulse fall time constant (0.4ns)",f"Event",plotDir,"h_pulse_FWHM")
        project(intree,h_pre_range,"pre_range","",basename,"pre_range (V)","Event",plotDir,"h_pre_range")
        project(intree,h_eff,"1","pulse_fall_range>0.1",basename,"Pulse detected","Event",plotDir,"h_eff")
        # Rebin
        # h_pulse_fall_range_rebin1 = rebin(h_pulse_fall_range,f'{basename}_rebin',"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range_rebin1",True)
        # h_pulse_fall_range_rebin2 = rebin(h_pulse_fall_range_rebin1,f'{basename}_rebin',"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range_rebin2",True)
        # h_pulse_fall_range_rebin3 = rebin(h_pulse_fall_range_rebin2,f'{basename}_rebin',"pulse_range (V)",f"Event/{(range_max-range_min)/nbin:.4f}V",plotDir,"h_pulse_fall_range_rebin3",True)
        # Fit
        mean, mean_error, std, std_error, const, const_error, integral, fitResults = fit_histo(h_pulse_fall_range, fit_range_min, fit_range_max, f"fit_pulse_fall_range_{basename}", 'pulse_range (V)', "", f"{plotDir}/fit_pulse_fall_range.png")
        # Graph
        if (args.initPlot):
            for i,entry in enumerate(intree):
                g_pulse_fall_range.SetPoint(i,i,entry.pulse_fall_range)
            c_graph = ROOT.TCanvas()
            g_pulse_fall_range.GetXaxis().SetTitle("Event")
            g_pulse_fall_range.GetYaxis().SetTitle("Pulse Fall Range (V)")
            g_pulse_fall_range.Draw("AP")
            c_graph.SaveAs(f"{plotDir}/fall_range_event.png")
        # Calculate
        if (const_error>10):
            eff = h_eff.Integral()/intree.GetEntries()
            pulse_range = h_pulse_fall_range.GetMean()
            pulse_range_error = h_pulse_fall_range.GetRMS()
            pulse_integrals[basename] = 0
        else:
            eff = h_eff.Integral()/intree.GetEntries()
            # eff = integral/intree.GetEntries()
            pulse_range = mean
            pulse_range_error = std
            pulse_integrals[basename] = integral
        pre_range = h_pre_range.GetMean()
        pulse_fall_tau = h_pulse_fall_tau.GetMean()
        pulse_FWHM = h_pulse_FWHM.GetMean()
        pre_range_err = h_pre_range.GetRMS()
        try:
            pulse_range_stderror = h_pulse_fall_range.GetRMS()/math.sqrt(h_pulse_fall_range.Integral())
        except ZeroDivisionError:
            pulse_range_stderror = 0
        # Append sweep variables
        if photon_number not in Photons:
            Photons.append(photon_number)
        if bias_voltage not in BVs:
            BVs.append(bias_voltage)
            BCs.append(bias_current)
        if polarization not in Polars:
            Polars.append(polarization)
        # if sample_temperature not in Temps:
        Temps.append(sample_temperature)
        # Fill stats dict
        effs[basename] = eff
        pulse_ranges[basename] = pulse_range
        pulse_range_errs[basename] = pulse_range_error
        pre_ranges[basename] = pre_range
        pre_range_errs[basename] = pre_range_err
        pulse_fall_taus[basename] = pulse_fall_tau
        pulse_FWHMs[basename] = pulse_FWHM
        spectrums[basename]=graph2list(Pulse_avg_display)
        avgMaxs[basename]=max(spectrums[basename])
        avgMins[basename]=min(spectrums[basename])
        avgs[basename]=max(spectrums[basename])-min(spectrums[basename])
        range_avgs[basename]=pulse_ranges[basename]-avgs[basename]
        # Histograms
        h_pulse_fall_taus[basename] = h_pulse_fall_tau.Clone()
        h_pulse_FWHMs[basename] = h_pulse_FWHM.Clone()
        h_pre_ranges[basename] = h_pre_range.Clone()
        h_pulse_fall_ranges[basename] = h_pulse_fall_range.Clone()
        h_pulse_fall_ranges_limitfit[basename] = h_pulse_fall_range_limitfit.Clone()
        h_pulse_fall_taus[basename].SetDirectory(0)
        h_pulse_FWHMs[basename].SetDirectory(0)
        h_pulse_fall_ranges[basename].SetDirectory(0)
        h_pre_ranges[basename].SetDirectory(0)
        h_pulse_fall_ranges_limitfit[basename].SetDirectory(0)
        print(f"{basename}: {eff*100:.1f}%, {pulse_range*1000:.1f}mV+-{pulse_range_error*1000:.2f}mV, status:{fitResults.Status()}")

def plots():
    print("\n==================== Start plotting ====================")
    # Sort sweep variables
    Photons.sort()
    print(Photons)
    BVs.sort()
    BCs.sort()
    Polars.sort()
    # Graphs of stats vs sweep. variables
    if (args.sweepBias):
        sweepBias_plotDir = args.in_filenames[0].rsplit('/',2)[0]
        Graph_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,effs,title="g_eff",ytit="Pulse Detection Efficiency",ymin=0,ymax=1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias_err(Photons[0],Polars[0],Temps[0],BVs,BCs,pulse_ranges,pulse_range_errs, title="g_pulse_range",ytit="Pulse range mean (V)",ymin=0,ymax=max(pulse_ranges.values())*1.6,plotDir=sweepBias_plotDir)
        Graph_sweep_bias_err(Photons[0],Polars[0],Temps[0],BVs,BCs,pre_ranges,pre_range_errs,title="g_pre_range",ytit="Pre range mean (V)",ymin=min(pre_ranges.values())*0.8,ymax=max(pre_ranges.values())*1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,pulse_integrals,title="g_integral",ytit="Signal integrals",ymin=min(pulse_integrals.values())*0.8,ymax=max(pulse_integrals.values())*1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,pulse_fall_taus,title="g_fall_tau",ytit="Fall time constant (0.4ns)",ymin=min(pulse_fall_taus.values())*0.8,ymax=max(pulse_fall_taus.values())*1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,pulse_FWHMs,title="g_FWHM",ytit="FWHM (0.4ns)",ymin=min(pulse_FWHMs.values())*0.8,ymax=max(pulse_FWHMs.values())*1.2,plotDir=sweepBias_plotDir)
        Graph_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,avgs,title="g_avg",ytit="Average Pulse Range (V)",ymin=min(avgMaxs.values())*0.8,ymax=max(avgMaxs.values())*1.2,plotDir=sweepBias_plotDir)
        spectrum_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,spectrums,sweepBias_plotDir)
        # multi_spectrum_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,spectrums,sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pulse_fall_ranges,"range",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pre_ranges,"pre_range",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pulse_fall_taus,"tau",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pulse_FWHMs,"FWHM",sweepBias_plotDir)
    if (args.sweepAll):
        sweepAll_plotDir = args.in_filenames[0].rsplit('/',5)[0]
        Graph_sweep(Polars[0],Temps[0],Photons,BVs,BCs,effs,title="g_eff",ytit="Pulse Detection Efficiency",plotDir=sweepAll_plotDir)
        Graph_sweep(Polars[0],Temps[0],Photons,BVs,BCs,pulse_ranges,title="g_pulse_range",ytit="Pulse range mean (V)",plotDir=sweepAll_plotDir)
        Graph_sweep(Polars[0],Temps[0],Photons,BVs,BCs,pre_ranges,title="g_pre_range",ytit="Pre range mean (V)",plotDir=sweepAll_plotDir)
        Graph_sweep(Polars[0],Temps[0],Photons,BVs,BCs,avgs,title="g_avg",ytit="Average Pulse Range (V)",plotDir=sweepAll_plotDir)
        Graph_sweep(Polars[0],Temps[0],Photons,BVs,BCs,avgMaxs,title="g_avg_Max",ytit="Average Pulse Max (V)",plotDir=sweepAll_plotDir)
        Graph_sweep(Polars[0],Temps[0],Photons,BVs,BCs,range_avgs,title="g_Range_Avg",ytit="Pulse Range Mean - Average (V)",plotDir=sweepAll_plotDir)
        multi_histo_canvas(Polars[0],Temps[0],Photons,BVs,BCs,h_pulse_fall_ranges,"Pulse range (V)",sweepAll_plotDir)
    if (args.sweepPolarization):
        sweepPolar_plotDir = args.in_filenames[0].split('degrees/')[0].split('/')[0]
        Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,effs,title="g_eff",ytit="Pulse Detection Efficiency",ymin=0,ymax=1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization_err(Photons[0],BVs[0],Temps[0],Polars,pulse_ranges,pulse_range_errs, title="g_pulse_range",ytit="Pulse range mean (V)",ymin=0,ymax=max(pulse_ranges.values())*1.6,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization_err(Photons[0],BVs[0],Temps[0],Polars,pre_ranges,pre_range_errs,title="g_pre_range",ytit="Pre range mean (V)",ymin=min(pre_ranges.values())*0.8,ymax=max(pre_ranges.values())*1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,pulse_integrals,title="g_integral",ytit="Signal integrals",ymin=min(pulse_integrals.values())*0.8,ymax=max(pulse_integrals.values())*1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,pulse_fall_taus,title="g_fall_tau",ytit="Fall time constant (0.4ns)",ymin=min(pulse_fall_taus.values())*0.8,ymax=max(pulse_fall_taus.values())*1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,pulse_FWHMs,title="g_FWHM",ytit="FWHM (0.4ns)",ymin=min(pulse_FWHMs.values())*0.8,ymax=max(pulse_FWHMs.values())*1.2,plotDir=sweepPolar_plotDir)

if __name__ == "__main__":
    ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
    # ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
    Photons,BVs,BCs,Polars,Temps = [],[],[],[],[] # List for sweep variables
    effs, pulse_ranges, pulse_range_errs, pulse_fall_taus, pulse_FWHMs, pulse_integrals, pre_ranges, pre_range_errs, noSignals, noSignals={},{},{},{},{},{},{},{},{},{} # List for stats
    avgMaxs,avgMins,avgs,range_avgs={},{},{},{}
    h_pre_ranges,h_pulse_fall_ranges,h_pulse_fall_ranges_limitfit,h_pulse_fall_taus,h_pulse_FWHMs={},{},{},{},{} # List of histos
    spectrums={} #List of graphs
    gettree() # loop over the input files
    plots() # Plot them together
