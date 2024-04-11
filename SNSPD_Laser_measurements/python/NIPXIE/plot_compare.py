#!/usr/bin/env python3

import ROOT
from array import array
import argparse
import matplotlib.pyplot as plt
import numpy as np

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

def getgraph(in_filename):
    plotDir= in_filename.rsplit("/",1)[0]
    infile = ROOT.TFile.Open(in_filename)
    g_eff = infile.Get('g_eff')
    g_pulse_range = infile.Get('g_pulse_range')
    g_pre_range = infile.Get('g_pre_range')
    return g_eff, g_pulse_range, g_pre_range

def plot_compare():
    g_effs, g_pulse_ranges, g_pres, powers=[],[],[],[]
    for in_filename in args.in_filenames:
        power = in_filename.split('nW/')[0].split('/')[-1]
        print(power)
        g_eff, g_pulse_range, g_pre = getgraph(in_filename)
        g_effs.append(g_eff)
        g_pres.append(g_pre)
        g_pulse_ranges.append(g_pulse_range)
        powers.append(power)

    c1 = ROOT.TCanvas()
    plot_Dir = args.in_filenames[0].rsplit("/",2)[0]
    leg_eff = ROOT.TLegend(0.2,0.8,0.8,0.9)
    leg_eff.SetNColumns(3)
    for i, (g_eff,power) in enumerate(zip(g_effs,powers)):
        if(i==0):
            g_eff.Draw("AP PMC")
            g_eff.GetYaxis().SetRangeUser(0,1.2)
        else: g_eff.Draw("PSame PMC")
        leg_eff.AddEntry(g_eff,f'{power}nW','p')
        # g_eff.SetMarkerColor(color(i))
        g_eff.SetMarkerStyle(i+20)
    leg_eff.Draw()
    c1.SaveAs("eff_compare.png")

    for i, g_pulse_range in enumerate(g_pulse_ranges):
        if(i==0):
            g_pulse_range.Draw("AP PMC")
            g_pulse_range.GetYaxis().SetRangeUser(0,0.35)
        else: g_pulse_range.Draw("PSame PMC")
        # g_pulse_range.SetMarkerColor(color(i))
        g_pulse_range.SetMarkerStyle(i+20)
    leg_eff.Draw()
    c1.SaveAs("amp_compare.png")

    leg_pre = ROOT.TLegend(0.2,0.8,0.8,0.9)
    leg_pre.SetNColumns(3)
    for i, (g_pre,power) in enumerate(zip(g_pres,powers)):
        if(i==0):
            g_pre.Draw("AP PMC")
            g_pre.GetYaxis().SetRangeUser(0.08,0.1)
        else: g_pre.Draw("PSame PMC")
        leg_pre.AddEntry(g_pre,f'{power}nW','p')
        # g_pre.SetMarkerColor(color(i))
        g_pre.SetMarkerStyle(i+20)
    leg_pre.Draw()
    c1.SaveAs("pre_compare.png")

if __name__ == "__main__":
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    plot_compare()
