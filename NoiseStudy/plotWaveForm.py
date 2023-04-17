#!/usr/bin/env python

import os
import errno
import ROOT
from array import array
import csv

import argparse
parser = argparse.ArgumentParser(description='draw stacked waveforms')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
args = parser.parse_args()

c1 = ROOT.TCanvas()
g_waveform_1={}
g_waveform_2={}
xaxis={}
ch1={}
ch2={}
for index,infile in enumerate(args.in_filenames):
    with open (infile) as f:
        Lines = f.readlines()

    Lines.pop()
    xaxis[index], ch1[index], ch2[index] = array('d'), array('d'), array('d')

    for i, Line in enumerate(Lines[2:],2):
        x, v1, v2 = Line.split(',')
        xaxis[index].append(float(x)*10e6)
        ch1[index].append(float(v1))
        ch2[index].append(float(v2))

    if index > 0:
        for sample, v in enumerate(ch1[0]):
            ch1[0][sample] += ch1[index][sample]
            # ch2[0][sample] += ch2[index][sample]

        g_waveform_1[index-1]=ROOT.TGraph(len(xaxis[0]), xaxis[0], ch1[0])
        g_waveform_2[index-1]=ROOT.TGraph(len(xaxis[0]), xaxis[0], ch2[0])

        g_waveform_1[index-1].SetTitle("%d File(s) Added" % index)
        g_waveform_1[index-1].GetXaxis().SetTitle("Time(us)")
        g_waveform_1[index-1].GetYaxis().SetTitle("ADC(V)")
        g_waveform_1[index-1].SetMarkerColor(2)
        g_waveform_1[index-1].Draw("APE")
        g_waveform_1[index-1].GetXaxis().SetRangeUser(-3,3)
        g_waveform_1[index-1].GetYaxis().SetRangeUser(-5,5)

        g_waveform_2[index-1].SetMarkerColor(4)
        g_waveform_2[index-1].Draw("PSame")
        ROOT.gPad.WaitPrimitive()
