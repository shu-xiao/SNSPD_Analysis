#!/usr/bin/env python3
# NI-5162 specs: https://www.ni.com/docs/en-US/bundle/pxie-5162-specs/page/specs.html

# python3 -m python.NIPXIE.plot_multiHistograms --sweepAll --res plots/SNSPD_5/Laser/Ch3/20240504/11p47K/Pulse/515/10000kHz/*nW/240degrees/*/*/*.root
# python3 -m python.NIPXIE.plot_multiHistograms --sweepAll plots/SNSPD_5/Laser/Ch3/20240503/4p68/Pulse/515/10000kHz/*nW/240degrees/*/*/*.root

import numpy as np
np.seterr(all='ignore')
import ROOT
from array import array
import argparse
import math
import json
import pandas as pd

from ..utils.plotUtils import markers, savefig
from ..utils.osUtils import createDir
from ..utils.fitfunctions import alt_expo
from ..utils.Results import *

import warnings

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--debug','-d',action="store_true",help='debug mode')
parser.add_argument('--initPlot','-i',action="store_true",help='save histograms')
parser.add_argument('--sweepBias',action="store_true",help='sweep bias mode')
parser.add_argument('--sweepAll',action="store_true",help='Run all plots for each sweep variable')
parser.add_argument('--sweepPolarization',action="store_true",help='Run plots for polarization')
parser.add_argument('--fit','-f',action="store_true",help='do fit')
parser.add_argument('--res',action="store_true",help='do resistivie mode analysis')
args = parser.parse_args()

def getkey(photon_number,bias_voltage,polarization,sample_temperature):
    key = str(photon_number) + '_' + str(bias_voltage) + 'mV_' + str(polarization) + 'degrees_' + str(sample_temperature) + 'K'
    return key

def Normalized_Graph_sweep(statname, ytit="", plotDir="./"):
    createDir(f"{plotDir}/Norm")
    fig, ax = plt.subplots()
    yvals,xvals=[],[]
    photon_subset = [699, 777, 855, 933, 1011, 1088, 1166, 1244, 1322, 1399, 1477, 1555, 1788, 2333] if not args.res else [1555, 3111, 4666, 6222, 7777, 9333, 11666, 15555, 23333]
    bias_subset = [800,850,900,950,1000,1030] if not args.res else [84, 94, 104, 124, 134, 144, 154, 164, 174]
    index=0
    for bv,n_bv in zip(Vars["bv"],Vars["n_bv"]):
        if (bv not in bias_subset): continue
        for (pho,n_pho) in zip(Vars["pho"],Vars["n_pho"]):
            key = getkey(pho,bv,Vars["pol"][0],Vars["temp"][0])
            try:
                value = Stats[key][statname]
                yvals.append(value)
                xvals.append(n_pho)
            except KeyError:
                print(key)
                continue
        ax.plot(xvals, yvals, marker=markers(index), markersize=8, fillstyle='full', linestyle='-', alpha=1, label=f'{n_bv:.2f}')
        index+=1
        xvals.clear()
        yvals.clear()
    # ax.grid(True)
    if(args.res): ax.set_xlabel(r'Laser Intensity ($P$ / $P_0^{cal}$)',fontsize=15)
    else: ax.set_xlabel(r'Laser Intensity ($P$ / $P_0^{sc})$',fontsize=15)
    ax.set_ylabel(ytit,fontsize=15)
    # ax.set_title(f'{title}_sweep_intensity',fontsize=15)
    handles, labels = ax.get_legend_handles_labels() # Get the current handles and labels from the axes
    handles = handles[::-1] # Reverse the order of handles and labels
    labels = labels[::-1]
    Ibias=r'$I_{bias}$'
    Isw=r'$I_{sw}$'
    if (Vars["temp"][0]=="11p47K"): legend = ax.legend(handles, labels, title=f'{Ibias} / {Isw}({SampleT_norm:.2f}Tc)')
    else: legend = ax.legend(handles, labels, title=f'{Ibias} / {Isw}({SampleT_norm:.2f}Tc)',loc='upper right', bbox_to_anchor=(0.93, 0.98))
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.tight_layout()
    savefig(fig,f"{plotDir}/Norm/{statname}_sweep_intensity")
    ax.set_yscale('log')
    savefig(fig,f"{plotDir}/Norm/{statname}_sweep_intensity_log")
    plt.close("all")

    fig1, ax1 = plt.subplots()
    for i, (pho,n_pho) in enumerate(zip(Vars["pho"],Vars["n_pho"])):
        for ibv, (bv,n_bv) in enumerate(zip(Vars["bv"],Vars["n_bv"])):
            if (bv not in bias_subset): continue
            key = getkey(pho,bv,Vars["pol"][0],Vars["temp"][0])
            try:
                value = Stats[key][statname]
                yvals.append(value)
                xvals.append(n_bv)
            except KeyError:
                print(key)
                continue
        ax1.plot(xvals, yvals, marker=markers(i), markersize=8, fillstyle='full', alpha=1, label=f'{n_pho}x')
        xvals.clear()
        yvals.clear()
    # ax1.grid(True)
    ax1.set_xlabel(r'Bias Current ($\mu$A)',fontsize=15)
    ax1.set_ylabel(ytit,fontsize=15)
    ax1.legend(title='Laser Intensity')
    fig1.tight_layout()
    savefig(fig,f"{plotDir}/Norm/{statname}_sweep_bias")
    plt.close()

def Normalized_Graph_sweep_singleIntensity(pho, statname, ytit="", plotDir="./"):
    createDir(f"{plotDir}/Norm")
    fig, ax = plt.subplots()
    yvals,xvals=[],[]
    index=0
    for bv,n_bv in zip(Vars["bv"],Vars["n_bv"]):
        key = getkey(pho,bv,Vars["pol"][0],Vars["temp"][0])
        try:
            value = Stats[key][statname]
            yvals.append(value)
            xvals.append(n_bv)
        except KeyError:
            print(key)
            continue
    ax.plot(xvals, yvals, marker="o", markersize=8, fillstyle='full', linestyle='-', alpha=1)
    ax.set_xlabel(r'Normalized Bias Current ($I_{bias}$ / $I_{sw}$)',fontsize=15)
    ax.set_ylabel(ytit,fontsize=15)
    # ax.set_title(f'{title}_sweep_intensity',fontsize=15)
    plt.tight_layout()
    savefig(fig,f"{plotDir}/Norm/{statname}_sweep_bias_{pho}")
    plt.close('all')

def Resolution(plotDir):
    createDir(f"{plotDir}/Reso/")
    fig, ax = plt.subplots()
    index=0
    for i, (bv,n_bv) in enumerate(zip(Vars["bv"],Vars["n_bv"])):
        if (bv != 174): continue
        yvals,xvals=[],[]
        for (pho,n_pho) in zip(Vars["pho"],Vars["n_pho"]):
            if (pho < 2500): continue
            key = getkey(pho,bv,Vars["pol"][0],Vars["temp"][0])
            try:
                mean = Stats[key]["avg"]
                sigma = Stats[key]["pulse_range_ptp_err"] / Stats[key]["avg"]
                xvals.append(mean)
                yvals.append(sigma)
            except KeyError:
                print(key)
                continue
        ax.plot(xvals, yvals, marker=markers(i), markersize=8, fillstyle='full', linestyle='-', alpha=1, label=n_bv)
    ax.set_xlabel('Signal Amplitude (V)',fontsize=15)
    ax.set_ylabel('Relative Resolution',fontsize=15)
    ax.legend()
    plt.tight_layout()
    savefig(fig,f"{plotDir}/Reso/resolution")
    plt.close('all')

def multi_histo_matplotlib(histname,xtit,plotDir):
    photon_subset = [1.0, 3.0, 5.0, 7.5, 10.0, 15.0] if args.res else [1.0, 1.3, 1.5, 2.0]
    for (bv,n_bv) in zip(Vars["bv"],Vars["n_bv"]):
        fig, ax = plt.subplots()
        for (pho,n_pho) in zip(Vars["pho"],Vars["n_pho"]):
            if (n_pho not in photon_subset): continue
            key = getkey(pho,bv,Vars["pol"][0],Vars["temp"][0])
            try:
                [lst1, lst2] = Stats[key][histname]
                bin_contents = np.array(lst1)
                bin_edges = np.array(lst2)
                non_zero_bins = np.where(bin_contents > 0)[0] # Find the rightmost bin with non-zero content
                if len(non_zero_bins) > 0:
                    min_non_zero_bin_edge = bin_edges[non_zero_bins[0] + 1]
                    max_non_zero_bin_edge = bin_edges[non_zero_bins[-1] + 1]
                else:
                    max_non_zero_bin_edge = bin_edges[-1]
                ax.hist(bin_edges[:-1], bins=bin_edges, weights=bin_contents, histtype='step', linewidth=2, label=f'{n_pho}x')
            except KeyError:
                print(key)
                continue
        # Customize the plot
        ax.set_yscale('log')
        if (args.res):
            ax.set_xlim(0, max_non_zero_bin_edge*1.2)
            ax.set_ylim(top=0.4)
        else: ax.set_xlim(-0.1, max_non_zero_bin_edge*1.2)
        ax.set_xlabel(xtit, fontsize=15)
        ax.set_ylabel('Probability', fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=10)
        if(args.res): legend = ax.legend(title=r'$P$ / $P_0^{cal}$')
        else: legend = ax.legend(title=r'$P$ / $P_0^{sc}$')
        legend.get_title().set_fontsize('large')  # Adjust the font size of the legend title
        Ibias=r'$I_{bias}$'
        Isw=r'$I_{sw}$'
        ax.text(0.3, 0.95, f'{Ibias} = {n_bv}{Isw}({SampleT_norm:.2f}Tc)', transform=ax.transAxes, verticalalignment='top', fontsize=15)
        fig.tight_layout()
        savefig(fig,f"{plotDir}/multi/histos_sweep_photon_{bv}mV")
        fig.clf()  # Clear the figure
        plt.close(fig)

def spectrum_sweep_bias(plotDir="./"):
    photon_subset = [699, 777, 855, 933, 1011, 1088, 1166, 1244, 1322, 1399, 1477, 1555, 1788, 2333] if not args.res else [1555, 3111, 4666, 6222, 7777, 9333, 11666, 15555, 23333]
    bias_subset = [800,850,900,950,1000,1030] if not args.res else [84, 94, 104, 124, 134, 144, 154, 164, 174]
    createDir(f"{plotDir}/spectrums")
    xvals=[]
    for i in range(1000):
        xvals.append(i*0.4-122)
    xmin,xmax = 0,18
    for (pho,n_pho) in zip(Vars["pho"],Vars["n_pho"]):
        fig, ax = plt.subplots()
        pad_xnum = 5
        # fig_multi, ax_multi = plt.subplots(int(len(Vars["bv"])/pad_xnum+1), pad_xnum, figsize=(36, 6*(int(len(Vars["bv"])/pad_xnum)+1)), sharex=True, sharey=True)
        index=0
        for (bv,n_bv) in zip(Vars["bv"],Vars["n_bv"]):
            if (bv not in bias_subset): continue
            key = getkey(pho,bv,Vars["pol"][0],Vars["temp"][0])
            try:
                graph = Stats[key]["spectrum"]
            except KeyError:
                print(key)
                continue
            ax.plot(xvals, graph, marker=markers(index), markersize=6, fillstyle='none', alpha=0.75, label=f'{n_bv}')
            # ax multi
            # ax1 = ax_multi.flatten()[index]
            # ax1.plot(xvals, graph, marker='o', markersize=6, fillstyle='none', alpha=0.75, label=f'{bv}mV')
            # ax1.tick_params(axis='both', which='major', labelsize=20) # Set larger ticks and labels
            # ax1.legend(fontsize=30)
            # ax1.grid()
            # if index // pad_xnum == int(len(Vars["bv"])/pad_xnum):  ax1.set_xlabel('Time (ns)', fontsize=30) # Add x and y titles only to the leftmost and bottom-most subplot
            # if index % pad_xnum == 0:  ax1.set_ylabel('Voltage (V)', fontsize=30)
            # ax1.set_xlim(xmin,xmax)
            # ax1.legend()
            index+=1
        # fig_multi.tight_layout()
        # fig_multi.savefig(f"{plotDir}/spectrums/multi_{int(n_pho)}.png", format='png')
        # print(f"{plotDir}/spectrums/multi_{int(n_pho)}.png")
        # # ax.grid(True)
        ax.set_xlabel('Time (ns)',fontsize=15)
        ax.set_ylabel('Voltage (V)',fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_xlim(xmin,xmax)
        handles, labels = ax.get_legend_handles_labels() # Get the current handles and labels from the axes
        handles = handles[::-1] # Reverse the order of handles and labels
        labels = labels[::-1]
        Ibias=r'$I_{bias}$'
        Isw=r'$I_{sw}$'
        legend = ax.legend(handles,labels,title=f'{Ibias} / {Isw}({SampleT_norm:.2f}Tc)')
        legend.get_title().set_fontsize('large')  # Adjust the font size of the legend title
        P0=r'$P_{0}^{cal}$' if (args.res) else r'$P_{0}^{sc}$'
        # ax.text(0.25, 0.95, f'P = {n_pho}{P0}', transform=ax.transAxes, verticalalignment='top', fontsize=15)
        fig.tight_layout()
        savefig(fig,f"{plotDir}/spectrums/comb_{pho}")
        plt.close("all")

def plots():
    print("\n==================== Start plotting ====================")
    if (args.sweepAll):
        sweepAll_plotDir = args.in_filenames[0].rsplit('/',5)[0]
        createDir(sweepAll_plotDir)
        createDir(f"{sweepAll_plotDir}/multi")
        Normalized_Graph_sweep("pulse_range_ptp",ytit="Signal Amplitude (V)",plotDir=sweepAll_plotDir)
        Normalized_Graph_sweep("pulse_fall_tau",ytit="Fall Time Constant (ns)",plotDir=sweepAll_plotDir)
        Normalized_Graph_sweep_singleIntensity(15555 if args.res else 1555, "pulse_fall_tau", ytit="Fall Time Constant (ns)", plotDir=sweepAll_plotDir)
        Normalized_Graph_sweep("avg",ytit="Signal Amplitude (V)",plotDir=sweepAll_plotDir)
        Normalized_Graph_sweep("resist",ytit=r"Resistance ($Omega$)",plotDir=sweepAll_plotDir)
        Resolution(sweepAll_plotDir)
        multi_histo_matplotlib("h_pulse_fall_range_ptp","Signal Amplitude (V)",sweepAll_plotDir)
        spectrum_sweep_bias(sweepAll_plotDir)
    if (args.sweepPolarization):
        sweepPolar_plotDir = args.in_filenames[0].split('degrees/')[0].split('/')[0]
        Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,effs,title="g_eff",ytit="Pulse Detection Efficiency",ymin=0,ymax=1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization_err(Photons[0],BVs[0],Temps[0],Polars,pulse_ranges,pulse_range_errs, title="g_pulse_range",ytit="Pulse range mean (V)",ymin=0,ymax=max(pulse_ranges.values())*1.6,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization_err(Photons[0],BVs[0],Temps[0],Polars,pre_ranges,pre_range_errs,title="g_pre_range",ytit="Pre range mean (V)",ymin=min(pre_ranges.values())*0.8,ymax=max(pre_ranges.values())*1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,pulse_integrals,title="g_integral",ytit="Signal integrals",ymin=min(pulse_integrals.values())*0.8,ymax=max(pulse_integrals.values())*1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,pulse_fall_taus,title="g_fall_tau",ytit="Fall time constant (0.4ns)",ymin=min(pulse_fall_taus.values())*0.8,ymax=max(pulse_fall_taus.values())*1.2,plotDir=sweepPolar_plotDir)
        # Graph_sweep_polarization(Photons[0],Temps[0],BVs,BCs,Polars,pulse_FWHMs,title="g_FWHM",ytit="FWHM (0.4ns)",ymin=min(pulse_FWHMs.values())*0.8,ymax=max(pulse_FWHMs.values())*1.2,plotDir=sweepPolar_plotDir)
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
        spectrum_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,spectrum,sweepBias_plotDir)
        # multi_spectrum_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,spectrum,sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pulse_fall_ranges,"range",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pre_ranges,"pre_range",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pulse_fall_taus,"tau",sweepBias_plotDir)
        multi_histo_canvas_sweep_bias(Photons[0],Polars[0],Temps[0],BVs,BCs,h_pulse_FWHMs,"FWHM",sweepBias_plotDir)

if __name__ == "__main__":
    # Suppress the specific warning
    warnings.filterwarnings("ignore", category=UserWarning, message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")

    # ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
    # ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
    Vc=130 if args.res else 1060
    Ic=13.6e-6 if args.res else 106e-6
    SampleTc = 12
    SampleT=11.5 if args.res else 4.7
    SampleT_norm = SampleT/SampleTc

    Vars, Stats, subsets = {}, {}, {}
    Vars["pho"], Vars["bv"], Vars["bc"], Vars["n_pho"], Vars["n_bv"], Vars["n_bc"], Vars["pol"], Vars["temp"] = [],[],[],[],[],[],[],[] # List for sweep variables

    # All bias: [74, 84, 94, 104, 114, 124, 134, 144, 154, 164, 174], [500, 600, 700, 800, 850, 900, 950, 1000, 1010, 1020, 1030, 1040, 1050]
    subsets["photon"] = [699,1011,1399,2333]
    subsets["bv"] = [84, 94, 104, 124, 134, 144, 154, 164, 174] if args.res else [800, 850, 900, 950, 1000, 1010, 1020, 1030, 1050]

    for i, in_filename in enumerate(args.in_filenames):
        print (f"{i}/{len(args.in_filenames)}: {in_filename}")
        runFlag, info = getvars(in_filename,Vars,subsets)
        if (runFlag==False): continue
        key = getkey(info['photon_number'],info['bias_voltage'],info['polarization'],info['sample_temperature'])
        print(key)
        if(args.initPlot): stat=getstat(in_filename,key) # loop over the input files
        else: stat=readstat(in_filename,key)
        Stats[key]=stat
    SortVars(Vars,Vc,Ic)
    plots() # Plot them together
