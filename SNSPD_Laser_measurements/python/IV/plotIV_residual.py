#!/usr/bin/env python3

import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd
import numpy as np

from ..utils.osUtils import createDir
from ..utils.plotUtils import savefig,markers

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
args = parser.parse_args()

def sort_temp(temp,var):
    temp_array = np.array(temp)
    var_array = np.array(var)
    sorted_temp_array = np.sort(temp_array)
    sorted_var_array= var_array[temp_array.argsort()]
    return sorted_temp_array, sorted_var_array

def read_file(in_filename,subset):
    print(f"Processing {in_filename}")
    Vars={}
    # Info
    temp_str = in_filename.split('/')[-1].split('K')[0].split('_')[-1]
    Vars["temp"] = float(temp_str.replace('p', '.'))
    Vars["Sample"] = in_filename.split('SNSPD_rawdata/')[1].split('/')[0]
    Vars["Ch"] = in_filename.split('SNSPD_rawdata/')[1].split('/')[2]
    Vars["Source"] = in_filename.split('IV_')[1].split('_Source/')[0]
    Vars["basename"] = in_filename.rsplit('/')[-1].split('.txt')[0]
    # Read the data from the file
    df = pd.read_csv(in_filename, header=None, delim_whitespace=True, names=['Volts', 'Currents', 'Resists'])
    split_index = len(df) // 4
    df_1 = df.iloc[:split_index].copy()
    df_1["Currents"] *= 1e6
    df_1['Volts'] *= 1e3
    df_1['Resists'] /= 1e3
    subset_df = df_1[df_1['Currents'] < subset]
    return subset_df, Vars

def bare():
    Bare_file="/wk_cms3/wuhsinyeh/SNSPD/SNSPD_rawdata/Bare/IV_Current_Source/Ch0/Bare_Ch0_300K_20240524_1310.txt"
    df_bare, Vars_bare = read_file(Bare_file, 1000)
    cubic_spline = CubicSpline(df_bare['Currents'], df_bare['Resists'])
    return cubic_spline

def residual(df,cubic_spline):
    fitted_y_data = cubic_spline(df['Currents'])
    residuals = df['Resists'] - fitted_y_data
    return [df['Currents'], residuals]

def find_ic(bc,resist):
    slopes = np.diff(resist) / np.diff(bc) # Calculate the slope between consecutive points
    slope_changes = np.abs(np.diff(slopes))
    # print(slope_changes)
    # turning_point_indices = np.where(slope_changes >= 0.5)[0] + 1
    # turning_points = [(bc[i], resist[i]) for i in turning_point_indices] # Coordinates of the turning points
    # print(f'Turning Point: Bias Current = {turning_points[0]} µA, Resistance = {turning_points[1]} kΩ')

def plot(temps,graphs,xtitle,ytitle,plotDir,savetitle,xmin=-1,xmax=-1,ymin=-1,ymax=-1):
    fig,ax = plt.subplots()
    for i,(graph,temp) in enumerate(zip(graphs,temps)):
        print(temp)
        ax.plot(graph[0], graph[1], marker=markers(i), fillstyle='none', linestyle='none', label=f'{temp/SampleTc:.2f}')
        find_ic(graph[0], graph[1])
    ax.set_ylabel(ytitle,fontsize=15)
    ax.set_xlabel(xtitle,fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(title=r'$T$ / $T_{C}$')
    fig.tight_layout()
    savefig(fig,f"{plotDir}/{savetitle}")
    plt.show()
    ax.set_xlim(13,14.5)
    ax.set_ylim(-0.15,0.3)
    savefig(fig,f"{plotDir}/{savetitle}_zoomin")

if __name__ == '__main__':
    SampleTc = 12.5
    resists, temps=[], []
    bare_spline = bare()
    for i, in_filename in enumerate(args.in_filenames):
        df, Vars = read_file(in_filename, 1000)
        resists.append(residual(df,bare_spline))
        temps.append(Vars["temp"])

    plotDir = 'plots/' + Vars["Sample"] + '/IV/' + Vars["Ch"]
    createDir(plotDir)
    # sort_temp(temps,resists)
    plot(temps, resists, r"Bias Current ($\mu$A)", r"Resistance (k$\Omega$)", plotDir, "IR_residual")
