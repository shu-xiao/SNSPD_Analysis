#!/usr/bin/env python3

from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# User defined functions
from pulseLaserCalib_tdms import SingleTDMS_analysis
from Timing_Analyzer import *
from tdmsUtils import *
from plotUtils import *

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-d',default="./",type=str,help='output directory')
parser.add_argument('--report','-r',default=100,type=int,help='report every x events')
parser.add_argument('--debug_report','-b',default=-1,type=int,help='report every x events')
parser.add_argument('--display_report','-p',default=-1,type=int,help='report every x events')
args = parser.parse_args()

if __name__ == "__main__":

    for in_filename in args.in_filenames:
        sweep_voltage_current = in_filename.rsplit('/',1)[1].split('mV')[0].rsplit('_',1)[1]
        SingleTDMS_analysis(in_filename)
