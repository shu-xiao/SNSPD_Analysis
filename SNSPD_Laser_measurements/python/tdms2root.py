#!/usr/bin/env python

from nptdms import TdmsFile
import pandas as pd
import numpy as np
from array import array
import ROOT
import errno
import os

# User defined functions
from tdmsUtils import *

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-d',default="./output_root",type=str,help='output directory')
parser.add_argument('--report','-r',default=100,type=int,help='report every x events')
args = parser.parse_args()

if __name__ == "__main__":

    # make outputDir
    try:
        os.makedirs(args.outputDir)
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('output directory exists.')
        else:
            raise

    for i, in_filename in enumerate(args.in_filenames):
        print('Processing %s (%d/%d)'% (in_filename, i, len(args.in_filenames)))
        if (in_filename.find('txt')!=-1): baseName = in_filename.rsplit('/',1)[1].split('.txt')[0]
        elif (in_filename.find('tdms')!=-1): baseName = in_filename.rsplit('/',1)[1].split('.tdms')[0]
        with TdmsFile.open(in_filename) as tdms_file:
            # Read Meta Data (Basic information)
            metadata = tdms_file.properties
            metadata_df = pd.DataFrame(metadata.items(), columns=['metaKey', 'metaValue'])
            print(metadata_df)
            totalEvents = int(metadata_df.loc[metadata_df['metaKey'] == 'Total Events', 'metaValue'].iloc[0])
            recordLength = int(metadata_df.loc[metadata_df['metaKey'] == 'record length', 'metaValue'].iloc[0])
            # Read Groups and Channels
            Read_Groups_and_Channels(tdms_file)

            # Write metadata to json file
            metadata_df.to_json(f'{args.outputDir}/{baseName}.json',orient="records",lines=True)
            # Write metadata to txt file
            metadata_df = metadata_df.reset_index()  # make sure indexes pair with number of rows
            with open(f'{args.outputDir}/{baseName}.txt', 'w') as txtFile:
                for index, row in metadata_df.iterrows():
                    txtFile.write(f"{row['metaKey']} : {row['metaValue']}\n")


            ch1 = array( 'f', [ 0 ] * recordLength)
            ch2 = array( 'f', [ 0 ] * recordLength)
            t = ROOT.TTree("SNSPD_data", "SNSPD_data")
            t.Branch('ch1',ch1,"ch1[%d]/F" %(recordLength))
            t.Branch('ch2',ch2,"ch2[%d]/F" %(recordLength))

            for event, chunk in enumerate(tdms_file.data_chunks()):
                # Loop progress
                if ((event+1)%args.report==0): print (f"==========Processing {event}/{totalEvents} event==========")
                # Skip chunk larger than totalEvents
                if (event > totalEvents-1): break
                # Read ch1 into np array
                ch1_chunk = chunk['ADC Readout Channels']['ch1']._data()
                ch2_chunk = chunk['ADC Readout Channels']['ch2']._data()

                for i in range(0, recordLength):
                    ch1[i] = ch1_chunk[i];
                    ch2[i] = ch2_chunk[i];
                t.Fill()

            # Output
            rootFileName = f"{args.outputDir}/{baseName}.root"
            rootFile = ROOT.TFile.Open(rootFileName,"RECREATE")
            rootFile.cd()
            t.Write()
            rootFile.Close()
