#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('in_filenames',nargs="+",help='input filenames')
parser.add_argument('--outputDir','-o',default="./",type=str,help='output directory')
parser.add_argument('--report','-r',default=10000,type=int,help='report every x events')
parser.add_argument('--debug','-d',action="store_true",help='debug mode')
args = parser.parse_args()

# Read the data from the provided text
df = pd.read_csv(args.in_filenames[0], delimiter='\t')

# Extract the first column and get the column names
x = df.iloc[:, 0]
column_names = df.columns[1:]

# Plot the value of each column against the first column
for column in column_names:
    values = df[column]
    plt.plot(x, values, label=column)

# Set the plot labels and title
plt.xlabel(r'Wavelength $\lambda$ [nm]',fontsize=15)
plt.ylabel('Value',fontsize=15)
# if (args.in_filenames[0].find("real")!=-1): plt.ylabel(r'$Re(\varepsilon)$',fontsize=15)
# else: plt.ylabel(r'$Im(\varepsilon)$',fontsize=15)
# plt.title("RF Power = 120W",fontsize=20)
# plt.title(r"Ar:$N_{2}$ = 12:0.2",fontsize=20)
# plt.title('Refractive index',fontsize=20)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
# plt.legend(title=r'      Ar:$N_{2}$')
plt.legend(title='n+ik',fontsize='x-large')
plt.grid(True)
plt.grid(True, which='minor', linestyle=':', linewidth='0.5')  # Minor grid
plt.minorticks_on()  # Enable minor ticks
plt.tight_layout()
plt.savefig(f"{args.outputDir}/{args.in_filenames[0].split('/')[-1].split('.txt')[0]}.png")
plt.show()
