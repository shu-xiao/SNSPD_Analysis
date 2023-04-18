#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import config

def event_display(np,title='Waveform'):
    # Create a line plot of the data
    plt.plot(range(len(np)), np)
    # Add labels to the plot
    plt.title(title)
    plt.xlabel('Index')
    plt.ylabel('Value')
    # Display the plot
    if (config.DISPLAY): plt.show()

def event_display_2ch(np1, np2, title):
    # Create a new figure
    fig, ax = plt.subplots()

    # Plot the three arrays
    ax.plot(range(len(np2)), np2-0.15, label='ch2', marker='o',fillstyle='none')
    ax.plot(range(len(np1)), np1, label='ch1', marker='o',fillstyle='none')
    # Add labels to the plot
    ax.legend()
    ax.set_title(title)
    ax.set_xlabel('Index')
    ax.set_ylabel('Value')
    # ax.set_ylim(-0.25,0.15)
    # ax.set_xlim(760,780)
    plt.grid()
    # Display the plot
    if (config.DISPLAY): plt.show()

def display_spline_fit(spline_func, x_index):
    x_spline_range = np.linspace(x_index.min(), x_index.max(), num=10000)
    y_spline = spline_func(x_spline_range)
    plt.plot(x_spline_range, y_spline, '-', label='Spline Fit')
    if (config.DISPLAY): plt.show()


def plot_histo(np1, nbin, rangemin, rangemax, xTitle, title, saveTitle):
    fig, ax = plt.subplots()
    ax.hist(np1, bins=nbin, range=(rangemin, rangemax), alpha=0.5, color='blue', edgecolor='black')
    mean, std = np.mean(np1), np.std(np1)
    textstr = f'$\mu={mean:.4f}$\n$\sigma={std:.4f}$'
    ax.text(0.73, 0.93, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top')
    ax.set_xlabel(xTitle)
    ax.set_ylabel('Events')
    ax.set_title(title)
    plt.savefig(saveTitle)
    if (config.DISPLAY): plt.show()
    return mean, std

def plot_2histo(np1, np2, nbin, rangemin, rangemax, label1, label2, title, saveTitle='test.png'):
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(8,5))
    n1, bins1, patches1 = ax1.hist(np1, bins=nbin, range=(rangemin, rangemax), alpha=0.6, color='red', label=label1)
    n2, bins2, patches2 = ax2.hist(np2, bins=nbin, range=(rangemin, rangemax), alpha=0.6, color='blue', label=label2)
    ax1.legend()
    ax2.legend()
    # Add mean and standard deviation to the plot
    mean1, std1 = np.mean(np1), np.std(np1)
    mean2, std2 = np.mean(np2), np.std(np2)
    textstr1 = f'$\mu_1={mean1:.4f}$\n$\sigma_1={std1:.4f}$'
    textstr2 = f'$\mu_2={mean2:.4f}$\n$\sigma_2={std2:.4f}$'
    ax1.text(0.8, 0.75, textstr1, transform=ax1.transAxes, fontsize=14, verticalalignment='top')
    ax2.text(0.8, 0.75, textstr2, transform=ax2.transAxes, fontsize=14, verticalalignment='top')
    fig.suptitle(title)
    plt.savefig(saveTitle)
    if (config.DISPLAY): plt.show()
    array = [mean1, std1, mean2, std2]
    return array
