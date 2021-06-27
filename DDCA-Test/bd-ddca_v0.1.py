# -*- coding: utf-8 -*-
"""
Created on Wed May 12 23:16:54 2021

Implementation of Deterministic Dendritic Cell Algorithm
  Based on original code by Dr Julie Greensmith (27/03/2008)
  Subsequently modified by Feng Gu on 11/07/2008

@author: Benjamin Donnachie  <benjamin.donnachie@open.ac.uk>

Change log:
    
    2021/06/27 initial implementation

"""

# Import packages
import os
import pandas
import matplotlib.pyplot

# Later - read files in from directory.
#  stages:
#    Import CSV data into pandas
#    Resample per 5mins
#    Graph timeline
#    Calculate safe and danger values
#    Graph safe and danger values
#    Run through DDCA to calculate MCAV and Ka
#    Graph MCAV and Ka
#   Repeat

# Import data from .csv file
#fname = os.path.join("C:\\Users\\benja\\Documents\\plaso",
#                     "AH_Case1-Webserver.E01.plaso.dyn")

fname = os.path.join("C:\\Users\\benja\\Documents\\plaso",
                     "ENISA_coloserver1337.myhosting.ex-disk1.plaso.dyn")

# Read csv
log2timeline = pandas.read_csv(fname)

# Convert string datetime to datetime, 
log2timeline['datetime']= pandas.to_datetime(log2timeline['datetime'], errors='coerce')

# log2timeline.info()

# Create a summary of activity in 5min chunks
summary = log2timeline.set_index("datetime").resample('5T').apply('count')

# Remove any empty summary chunks
summary = summary.loc[(summary != 0).any(axis=1)]

# Plot summary

matplotlib.pyplot.yscale("Log")
matplotlib.pyplot.scatter(summary.index.values, summary["tag"])
matplotlib.pyplot.show()
