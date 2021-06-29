# -*- coding: utf-8 -*-
"""
Created on Wed May 12 23:16:54 2021

Implementation of Deterministic Dendritic Cell Algorithm
  Based on original code by Dr Julie Greensmith (27/03/2008)
  Subsequently modified by Feng Gu on 11/07/2008

This code take exports from plaso in dynamic 

@author: Benjamin Donnachie  <benjamin.donnachie@open.ac.uk>

Change log:
    
    2021/06/27 initial implementation

"""

# Import packages
import os
import pandas as pd
import matplotlib.pyplot as plt

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

directory = "C:\\Users\\benja\\Documents\\GitHub\\Machine-learning-research\\Data sets"
for entry in os.scandir(directory):
    if (entry.path.endswith(".dyn") and entry.is_file()):
        print(entry.path)
        # Read csv
        log2timeline = pd.read_csv(entry.path)

        # Convert string datetime to datetime, 
        log2timeline['datetime']= pd.to_datetime(log2timeline['datetime'], errors='coerce')

        # Create a summary of activity in min chunks (Replace 1T with variable)
        summary_timeline = log2timeline.set_index("datetime").resample('1T').apply('count')

         # Remove any empty summary chunks
        summary_timeline = summary_timeline.loc[(summary_timeline != 0).any(axis=1)]
        
        # Remove unneeded columns and rename to tidy up
        summary_timeline.rename(columns={'tag': 'count'}, inplace=True)
        summary_timeline.drop(summary_timeline.columns.difference(['count']), axis=1, inplace=True)
        
        
        # Calculate safe and danger
        
        #  Danger - percentage of how close to half maximum (ie 2 * 100 = 200)                   
        summary_timeline['danger'] = summary_timeline['count'] * 200 / summary_timeline.max()['count']
        # Cap at 100
        summary_timeline['danger'].loc[(summary_timeline['danger'] > 100)] = 100
        summary_timeline['danger'].loc[(summary_timeline['danger'] < 0)] = 0
        
        #  Safe - simple percentage ratio of how close current value is to previous value
        
        summary_timeline['safe'] = 100 - (100 * summary_timeline['count'].pct_change().abs())
        # Cap at 100
        summary_timeline['safe'].loc[(summary_timeline['safe'] > 100)] = 100
        summary_timeline['safe'].loc[(summary_timeline['safe'] < 0)] = 0

        # Plot summary

        plt.title(entry.path)
        plt.yscale("Log")
        plt.scatter(summary_timeline.index.values, summary_timeline['count'])
        plt.show()
        
        plt.title("Safe / Danger " + entry.path)
        plt.yscale("Linear")
        plt.scatter(summary_timeline.index.values, summary_timeline['danger'])
        plt.scatter(summary_timeline.index.values, summary_timeline['safe'])
        plt.show()


        # initialise DDCA
        #initDDCA()
        
        # Iterate over DataFrame - passing index as antigen together with
        #   safe and danger values
        for i in range(len(summary_timeline)) :
           print(i, summary_timeline['danger'][i], summary_timeline['safe'][i])
        
        #signalantigenDDCA()




