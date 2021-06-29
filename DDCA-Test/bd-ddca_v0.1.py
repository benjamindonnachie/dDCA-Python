# -*- coding: utf-8 -*-
"""
Created on Wed May 12 23:16:54 2021

Implementation of Deterministic Dendritic Cell Algorithm
  Based on original code by Dr Julie Greensmith (27/03/2008)
  Subsequently modified by Feng Gu on 11/07/2008

This code take exports from plaso in dynamic format (essentially CSV)

@author: Benjamin Donnachie  <benjamin.donnachie@open.ac.uk>

Change log:
    
    2021/06/27 initial implementation

"""

# Import packages
import os
import pandas as pd
import numpy as np
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

def percent_cap(a):
    if (a > 100) :
        return 100
    if (a < 0) :
        return 0
    return a

directory = "C:\\Users\\benja\\Documents\\plaso"
for entry in os.scandir(directory):
    if (entry.path.endswith(".dyn") and entry.is_file()):
        print("Loading data from ", entry.path)
        
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
        
        # Correct first safe value - calculated as NAN
        summary_timeline['safe'][0] = 0


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
        CELLS = 2
        MAXmigration = 100
        ANTIGEN = 99999
        # calculate lifespan and then allocated
        tr_interval = MAXmigration / (CELLS - 1)
        # Comments from original code included.
        DDCAcell = {'lifespan': np.linspace(0, MAXmigration, CELLS), # migration threshold countdown
                                  'k': np.zeros(CELLS), # anomaly output variable
                                  'antigen' : np.zeros((CELLS, ANTIGEN), dtype=np.int32),
                                  'iter' : np.zeros(CELLS, dtype=np.int32), # iterations of signal updates received
                                  'incarnations' : np.zeros(CELLS, dtype=np.int32),
                                  'id' : range(CELLS),
                                  'totIter' : np.zeros(CELLS, dtype=np.int32),
                                  'totAg' : np.zeros(CELLS, dtype=np.int32) # Total amount of antigen collected per incarnation
                                  }
        
        agtype = { 's' : np.zeros(ANTIGEN), # defaults to float
                  'm' : np.zeros(ANTIGEN),
                  'k' : np.zeros(ANTIGEN)
                  }
        
        cellIndex = 0;
        
        # Iterate over summary timeline - use index as antigen value to allow linking back
        for i in range(len(summary_timeline)) :
        #for i in range(10) :
           #print(i, 
           #      summary_timeline['danger'][i], 
           #      summary_timeline['safe'][i])
           
           # do_signals()
           
           # csm = safe + danger
           csm = summary_timeline['danger'][i] + summary_timeline['safe'][i]
           # k = (mature - semi) = (D - S) - S = danger - 2(safe)
           k = summary_timeline['danger'][i] - (2 * summary_timeline['safe'][i])
           
           #print ("slice (antigen) ", i, "csm is ", csm, " k is ", k)
           
           # For each DDCAcell, signal update.
           for counter in range (CELLS):
               # update_DC()
               DDCAcell['lifespan'][counter] -= csm
               DDCAcell['k'][counter] += k
               DDCAcell['iter'][counter] += 1
               #print ("Updated lifespan ",DDCAcell['lifespan'][counter]," k ",DDCAcell['k'][counter], 
               #       " iter ", DDCAcell['iter'][counter])
               if (DDCAcell['lifespan'][counter] <= 0) :
                   #print ("*** Lifespan below zero, logging antigen ***")
                   # log antigen()
                   for antigenCounter in range (ANTIGEN) :
                       if (DDCAcell['antigen'][counter][antigenCounter] > 0) :
                           DDCAcell['totAg'][counter] += DDCAcell['antigen'][counter][antigenCounter]
                           for yyCounter in range (DDCAcell['antigen'][counter][antigenCounter]) :
                               #print ("yyCounter entered - cell ", counter, " antigen ", antigenCounter)
                               agtype['k'][antigenCounter] += DDCAcell['k'][counter]
                               if (DDCAcell['k'][counter] > 0) :
                                   agtype['m'][antigenCounter] += 1
                               else:
                                   agtype['s'][antigenCounter] += 1
                                   
                           DDCAcell['antigen'][counter][antigenCounter] = 0
                   
                   # update_DC() resumes
                   DDCAcell['lifespan'][counter] = counter * tr_interval
                   DDCAcell['k'][counter] = 0
                   DDCAcell['totIter'][counter] += DDCAcell['iter'][counter]
                   DDCAcell['iter'][counter] = 0
                   DDCAcell['totAg'][counter] = 0
                   DDCAcell['incarnations'][counter] += 1
               
           # Antigen update - do_antigen()
           # Distribute between cells
           cellIndex += 1
           cellIndex %= CELLS
           #
           DDCAcell['antigen'][cellIndex][i] += 1
          
        # Flush cells
        for counter in range (CELLS):
            # Repeat log_antigen() - make function
            for antigenCounter in range (ANTIGEN) :
                if (DDCAcell['antigen'][counter][antigenCounter] > 0) :
                    DDCAcell['totAg'][counter] += DDCAcell['antigen'][counter][antigenCounter]
                    for yyCounter in range (DDCAcell['antigen'][counter][antigenCounter]) :
                        #print ("yyCounter2 entered - cell ", counter, " antigen ", antigenCounter)
                        agtype['k'][antigenCounter] += DDCAcell['k'][counter]
                        if (DDCAcell['k'][counter] > 0) :
                            agtype['m'][antigenCounter] += 1
                        else:
                            agtype['s'][antigenCounter] += 1
                                   
                    DDCAcell['antigen'][counter][antigenCounter] = 0
            
           # dc_stats
           #for counter in range (CELLS): 
           #    if DDCAcell['incarnations'][counter] > 0 :
           #        iterIncarn = DDCAcell['totIter'][counter] / DDCAcell['incarnations'][counter]
           #    else :
           #        iterIncarn = DDCAcell['totIter'][counter]
           
           
           
        # Now result() 
        summary_timeline['mcav'] = np.zeros(len(summary_timeline))
        summary_timeline['ka'] = np.zeros(len(summary_timeline))
           
        for antigenCounter in range (ANTIGEN) :
            if (agtype['m'][antigenCounter] + agtype['s'][antigenCounter] != 0) :
                #print ("Updating MCAV and k at ", antigenCounter)
                summary_timeline['mcav'][antigenCounter] = agtype['m'][antigenCounter] / (agtype['m'][antigenCounter] + agtype['s'][antigenCounter])
                summary_timeline['ka'][antigenCounter] = agtype['k'][antigenCounter] / (agtype['m'][antigenCounter] + agtype['s'][antigenCounter])
                   
        plt.title("MCAV " + entry.path)
        plt.yscale("Linear")
        plt.scatter(summary_timeline.index.values, summary_timeline['mcav'])
        plt.show()

        plt.title("Ka " + entry.path)
        plt.yscale("Linear")
        plt.scatter(summary_timeline.index.values, summary_timeline['ka'])
        plt.show()
                    