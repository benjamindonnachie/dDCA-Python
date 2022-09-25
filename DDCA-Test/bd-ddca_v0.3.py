# -*- coding: utf-8 -*-
"""

Implementation of Deterministic Dendritic Cell Algorithm
  Based on original code by Dr Julie Greensmith (27/03/2008)
    Which was subsequently modified by Feng Gu on 11/07/2008

@author: Benjamin Donnachie  <benjamin.donnachie@open.ac.uk>

Licence: TBC


This code take exports from plaso in dynamic format (essentially CSV)

Basic structure:
   Read *.dyn files in from directory.
    Import CSV data into pandas
    Resample per Xmins
    Graph timeline
    Calculate safe and danger values
    Graph safe and danger values
    Run through DDCA to calculate MCAV and Ka
    Graph MCAV and Ka
   Repeat

To do:
    Runs slowly compared with original C code : replace for loops with vectors.
    Implement functions


Change log:
    
    2021/06/30 initial implementation
    2021/07/29 Improved graphing
    2021/07/30 Use pickle files to avoid recalculating

"""

# Import packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from math import log, floor

# Data file directory
#directory = "C:\\Users\\benja\\Documents\\plaso-v20210606-Nov_run"
directory = "/home/benjamindonnachie/Documents/plaso-v20210606-Nov_run/"

# Set up
SAMPLESIZE='12H'

# Extend pickle filename - increment if DDCA code changes
PICKLE='a'

# initialise DDCA
CELLS = 2
MAXmigration = 100
#ANTIGEN = 99999
# calculate lifespan and then allocated
tr_interval = MAXmigration / (CELLS - 1)

for entry in os.scandir(directory):
    if (entry.path.endswith(".dyn") and entry.is_file()):
        try :
            # If results file exists, reload and skip potentially lengthly DDCA calculations
            summary_timeline = pd.read_pickle(entry.path + "-" + SAMPLESIZE + PICKLE + ".pkl")
            print("Loaded results from " + entry.path + "-" + SAMPLESIZE + PICKLE + ".pkl")
        except:
            print("Unable to load results pickle file, loading data from ", entry.path)
        
            # Read csv
            log2timeline = pd.read_csv(entry.path)

            # Convert string datetime to datetime, 
            log2timeline['datetime']= pd.to_datetime(log2timeline['datetime'], errors='coerce')

            # Create a summary of activity in min chunks (Replace 1T with variable)
            summary_timeline = log2timeline.set_index("datetime").resample(SAMPLESIZE).apply('count')
        
#            # Remove any empty summary chunks
#            summary_timeline = summary_timeline.loc[(summary_timeline != 0).any(axis=1)]
        
            # Remove unneeded columns and rename to tidy up
            summary_timeline.rename(columns={'tag': 'count'}, inplace=True)
            summary_timeline.drop(summary_timeline.columns.difference(['count']), axis=1, inplace=True)
        
            # Use number of slices to get ANTIGENs
            ANTIGEN = len(summary_timeline)
            
            print("Using ANTIGEN count of", ANTIGEN)
            
            # Calculate safe and danger
            #  Danger - percentage of how close to half maximum (ie 2 * 100 = 200)                   
            summary_timeline['danger'] = summary_timeline['count'] * 200 / summary_timeline.max()['count']
             # Cap at 0 - 100
            summary_timeline['danger'].loc[(summary_timeline['danger'] > 100)] = 100
            summary_timeline['danger'].loc[(summary_timeline['danger'] < 0)] = 0
        
            #  Safe - simple percentage ratio of how close current value is to previous value
            summary_timeline['safe'] = 100 - (100 * summary_timeline['count'].pct_change().abs())
            # Cap at 0 - 100
            summary_timeline['safe'].loc[(summary_timeline['safe'] > 100)] = 100
            summary_timeline['safe'].loc[(summary_timeline['safe'] < 0)] = 0
        
            # Correct first safe value - calculated as NAN
            summary_timeline['safe'][0] = 0
            
            # Temporary work around for pct_change() returning NaN with repeated zero counts.
            summary_timeline = summary_timeline.fillna(100)
        
       # For comparison with original code, write output to text file
#        outputFile = open (entry.path + ".log", 'w')
#       print ("Writing antigen, safe and danger signals to ", entry.path + ".log")
#        for i in range(len(summary_timeline)) :
#            outputFile.write(str(i) + " signal " + str(summary_timeline['danger'][i]) +
#                             " " + str(summary_timeline['safe'][i]) + "\n")
#            outputFile.write(str(i) + " antigen " + str(i) + "\n")
#        outputFile.close()
        # initialise DDCA
        
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
                # Repeat log_antigen() - should make this a function
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
            
           # dc_stats (Only concerned with results here)
           #for counter in range (CELLS): 
           #    if DDCAcell['incarnations'][counter] > 0 :
           #        iterIncarn = DDCAcell['totIter'][counter] / DDCAcell['incarnations'][counter]
           #    else :
           #        iterIncarn = DDCAcell['totIter'][counter]
           
        # Now result() 
        #  Add mcav and ka fields.
            summary_timeline['mcav'] = np.zeros(len(summary_timeline))
            summary_timeline['ka'] = np.zeros(len(summary_timeline))

            for antigenCounter in range (ANTIGEN) :
                if (agtype['m'][antigenCounter] + agtype['s'][antigenCounter] != 0) :
                    summary_timeline['mcav'][antigenCounter] = agtype['m'][antigenCounter] / (agtype['m'][antigenCounter] + agtype['s'][antigenCounter])
                    summary_timeline['ka'][antigenCounter] = agtype['k'][antigenCounter] / (agtype['m'][antigenCounter] + agtype['s'][antigenCounter])
           
            try :
                # For comparison with original code, write output to text file
                outputFile = open (entry.path + "-" + SAMPLESIZE + PICKLE + ".output.csv", 'w')
                print ("Writing antigen profiles to " + entry.path + "-" + SAMPLESIZE + PICKLE + ".output.csv")
        
                # write complete results to csv
                outputFile.write("antigen,datatime,count,danger,safe,mcav,ka\n")
                
                for i in range(len(summary_timeline)) :
                    outputFile.write(str(i) + "," + 
                                 str(summary_timeline.index[i]) + "," +
                                 str(summary_timeline['count'][i]) + "," +
                                 str(summary_timeline['danger'][i]) + "," +
                                 str(summary_timeline['safe'][i]) + "," +
                                 str(summary_timeline['mcav'][i]) + "," +
                                 str(summary_timeline['ka'][i]) + "\n")
                outputFile.close()
            except:
                print("*** Failed writing to output file ***")
                
            try :
                # Also save as a pickle
                summary_timeline.to_pickle(entry.path + "-" + SAMPLESIZE + PICKLE + ".pkl")
                print("Successful pickle file")
                #print("Did not write pickle")
            except :
                print("*** Pickle save failed ***")
                
        # Graph results
         # Stacked vertical plots
        fig, axs = plt.subplots(3, sharex=True)
        fig.suptitle(entry.path)
        plt.xlabel("Time segment - per " + SAMPLESIZE)
        
        # Plot count
        axs[0].set_yscale('log')
        axs[0].set_ylim(0, 10 ** (1 + floor(log(summary_timeline.max()['count'],10))))
        axs[0].scatter(summary_timeline.index.values, summary_timeline['count'], 
                       marker="*", label='Count')
        axs[0].legend(loc='best')
        axs[0].set_ylabel("Activity count")
        
        # Plot danger / safe
        axs[1].scatter(summary_timeline.index.values, summary_timeline['safe'] + 0.33, 
                       color='g', marker='x', label='Safe')
        axs[1].scatter(summary_timeline.index.values, summary_timeline['danger'] - 0.33, 
                       color='r', marker='+', label='Danger')
        axs[1].legend(loc='best')
        axs[1].set_ylabel("Calculated inputs")
        
        # Plot Ka using MCAV to change colour of plot
        axs[2].scatter(summary_timeline.index.values, summary_timeline['ka'], 
                       c=summary_timeline['mcav'], cmap=plt.cm.seismic, label='ka')
        axs[2].axhline(color="red", linestyle=":")
        axs[2].legend(loc='best')
        axs[2].set_ylabel("Anomaly output")
        
        fig.show()
              
        #break