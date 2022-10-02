#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 19:49:25 2022

@author: benjamindonnachie
"""

import ddca
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import log, floor

def main (x):
    # Set up common values [in one place]
    SAMPLESIZE='1min'  # Adjust according to samplesize wanted

    # Extend pickle filename - increment if DDCA code changes
    PICKLE='fa' # Added k and csm values to pickle.
    PICKLE='zz' # force recalc.

    # initialise DDCA
    CELLS = 2
    MAXmigration = 100
    
    MAXvalue = 100 # Original DDCA used 50?
    
    STARTdate = '2015-01-01 00:00:00'
    ENDdate = '2018-12-31 23:59:59'
    SHORTdate = '2015-2018'

    try :
            # If results file exists, reload and skip potentially lengthly DDCA calculations
            summary_timeline = pd.read_pickle(x + "-" + SAMPLESIZE + PICKLE + SHORTdate + ".pkl")
            print("Loaded results from " + x + "-" + SAMPLESIZE + PICKLE + SHORTdate +".pkl")
    except:
            print("Unable to load results pickle file, loading data from " + x)
        
            # Read csv
            log2timeline = pd.read_csv(x)

            # Convert string datetime to datetime, 
            log2timeline['datetime']= pd.to_datetime(log2timeline['datetime'], errors='coerce')

           # count of number of fields in yara_match field before stripping invalid dates
            log2timeline['yara'] = log2timeline['yara_match'].str.split().str.len()            

            # Create a summary of activity in min chunks (Replace 1T with variable)
            summary_timeline = log2timeline.set_index("datetime").resample(SAMPLESIZE).apply('count')
             
            # Remove unneeded columns and rename to tidy up
            summary_timeline.rename(columns={'tag': 'count'}, inplace=True)
            summary_timeline.drop(summary_timeline.columns.difference(['count', 'yara']), axis=1, inplace=True)

            # Restrict timeline
            summary_timeline = summary_timeline[(summary_timeline.index >= STARTdate) & (summary_timeline.index <= ENDdate)]
            
            # Use number of slices to get ANTIGENs
            ANTIGEN = len(summary_timeline)
            
            print("Using ANTIGEN count of", ANTIGEN)
                        
            # Calculate safe and danger
            #   Danger -  count of yara matches calculated above and then resampled.
            
            # Now convert yara count to percentage
            summary_timeline['danger'] = np.log(summary_timeline['yara']) * (MAXvalue / log(summary_timeline['yara'].max()))
            np.nan_to_num(summary_timeline['danger'], copy=False, nan=0.0, posinf=MAXvalue, neginf=0)
       
            ##  Safe - simple percentage ratio of how close current value is to previous value
            # np.diff calculates using out[i] = a[i+1] - a[i]
            #summary_timeline['safe'] = np.diff(summary_timeline['count'], prepend=0)
            # change to out[i] = a[i] - a[i-1]
            summary_timeline['safe'] = np.diff(summary_timeline['count'], prepend=0)
            
            summary_timeline['safe'] = summary_timeline['safe'].abs()
            summary_timeline['safe'] = MAXvalue - (np.log(summary_timeline['safe']) * MAXvalue / log(summary_timeline['safe'].max()))
            
            np.nan_to_num(summary_timeline['safe'], copy=False, nan=100, posinf=MAXvalue, neginf=0)
        
                    
        
            # Use matricies to calculate csm and k across entire dataset
            # csm = safe + danger
            summary_timeline['csm'] = summary_timeline['danger'] + summary_timeline['safe']
            # k = (mature - semi) = (D - S) - S = danger - 2(safe)
            summary_timeline['k'] = summary_timeline['danger'] - (2 * summary_timeline['safe']) 
            
            testdDCA = ddca.dDCA(CELLS, ANTIGEN, MAXmigration)
            
            for i in range (ANTIGEN) :
                if (i % 5000 == 0) : print (".", end='')
                testdDCA.doAntigen(i)
                testdDCA.doSignals(summary_timeline['danger'][i], summary_timeline['safe'][i])
         
            testdDCA.results()
            
            summary_timeline['mcav'] = testdDCA.mcav
            summary_timeline['ka'] = testdDCA.ka
                       
#            try :
#                # Save results as a pickle
#                summary_timeline.to_pickle(x + "-" + SAMPLESIZE + PICKLE + SHORTdate + ".pkl")
#                print("Successful pickle file for " + x)
#            except :
#                print("*** Pickle save failed for " + x + "***")
                
    
    # Either calculated or loaded results from file.
    
    
    # calculate clusters over k (i.e. combined safe and danger)
    
    #print("Calculating k-means...")
    
    #x = summary_timeline['k'].values.reshape(-1,1)
    #kmeans = KMeans(n_clusters=2)
    #kmeans.fit(x)
    #summary_timeline['k_means'] = kmeans.labels_

    #print("Calculating svm...")
    
    #clf = OneClassSVM(gamma='auto').fit(x)
    #summary_timeline['svm'] = clf.predict(summary_timeline['k'])
    
    # Graph results
    # Stacked vertical plots
    fig, axs = plt.subplots(4, sharex=True)
    
    fig.set_size_inches(25.6, 13.43)
    fig.set_dpi(100)
    
#    plt.tight_layout()
    
    fig.suptitle(x)
    plt.xlabel("Time segment - per " + SAMPLESIZE)
        
    # Plot count
    axs[0].set_yscale('log')
    axs[0].set_ylim(0, 10 ** (1 + floor(log(summary_timeline.max()['count'],10))))
    axs[0].scatter(summary_timeline.index.values, summary_timeline['count'], 
                       marker="*", label='Count')
    axs[0].scatter(summary_timeline.index.values, summary_timeline['yara'], 
                       marker='.', label='Yara')
    axs[0].legend(loc='best')
    axs[0].set_ylabel("Activity count")
        
    # Plot danger / safe
    axs[1].scatter(summary_timeline.index.values, summary_timeline['safe'] + 0.33, 
                       color='g', marker='x', label='Safe')
    axs[1].scatter(summary_timeline.index.values, summary_timeline['danger'] - 0.33, 
                       color='r', marker='+', label='Danger')
    axs[1].legend(loc='best')
    axs[1].set_ylabel("Calculated inputs")
        
    # Plot danger / safe
    axs[2].scatter(summary_timeline.index.values, summary_timeline['csm'] + 0.33, 
                       color='orange', marker='x', label='csm')
    axs[2].scatter(summary_timeline.index.values, summary_timeline['k'] - 0.33, 
                       color='b', marker='+', label='k')
    axs[2].legend(loc='best')
    axs[2].set_ylabel("Derived inputs")
        
    # Plot Ka using MCAV to change colour of plot
    axs[3].scatter(summary_timeline.index.values, summary_timeline['ka'], 
                       c=summary_timeline['mcav'], cmap=plt.cm.seismic, label='ka')
    axs[3].axhline(color="red", linestyle=":")
    axs[3].legend(loc='best')
    axs[3].set_ylabel("Anomaly output")
    
    #plt.xlim([datetime.date(2015, 8, 1), datetime.date(2015, 9, 30)])  # If plotting whole file

    #fig.savefig(x + '-' + SAMPLESIZE + PICKLE + SHORTdate + '.png')

    #plt.close(fig)    
    fig.show()   # Processes don't have GUI access
              
    # Testing - stop after one
        
    # break

if __name__ == '__main__' :
    # Data file directory
    #directory = "/home/benjamindonnachie/Documents/plaso-v20210606-Nov_run/"
    #directory = "/home/benjamindonnachie/Documents/plaso-v20220129-March_run/"
    #directory = "/run/media/benjamindonnachie/PhD Data/202203 Run_Yara/MT_Test/"
    directory = "/Volumes/EXTERNAL/PhD/Upgrade/202203 Run_Yara/MT_test/"
 
    file_list = []
    STARTdate = '2015-01-01 00:00:00'
    ENDdate = '2018-12-31 23:59:59'
    SHORTdate = '2015-2018'
    SAMPLESIZE='1min'  # Adjust according to samplesize wanted

    # Extend pickle filename - increment if DDCA code changes
    PICKLE='ea' # Added k and csm values to pickle.

    for entry in os.scandir(directory):
        if (entry.path.endswith(".dyn") and entry.is_file()):
            file_list.append(entry.path)
            main(entry.path) # remove if multi-threaded and use below instead

#    with Pool() as p:
#      p.map(main, file_list)

#    main ("/Volumes/EXTERNAL/PhD/Upgrade/202203 Run_Yara/MT_test/20220309T203701-Case1-Webserver.E01.plaso_yara.dyn")
    #main ("/Volumes/PhD Data/202203 Run_Yara/MT_Test/20220309T203701-Case1-Webserver.E01.plaso_yara.dyn")

#main()

