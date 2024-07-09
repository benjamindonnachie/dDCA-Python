#!/usr/bin/env python3
"""

Implementation of Deterministic Dendritic Cell Algorithm
  Based on original code by Dr Julie Greensmith (27/03/2008)
    Which was subsequently modified by Feng Gu on 11/07/2008

@author: Benjamin Donnachie  <benjamin.donnachie@open.ac.uk>

Licence: CC-BY-SA-NC

Change log:
    
    2021/06/30 initial implementation
    2021/07/29 Improved graphing
    2021/07/30 Use pickle files to avoid recalculating
    2022/03/10 Added support for yara to calculate danger signal 
      and "Not set" or "error" times (needs additional psort options)
    2022/03/24 Improved loop handling for do_signals
    2022/04/10 Reverse antigen and signal.  Introduce date ranges for graph
    2022/04/11 Change danger signal to use logarithmic ratios.  Use date ranges for calculations.
    2022/09/25 Implement object based code
    2024/06/30 Replaced for loops with vectors and general tidy
    2024/07/06 Tuning.  Implemented 'single antigen' variant
    2024/07/08 Replace lists with np arrays.  Restored previous multi-antigen 
                functionality using lists for local cell antigen profile. 
                Removed 'legacy' historicl variables used for metrics.
                Significant performance increase - from 9.66 sec on test dataset 
                of 13,128 entries down to 0.64 sec (Apple Silicon M3)
                
Using docker for reproducibility.  Invoke plaso / log2timeline with:
    
    
    docker run -v /mnt/evidence/:/data/ log2timeline/plaso:latest log2timeline --vss_stores all --partitions all --volumes all --hashers md5,sha1,sha256 --parsers win7_slow --yara-rules /data/yara-rules-full_20240602.yar --storage_file /data/evidence.plaso /data/evidence.E01
    docker run -v /mnt/evidence/:/data/ log2timeline/plaso:latest psort --analysis browser_search,chrome_extension,sessionize,unique_domains_visited -o null /data/evidence.plaso
    docker run -v /mnt/evidence/:/data/ log2timeline/plaso:latest psort -o dynamic --fields datetime,MACB,source,sourcetype,type,user,host,short,desc,version,filename,inode,notes,format,extra,yara_match,file_entropy,md5_hash,sha1_hash,sha256_hash,tag,message,message_short -w /data/evidence.plaso.csv /data/evidence.plaso
    
    Then use evidence.plaso.csv to derive antigen, safe and danger signals
      for this code.

"""
import numpy as np
import itertools

class dDCA_cells:
    def __init__(self, cellid, TotCells, maxMigration, antigenCount):
            # Initialise dDCA cell using variables from original code
            self.maxLifespan = cellid * maxMigration / (TotCells - 1) # evenly distribute migration count across cells
            self.lifespan = self.maxLifespan
            self.k = float(0.0)
            self.antigen = []

    def reset(self):
            self.lifespan = self.maxLifespan
            self.k = float(0.0)
        

class dDCA:
    def __init__(self, cells, antigenCount, maxMigration):
        # Store set up variables and initialise
        self.cells = cells
        self.maxMigration = maxMigration
        self.antigen = antigenCount
                
        # Overall master agtype antigen profile
        self.s = np.zeros(self.antigen, 'i')
        self.m = np.zeros(self.antigen, 'i')
        self.k = np.zeros(self.antigen, 'f')

        # Create dendritic cells.
        self.DCs = [ dDCA_cells(eachCell, cells, maxMigration, antigenCount) 
               for eachCell in range (cells) ]
        self.activeCell = itertools.cycle(self.DCs)
        # Start with first cell
        next(self.activeCell)

    def doSignals (self, danger, safe):
        csm = danger + safe
        k = danger - (2 * safe)
        
        for eachCell in self.DCs:  
            eachCell.lifespan -= csm
            eachCell.k += k
            if (eachCell.lifespan <= 0): # If lifespan exhausted, pass antigen to master antigen profile
                if (eachCell.antigen): self.logAntigen(eachCell)
                eachCell.reset() # reinitialise

    def doAntigen (self, antigen):    
        # Adcance cell and record antigen profile
        next(self.activeCell).antigen.append(antigen)        

    def logAntigen(self, thisCell):
        for eachAntigen in thisCell.antigen:
            self.k[eachAntigen] += thisCell.k
            if (thisCell.k > 0): # gt zero, log as anomaly (mature)
                self.m[eachAntigen] += 1
            else: # otherwise log as safe - semi-mature
                self.s[eachAntigen] += 1
        thisCell.antigen = []

    def results(self):
        # First flush antigen
        [ self.logAntigen(eachCell) for eachCell in self.DCs]
        
        # Calculate anomaly metrics (mcav and ka)
        self.mcav = self.m / (self.m + self.s)
        self.ka = self.k / (self.m + self.s)


if __name__ == '__main__':

    # Tested against original dDCA (with correction for time interval typed as
    #  a float and then an int) and only difference relates to number of
    #  significant digits used.
    
    import pandas as pd
    import time
    import sys
    
    # values to initialise DDCA
    CELLS = 100
    MAXmigration = 100
    MAXvalue = 100
    
    x = '20240212-decrypted-Windows_Server_2022.dd.plaso-with_messsage.csv-preprocessed_danger.pkl-preprocessed_summary-60min--100_cells-100_maxmig.pkl'

    print ("Loading pickles...")
    
    spinner = itertools.cycle(['-', '/', '|', '\\'])
    
    summary_timeline = pd.read_pickle(x)
    
    STARTdate = '2023-01-01 00:00:00'
    ENDdate = '2024-06-30 23:59:59'
    
    summary_timeline_test = summary_timeline [(summary_timeline.index >= STARTdate) & (summary_timeline.index <= ENDdate)]
    
    ANTIGEN = len(summary_timeline_test)
    print("Using ANTIGEN count of", ANTIGEN)
    testdDCA = dDCA(CELLS, ANTIGEN, MAXmigration)
    
    start = time.time()
    
    for curAntigen in range (ANTIGEN) :
        if (curAntigen % 5000 == 0) : 
            sys.stdout.write(next(spinner))   # write the next character
            sys.stdout.flush()                # flush stdout buffer (actual character display)
            sys.stdout.write('\b')            # erase the last written char
            
        # As per paper, antigen first
        # Gu, Feng, Julie Greensmith, and Uwe Aickelin. ‘Integrating Real-Time Analysis with the Dendritic Cell Algorithm through Segmentation’. In Proceedings of the 11th Annual Conference on Genetic and Evolutionary Computation, 1203–10. Montreal Québec Canada: ACM, 2009. https://doi.org/10.1145/1569901.1570063.
    
        testdDCA.doAntigen(curAntigen)
        
        safe = summary_timeline_test.iloc[curAntigen]['safe']
        danger = summary_timeline_test.iloc[curAntigen]['danger_max']
    
        testdDCA.doSignals(danger, safe)
     
    testdDCA.results()
    
    end = time.time()
    
    summary_timeline_test['mcav_newdDCA'] = testdDCA.mcav
    summary_timeline_test['ka_newdDCA'] = testdDCA.ka
        
    print("Time taken: ", end - start)