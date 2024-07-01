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

Using docker for reproducibility.  Invoke plaso / log2timeline with:
    
    
    docker run -v /mnt/evidence/:/data/ log2timeline/plaso:latest log2timeline --vss_stores all --partitions all --volumes all --hashers md5,sha1,sha256 --parsers win7_slow --yara-rules /data/yara-rules-full_20240602.yar --storage_file /data/evidence.plaso /data/evidence.E01
    docker run -v /mnt/evidence/:/data/ log2timeline/plaso:latest psort --analysis browser_search,chrome_extension,sessionize,unique_domains_visited -o null /data/evidence.plaso
    docker run -v /mnt/evidence/:/data/ log2timeline/plaso:latest psort -o dynamic --fields datetime,MACB,source,sourcetype,type,user,host,short,desc,version,filename,inode,notes,format,extra,yara_match,file_entropy,md5_hash,sha1_hash,sha256_hash,tag,message,message_short -w /data/evidence.plaso.csv /data/evidence.plaso
    
    Then use evidence.plaso.csv to derive antigen, safe and danger signals
      for this code.

"""
import numpy as np


class dDCA_cells:
    def __init__(self, cellid, maxCells, maxMigration, antigenCount):
            # Initialise dDCA cell using variables from original code
            self.id = cellid  # id of this cell
            # evenly distribute migration count across cells
            self.maxLifespan = cellid * maxMigration / (maxCells - 1)
            self.lifespan = self.maxLifespan
            self.k = float(0.0)
            self.TotIter = 0
            self.iter = 0
            self.TotAntigen = 0
            self.incarnations = 0
            self.antigen = np.zeros(antigenCount, dtype=np.int32)

    def reset(self):
            self.lifespan = self.maxLifespan
            self.k = float(0.0)
            self.TotIter += self.iter
            self.iter = 0
            self.TotAntigen = 0 # Total antigen a DC collected per incarnation
            self.incarnations += 1
        

class dDCA:
    def __init__(self, cells, antigenCount, maxMigration):
        # Store set up variables and initialise
        self.cells = cells
        self.maxMigration = maxMigration
        self.antigen = antigenCount
        self.cellIndex = 0  # counter to distribute antigen between cells
        
        # Overall agtype
        self.s = np.zeros(antigenCount)
        self.m = np.zeros(antigenCount)
        self.k = np.zeros(antigenCount)
        self.mcav = np.zeros(antigenCount)
        self.ka = np.zeros(antigenCount)
        
        # Create dendritic cells.
        self.DCs = [ dDCA_cells(thisCell, cells, maxMigration, antigenCount) 
               for thisCell in range (cells) ]



    def doSignals (self, danger, safe):
        # Combined with updateDC for cleaner code
        for eachCell in (self.DCs):
            eachCell.lifespan -= (danger + safe)
            eachCell.k += (danger) - (2 * safe)
            eachCell.iter += 1
            if (eachCell.lifespan <= 0):
                # If lifespan exhausted, pass antigen to master antigen profile
                self.logAntigen(eachCell.id)
                # reinitialise
                eachCell.reset()
    

    def doAntigen (self, antigen):
        # Increment counter to distribute between cells.
        self.cellIndex += 1
        self.cellIndex %= self.cells
        
        self.DCs[self.cellIndex].antigen[antigen] += 1
        


    def logAntigen(self, cellid):
        # First check if anything to log...
        if (np.count_nonzero(self.DCs[cellid].antigen) > 0) :
            self.DCs[cellid].TotAntigen += self.DCs[cellid].antigen.sum()
            self.k += (self.DCs[cellid].antigen) * self.DCs[cellid].k
            if (self.DCs[cellid].k > 0): # gt zero, log as anomaly (mature)
                self.m += self.DCs[cellid].antigen
            else: # otherwise log as safe - semi-mature
                self.s +=  self.DCs[cellid].antigen
        
            self.DCs[cellid].antigen = np.zeros(self.antigen, dtype=np.int32)

          
    def results(self):
        # First flush antigen
        [ self.logAntigen(p) for p in range (self.cells)]
        # Calculate anomaly metrics (mcav and ka)
        self.mcav = self.m / (self.m + self.s)
        self.ka = self.k / (self.m + self.s)


if __name__ == '__main__':

    # Tested against original dDCA (with correction for time interval typed as
    #  a float and then an int) and only difference relates to number of
    #  significant digits used.
    
    import pandas as pd
    import time
    
    # values to initialise DDCA
    CELLS = 100
    MAXmigration = 100
    MAXvalue = 100
    
    x = 'C:\\Users\\benja\\Documents\\Test images\\Initial analysis_Win network_20240601\\20240212-decrypted-Windows_Server_2022.dd.plaso-with_messsage.csv-preprocessed_danger.pkl-preprocessed_summary-60min--100_cells-100_maxmig.pkl'

    print ("Loading pickles...")
    
    summary_timeline = pd.read_pickle(x)
    
    STARTdate = '2023-01-01 00:00:00'
    ENDdate = '2024-06-30 23:59:59'
    
    summary_timeline_test = summary_timeline [(summary_timeline.index >= STARTdate) & (summary_timeline.index <= ENDdate)]
    
    ANTIGEN = len(summary_timeline_test)
    print("Using ANTIGEN count of", ANTIGEN)
    testdDCA = dDCA(CELLS, ANTIGEN, MAXmigration)
    
    start = time.time()
    
    for curAntigen in range (ANTIGEN) :
        if (curAntigen % 5000 == 0) : print ("*", end='')
    
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
    