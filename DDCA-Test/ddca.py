#!/usr/bin/env python3
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

Notes:
    Check do_antigen() and do_signals() are in the correct order.

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

Invoke plaso / log2timeline with:
    
    log2timeline.py --vss-stores all --partitions all --volumes all --yara-rules index.yar evidence.E01 
    
    psort.py -o dynamic --additional_fields yara_match --dynamic-time --status-view window -w outfile_yara.dyn infile.plaso 
    
    Then use outfile_yara.dyn as input to this code.

"""
import numpy as np
import pandas as pd

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
        self.DCs = [ dDCA_cells(i, cells, maxMigration, antigenCount) 
               for i in range (cells) ]


    def doSignals (self, danger, safe):
        csm = danger + safe
        k = danger - (2 * safe)
        
        [ self.updateDC(k, csm, j) for j in range (self.cells) ]
    

    def doAntigen (self, antigen):
        # Increment counter to distribute between cells.
        self.cellIndex += 1
        self.cellIndex %= self.cells
        
        self.DCs[self.cellIndex].antigen[antigen] += 1
        

    def updateDC(self, K, csm, j):
        # consider moving to a function of dDCA_cells - although needs to call
        # logAntigen function, could return lifespan.
        self.DCs[j].lifespan -= csm
        self.DCs[j].k += K
        self.DCs[j].iter += 1
        
        if (self.DCs[j].lifespan <= 0):
            # If lifespan exhausted, pass antigen to master antigen profile
            self.logAntigen(j)
            # reinitialise
            self.DCs[j].reset()
           

    def logAntigen(self, cellid):
        if (np.count_nonzero(self.DCs[cellid].antigen) > 0) :
            for antigenCounter in (self.DCs[cellid].antigen.nonzero()) :
                self.DCs[cellid].TotAntigen += self.DCs[cellid].antigen[antigenCounter]
                #  Note, antigenCounter starts at zero so need to +1 below.
                self.k[antigenCounter] += self.DCs[cellid].k * (1 + self.DCs[cellid].antigen[antigenCounter])
                if (self.DCs[cellid].k > 0) :
                    self.m[antigenCounter] += (1 + self.DCs[cellid].antigen[antigenCounter])
                else:
                    self.s[antigenCounter] += (1 + self.DCs[cellid].antigen[antigenCounter])                                
            self.DCs[cellid].antigen = np.zeros(self.antigen, dtype=np.int32)

          
    def results(self):
        # First flush antigen
        [ self.logAntigen(p) for p in range (self.cells)]
        # Calculate anomaly metrics (mcav and ka)
        self.mcav = self.m / (self.m + self.s)
        self.ka = self.k / (self.m + self.s)


if __name__ == '__main__':
    # if running as a script, create test object
#    mydDCAtest = dDCA (2, 8, 100)
#    mydDCAtest.doSignals(25,75)
#    mydDCAtest.doAntigen(1)
#    mydDCAtest.doSignals(50,50)
#    mydDCAtest.doAntigen(2)
#    mydDCAtest.doSignals(100,0)
#    mydDCAtest.doAntigen(3)
#    mydDCAtest.results()
    
    summary_timeline = pd.read_pickle("/Volumes/EXTERNAL/PhD/Upgrade/202203 Run_Yara/MT_test/20220309T203701-Case1-Webserver.E01.plaso_yara.dyn-1minea2015-2018.pkl")
    
    agCount = len(summary_timeline)
    
    #if (agCount >200000): agCount=200000
    
    mydDCAtest = dDCA(2, agCount, 100)
    #mydDCAtest = dDCA(100, agCount, 1000)
        
    print ("Using an ANTIGEN count of ", agCount)
    
    for i in range (agCount) :
        if (i % 5000 == 0) : print (".", end='')
        mydDCAtest.doAntigen(i)
        mydDCAtest.doSignals(summary_timeline['danger'][i], summary_timeline['safe'][i])
    print("")
    
    mydDCAtest.results()
        
    if (agCount < len(summary_timeline)): summary_timeline = summary_timeline.head(agCount)
    
    summary_timeline['new mcav'] = mydDCAtest.mcav
    summary_timeline['new ka'] = mydDCAtest.ka
    summary_timeline['mcav difference'] = summary_timeline['new mcav'] - summary_timeline['mcav']
    summary_timeline['ka difference'] = summary_timeline['new ka'] - summary_timeline['ka']
    
    print ("Aggregated difference in results is.... ", 
           (summary_timeline['mcav difference'] + summary_timeline['ka difference']).sum())
