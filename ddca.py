#!/usr/bin/env python3
"""Authoritative Python implementation of the deterministic DCA (dDCA).

This research implementation was written by Benjamin Donnachie. It is based
on a C implementation shared by Dr Julie Greensmith (27 March 2008), which
records later modification by Feng Gu (11 July 2008).

The module is one component of a wider research pipeline. Upstream stages must
derive integer antigen identifiers and floating-point danger/safe signals.
Present each antigen before its corresponding signal pair, then call
``results()`` after the stream has been consumed.

Copyright (c) 2021-2026 Benjamin Donnachie
SPDX-License-Identifier: CC-BY-NC-SA-4.0
"""

import itertools

import numpy as np

__version__ = "9.0.0"


class dDCA_cells:
    """State for one cell in the deterministic DCA population.

    The legacy class and attribute names are retained for pipeline
    compatibility. ``antigenCount`` is retained in the constructor signature
    for compatibility but is not used by an individual cell.
    """

    def __init__(self, cellid, TotCells, maxMigration, antigenCount):
        self.cellid = cellid

        # A single cell cannot use the original evenly distributed expression,
        # because TotCells - 1 would be zero.
        if TotCells == 1:
            self.maxLifespan = maxMigration
        else:
            self.maxLifespan = cellid * maxMigration / (TotCells - 1)

        self.lifespan = self.maxLifespan
        self.k = float(0.0)
        # Dynamic per-cell collections remain lists: this is deliberate and
        # avoids repeatedly growing dense arrays in the antigen hot path.
        self.antigen = []
        self.TotIter = 0
        self.iter = 0
        self.TotAntigen = 0
        self.incarnations = 0

    def reset(self):
        """Reset the cell after migration while retaining lifetime statistics."""
        # Deliberate divergence from the shared C implementation: its reset
        # interval is declared as int, although initialisation uses float. That
        # truncation is believed unintentional, so the original floating-point
        # maximum lifespan is preserved across every cell incarnation here.
        self.lifespan = self.maxLifespan
        self.k = float(0.0)
        self.TotIter += self.iter
        self.TotAntigen = 0
        self.incarnations += 1
        self.iter = 0


class dDCA:
    """Deterministic DCA population and aggregate antigen profile."""

    def __init__(self, cells, antigenCount, maxMigration):
        self.cells = cells
        self.maxMigration = maxMigration
        self.antigen = antigenCount
        self.TotDanger = 0
        self.TotSafe = 0
        self.Is = 0
        self.Tk = 0
        self.Tmcav = 0

        # Fixed-size aggregate profiles remain NumPy arrays. Dynamic per-cell
        # antigen collections are lists (see dDCA_cells).
        self.s = np.zeros(self.antigen, "i")
        self.m = np.zeros(self.antigen, "i")
        self.k = np.zeros(self.antigen, "f")

        self.DCs = [
            dDCA_cells(eachCell, cells, maxMigration, antigenCount)
            for eachCell in range(cells)
        ]
        self.activeCell = itertools.cycle(self.DCs)

        # The shared C implementation increments cell_index before assigning
        # antigen. Pre-advancing the cycle reproduces that behaviour: with more
        # than one cell, the first antigen is assigned to cell 1.
        next(self.activeCell)

    def doSignals(self, danger, safe):
        """Apply one danger/safe signal pair to every cell."""
        csm = danger + safe
        k = danger - (2 * safe)
        self.TotDanger += danger
        self.TotSafe += safe
        self.Is += 1

        for eachCell in self.DCs:
            eachCell.lifespan -= csm
            eachCell.k += k
            eachCell.iter += 1
            if eachCell.lifespan <= 0:
                if eachCell.antigen:
                    self.logAntigen(eachCell)
                eachCell.reset()

    def doAntigen(self, antigen):
        """Assign one integer antigen identifier to the next active cell."""
        next(self.activeCell).antigen.append(antigen)

    def logAntigen(self, thisCell):
        """Merge a migrated cell's antigen into the aggregate profiles."""
        if not thisCell.antigen:
            return

        # The shared C implementation traverses antigen identifiers in numeric
        # order. Sorting the list before accumulation reproduces that order.
        thisCell.antigen.sort()
        for eachAntigen in thisCell.antigen:
            self.k[eachAntigen] += thisCell.k
            thisCell.TotAntigen += 1
            if thisCell.k > 0:
                self.m[eachAntigen] += 1
            else:
                self.s[eachAntigen] += 1
        thisCell.antigen = []

    def results(self):
        """Flush pending antigen and calculate MCAV, Kalpha and thresholds."""
        for eachCell in self.DCs:
            self.logAntigen(eachCell)

        self.mcav = self.m / (self.m + self.s)
        self.ka = self.k / (self.m + self.s)

        totalIterIncarn = 0
        for eachCell in self.DCs:
            if eachCell.incarnations > 0:
                totalIterIncarn += eachCell.TotIter / eachCell.incarnations
            else:
                totalIterIncarn += eachCell.TotIter

        aveIterIncarn = totalIterIncarn / self.cells

        # Kalpha threshold derived from equations 7 and 8 in Greensmith and
        # Aickelin (2008), using danger - 2*safe as the weighted signal sum.
        self.Tk = (
            (self.TotDanger - (2 * self.TotSafe)) / self.Is
        ) * aveIterIncarn

        # MCAV threshold described by Greensmith and Aickelin (2008) as the
        # ratio of total danger signals to total safe signals in the dataset.
        self.Tmcav = self.TotDanger / self.TotSafe

    def dc_stats(self):
        """Return per-cell lifetime statistics as a pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame(
            {
                "cell_id": [eachCell.cellid for eachCell in self.DCs],
                "incarnations": [eachCell.incarnations for eachCell in self.DCs],
                "TotIter": [eachCell.TotIter for eachCell in self.DCs],
            }
        )
