#!/usr/bin/env python
# encoding: utf-8
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    This module configures the look and feel of pconpy via matplotlib
    configurations.

Usage:
    mplutils.py

"""

import os
import matplotlib
import matplotlib as mpl
mpl.use('Agg')
import pylab
import numpy as np


PWD = os.path.abspath(os.path.dirname(__file__))

def px2pt(p):
    """Convert pixels to points.

    """
    return p * 72. / 96.


def init_spines(hidden=[]):
    """Initialise the plot frame, hiding the selected spines.

    Arguments:
        hidden -- list of spines to hide (default=[]). For example, set hidden
        to ["top", "right"] to hide both the top and right axes borders from
        the plot.

    Returns:
        None

    """

    ax = pylab.gca()

    all_spines = ["top", "bottom", "right", "left", "polar"]

    for spine in all_spines:
        if spine in hidden:
            ax.spines[spine].set_visible(False)
        else:
            try:
                ax.spines[spine].set_visible(True)
                ax.spines[spine].set_linewidth(px2pt(0.75))
            except KeyError:
                pass

    return

def init_pylab(family=None):
    """Initialise and clean up the look and feel of the plotting area.

    """

    mpl.rc("lines", linewidth=px2pt(1))
    mpl.rc("xtick", **{"direction" : "out" })
    mpl.rc("ytick", **{"direction" : "out" })
    mpl.rc("legend", frameon=False, fontsize=10., numpoints=1)

    if family is not None:
        mpl.rc("font", **{ "family" : family })

    pylab.tick_params(axis="x", which="both", top="off")
    pylab.tick_params(axis="y", which="both", right="off")

    init_spines()

    return


if __name__ == '__main__':
    pass
