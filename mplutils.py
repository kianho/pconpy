#!/usr/bin/env python
# encoding: utf-8
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    TODO

Usage:
    pplutils.py

Options:

"""

import os
import matplotlib
import matplotlib as mpl
mpl.use('Agg')
import pylab
import numpy as np


PWD = os.path.abspath(os.path.dirname(__file__))

def pt2px(p):
    return p * 96. / 72.

def px2pt(p):
    return p * 72. / 96.

def init_spines(hidden=[]):
    """
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

def init_pylab():
    """Initialise the pylab look and feel like prettyplotlib

    """

    mpl.rc("lines", linewidth=px2pt(1))
    mpl.rc("xtick", **{"direction" : "out" })
    mpl.rc("ytick", **{"direction" : "out" })
    mpl.rc("font", **{"family" : "FreeSans", "size" : 10.})
    mpl.rc("legend", frameon=False, fontsize=10., numpoints=1)

    pylab.tick_params(axis="x", which="both", top="off")
    pylab.tick_params(axis="y", which="both", right="off")

    init_spines()

    return


init_pylab()


if __name__ == '__main__':


    for i in xrange(3):
        x = np.arange(1000)
        y = np.random.normal(size=1000).cumsum()
        pylab.plot(x, y)

    pylab.title("This is a title")
    pylab.xlabel("This is the x-title")
    pylab.ylabel("This is the y-title")
    pylab.semilogy()
    pylab.savefig("./test.pdf")
