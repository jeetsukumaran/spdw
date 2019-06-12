#! /usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
##
##  Copyright 2019 Jeet Sukumaran.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Species Delimitation Workshop: Plot inverse-gamma distribution.
"""

import sys
import os
import random
import argparse
import dendropy

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import invgamma

import spdw

def fit_exponential(df):
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    x = df["coalescent_event_idx"]
    y = df["waiting_time"]
    yn = y + 0.2*np.random.normal(size=len(x))

    popt, pcov = curve_fit(func, x, yn)
    plt.figure()
    plt.plot(x, yn, 'ko', label="Original Noised Data")
    plt.plot(x, func(x, *popt), 'r-', label="Fitted Curve")
    plt.legend()
    plt.show()

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran'
__copyright__ = 'Copyright (C) 2019 Jeet Sukumaran.'

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("alpha",
            type=float,
            default=3.00,
            help="Alpha parameter value [default=%(default)s].")
    parser.add_argument("beta",
            type=float,
            default=0.001,
            help="Beta parameter value [default=%(default)s].")
    parser.add_argument("-f", "--input-format",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="Input data format (default='%(default)s')")

    args = parser.parse_args()
    args.output_prefix = None
    args.show_plot_on_screen = True
    a = args.alpha
    b = args.beta
    rv = invgamma(a, loc=0, scale=b)
    print("Mean: {}".format(rv.mean()))
    print("Variance: {}".format(rv.var()))
    fig, ax = plt.subplots(1, 1)
    x = np.linspace(
            invgamma.ppf(0.01, a, scale=b),
            invgamma.ppf(0.99, a, scale=b), 100)
    ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
    ax.plot(x, invgamma.pdf(x, a, scale=b), 'r-', lw=5, alpha=0.6, label='invgamma pdf')
    spdw.render_output(args, "InverseGamma")
    # print("Mean: {}".format(b / (a-1)))
    # print("Variance: {}".format( pow(b,2) / (pow(a-1,2) * (a-2))))

if __name__ == '__main__':
    main()



