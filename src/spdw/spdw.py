#! /usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def render_output(args, name):
    if args.output_prefix:
        # ax.get_figure().savefig(
        # ax.savefig(
        # ax.figure.savefig(
        plt.savefig(
                compose_output_path(args.output_prefix, name, args.output_format),
                bbox_inches="tight")
    if args.show_plot_on_screen is not False:
        plt.show()

