#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

for arg in sys.argv[1:]:
    mean = float(arg)
    a = 3
    b = mean * ( a- 1 )
    print("{}:  {}, {}".format(mean, a, b))
