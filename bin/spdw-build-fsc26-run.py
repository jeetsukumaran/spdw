#! /usr/bin/env python
# -*- coding: utf-8 -*-

fsc_parfile_template = """\
//Number of population samples (demes)
{num_samples} samples to simulate :
//Population effective sizes (number of genes)
{population_effective_sizes}
//Samples sizes
{sample_sizes}
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
{num_historical_events} historical event
{historical_events}
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA {num_loci} 0.0000 0.0005 0
"""

historical_events_template="""\
{div_time} {child} {parent} 1 1 0 0"""

h = historical_events_template.format(div_time=1000, child=1, parent=0)
fsc_parfile_str = fsc_parfile_template.format(
        num_samples=2,
        population_effective_sizes="1000\n1000",
        sample_sizes="5\n5",
        num_historical_events=1,
        historical_events=h,
        num_loci=10)
out = open("x.par", "w")
with out:
    out.write(fsc_parfile_str)

