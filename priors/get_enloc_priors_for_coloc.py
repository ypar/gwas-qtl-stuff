#!/usr/bin/env python3

# ypar

import pandas
import numpy
import sys
import os

# example files here
# filename = 'gtex/enloc_v8/enloc/Artery_Tibial_LDL_Cholesterol/output/Artery_Tibial_LDL_Cholesterol.enrich.est'

tissue = sys.argv[1]
trait = sys.argv[2]
tissuetrait = tissue + '_' + trait

# tissue = 'Artery_Tibial'
# trait = 'LDL_Cholesterol'

rootdir = '~/gtex/enloc_v8/'
estimatefile = rootdir + 'enloc/' + tissuetrait + '/output/' + tissuetrait + '.enrich.est'
qtlpipfile = rootdir + 'dapg/' + tissue + '/eqtl.avg.pip'
outfile = rootdir + 'enloc_enrichment_and_priors.txt'

colheader = ['tissue','trait','p1_trait_prior','p2_tissue_prior','p12','pm1','pm2', 'enloc_enrich_intercept_effect(inteffect)', 'inteffect_95ci_lower', 'inteffect_95ci_upper', 'enloc_enrich_log_odds_ratio_effect(logor)', 'logor_95ci_lower', 'logor_95ci_upper']

# outdata = pandas.DataFrame(columns=colheader)

if os.path.isfile(outfile):
  outdata = pandas.read_csv(outfile, sep='\t', header=0)
else:
  outdata = pandas.DataFrame(columns=colheader)


estimate = pandas.read_csv(
  estimatefile,
  delim_whitespace=True,
  header=None,
  names=['effect', 'effect_95l', 'effect_95u']
)

a0 = estimate['effect'][0]
a1 = estimate['effect'][1]

inteffect = estimate['effect'][0]
int95l = estimate['effect_95l'][0]
int95u = estimate['effect_95u'][0]
logor = estimate['effect'][1]
logor95l = estimate['effect_95l'][1]
logor95u = estimate['effect_95u'][1]

qtlpip = pandas.read_csv(qtlpipfile, header=None, names=['prd1'])
prd1 = qtlpip['prd1'][0]

print(a0, a1, prd1)

p1 = (numpy.exp(a0)/(1 + numpy.exp(a0))) * (1 - prd1)
p2 = 1 / (1 + numpy.exp(a0 + a1)) * prd1
p12 = numpy.exp(a0 + a1)/(1 + numpy.exp(a0 + a1)) * prd1
pm1 = p1 + p12
pm2 = p2 + p12

outline = pandas.Series(
  [tissue, trait, p1, p2, p12, pm1, pm2, inteffect, int95l, int95u, logor, logor95l, logor95u],
  index=colheader
)

outdata = outdata.append(outline, ignore_index=True)

outdata.to_csv(outfile, sep='\t', index=False, header=True, float_format='%.5E')


