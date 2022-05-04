##-----------------------------------------------------------------------------
## Facets Re-Analysis of solid and CSF samples
##-----------------------------------------------------------------------------
##
## start: 05/04/2022
## chris-kreitzer
##

clean()
gc()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')

library(facets)
library(patchwork)

x = facets::readSnpMatrix(filename = '06287_BJ/FROZENPOOLEDNORMAL_IMPACT505_V2.rg.md.abra.printreads__s_C_001327_R018_d.rg.md.abra.printreads.dat.gz')
x1 = facets::preProcSample(x, gbuild = 'hg19')
x2 = facets::procSample(x1, cval = 100)
x3 = facets::emcncf(x2)


y = facetsSuite::run_facets(x, cval = 100, genome = 'hg19', snp_nbhd = 300)
a = facetsSuite::cnlr_plot(y)
b = facetsSuite::valor_plot(y)
c = facetsSuite::icn_plot(y)

a/b/c

facets_fit_qc(y)
