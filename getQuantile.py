from ROOT import *
from math import sqrt, log

import sys
RF = RooFit

CLbound = 68
CLbound = 90

verbose = False
maxtrys = 50
tolerance = 0.0001
massLo, massHi = 105., 160.
massMu = 125.09

tf = TFile('workspace/WS-HGam-STXS.root')
ws = tf.Get('combWS')
mc = ws.obj('ModelConfig')

procs = ["ggH","VBF","WH","ZH","ttH","bbH","tHjb","WtH"]

obsPrefix = "atlas_invMass_"
sigModPrefix = "PDF__ggH_"
sbModPrefix = "modelSB_"

output = []

pdf = mc.GetPdf()
cat = pdf.indexCat()
dat = ws.obj('combData')

datList = dat.split( cat, True )
nCats   = datList.GetEntries()

sigmas = {}

for icat in xrange( nCats ):
  cat.setBin(icat)
  catName = cat.getLabel()

  ws.var('NFcontBkg').setVal( 0. )
  pdfi = pdf.getPdf( catName )
  dati = datList.At( icat )

  obs = pdfi.getObservables( dati ).first()
  obsName = obs.GetName()

  print "%s - %s" % (catName,obsName)
  ws.var('NFcontBkg').setVal(0.0)

  obs.setRange("PartRange", massLo, massHi )
  integ = pdf.createIntegral(RooArgSet(obs),RooArgSet(obs),"PartRange")

  # try each lower edge and scan for upper edge
  alpha = 0.01
  bestRange, bestLo, bestHi = 999, 105, 160
  for iLo in xrange( int((massMu-massLo)/alpha) ):
    edgeLo = massLo + alpha * iLo
    itry, xmin, xmax, diff = 0, massMu, massHi, 999
    while (True):
      edgeHi = (xmax+xmin)*0.5
      obs.setRange( "PartRange", edgeLo, edgeHi )
      quant = integ.getVal()
      prevdiff = diff
      diff =  ( quant - 0.01*CLbound )
      if (verbose):
        print '    edgeHi = %f    quantile = %f    (%f)' % ( edgeHi, quant, diff )
      if ( abs(diff) < tolerance ): break # found range within tolerance, break
      if (diff < 0): xmin = edgeHi
      else: xmax = edgeHi
      itry += 1
      if (itry > maxtrys):
        if (verbose): print ":: ERROR!! :: Did not converge in", maxtrys, "attempts"
        break

    if (edgeHi == massHi): break
    range = (edgeHi-edgeLo)
    if (range < bestRange): bestRange, bestLo, bestHi = range, edgeLo, edgeHi

  #print bestRange, bestLo, bestHi
  #obs.setRange( "PartRange", bestLo, bestHi )
  #intSig = pdf.createIntegral(RooArgSet(obs),"PartRange")
  #print intSig.getVal()
  sigma = bestRange/2.0
  sigmas[catName] = sigma

print '%15s : %6s' % ('Category', '%d%% SR' % CLbound )
for catName in sigmas:
  print '%15s : %6.2f' % (catName, sigmas[catName])



import sys; sys.exit()
if False:
  sbModName = sbModPrefix+catName
  sbpdf = ws.obj(sbModName)

  #intsb = sbpdf.createIntegral(RooArgSet(obs),"QTRange")
  intsb = sbpdf.createIntegral(RooArgSet(obs),RooArgSet(obs),"QTRange")

  ws.loadSnapshot("ucmles_0") # background only, post-fit snapshot
  obs.setRange( "QTRange", edgeLo, edgeHi )
  ws.var("mu").setVal(0.0)
  b = intsb.getVal() * ws.obj('yield_background_'+catName).getVal()

  ws.var("mu").setVal(1.0)
  Ns = 0.
  for proc in procs: Ns += ws.obj("yield_"+proc+"_"+catName).getVal()
  s = CLbound*0.01*Ns
  s = Ns

  #ws.loadSnapshot("ucmles_1") # signal+background, post-fit snapshot
  ws.loadSnapshot("ucmles_prod5") # signal+background, post-fit snapshot
  Ns2 = 0.
  for proc in procs: Ns2 += ws.obj("yield_"+proc+"_"+catName).getVal()
  s2 = CLbound*0.01*Ns2

  f = s/float(s+b)
  Z = sqrt( 2*((s+b)*log(1+s/b)-s) )

  CL = CLbound
  #line = "  %20s :    B%d = %8.2f    S%d = %8.2f" % (catName, CL, b, CL, s )
  #line = " %25s  &  %8.2f  &  %6.2f  &  %4.2f  &  %4.2f " % (catName, b, s, f, Z )

  labels = ( "Category", "$B_{%d}$"%CL, "$S_{%d}$"%CL, "$f_{%d}$"%CL, "$Z_{%d}$"%CL, "$S^{fit}_{%d}$"%CL )

  if (len(output) < 1):
    line = " %25s  &  %8s  &  %8s  &  %8s  &  %8s  &  %14s  \\\\" % labels
    output.append(line)
    output.append("\\hline")

  #line = " %25s  :  %8.2f  " % (cattitles[icat], s/float(b))
  line = " %25s  &  %8.2f  &  %8.2f  &  %8.2f  &  %8.2f  &  %14.2f  \\\\" % (cattitles[icat], b, s, f, Z, s2 )
  output.append(line)
