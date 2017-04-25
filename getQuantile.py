import yaml
from ROOT import *
from math import sqrt, log

import sys
RF = RooFit

CLbound = 68
CLbound = 90

verbose = False
maxtrys = 50
tolerance = 0.003
massLo, massHi = 120., 130.
massMu = 125.09

tf = TFile('ucmles/WS-HGam-STXSMerged.root')
ws = tf.Get('combWS')
mc = ws.obj('ModelConfig')

procs = ["ggH","VBF","WH","ZH","ttH","bbH","tHjb","WtH"]
sigYields = yaml.safe_load(open('sigYields.yml'))['SignalYields']

obsPrefix = "atlas_invMass_"
sigModPrefix = "PDF__ttH_"
sbModPrefix = "_modelSB_"

output = []

pdf = mc.GetPdf()
cat = pdf.indexCat()
dat = ws.obj('combData')

datList = dat.split( cat, True )
nCats   = datList.GetEntries()

sigmas = {}
table = []

for icat in xrange( nCats ):
  cat.setBin(icat)
  catName = cat.getLabel()

  ws.var('mu').setVal( 0. )
  pdfi = pdf.getPdf( catName )
  dati = datList.At( icat )

  obs = pdfi.getObservables( dati ).first()
  obsName = obs.GetName()

  print "%s - %s" % (catName,obsName)
  ws.var('NFcontBkg').setVal(0.0)

  obs.setRange("PartRange", massLo, massHi )
  integ = pdf.createIntegral(RooArgSet(obs),RooArgSet(obs),"PartRange")

  # try each lower edge and scan for upper edge
  alpha = 0.1
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
        print '    edgeLo = %f    edgeHi = %f    quantile = %f    (%f)' % ( edgeLo, edgeHi, quant, diff )
      if ( abs(diff) < tolerance ): break # found range within tolerance, break
      if (diff < 0): xmin = edgeHi
      else: xmax = edgeHi
      itry += 1
      if (itry > maxtrys):
        if (verbose): print ":: ERROR!! :: Did not converge in", maxtrys, "attempts"
        break

    if (edgeHi == massHi): break
    rangeX = (edgeHi-edgeLo)
    if (rangeX < bestRange): bestRange, bestLo, bestHi = rangeX, edgeLo, edgeHi

  obs.setRange( "PartRange", bestLo, bestHi )
  intSig = pdf.createIntegral(RooArgSet(obs),"PartRange")
  print intSig.getVal()
  sigma = bestRange/2.0
  sigmas[catName] = sigma
  table.append( [ catName, round(bestLo,2), round(bestHi,2) ] )

  print ''
  print ''
  print 'BestRange:', [bestLo, bestHi]
  print 'Sigma68:', (bestRange/2.0)
  print ''
  print ''

  sbModName = sbModPrefix+catName
  sbpdf = ws.obj(sbModName)

  #intsb = sbpdf.createIntegral(RooArgSet(obs),"QTRange")

  #ws.loadSnapshot("ucmles_0") # background only, post-fit snapshot
  ws.var('NFcontBkg').setVal(1.0)
  obs.setRange( "QTRange", bestLo, bestHi )
  obs.setRange( "FullRange", 105., 160.)
  ws.var("mu").setVal(0.0)

  #intsb  = sbpdf.createIntegral(RooArgSet(obs),"QTRange")
  #intsbf = sbpdf.createIntegral(RooArgSet(obs),"FullRange")
  intsb  = sbpdf.createIntegral(RooArgSet(obs),RooArgSet(obs),'QTRange')
  intsbf = sbpdf.createIntegral(RooArgSet(obs),RooArgSet(obs))

  b = intsb.getVal() * ws.obj('nbkg_Hgg_'+catName).getVal()

  #ws.var("mu").setVal(1.0)
  #Ns = 0.
  #for proc in procs: Ns += ws.obj("yield_"+proc+"_"+catName).getVal()

  Ns = sigYields[catName]
  s = CLbound*0.01*Ns
  s = Ns

  #ws.loadSnapshot("ucmles_1") # signal+background, post-fit snapshot
  #ws.loadSnapshot("ucmles") # signal+background, post-fit snapshot
  #Ns2 = 0.
  #for proc in procs: Ns2 += ws.obj("yield_"+proc+"_"+catName).getVal()
  #s2 = CLbound*0.01*Ns2
  #s2 = s

  f = s/float(s+b)
  Z = sqrt( 2*((s+b)*log(1+s/b)-s) )

  CL = CLbound
  #line = "  %20s :    B%d = %8.2f    S%d = %8.2f" % (catName, CL, b, CL, s )
  #line = " %25s  &  %8.2f  &  %6.2f  &  %4.2f  &  %4.2f " % (catName, b, s, f, Z )

  if (len(output) < 1):
    labels = ( "Category", "$B_{%d}$"%CL, "$S_{%d}$"%CL, "$f_{%d}$"%CL, "$Z_{%d}$"%CL )
    line = " %25s  &  %8s  &  %8s  &  %8s  &  %8s  \\\\" % labels
    output.append(line)
    output.append("\\hline")

  line = " %25s  &  %8.2f  &  %8.2f  &  %8.2f  &  %8.2f  \\\\" % (catName, b, s, f, Z )
  output.append(line)

  print '\n\n'
  print line
  print 'Bg yield:', ws.obj('nbkg_Hgg_'+catName).getVal()
  print 'Int Part:', intsb.getVal()
  print 'Int Full:', intsbf.getVal()
  print 'Int Range:', [bestLo, bestHi]
  print '\n\n'


print '\n'
for row in table: print row
print '\n'

print '%15s : %6s' % ('Category', '%d%% SR' % CLbound )
for catName in sigmas:
  print '%15s : %6.2f' % (catName, sigmas[catName])

print '\n'
for line in output:
  print line
