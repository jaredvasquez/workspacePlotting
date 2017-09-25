import yaml
from ROOT import *
from math import sqrt, log

import sys
RF = RooFit

CLbound = 68
CLbound = 90

verbose = False
maxtrys = 50
tolerance = 0.005
massLo, massHi = 120., 130.
massMu = 125.09

tf = TFile('ucmles/WS-HGam-STXSMerged.root')
ws = tf.Get('combWS')
mc = ws.obj('ModelConfig')

procs = ["ggH","VBF","WH","ZH","ttH","bbH","tHjb","WtH"]
sigYields = yaml.safe_load(open('sigYields.yml'))['SignalYields']

obsPrefix = "atlas_invMass_"
sigModPrefix = "PDF__ttH_"
sbModPrefix = "modelSB_"

output = []

pdf = mc.GetPdf()
cat = pdf.indexCat()
dat = ws.obj('combData')

datList = dat.split( cat, True )
nCats   = datList.GetEntries()

import yaml
ranges = yaml.safe_load(open('ranges.yml'))

if True:
  sbModName = sbModPrefix+catName
  sbpdf = ws.obj(sbModName)

  #intsb = sbpdf.createIntegral(RooArgSet(obs),"QTRange")
  intsb = sbpdf.createIntegral(RooArgSet(obs),RooArgSet(obs),"QTRange")

  ws.loadSnapshot("ucmles_0") # background only, post-fit snapshot
  obs.setRange( "QTRange", edgeLo, edgeHi )
  ws.var("mu").setVal(0.0)
  b = intsb.getVal() * ws.obj('yield_background_'+catName).getVal()

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
  s2 = s

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
