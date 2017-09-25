#!/usr/bin/env python
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import *
import yaml

RF = RooFit
gROOT.SetBatch(ROOT.kTRUE)

rename = yaml.load(open('rename_chan.yml'))

# fix TLatex from making everything bold
# ---------------------------------------------------------------
def DrawNDC(self, x, y, text): self.DrawLatexNDC( x, y, '#bf{ %s }' % text )
TLatex.DrawNDC = DrawNDC

def iterate( args ):
  iter = args.createIterator()
  var = iter.Next()
  while var:
    yield var
    var = iter.Next()

def drawLegend(x,y):
  tl = TLine()
  tt = TLatex()
  tt.SetTextSize(0.038)
  tl.SetLineWidth(3)
  tl.SetLineColor(kRed)
  tl.DrawLineNDC(x,y,x+0.04,y)
  tt.DrawNDC( x+0.05, y-0.015, 'Fitted Signal')
  y -= 0.05
  tl.SetLineColor(kGreen+1)
  tl.DrawLineNDC(x,y,x+0.04,y)
  tt.DrawNDC( x+0.05, y-0.015, 'SM Signal')
  y -= 0.05
  tl.SetLineStyle(2)
  tl.SetLineColor(kBlue)
  tl.DrawLineNDC(x,y,x+0.04,y)
  tt.DrawNDC( x+0.05, y-0.015, 'BG Only')

import os, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option( '-f', '--file', type=str, default='workspaces/WS-HGam-STXS.root' )

parser.add_option( '-p', '--poi', type=str, default='mu' )
parser.add_option( '-d', '--dataset', type=str, default='AsimovSB' )
parser.add_option( '-w', '--workspace', type=str, default='combWS' )
parser.add_option( '--modelconfig', type=str, default='ModelConfig' )
parser.add_option( '--nbins', type=int, default=55 )

opt, args = parser.parse_args()

tf = TFile( opt.file )
ws = tf.Get( opt.workspace )
mc = ws.obj( opt.modelconfig )

#ws.saveSnapshot( 'orig', mc.GetParametersOfInterest() )

dat = ws.obj( opt.dataset )
mu  = ws.obj( opt.poi )

pdf = mc.GetPdf()
cat = pdf.indexCat()

datList = dat.split( cat, True )
nCats   = datList.GetEntries()
#bins = RooBinning( 55, 105.,160. )

can = TCanvas()
can.SetMargin( 0.12, 0.04, 0.14, 0.04 )
#can.SetTopMargin( 0.08 )
#can.SetRightMargin( 0.05 )

outPATH = 'plots/dataPdf/%s' % opt.dataset
if not os.path.exists( outPATH ):
  os.makedirs( outPATH )


for icat in xrange( nCats ):
  cat.setBin(icat)
  chanName = cat.getLabel()
  print '\n--> Category : {}'.format( chanName )

  pdfi = pdf.getPdf( chanName )
  dati = datList.At( icat )

  obs = pdfi.getObservables( dati ).first()
  obs.Print()

  RFNorm = RF.Normalization(1.0,RooAbsReal.RelativeExpected)

  frame = obs.frame( 105, 160, opt.nbins )
  frame.SetTitle("")
  dati.plotOn( frame, RF.XErrorSize(0), RF.DataError( RooAbsData.Poisson ) )

  # Draw BG only
  ws.loadSnapshot('original')
  mu.setVal( 0. )
  pdfi.plotOn( frame, RF.LineStyle(2), RF.LineColor( kBlue ), RFNorm )

  # Draw SM expectation
  mu.setVal( 1.0 )
  for poi in iterate(mc.GetParametersOfInterest()): poi.setVal( 1. )
  pdfi.plotOn( frame, RF.LineStyle(1), RF.LineColor( kGreen+1 ), RFNorm )

  # Draw post-fit
  mu.setVal( 1. )
  ws.loadSnapshot('ucmles')
  pdfi.plotOn( frame, RF.LineStyle(1), RF.LineColor( kRed ), RFNorm )

  # Re-draw data to go over pdfs
  dati.plotOn( frame, RF.XErrorSize(0), RF.DataError( RooAbsData.Poisson ) )

  frame.Draw()

  # Aesthetics
  can.cd()
  frame.SetMinimum( 1.0E-03 )
  frame.SetMaximum( frame.GetMaximum() * 1.2 )
  frame.GetYaxis().SetTitleSize(0.042)
  frame.GetXaxis().SetTitleSize(0.042)
  frame.GetYaxis().SetTitleOffset(1.25)
  frame.GetXaxis().SetTitleOffset(1.25)
  frame.GetYaxis().SetTitle('Events / GeV')
  frame.GetXaxis().SetTitle('m_{#gamma#gamma} [GeV]')
  frame.Draw('SAME')

  drawLegend(0.40,0.92)

  l = TLatex()
  l.SetTextSize(0.05)
  l.DrawNDC( 0.69, 0.900, '#bf{#it{ATLAS}} Internal')
  l.SetTextSize(0.042)
  l.DrawNDC( 0.66, 0.845, '#sqrt{s} = 13 TeV, 36.1 fb^{-1}')
  l.DrawNDC( 0.62, 0.795, 'H#rightarrow#gamma#gamma, m_{H} = 125.09 GeV')
  l.DrawNDC( 0.15, 0.900, rename[chanName])

  can.Update()
  can.SaveAs( os.path.join( outPATH, 'Cat%d_%s.pdf' % (icat,chanName) ) )
  can.Clear()

print ''
print '*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'
print '   Plots written to %s' % outPATH
print '*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'
print ''
