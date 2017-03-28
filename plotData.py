import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import *

RF = RooFit
gROOT.SetBatch(ROOT.kTRUE)

# fix TLatex from making everything bold
# ---------------------------------------------------------------
def DrawNDC(self, x, y, text): self.DrawLatexNDC( x, y, '#bf{ %s }' % text )
TLatex.DrawNDC = DrawNDC

import os, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option( '-f', '--file', type=str, default='workspace/WS-HGam-STXS.root' )

parser.add_option( '-p', '--poi', type=str, default='mu' )
parser.add_option( '-d', '--dataset', type=str, default='AsimovSB' )
parser.add_option( '-w', '--workspace', type=str, default='combWS' )
parser.add_option( '--modelconfig', type=str, default='ModelConfig' )
parser.add_option( '--nbins', type=int, default=55 )

opt, args = parser.parse_args()

tf = TFile( opt.file )
ws = tf.Get( opt.workspace )
mc = ws.obj( opt.modelconfig )

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
  mu.setVal( 0. )
  pdfi.plotOn( frame, RF.LineStyle(2), RF.LineColor( kBlue ), RFNorm )

  # Draw SM expectation
  mu.setVal( 1. )
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

  l = TLatex()
  l.SetTextSize(0.05)
  l.DrawNDC( 0.69, 0.900, '#bf{#it{ATLAS}} Internal')
  l.SetTextSize(0.042)
  l.DrawNDC( 0.66, 0.845, '#sqrt{s} = 13 TeV, 36.1 fb^{-1}')
  l.DrawNDC( 0.62, 0.795, 'H#rightarrow#gamma#gamma, m_{H} = 125.09 GeV')
  l.DrawNDC( 0.15, 0.900, chanName)

  can.Update()
  can.SaveAs( os.path.join( outPATH, 'Cat%d_%s.pdf' % (icat,chanName) ) )
  can.Clear()

print ''
print '*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'
print '   Plots written to %s' % outPATH
print '*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'
print ''
