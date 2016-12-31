import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import *

RF = RooFit
gROOT.SetBatch(ROOT.kTRUE)

import os, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option( '-f', '--file', type=str, default='workspace/WS-HGam-STXS.root' )

parser.add_option( '-p', '--poi', type=str, default='mu' )
parser.add_option( '-d', '--dataset', type=str, default='combData' )
parser.add_option( '-w', '--workspace', type=str, default='combWS' )
parser.add_option( '--modelconfig', type=str, default='ModelConfig' )
parser.add_option( '--nbins', type=int, default=55 )

opt, args = parser.parse_args()

outPATH = 'plots/sidebands/%s' % opt.dataset
if not os.path.exists( outPATH ):
  os.makedirs( outPATH )

def iterate( args ):
  iter = args.createIterator()
  var = iter.Next()
  while var:
    yield var
    var = iter.Next()

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
can.SetTopMargin( 0.08 )
can.SetRightMargin( 0.05 )

for icat in xrange( nCats ):
  cat.setBin(icat)
  chanName = cat.getLabel()
  print '\n--> Category : {}'.format( chanName )

  pdfi = pdf.getPdf( chanName )
  dati = datList.At( icat )

  obs = pdfi.getObservables( dati ).first()
  obs.Print()
  mu.setVal( 0. )

  obs.setRange( 'SBLo', 105, 120 )
  obs.setRange( 'SBHi', 130, 160 )

  nps = mc.GetNuisanceParameters()
  for np in iterate(nps): np.setConstant()
  pdfi.fitTo( dati, RF.Range('SBLo,SBHi'), RF.Extended(kTRUE) )

  frame = obs.frame( 105, 160, opt.nbins )
  dati.plotOn( frame, RF.XErrorSize(0) )
  #dati.plotOn( frame, RF.DataError( RooAbsData.Poisson ), RF.XErrorSize(0) )
  frame.Draw()

  # Draw BG only
  pdfi.plotOn( frame, RF.Range('SBLo,SBHi'), RF.LineStyle(2), RF.LineColor( kBlue ) )

  # Aesthetics
  can.cd()
  frame.SetMinimum( 1.0E-03 )
  frame.SetMaximum( frame.GetMaximum() * 1.6 )
  frame.GetYaxis().SetTitle('Events / GeV')
  frame.GetXaxis().SetTitle('m_{#gamma#gamma} [GeV]')
  frame.Draw('SAME')

  l = TLatex()
  l.SetNDC()
  l.SetTextSize(0.042)
  l.DrawLatex( 0.650, 0.785, '#sqrt{s} = 13 TeV, 13.3 fb^{-1}')
  l.DrawLatex( 0.610, 0.740, 'H#rightarrow#gamma#gamma, m_{H} = 125.09 GeV')

  can.Update()
  can.SaveAs( os.path.join( outPATH, 'Cat%d_%s.png' % (icat,chanName) ) )
  can.Clear()

print ''
print '*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'
print '   Plots written to %s' % outPATH
print '*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'
print ''
