from ROOT import *
RF = RooFit

import os
outPATH = 'plots/asimov/'
if not os.path.exists( outPATH ):
  os.makedirs( outPATH )

tf = TFile('workspace/WS-HGam-STXS.root')
ws = tf.Get('combWS')
mc = ws.obj('ModelConfig')

dat = ws.obj('AsimovSB')
mu  = ws.obj('mu')

pdf = mc.GetPdf()
cat = pdf.indexCat()

datList = dat.split( cat, True )
nCats   = datList.GetEntries()

bins = RooBinning( 55, 105.,160. );

can = TCanvas()
can.SetTopMargin( 0.08 )
can.SetRightMargin( 0.05 )

for icat in xrange(nCats):
  cat.setBin(icat)
  chanName = cat.getLabel()
  print '  {}'.format( chanName )

  pdfi = pdf.getPdf( chanName )
  dati = datList.At( icat )

  obs = pdfi.getObservables( dati ).first()
  obs.Print()

  frame = obs.frame()
  dati.plotOn( frame, RF.XErrorSize(0), RF.Binning(bins) )
  #dati.plotOn( frame, RF.DataError( RooAbsData.Poisson ), RF.XErrorSize(0), RF.Binning(bins) )
  frame.Draw()

  # Draw BG only
  mu.setVal( 0. )
  pdfi.plotOn( frame, RF.LineStyle(2), RF.LineColor( kBlue ) )

  # Draw SM expectation
  mu.setVal( 1. )
  pdfi.plotOn( frame, RF.LineStyle(1), RF.LineColor( kRed ) )

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

#nCats = cat.numBins()
#dataList = data.split( cat, True )
#print 'nCats = %d' % nCats
#print mc.Print()
