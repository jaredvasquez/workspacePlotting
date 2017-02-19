from math import sqrt
from ROOT import *

vev = 246. # GeV
kappas = {
       # class, mass, kappa, error
  't' : ('F', 173.34, 1.00, 0.94, 0.19 ),
  'W' : ('V', 80.385, 1.00, 0.10, 0.12 ),
  'Z' : ('V', 91.188, 1.00, 2.00, 4.02 ),
}

labels = []

gROOT.LoadMacro("AtlasStyle/AtlasStyle.C")
SetAtlasStyle()

tc = TCanvas(); tc.cd()
tc.SetLogx(); tc.SetLogy()

tg = TGraphAsymmErrors()
tg.SetLineWidth(2)
tg.SetMarkerSize(1)
tg.SetMarkerStyle(20)


for i, Ci in enumerate(list(kappas)):
  ptype, mass, kappa, errHI, errLO = kappas[Ci]

  y = kappa*mass/vev
  yHI = errHI*mass/vev
  yLO = errLO*mass/vev
  if (ptype=='V'):
    y = sqrt(kappa)*mass/vev
    yHI *= (0.5/sqrt(kappa))
    yLO *= (0.5/sqrt(kappa))

  tg.SetPoint( i, mass, y )
  tg.SetPointError( i, 0, 0, yHI, yLO )

  labels.append( ( Ci, mass, y-yLO ) )
  print '  %s = %.2f +/- (%+.2f, %+.2f)' % (Ci,y,yHI,-1*yLO), ptype
print ''

tl = TF1( 'tl', 'x/[0]', 0.5, 500 );
tl.SetParameter(0, vev)
tl.SetLineColor(kBlue)
tl.SetLineStyle(2)

#tg.GetXaxis().SetTitle('Particle Mass [GeV]')
#tg.GetYaxis().SetTitle('#kappa_{F} #frac{m_{F}}{v} or #sqrt{#kappa_{V}} #frac{m_{V}}{v}')

h = TH1F( 'h', '', 1000, 5, 500 )
h.GetXaxis().SetRangeUser(5,500)
h.GetYaxis().SetRangeUser(0.01,3.5)
h.GetXaxis().SetTitle('Particle Mass [GeV]')
h.GetYaxis().SetTitle('#kappa_{F} #frac{m_{F}}{v} or #sqrt{#kappa_{V}} #frac{m_{V}}{v}')
h.Draw('HIST')

#tg.Draw('APE')
#tg.GetXaxis().SetRangeUser(0.5,500)
#tg.GetYaxis().SetRangeUser(0.003,3.5)
tl.Draw('SAME')
tg.Draw('PE SAME')

gPad.RedrawAxis()
gPad.SetTicks(1,1)

# Draw labels for each point
tt = TLatex();
tt.DrawLatex(  56., 0.37, 'W' )
tt.DrawLatex( 100., 0.22, 'Z' )
tt.DrawLatex( 145., 0.90, 't' )

# Draw Plot Labels
tt.SetNDC()
tt.DrawLatex( 0.21, 0.87, '#it{#bf{ATLAS}} Internal' )
tt.SetTextSize( 0.04 )
tt.DrawLatex( 0.20, 0.815, '#sqrt{s} = 13 TeV, 36.5 fb^{-1}' )

ll = TLine(); ll.SetLineWidth(2)
ll.DrawLineNDC( 0.22, 0.76, 0.28, 0.76 )
tt.DrawLatex( 0.3, 0.75, 'Observed' )
ll.SetLineColor(kBlue); ll.SetLineStyle(2)
ll.DrawLineNDC( 0.22, 0.70, 0.28, 0.70 )
tt.DrawLatex( 0.3, 0.69, 'SM Expected' )


tc.SaveAs('plots/kappaMass.pdf')
tc.SaveAs('plots/kappaMass.png')

#raw_input('done?')
