import ROOT as r
import numpy as np
from scipy.optimize import root, fsolve

r.gStyle.SetOptStat(0)

# Fake profile LH curve
def f(x):
  if (x<1): return 1.9*(x-1)**2
  return 1.2*(x-1)**2

# Prepare TGraphs
colors = [None, r.kGreen-9, r.kYellow-9]
tg = [ r.TGraph() for i in xrange(3) ]
for i in xrange(3):
  tg[i].SetLineColor(r.kBlack)
  tg[i].SetLineWidth(1)
  if colors[i]: tg[i].SetFillColor( colors[i] )

# Get points
npts = 121
xmin, xmax, ymax = -0.5, 3.5, 6.0
pts = [ (x, f(x)) for x in np.linspace( xmin, xmax, npts ) ]

# Fill TGraphs
for i in xrange(npts): tg[0].SetPoint( i, pts[i][0], pts[i][1] )

# Get spline and find 1 sigma and 2 sigma intercepts
sp = r.TSpline3('s',tg[0])
x0  = root(lambda x : sp.Eval(x), x0=1.0).x[0]
x2p = root(lambda x: np.abs(4 - sp.Eval(x)), x0=xmax).x[0]
x2m = root(lambda x: np.abs(4 - sp.Eval(x)), x0=xmin).x[0]
x1p = root(lambda x: np.abs(1 - sp.Eval(x)), x0=xmax).x[0]
x1m = root(lambda x: np.abs(1 - sp.Eval(x)), x0=xmin).x[0]
xbs = [ None, (x1m,x1p), (x2m,x2p) ]
errors = [ ( abs(x0-x1p), -abs(x0-x1m) ), ( abs(x0-x2p), -abs(x0-x2m) ) ]

print ' mu = %.3f +/- (%.3f,%.3f) ++/-- (%.3f,%.3f)' % ( x0, errors[0][0], errors[0][1], errors[1][0], errors[1][1] )

# Make 1 sigma and 2 sigma bands
n = 501
xs = np.linspace(x1m,x1p,n)
for i in xrange(n):
  tg[1].SetPoint( i,   xs[i],  sp.Eval(xs[i]) )
  tg[1].SetPoint( n+i, xs[n-i-1], ymax )

n = 501
xs = np.linspace(x2m,x2p,n)
for i in xrange(n):
  tg[2].SetPoint( i,   xs[i],  sp.Eval(xs[i]) )
  tg[2].SetPoint( n+i, xs[n-i-1], ymax )

# Prep the canvas
can = r.TCanvas()
can.cd()
can.SetMargin( 0.10, 0.05, 0.15, 0.05 )
h = r.TH1F('hist','',npts,xmin,xmax)
h.SetMaximum(ymax+0.01)
h.SetMinimum(1e-06)
h.GetXaxis().SetTitle('#mu')
h.GetYaxis().SetTitle('#lambda(#mu)')
h.GetYaxis().SetTitleOffset(0.7)
h.GetXaxis().SetTitleSize(0.05)
h.GetYaxis().SetTitleSize(0.05)
h.Draw('HIST')

# Draw guides
tl = r.TLine()
tt = r.TLatex()
tl.SetLineStyle(2)

colors = [None, r.kGreen-3, r.kOrange-2]
for i in xrange(1,3):
  tt.SetTextColor( colors[i] )
  tl.SetLineColor( colors[i] )
  tt.DrawLatex( 0.91*xmax, i**2-0.35, '%d #sigma' % i)
  tl.DrawLine( xmin, i**2, xmax, i**2 )
  if (i==1):
    tl.DrawLine( xbs[i][0], 0, xbs[i][0], ymax )
    tl.DrawLine( xbs[i][1], 0, xbs[i][1], ymax )

# Draw contours
tg[2].Draw('F SAME')
tg[1].Draw('F SAME')
tg[0].Draw('C SAME')

# Draw more guides
tl.SetLineColor(r.kBlack)
r.gStyle.SetLineStyleString(11,'15 45');
for i in xrange( int(xmin), int(xmax)+1 ):
  if (i==1):
    tl.SetLineWidth(1)
    tl.SetLineStyle(3)
  else:
    tl.SetLineWidth(1)
    tl.SetLineStyle(11)
  tl.DrawLine( i, 0, i, ymax )

tl.SetLineWidth(1)
tl.SetLineStyle(1)
tl.DrawLine( x0, 0, x0, ymax )

# Tidy up
can.SaveAs('profile.png')
can.SaveAs('profile.pdf')
h.Draw('AXIS SAME')
raw_input('Done?')

