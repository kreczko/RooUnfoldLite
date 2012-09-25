//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldResponse.cxx 309 2011-10-10 20:40:14Z T.J.Adye $
//
// Description:
//      Response Matrix
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
 <p> Class to create response object as used in RooUnfold </p>
 <p> Contains measured and truth distributions as TH1s and the response matrix as a TH2. Also contains methods for handling these data</p>
<p> Can handle 1,2 or 3 dimensional histograms and return vectors and matrices of their bin content and error (1 and 2D distributions respectively).
 Conversely can also convert these vectors and matrices into TH1Ds and TH2Ds. </p>
<p> Can also take a variety of parameters as inputs. This includes maximum and minimum values, distributions and vectors/matrices of values. </p>
<p> This class does the numerical modifications needed to allow unfolding techniques to work in the unfolding routines used in RooUnfold. </p>
END_HTML */

/////////////////////////////////////////////////////////////

#include "../interface/RooUnfoldResponse.h"

#include <iostream>
#include <assert.h>
#include <cmath>

#include "TClass.h"
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TMatrixD.h"

using std::cout;
using std::cerr;
using std::endl;
using std::pow;
using std::sqrt;

ClassImp (RooUnfoldResponse);

RooUnfoldResponse::RooUnfoldResponse (const RooUnfoldResponse& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  // RooUnfoldResponse copy constructor
  Init();
  Setup (rhs);
}

RooUnfoldResponse::RooUnfoldResponse (int nb, double xlo, double xhi,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponse constructor - simple 1D case with same binning, measured vs truth
  Init();
  Setup (nb, xlo, xhi);
}

RooUnfoldResponse::RooUnfoldResponse (int nm, double mlo, double mhi, int nt, double tlo, double thi,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponse constructor - simple 1D case
  Init();
  Setup (nm, mlo, mhi, nt, tlo, thi);
}

RooUnfoldResponse::RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2* response,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponse constructor - create from already-filled histograms
  // "response" gives the response matrix, measured X truth.
  // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  // in "truth" for unmeasured events (inefficiency).
  // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  Init();
  Setup (measured, truth, response);
}

RooUnfoldResponse::RooUnfoldResponse (const TH1* measured, const TH1* truth,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponse constructor - measured and truth only used for shape
  Init();
  Setup (measured, truth);
}

RooUnfoldResponse&
RooUnfoldResponse::operator= (const RooUnfoldResponse& rhs)
{
  // RooUnfoldResponse assignment operator
  if (this == &rhs) return *this;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  return Setup (rhs);
}

void
RooUnfoldResponse::Add (const RooUnfoldResponse& rhs)
{
  // Add another RooUnfoldResponse, accumulating contents
  if (_res == 0) {
    Setup (rhs);
    return;
  }
  assert (_mdim==rhs._mdim);
  assert (_tdim==rhs._tdim);
  assert (_mes != 0 && rhs._mes != 0);
  assert (_fak != 0 && rhs._fak != 0);
  assert (_tru != 0 && rhs._tru != 0);
  assert (_res != 0 && rhs._res != 0);
  if (_cached) ClearCache();
  _mes->Add (rhs._mes);
  _fak->Add (rhs._fak);
  _tru->Add (rhs._tru);
  _res->Add (rhs._res);
}

RooUnfoldResponse&
RooUnfoldResponse::Reset()
{
  // Resets object to initial state.
  ClearCache();
  delete _mes;
  delete _fak;
  delete _tru;
  delete _res;
  return Setup();
}

RooUnfoldResponse&
RooUnfoldResponse::Init()
{
  _overflow= 0;
  return Setup();
}

RooUnfoldResponse&
RooUnfoldResponse::Setup()
{
  _tru= _mes= _fak= 0;
  _res= 0;
  _vMes= _eMes= _vFak= _vTru= _eTru= 0;
  _mRes= _eRes= 0;
  _nm= _nt= _mdim= _tdim= 0;
  _cached= false;
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const RooUnfoldResponse& rhs)
{
  // Copy data from another RooUnfoldResponse
  _overflow= rhs._overflow;
  return Setup (rhs.Hmeasured(), rhs.Htruth(), rhs.Hresponse());
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (int nm, double mlo, double mhi, int nt, double tlo, double thi)
{
  // set up simple 1D case
  Reset();
  bool oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= new TH1D ("measured", "Measured", nm, mlo, mhi);
  _fak= new TH1D ("fakes",    "Fakes",    nm, mlo, mhi);
  _tru= new TH1D ("truth",    "Truth",    nt, tlo, thi);
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  SetNameTitleDefault ("response", "Response");
  _res= new TH2D (GetName(), GetTitle(), nm, mlo, mhi, nt, tlo, thi);
  TH1::AddDirectory (oldstat);
  return *this;
}

void
RooUnfoldResponse::ReplaceAxis (TAxis* dest, const TAxis* source)
{
  // Replaces an axis with that of a different histogram
  TObject* hist= dest->GetParent();
  source->Copy (*dest);
  dest->SetParent (hist);
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const TH1* measured, const TH1* truth)
{
  // set up - measured and truth only used for shape
  Reset();
  bool oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= (TH1*) measured ->Clone();
  _mes->Reset();
  _fak= (TH1*) _mes     ->Clone("fakes");
  _fak->SetTitle("Fakes");
  _tru= (TH1*) truth    ->Clone();
  _tru->Reset();
  _mdim= _mes->GetDimension();
  _tdim= _tru->GetDimension();
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  SetNameTitleDefault();
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  _res= new TH2D (GetName(), GetTitle(), _nm, 0.0, double(_nm), _nt, 0.0, double(_nt));
  if (_mdim==1) ReplaceAxis (_res->GetXaxis(), _mes->GetXaxis());
  if (_tdim==1) ReplaceAxis (_res->GetYaxis(), _tru->GetXaxis());
  TH1::AddDirectory (oldstat);
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const TH1* measured, const TH1* truth, const TH2* response)
{
  // Set up from already-filled histograms.
  // "response" gives the response matrix, measured X truth.
  // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  // in "truth" for unmeasured events (inefficiency).
  // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  Reset();
  bool oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _res= (TH2*) response->Clone();
  if (measured) {
    _mes= (TH1*) measured->Clone();
    _fak= (TH1*) measured->Clone("fakes");
    _fak->Reset();
    _fak->SetTitle("Fakes");
    _mdim= _mes->GetDimension();
  } else {
    _mes= new TH1D ("measured", "Measured", response->GetNbinsX(), 0.0, 1.0);
    ReplaceAxis (_mes->GetXaxis(), _res->GetXaxis());
    _fak= (TH1*) _mes->Clone("fakes");
    _fak->SetTitle("Fakes");
    _mdim= 1;
  }
  if (truth) {
    _tru= (TH1*)  truth   ->Clone();
    _tdim= _tru->GetDimension();
  } else {
    _tru= new TH1D ("truth",    "Truth",    response->GetNbinsY(), 0.0, 1.0);
    ReplaceAxis (_tru->GetXaxis(), _res->GetYaxis());
    _tdim= 1;
  }
  TH1::AddDirectory (oldstat);
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  if (_nm != _res->GetNbinsX() || _nt != _res->GetNbinsY()) {
    cerr << "Warning: RooUnfoldResponse measured X truth is " << _nm << " X " << _nt
         << ", but matrix is " << _res->GetNbinsX()<< " X " << _res->GetNbinsY() << endl;
  }

  int first=1, nm= _nm, nt= _nt, s= _res->GetSumw2N();
  if (_overflow) {
    first= 0;
    nm += 2;
    nt += 2;
  }

  if (!measured || _mes->GetEntries() == 0.0) {
    // Similar to _res->ProjectionX() but without stupid reset of existing histograms
    // Always include under/overflows in sum of truth.
    for (int i= 0; i<nm; i++) {
      double nmes= 0.0, wmes= 0.0;
      for (int j= 0; j<_nt+2; j++) {
               nmes +=      _res->GetBinContent (i+first, j);
        if (s) wmes += pow (_res->GetBinError   (i+first, j), 2);
      }
      int bin= GetBin (_mes, i, _overflow);
             _mes->SetBinContent (bin,      nmes );
      if (s) _mes->SetBinError   (bin, sqrt(wmes));
    }
  } else {
    // Fill fakes from the difference of _mes - _res->ProjectionX()
    // Always include under/overflows in sum of truth.
    int sm= _mes->GetSumw2N();
    for (int i= 0; i<nm; i++) {
      double nmes= 0.0, wmes= 0.0;
      for (int j= 0; j<_nt+2; j++) {
               nmes +=      _res->GetBinContent (i+first, j);
        if (s) wmes += pow (_res->GetBinError   (i+first, j), 2);
      }
      int bin= GetBin (_mes, i, _overflow);
      double fake= _mes->GetBinContent (bin) - nmes;
      if (!s) wmes= nmes;
      _fak->SetBinContent (bin, fake);
      _fak->SetBinError   (bin, sqrt (wmes + (sm ? pow(_mes->GetBinError(bin),2) : _mes->GetBinContent(bin))));
    }
    _fak->SetEntries (_fak->GetEffectiveEntries());  // 0 entries if 0 fakes
  }

  if (!truth || _tru->GetEntries() == 0.0) {
    // similar to _res->ProjectionY() but without stupid reset of existing histograms
    // Always include under/overflows in sum of measurements.
    for (int j= 0; j<nt; j++) {
      double ntru= 0.0, wtru= 0.0;
      for (int i= 0; i<_nm+2; i++) {
               ntru +=      _res->GetBinContent (i, j+first);
        if (s) wtru += pow (_res->GetBinError   (i, j+first), 2);
      }
      int bin= GetBin (_tru, j, _overflow);
             _tru->SetBinContent (bin,      ntru);
      if (s) _tru->SetBinError   (bin, sqrt(wtru));
    }
  }

  SetNameTitleDefault();
  return *this;
}

void
RooUnfoldResponse::ClearCache()
{
  delete _vMes; _vMes= 0;
  delete _eMes; _eMes= 0;
  delete _vFak; _vFak= 0;
  delete _vTru; _vTru= 0;
  delete _eTru; _eTru= 0;
  delete _mRes; _mRes= 0;
  delete _eRes; _eRes= 0;
  _cached= false;
}

int
RooUnfoldResponse::Fill (double xr, double xt, double w)
{
  // Fill 1D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  _mes->Fill (xr, w);
  _tru->Fill (xt, w);
  return _res->Fill (xr, xt, w);
}

int
RooUnfoldResponse::Fill (double xr, double yr, double xt, double yt, double w)
{
  // Fill 2D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  ((TH2*)_mes)->Fill (xr, yr, w);
  ((TH2*)_tru)->Fill (xt, yt, w);
  return _res->Fill (double(FindBin (_mes, xr, yr))+.5, double(FindBin (_tru, xt, yt))+.5, w);
}

int
RooUnfoldResponse::Fill (double xr, double yr, double zr, double xt, double yt, double zt, double w)
{
  // Fill 3D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  ((TH3*)_mes)->Fill (xr, yr, zr, w);
  ((TH3*)_tru)->Fill (xt, yt, zt, w);
  return _res->Fill (double(FindBin (_mes, xr, yr, zt))+.5, double(FindBin (_tru, xt, yt, zt))+.5, w);
}

int
RooUnfoldResponse::FindBin(const TH1* h, double x, double y)
{
  // Get vector index (0..nx*ny-1) for (x,y) coordinates
  int nx=   h->GetNbinsX();
  int ny=   h->GetNbinsY();
  int binx= h->GetXaxis()->FindBin(x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny;
  int biny= h->GetYaxis()->FindBin(y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny;
  return binx + nx*biny;
}

int
RooUnfoldResponse::FindBin(const TH1* h, double x, double y, double z)
{
  // Get vector index (0..nx*ny*nz-1) for (x,y,z) coordinates
  int nx=   h->GetNbinsX();
  int ny=   h->GetNbinsY();
  int nz=   h->GetNbinsZ();
  int binx= h->GetXaxis()->FindBin(x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny*nz;
  int biny= h->GetYaxis()->FindBin(y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny*nz;
  int binz= h->GetZaxis()->FindBin(z) - 1;
  if (binz <  0)  return -1;
  if (binz >= nz) return nx*ny*nz;
  return binx + nx*(biny + ny*binz);
}

int
RooUnfoldResponse::GetBinDim (const TH1* h, int i)
{
  // Converts from vector index (0..nx*ny-1) or (0..nx*ny*nz-1) to multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) or (0..(nx+2)*(ny+2)*(nz+2)-1), skipping under/overflow bins.
  int ndim= h->GetDimension(), nx= h->GetNbinsX();
  if        (ndim == 2) {
//  cout << i << " -> " << "(" << i%nx+1 << "," << i/nx+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*(i/nx+1);
  } else if (ndim == 3) {
    int ny= h->GetNbinsY();
//  cout << i << " -> " << "(" << i%nx+1 << "," << (i/nx)%ny+1 << "," << i/(nx*ny)+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*((i/nx)%ny+1 + (ny+2)*(i/(nx*ny)+1));
  }
  return i+1;   // not used: 1D handled by inline GetBin() (and handling UseOverflow), don't support >3D.
}

int
RooUnfoldResponse::Miss1D (double xt, double w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 1D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==1);
  if (_cached) ClearCache();
  return _tru->Fill (xt, w);
}

int
RooUnfoldResponse::Miss2D (double xt, double yt, double w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 2D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==2);
  if (_cached) ClearCache();
  return ((TH2*)_tru)->Fill (xt, yt, w);
}

int
RooUnfoldResponse::Miss (double xt, double yt, double zt, double w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 3D Response Matrix
  assert (_tru != 0);
  assert (_tdim==3);
  if (_cached) ClearCache();
  return ((TH3*)_tru)->Fill (xt, yt, zt, w);
}

int
RooUnfoldResponse::Fake1D (double xr, double w)
{
  // Fill fake event (reconstructed event with no truth) into 1D Response Matrix (with weight)
  assert (_fak != 0 && _mes != 0);
  assert (_mdim==1);
  if (_cached) ClearCache();
         _mes->Fill (xr, w);
  return _fak->Fill (xr, w);
}

int
RooUnfoldResponse::Fake2D (double xr, double yr, double w)
{
  // Fill fake event (reconstructed event with no truth) into 2D Response Matrix (with weight)
  assert (_mes != 0);
  assert (_mdim==2);
  if (_cached) ClearCache();
         ((TH2*)_fak)->Fill (xr, yr, w);
  return ((TH2*)_mes)->Fill (xr, yr, w);
}

int
RooUnfoldResponse::Fake (double xr, double yr, double zr, double w)
{
  // Fill fake event (reconstructed event with no truth) into 3D Response Matrix
  assert (_mes != 0);
  assert (_mdim==3);
  if (_cached) ClearCache();
         ((TH3*)_mes)->Fill (xr, yr, zr, w);
  return ((TH3*)_fak)->Fill (xr, yr, zr, w);
}

TH1D*
RooUnfoldResponse::H2H1D(const TH1* h, int nb)
{
  if (dynamic_cast<const TH1D*>(h)) return dynamic_cast<TH1D*>(h->Clone());
  TH1D* h1d= new TH1D(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
  for (int i= 0; i < nb; i++) {
    int j= GetBin (h, i);  // don't bother with under/overflow bins (not supported for >1D)
    h1d->SetBinContent (i+1, h->GetBinContent (j));
    h1d->SetBinError   (i+1, h->GetBinError   (j));
  }
  return h1d;
}

TH2D*
RooUnfoldResponse::HresponseNoOverflow() const
{
  const TH2* h= Hresponse();
  int nx= h->GetNbinsX(), ny= h->GetNbinsY();
  if (_overflow) {  // implies truth/measured both 1D
    double xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nx;
    double ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax(), yb= (yhi-ylo)/ny;
    nx += 2; ny += 2;
    TH2D* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo-xb, xhi+xb, ny, ylo-yb, yhi+yb);
    for (int i= 0; i < nx; i++) {
      for (int j= 0; j < ny; j++) {
        hx->SetBinContent (i+1, j+1, h->GetBinContent (i, j));
        hx->SetBinError   (i+1, j+1, h->GetBinError   (i, j));
      }
    }
    return hx;
  } else if (dynamic_cast<const TH2D*>(h)) {
    TH2D* hx= dynamic_cast<TH2D*>(h->Clone());
    // clear under/overflows
    for (int i= 0; i <= nx+1; i++) {
      hx->SetBinContent (i, 0,    0.0);
      hx->SetBinContent (i, ny+1, 0.0);
    }
    for (int i= 1; i <= ny;   i++) {
      hx->SetBinContent (0,    i, 0.0);
      hx->SetBinContent (nx+1, i, 0.0);
    }
    return hx;
  } else {
    double xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax();
    double ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax();
    TH2D* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo, xhi, ny, ylo, yhi);
    for (int i= 0; i < nx+2; i++) {
      for (int j= 0; j < ny+2; j++) {
        hx->SetBinContent (i, j, h->GetBinContent (i, j));
        hx->SetBinError   (i, j, h->GetBinError   (i, j));
      }
    }
    return hx;
  }
}

TVectorD*
RooUnfoldResponse::H2V  (const TH1* h, int nb, bool overflow)
{
  // Returns TVectorD of the bin contents of the input histogram
  if (overflow) nb += 2;
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (int i= 0; i < nb; i++) {
    (*v)(i)= GetBinContent (h, i, overflow);
  }
  return v;
}

void
RooUnfoldResponse::V2H (const TVectorD& v, TH1* h, int nb, bool overflow)
{
  // Sets the bin content of the histogram as that element of the input vector
  h->Reset();  // in particular, ensure under/overflows are reset
  if (overflow) nb += 2;
  for (int i= 0; i < nb; i++) {
    int j= GetBin (h, i, overflow);
    h->SetBinContent (j, v(i));
  }
}

TVectorD*
RooUnfoldResponse::H2VE (const TH1* h, int nb, bool overflow)
{
  // Returns TVectorD of bin errors for input histogram
  if (overflow) nb += 2;
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (int i= 0; i < nb; i++) {
    (*v)(i)= GetBinError (h, i, overflow);
  }
  return v;
}

TMatrixD*
RooUnfoldResponse::H2M  (const TH2* h, int nx, int ny, const TH1* norm, bool overflow)
{
  // Returns Matrix of values of bins in a 2D input histogram
  int first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (int j= 0; j < ny; j++) {
    double nTrue= norm ? GetBinContent (norm, j, overflow) : 1.0;
    if (nTrue == 0.0) {
      for (int i= 0; i < nx; i++) {
        (*m)(i,j)= 0.0;
      }
    } else {
      for (int i= 0; i < nx; i++) {
        (*m)(i,j)= h->GetBinContent(i+first,j+first) / nTrue;
      }
    }
  }
  return m;
}

TMatrixD*
RooUnfoldResponse::H2ME (const TH2* h, int nx, int ny, const TH1* norm, bool overflow)
{
  // Returns matrix of bin errors for a 2D histogram.
  int first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (int j= 0; j < ny; j++) {
    double nTrue= norm ? GetBinContent (norm, j, overflow) : 1.0;
    if (nTrue == 0.0) {
      for (int i= 0; i < nx; i++) {
        (*m)(i,j)= 0.0;
      }
    } else {
      for (int i= 0; i < nx; i++) {
        // Assume Poisson nTrue, Multinomial P(mes|tru)
        (*m)(i,j)= h->GetBinError(i+first,j+first) / nTrue;
      }
    }
  }
  return m;
}

TH1*
RooUnfoldResponse::ApplyToTruth (const TH1* truth, const char* name) const
{
  // Apply the response matrix to the truth
  // Errors not set, since we assume original truth has no errors
  if (!Htruth()) return 0;  // Needed for checking binning if nothing else

  // If no truth histogram input, use training truth
  // If truth histogram input, make sure its binning is correct
  TVectorD* resultvect;
  if (truth) {
    if (truth->GetNbinsX() != _tru->GetNbinsX() ||
        truth->GetNbinsY() != _tru->GetNbinsY() ||
        truth->GetNbinsZ() != _tru->GetNbinsZ())
      cerr << "Warning: RooUnfoldResponse::ApplyToTruth truth histogram is a different size ("
           << (truth->GetNbinsX() * truth->GetNbinsY() * truth->GetNbinsZ()) << " bins) or shape from response matrix truth ("
           << ( _tru->GetNbinsX() *  _tru->GetNbinsY() *  _tru->GetNbinsZ()) << " bins)" << endl;
    resultvect= H2V (truth, GetNbinsTruth(), _overflow);
    if (!resultvect) return 0;
  } else {
    resultvect= new TVectorD (Vtruth());
  }

  (*resultvect) *= Mresponse();   // v= A*v

  // Turn results vector into properly binned histogram
  TH1* result= (TH1*) Hmeasured()->Clone (name);
  result->SetTitle (name);
  V2H (*resultvect, result, GetNbinsMeasured(), _overflow);
  delete resultvect;
  return result;
}

void
RooUnfoldResponse::SetNameTitleDefault (const char* defname, const char* deftitle)
{
  // Set object name and title
  const char* s= GetName();
  if (s[0] == '\0') {
    if (_res) s= _res->GetName();
    if (s[0] == '\0') {
      if (defname) SetName (defname);
      else if (_mes && _tru) {
        TString n= _mes->GetName();
        if (n.Length()) n.Append ("_");
        n.Append (_tru->GetName());
        if (!n.Length()) n= "response";
        SetName (n);
      }
    } else
      SetName (s);
  }
  s= GetTitle();
  if (s[0] == '\0') {
    if (_res) s= _res->GetTitle();
    if (s[0] == '\0') {
      if (deftitle) SetTitle (deftitle);
      else if (_mes && _tru) {
        TString n= _tru->GetTitle();
        if (n.Length()) n.Append (" #rightarrow ");
        n.Append (_mes->GetTitle());
        if (n.Length())
          n.Prepend ("Response ");
        else
          n= "Response";
        SetTitle (n);
      }
    } else
      SetTitle (s);
  }
}

void
RooUnfoldResponse::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    bool oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldResponse::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldResponse::Class()->WriteBuffer (R__b, this);
  }
}
