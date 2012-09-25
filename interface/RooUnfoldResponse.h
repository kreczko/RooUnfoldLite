//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldResponse.h 308 2011-10-10 18:13:54Z T.J.Adye $
//
// Description:
//      Response Matrix
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDRESPONSE_HH
#define ROOUNFOLDRESPONSE_HH

#include "TNamed.h"
#include "TMatrixD.h"
#include "TH1.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
#include "TVectorDfwd.h"
#else
class TVectorD;
#endif
class TH2;
class TH2D;
class TAxis;

class RooUnfoldResponse : public TNamed {

public:

  // Standard methods

  RooUnfoldResponse(); // default constructor
  RooUnfoldResponse (const char* name, const char* title); // named constructor
  RooUnfoldResponse (const TString& name, const TString& title); // named constructor
  RooUnfoldResponse (const RooUnfoldResponse& rhs); // copy constructor
  virtual ~RooUnfoldResponse(); // destructor
  virtual RooUnfoldResponse& operator= (const RooUnfoldResponse& rhs); // assignment operator

  // Special constructors

  RooUnfoldResponse (int nb, double xlo, double xhi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case with same binning, measured vs truth
  RooUnfoldResponse (int nm, double mlo, double mhi, int nt, double tlo, double thi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case
  RooUnfoldResponse (const TH1* measured, const TH1* truth, const char* name= 0, const char* title= 0);  // constructor - measured and truth only used for shape
  RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms

  // Set up an existing object

  virtual RooUnfoldResponse& Reset ();  // clear an existing object
  virtual RooUnfoldResponse& Setup (const RooUnfoldResponse& rhs);  // set up based on another instance
  virtual RooUnfoldResponse& Setup (int nb, double xlo, double xhi);  // set up simple 1D case with same binning, measured vs truth
  virtual RooUnfoldResponse& Setup (int nm, double mlo, double mhi, int nt, double tlo, double thi);  // set up simple 1D case
  virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth);  // set up - measured and truth only used for shape
  virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth, const TH2* response);  // set up from already-filled histograms

  // Fill with training data

  virtual int Fill (double xr, double xt, double w= 1.0);  // Fill 1D Response Matrix
  virtual int Fill (double xr, double yr, double xt, double yt, double w= 1.0);  // Fill 2D Response Matrix
  virtual int Fill (double xr, double yr, double zr, double xt, double yt, double zt, double w= 1.0);  // Fill 3D Response Matrix

          int Miss (double xt);  // Fill missed event into 1D Response Matrix
          int Miss (double xt, double w);  // Fill missed event into 1D (with weight) or 2D Response Matrix
          int Miss (double xt, double yt, double w);  // Fill missed event into 2D (with weight) or 3D Response Matrix
  virtual int Miss (double xt, double yt, double zt, double w);  // Fill missed event into 3D Response Matrix

          int Fake (double xr);  // Fill fake event into 1D Response Matrix
          int Fake (double xr, double w);  // Fill fake event into 1D (with weight) or 2D Response Matrix
          int Fake (double xr, double yr, double w);  // Fill fake event into 2D (with weight) or 3D Response Matrix
  virtual int Fake (double xr, double yr, double zr, double w);  // Fill fake event into 3D Response Matrix

  virtual void Add (const RooUnfoldResponse& rhs);

  // Accessors

  int        GetDimensionMeasured() const;   // Dimensionality of the measured distribution
  int        GetDimensionTruth()    const;   // Dimensionality of the truth distribution
  int        GetNbinsMeasured()     const;   // Total number of bins in the measured distribution
  int        GetNbinsTruth()        const;   // Total number of bins in the truth distribution

  const TH1*   Hmeasured()            const;   // Measured distribution, including fakes
  TH1*         Hmeasured();                    // Measured distribution, including fakes
  const TH1*   Hfakes()               const;   // Fakes distribution
  TH1*         Hfakes();                       // Fakes distribution
  const TH1*   Htruth()               const;   // Truth distribution, used for normalisation
  TH1*         Htruth();                       // Truth distribution, used for normalisation
  const TH2*   Hresponse()            const;   // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  TH2*         Hresponse();                    // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  TH2D*        HresponseNoOverflow()  const;   // Response matrix with under/overflow bins moved into histogram body

  const TVectorD& Vmeasured()         const;   // Measured distribution as a TVectorD
  const TVectorD& Emeasured()         const;   // Measured distribution errors as a TVectorD
  const TVectorD& Vfakes()            const;   // Fakes distribution as a TVectorD
  const TVectorD& Vtruth()            const;   // Truth distribution as a TVectorD
  const TVectorD& Etruth()            const;   // Truth distribution errors as a TVectorD
  const TMatrixD& Mresponse()         const;   // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  const TMatrixD& Eresponse()         const;   // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)

  double operator() (int r, int t) const;// Response matrix element (measured,truth)

  void   UseOverflow (bool set= kTRUE);      // Specify to use overflow bins
  bool UseOverflowStatus() const;            // Get UseOverflow setting
  double FakeEntries() const;                // Return number of bins with fakes

  static TH1D*     H2H1D(const TH1*  h, int nb);
  static TVectorD* H2V  (const TH1*  h, int nb, bool overflow= kFALSE);
  static TVectorD* H2VE (const TH1*  h, int nb, bool overflow= kFALSE);
  static TMatrixD* H2M  (const TH2* h, int nx, int ny, const TH1* norm= 0, bool overflow= kFALSE);
  static TMatrixD* H2ME (const TH2* h, int nx, int ny, const TH1* norm= 0, bool overflow= kFALSE);
  static void      V2H  (const TVectorD& v, TH1* h, int nb, bool overflow= kFALSE);
  static int   GetBin (const TH1*  h, int i, bool overflow= kFALSE);  // vector index (0..nx*ny-1) -> multi-dimensional histogram global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
  static double GetBinContent (const TH1* h, int i, bool overflow= kFALSE); // Bin content by vector index
  static double GetBinError   (const TH1* h, int i, bool overflow= kFALSE); // Bin error   by vector index

  TH1* ApplyToTruth (const TH1* truth= 0, const char* name= "AppliedResponse") const; // If argument is 0, applies itself to its own truth

private:

  virtual RooUnfoldResponse& Init();
  virtual RooUnfoldResponse& Setup();
  virtual void ClearCache();
  virtual void SetNameTitleDefault (const char* defname= 0, const char* deftitle= 0);
  virtual int Miss1D (double xt, double w= 1.0);  // Fill missed event into 1D Response Matrix (with weight)
  virtual int Miss2D (double xt, double yt, double w= 1.0);  // Fill missed event into 2D Response Matrix (with weight)
  virtual int Fake1D (double xr, double w= 1.0);  // Fill fake event into 1D Response Matrix (with weight)
  virtual int Fake2D (double xr, double yr, double w= 1.0);  // Fill fake event into 2D Response Matrix (with weight)

  static int FindBin   (const TH1* h, double x, double y);
  static int FindBin   (const TH1* h, double x, double y, double z);
  static int GetBinDim (const TH1* h, int i);
  static void ReplaceAxis(TAxis* axis, const TAxis* source);

  // instance variables

  int _mdim;     // Number of measured  dimensions
  int _tdim;     // Number of truth     dimensions
  int _nm;       // Total number of measured  bins (not counting under/overflows)
  int _nt;       // Total number of truth     bins (not counting under/overflows)
  TH1*  _mes;      // Measured histogram
  TH1*  _fak;      // Fakes    histogram
  TH1*  _tru;      // Truth    histogram
  TH2*  _res;      // Response histogram
  int _overflow; // Use histogram under/overflows if 1

  mutable TVectorD* _vMes;   //! Cached measured vector
  mutable TVectorD* _eMes;   //! Cached measured error
  mutable TVectorD* _vFak;   //! Cached fakes    vector
  mutable TVectorD* _vTru;   //! Cached truth    vector
  mutable TVectorD* _eTru;   //! Cached truth    error
  mutable TMatrixD* _mRes;   //! Cached response matrix
  mutable TMatrixD* _eRes;   //! Cached response error
  mutable bool    _cached; //! We are using cached vectors/matrices

public:

  ClassDef (RooUnfoldResponse, 1) // Respose Matrix
};

// Inline method definitions

inline
RooUnfoldResponse::RooUnfoldResponse()
  : TNamed()
{
  // RooUnfoldResponse default constructor. Use Setup() to set values.
  Init();
}

inline
RooUnfoldResponse::RooUnfoldResponse (const char*    name, const char*    title)
  : TNamed(name,title)
{
  // RooUnfoldResponse default named constructor. Use Setup() to set values.
  Init();
}

inline
RooUnfoldResponse::RooUnfoldResponse (const TString& name, const TString& title)
  : TNamed(name,title)
{
  // RooUnfoldResponse default named constructor. Use Setup() to set values.
  Init();
}

inline
RooUnfoldResponse::~RooUnfoldResponse()
{
  // RooUnfoldResponse destructor
  Reset();
}


inline
RooUnfoldResponse& RooUnfoldResponse::Setup (int nb, double xlo, double xhi)
{
  // constructor -  simple 1D case with same binning, measured vs truth
  return Setup (nb, xlo, xhi, nb, xlo, xhi);
}


inline
int RooUnfoldResponse::GetDimensionMeasured() const
{
  // Dimensionality of the measured distribution (1=1D, 2=2D, 3=3D)
  return _mdim;
}

inline
int RooUnfoldResponse::GetDimensionTruth() const
{
  // Dimensionality of the truth distribution (1=1D, 2=2D, 3=3D)
  return _tdim;
}

inline
int RooUnfoldResponse::GetNbinsMeasured() const
{
  // Total number of bins in the measured distribution
  return _nm;
}

inline
int RooUnfoldResponse::GetNbinsTruth() const
{
  // Total number of bins in the truth distribution
  return _nt;
}


inline
const TH1* RooUnfoldResponse::Hmeasured() const
{
  // Measured distribution, including fakes
  return _mes;
}


inline
TH1*         RooUnfoldResponse::Hmeasured()
{
  return _mes;
}


inline
const TH1* RooUnfoldResponse::Hfakes() const
{
  // Fakes distribution
  return _fak;
}


inline
TH1*         RooUnfoldResponse::Hfakes()
{
  return _fak;
}

inline
const TH1*   RooUnfoldResponse::Htruth() const
{
  // Truth distribution, used for normalisation
  return _tru;
}

inline
TH1*         RooUnfoldResponse::Htruth()
{
  return _tru;
}

inline
const TH2*   RooUnfoldResponse::Hresponse() const
{
  // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  return _res;
}

inline
TH2*         RooUnfoldResponse::Hresponse()
{
  return _res;
}


inline
const TVectorD& RooUnfoldResponse::Vmeasured() const
{
  // Measured distribution as a TVectorD
  if (!_vMes) _cached= (_vMes= H2V  (_mes, _nm, _overflow));
  return *_vMes;
}

inline
const TVectorD& RooUnfoldResponse::Vfakes() const
{
  // Fakes distribution as a TVectorD
  if (!_vFak) _cached= (_vFak= H2V  (_fak, _nm, _overflow));
  return *_vFak;
}

inline
const TVectorD& RooUnfoldResponse::Emeasured() const
{
  // Measured distribution errors as a TVectorD
  if (!_eMes) _cached= (_eMes= H2VE (_mes, _nm, _overflow));
  return *_eMes;
}

inline
const TVectorD& RooUnfoldResponse::Vtruth() const
{
  // Truth distribution as a TVectorD
  if (!_vTru) _cached= (_vTru= H2V  (_tru, _nt, _overflow)); 
  return *_vTru;
}

inline
const TVectorD& RooUnfoldResponse::Etruth() const
{
  // Truth distribution errors as a TVectorD
  if (!_eTru) _cached= (_eTru= H2VE (_tru, _nt, _overflow)); 
  return *_eTru;
}

inline
const TMatrixD& RooUnfoldResponse::Mresponse() const
{
  // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  if (!_mRes) _cached= (_mRes= H2M  (_res, _nm, _nt, _tru, _overflow)); 
  return *_mRes;
}

inline
const TMatrixD& RooUnfoldResponse::Eresponse() const
{
  // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)
  if (!_eRes) _cached= (_eRes= H2ME (_res, _nm, _nt, _tru, _overflow)); 
  return *_eRes;
}


inline
double RooUnfoldResponse::operator() (int r, int t) const
{
  // Response matrix element (measured,truth)
  return Mresponse()(r,t);
}

inline
int    RooUnfoldResponse::GetBin (const TH1* h, int i, bool overflow)
{
  // vector index (0..nx*ny-1) -> multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
  return (h->GetDimension()<2) ? i+(overflow ? 0 : 1) : GetBinDim(h,i);
}

inline
double RooUnfoldResponse::GetBinContent (const TH1* h, int i, bool overflow)
{
  // Bin content by vector index
  return h->GetBinContent (GetBin (h, i, overflow));
}

inline
double RooUnfoldResponse::GetBinError   (const TH1* h, int i, bool overflow)
{
  // Bin error   by vector index
  return h->GetBinError   (GetBin (h, i, overflow));
}


inline
int RooUnfoldResponse::Miss (double xt)
{
  // Fill missed event into 1D Response Matrix
  return Miss1D(xt);
}

inline
int RooUnfoldResponse::Miss (double xt, double w)
{
  // Fill missed event into 1D (with weight) or 2D Response Matrix
  return _tdim==2 ? Miss2D(xt,w) : Miss1D(xt,w);
}

inline
int RooUnfoldResponse::Miss (double xt, double yt, double w)
{
  // Fill missed event into 2D (with weight) or 3D Response Matrix
  return _tdim==3 ? Miss(xt,yt,w,1.0) : Miss2D(xt,yt,w);
}


inline
int RooUnfoldResponse::Fake (double xr)
{
  // Fill fake event into 1D Response Matrix
  return Fake1D(xr);
}

inline
int RooUnfoldResponse::Fake (double xr, double w)
{
  // Fill fake event into 1D (with weight) or 2D Response Matrix
  return _mdim==2 ? Fake2D(xr,w) : Fake1D(xr,w);
}

inline
int RooUnfoldResponse::Fake (double xr, double yr, double w)
{
  // Fill fake event into 2D (with weight) or 3D Response Matrix
  return _mdim==3 ? Fake(xr,yr,w,1.0) : Fake2D(xr,yr,w);
}


inline
void RooUnfoldResponse::UseOverflow (bool set)
{
  // Specify to use overflow bins
  _overflow= (set ? 1 : 0);
}

inline
bool RooUnfoldResponse::UseOverflowStatus() const
{
  // Get UseOverflow setting
  return _overflow;
}

inline
double RooUnfoldResponse::FakeEntries() const
{
  // Return number of fake entries
  return _fak->GetEntries();
}

#endif
