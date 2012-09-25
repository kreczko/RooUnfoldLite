//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBayes.h 291 2011-06-02 13:40:37Z T.J.Adye $
//
// Description:
//      Bayesian unfolding. Just an interface to RooUnfoldBayesImpl.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDBAYES_HH
#define ROOUNFOLDBAYES_HH

#include "RooUnfold.h"

#include "TVectorD.h"
#include "TMatrixD.h"

class TH1;
class TH2;
class RooUnfoldResponse;

class RooUnfoldBayes: public RooUnfold {

public:

	// Standard methods

	RooUnfoldBayes(); // default constructor
	RooUnfoldBayes(const char* name, const char* title); // named constructor
	RooUnfoldBayes(const TString& name, const TString& title); // named constructor
	RooUnfoldBayes(const RooUnfoldBayes& rhs); // copy constructor
	RooUnfoldBayes& operator=(const RooUnfoldBayes& rhs); // assignment operator
	virtual RooUnfoldBayes* Clone(const char* newname = 0) const;

	// Special constructors

	RooUnfoldBayes(const RooUnfoldResponse* res, const TH1* meas, int niter = 4, bool smoothit = false,
			const char* name = 0, const char* title = 0);

	void SetIterations(int niter = 4);
	void SetSmoothing(bool smoothit = false);
	int GetIterations() const;
	int GetSmoothing() const;
	const TMatrixD& UnfoldingMatrix() const;

	virtual void SetRegParm(double parm);
	virtual double GetRegParm() const;
	virtual void Reset();
	virtual void Print(Option_t* option = "") const;

	static TMatrixD& H2M(const TH2* h, TMatrixD& m, bool overflow);

protected:
	void Assign(const RooUnfoldBayes& rhs); // implementation of assignment operator
	virtual void Unfold();
	virtual void GetCov();
	virtual void GetSettings();

	void setup();
	void train();
	void getCovariance(bool doUnfoldSystematic = kFALSE);

	void smooth(TVectorD& PbarCi) const;
	double getChi2(const TVectorD& prob1, const TVectorD& prob2, double nevents) const;

private:
	void Init();
	void CopyData(const RooUnfoldBayes& rhs);

protected:
	// instance variables
	int _niter;
	int _smoothit;

	int _nc; // number of causes  (same as _nt)
	int _ne; // number of effects (same as _nm)
	double _nbartrue; // best estimate of number of true events

	TVectorD _nEstj; // Number of measured events from Effect E_j
	TVectorD _nCi; // Number of true events from cause C_i
	TVectorD _nbarCi; // Estimated number of true events from cause C_i
	TVectorD _efficiencyCi; // efficiency for detecting cause C_i

	TMatrixD _Nji; // mapping of causes to effects
	TMatrixD _Mij; // unfolding matrix
	TMatrixD _Vij; // covariance matrix
	TMatrixD _VnEstij; // covariance matrix of effects
	TMatrixD _dnCidnEj; // error propagation matrix

public:
ClassDef (RooUnfoldBayes, 1) // Bayesian Unfolding
};

// Inline method definitions

inline RooUnfoldBayes::RooUnfoldBayes() :
		RooUnfold() {
	// Default constructor. Use Setup() to prepare for unfolding.
	Init();
}

inline RooUnfoldBayes::RooUnfoldBayes(const char* name, const char* title) :
		RooUnfold(name, title) {
	// Basic named constructor. Use Setup() to prepare for unfolding.
	Init();
}

inline RooUnfoldBayes::RooUnfoldBayes(const TString& name, const TString& title) :
		RooUnfold(name, title) {
	// Basic named constructor. Use Setup() to prepare for unfolding.
	Init();
}

inline RooUnfoldBayes& RooUnfoldBayes::operator=(const RooUnfoldBayes& rhs) {
	// Assignment operator for copying RooUnfoldBayes settings.
	Assign(rhs);
	return *this;
}

inline
void RooUnfoldBayes::SetIterations(int niter) {
	// Set regularisation parameter (number of iterations)
	_niter = niter;
}

inline
void RooUnfoldBayes::SetSmoothing(bool smoothit) {
	// Enable smoothing
	_smoothit = smoothit;
}

inline
int RooUnfoldBayes::GetIterations() const {
	// Return regularisation parameter (number of iterations)
	return _niter;
}

inline
int RooUnfoldBayes::GetSmoothing() const {
	// Return smoothing setting
	return _smoothit;
}

inline const TMatrixD& RooUnfoldBayes::UnfoldingMatrix() const {
	// Access unfolding matrix (Mij)
	return _Mij;
}

inline
void RooUnfoldBayes::SetRegParm(double parm) {
	// Set regularisation parameter (number of iterations)
	SetIterations(int(parm + 0.5));
}

inline
double RooUnfoldBayes::GetRegParm() const {
	// Return regularisation parameter (number of iterations)
	return GetIterations();
}

#endif
