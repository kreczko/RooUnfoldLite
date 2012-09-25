// Author: Kerstin Tackmann, Andreas Hoecker, Heiko Lacker

/**********************************************************************************
 *                                                                                *
 * Project: TSVDUnfold - data unfolding based on Singular Value Decomposition     *
 * Package: ROOT                                                                  *
 * Class  : TSVDUnfold                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      Single class implementation of SVD data unfolding based on:               *
 *          A. Hoecker, V. Kartvelishvili,                                        *
 *          "SVD approach to data unfolding"                                      *
 *          NIM A372, 469 (1996) [hep-ph/9509307]                                 *
 *                                                                                *
 * Authors:                                                                       *
 *      Kerstin Tackmann <Kerstin.Tackmann@cern.ch>   - CERN, Switzerland         *
 *      Andreas Hoecker  <Andreas.Hoecker@cern.ch>    - CERN, Switzerland         *
 *      Heiko Lacker     <lacker@physik.hu-berlin.de> - Humboldt U, Germany       *
 *                                                                                *
 * Copyright (c) 2010:                                                            *
 *      CERN, Switzerland                                                         *
 *      Humboldt University, Germany                                              *
 *                                                                                *
 **********************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSVDUnfold                                                           //
//                                                                      //
// Data unfolding using Singular Value Decomposition (hep-ph/9509307)   //
// Authors: Kerstin Tackmann, Andreas Hoecker, Heiko Lacker             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef TSVDUNFOLD_H
#define TSVDUNFOLD_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TMatrixD
#include "TMatrixD.h"
#endif
#ifndef ROOT_TVectorD
#include "TVectorD.h"
#endif
#ifndef ROOT_TMatrixDSym
#include "TMatrixDSym.h"
#endif

class TH1D;
class TH2D;

class TSVDUnfold : public TObject {

public:

   // Constructor
   // Initialisation of unfolding
   // "bdat" - measured data distribution (number of events)
   // "bini" - reconstructed MC distribution (number of events)
   // "xini" - truth MC distribution (number of events)
   // "Adet" - detector response matrix (number of events)
   TSVDUnfold( const TH1D* bdat, const TH1D* bini, const TH1D* xini, const TH2D* Adet );
   TSVDUnfold( const TSVDUnfold& other );

   // Destructor
   virtual ~TSVDUnfold(); 

   // Set option to normalize unfolded spectrum to unit area
   // "normalize" - switch 
   void     SetNormalize ( bool normalize ) { fNormalize = normalize; }

   // Do the unfolding
   // "kreg"   - number of singular values used (regularisation)
   TH1D*    Unfold       ( int kreg );

   // Determine for given input error matrix covariance matrix of unfolded 
   // spectrum from toy simulation
   // "cov"    - covariance matrix on the measured spectrum, to be propagated
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   TH2D*    GetUnfoldCovMatrix( const TH2D* cov, int ntoys, int seed = 1 );

   // Determine covariance matrix of unfolded spectrum from finite statistics in 
   // response matrix
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   TH2D*    GetAdetCovMatrix( int ntoys, int seed=1 );

   // Regularisation parameter
   int    GetKReg() const { return fKReg; }

   // Obtain the distribution of |d| (for determining the regularization)
   TH1D*    GetD() const;

   // Obtain the distribution of singular values
   TH1D*    GetSV() const;
   
   // Helper functions
   static double ComputeChiSquared( const TH1D& truspec, const TH1D& unfspec, const TH2D& covmat, double regpar = 0.01 );

private: 
   
   // Helper functions for vector and matrix operations
   void            FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC ) const;
   static double GetCurvature       ( const TVectorD& vec, const TMatrixD& curv );

   void            InitHistos  ( );

   // Helper functions
   static void     H2V      ( const TH1D* histo, TVectorD& vec   );
   static void     H2Verr   ( const TH1D* histo, TVectorD& vec   );
   static void     V2H      ( const TVectorD& vec, TH1D& histo   );
   static void     H2M      ( const TH2D* histo, TMatrixD& mat   );
   static TMatrixD MatDivVec( const TMatrixD& mat, const TVectorD& vec, int zero=0 );
   static TVectorD CompProd ( const TVectorD& vec1, const TVectorD& vec2 );

   static TVectorD VecDiv                 ( const TVectorD& vec1, const TVectorD& vec2, int zero = 0 );
   static void     RegularisedSymMatInvert( TMatrixDSym& mat, double eps = 1e-3 );
   
   // Class members
   int       fNdim;        //! Truth and reconstructed dimensions
   int       fDdim;        //! Derivative for curvature matrix
   bool      fNormalize;   //! Normalize unfolded spectrum to 1
   int       fKReg;        //! Regularisation parameter
   TH1D*       fDHist;       // Distribution of d (for checking regularization)
   TH1D*       fSVHist;      // Distribution of singular values

   // Input histos
   const TH1D* fBdat;        // measured distribution (data)
   const TH1D* fBini;        // reconstructed distribution (MC)
   const TH1D* fXini;        // truth distribution (MC)
   const TH2D* fAdet;        // Detector response matrix

   // Evaluation of covariance matrices
   TH1D*       fToyhisto;    //! Toy MC histogram
   TH2D*       fToymat;      //! Toy MC detector response matrix
   bool      fToyMode;     //! Internal switch for covariance matrix propagation
   bool      fMatToyMode;  //! Internal switch for evaluation of statistical uncertainties from response matrix

   
   ClassDef( TSVDUnfold, 0 ) // Data unfolding using Singular Value Decomposition (hep-ph/9509307)   
};

#endif
