#ifndef AMESOS_CONTROL_H
#define AMESOS_CONTROL_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Teuchos_ParameterList.hpp"
using namespace Teuchos;
/*!
 \class Amesos_Control
 
 \brief Amesos_Control: Container for some control variables.

 \author Marzio Sala, SNL 9214

 \date Last updated on 24-May-05 (Champions' League Final day)
*/

class Amesos_Control
{
public:
  //! Default constructor.
  Amesos_Control()
  {
    AddToDiag_ = 0.0;
    AddZeroToDiag_ = false;
    rcond_threshold_ = 1e-12;
    refactorize_ = false;
    MaxProcesses_ = -1;
    ScaleMethod_ = 0;
    Reindex_ = 0;
  }

  //! Default destructor.
  ~Amesos_Control() {};

  void SetControlParameters( const Teuchos::ParameterList &ParameterList )  ;

  //! Add \c this value to the diagonal.
  double AddToDiag_;


  bool refactorize_;	    // if true, and if the Symbolic and Numeric
			    // objects have already been created, then
			    // attempt to "refactorize" (factor the matrix
			    // with no changes to the pivot order since the
			    // last call the klu_btf_factor).

  //! If error is greater than \c this value, perform symbolic and 
  //!  numeric factorization with full partial pivoting 
  double rcond_threshold_;  // if we refactorize, the factorization may suffer
			    // in numeric quality.  We compute rcond =
			    // min (abs (diag (U))) / max (abs (diag (U))).
			    // If this ratio is <= rcond_threshold_, then
			    // the "refactorization" is scrapped, and we factor
			    // with full partial pivoting instead.

  int ScaleMethod_;	    // most methods (KLU, UMFPACK, Mumps, ...) can scale
			    // the input matrix prior to factorization.  This can
			    // improve pivoting, reduce fill-in, and lead to a
			    // better quality factorization.  The options are:
			    // 0: no scaling
			    // 1: use the default method for the specific package
			    // 2: use the method's 1st alternative (if it has one)
			    // 3: use the method's 2nd alternative, and so on.
                            //
                            // Amesos_Klu is, at present, the only code which implements this

  //! Adds zero to diagonal of redistributed matrix (some solvers choke on a matrix with a partly empty diag)
  bool AddZeroToDiag_; 

  //! Set the matrix property.
  /*! Matrix property can be 
    - 0 : general unsymmetric matrix;
    - 1 : SPD;
    - 2 : general symmetric matrix.
    UNUSED - See bug #2331 and bug #2332
    */
  int MatrixProperty_;        

  int MaxProcesses_;                     // default is -1 ; If positive, distribute 
                                         // problem over MaxProcesses
  //! If true, the Amesos class should reindex the matrix to standard indexing (i.e. 0-n-1)
  //! at present, only Amesos_Klu supports this option.
  bool Reindex_ ; 


};

#endif
