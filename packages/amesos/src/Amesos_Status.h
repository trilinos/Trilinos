#ifndef AMESOS_STATUS_H
#define AMESOS_STATUS_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Teuchos_ParameterList.hpp"
using namespace Teuchos;



/*!
 \class Amesos_Status
 
 \brief Amesos_Status: Container for some status variables.

 \author Marzio Sala, SNL 9214

 \date Last updated on 24-May-05 (Champions' League Final day)
*/


class Amesos_Status
{
public:
  //! Default constructor.
  Amesos_Status()
  {
    IsSymbolicFactorizationOK_ = false;
    IsNumericFactorizationOK_ = false;
    PrintTiming_ = false;
    PrintStatus_ = false;
    ComputeVectorNorms_ = false;
    ComputeTrueResidual_ = false;
    verbose_ = 1;
    debug_ = 0;
    NumSymbolicFact_ = 0;
    NumNumericFact_ = 0;
    NumSolve_ = 0;  
    Threshold_ = 0.0;
    MyPID_ = 0;
    NumProcs_ = 1;
  }

  //! Default destructor.
  ~Amesos_Status() {};

  void SetStatusParameters( const Teuchos::ParameterList &ParameterList )  ;

  //! If \c true, SymbolicFactorization() has been successfully called.
  bool IsSymbolicFactorizationOK_;
  //! If \c true, NumericFactorization() has been successfully called.
  bool IsNumericFactorizationOK_;
  //! If \c true, prints timing information in the destructor.
  bool PrintTiming_;
  //! If \c true, print additional information in the destructor.
  bool PrintStatus_;
  //! If \c true, prints the norms of X and B in Solve().
  bool ComputeVectorNorms_;
  //! If \c true, computes the true residual in Solve().
  bool ComputeTrueResidual_;
  
  //! Toggles the output level.
  int verbose_;

  //! Sets the level of debug_ output
  int debug_;

  //! Number of symbolic factorization phases.
  int NumSymbolicFact_;
  //! Number of numeric factorization phases.
  int NumNumericFact_;
  //! Number of solves.
  int NumSolve_;  

  double Threshold_;

  int MyPID_;
  int NumProcs_;
};

#endif
