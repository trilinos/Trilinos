#ifndef AMESOS_STATUS_H
#define AMESOS_STATUS_H

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
    AddToDiag_ = 0.0;
    ComputeVectorNorms_ = false;
    ComputeTrueResidual_ = false;
    verbose_ = 0;
    NumSymbolicFact_ = 0;
    NumNumericFact_ = 0;
    NumSolve_ = 0;  
    Threshold_ = 0.0;
  }

  //! Default destructor.
  ~Amesos_Status() {};

  //! If \c true, SymbolicFactorization() has been successfully called.
  bool IsSymbolicFactorizationOK_;
  //! If \c true, NumericFactorization() has been successfully called.
  bool IsNumericFactorizationOK_;
  //! If \c true, prints timing information in the destructor.
  bool PrintTiming_;
  //! If \c true, print additional information in the destructor.
  bool PrintStatus_;
  //! Add \c this value to the diagonal.
  double AddToDiag_;
  //! If \c true, prints the norms of X and B in Solve().
  bool ComputeVectorNorms_;
  //! If \c true, computes the true residual in Solve().
  bool ComputeTrueResidual_;
  
  //! Toggles the output level.
  int verbose_;

  //! Number of symbolic factorization phases.
  int NumSymbolicFact_;
  //! Number of numeric factorization phases.
  int NumNumericFact_;
  //! Number of solves.
  int NumSolve_;  

  double Threshold_;
};

#endif
