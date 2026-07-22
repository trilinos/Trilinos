// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NOX_COARSEPROBLEM_INTERFACE_H
#define NOX_COARSEPROBLEM_INTERFACE_H

#include <iostream>

// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX headers
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"

// ML-headers
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"

// User Defined
#include "ml_nox_fineinterface.H"

// forward declarations
class Nox_Problem_Interface;
class Epetra_CrsMatrix;


namespace ML_NOX
{
class  Nox_CoarseProblem_Interface : public NOX::EpetraNew::Interface::Required,
                                     public NOX::EpetraNew::Interface::Jacobian,
                                     public NOX::EpetraNew::Interface::Preconditioner
{
public:
  Nox_CoarseProblem_Interface(ML_NOX::Ml_Nox_Fineinterface& interface,
                              int level,
                              int plevel,
                              Epetra_CrsMatrix** P,
                              const Epetra_BlockMap* this_RowMap,
                              int nlevel);

  ~Nox_CoarseProblem_Interface();

  //! recreate this interface. This is a constructor reusing previous data
  bool recreate(int plevel,Epetra_CrsMatrix** P, const Epetra_BlockMap* this_RowMap);
  
  //! Compute and return F (derived)
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

  //! Compute an explicit Jacobian (derived)
  bool computeJacobian(const Epetra_Vector& x);

  //! Computes a user supplied preconditioner based on input vector x.  
  //! Returns true if computation was successful. (derived)
  bool computePreconditioner(const Epetra_Vector& x,
		             NOX::Parameter::List* precParams = 0);
  
  // restrict a fine level vector to this coarse level
  Epetra_Vector* restrict_fine_to_this(const Epetra_Vector& xfine);
  
  // prolong a this coarse level vector to the fine level
  Epetra_Vector* prolong_this_to_fine(const Epetra_Vector& xcoarse);

  // return number of calls to the computeF-member
  int numcallscomputeF() { return nFcalls_; }
  
  // reset number of calls to the computeF-member to zero
  bool resetnumcallscomputeF() { nFcalls_=0; return true;}
  
  //! set the ptr to the prolonagation hierarchy
  bool setP(Epetra_CrsMatrix** P) { P_ = P; return true;}

  //! destroy the P-hierarchy WARNING: this may only be called once for all coarse interfaces!
  bool destroyP() 
  { 
     if (P_)
     {
        int i;
        for (i=0; i<nlevel_; i++)
           if (P_[i])
           {
              delete P_[i];
              P_[i] = 0;
           }
        delete [] (P_);
        P_ = 0;
     }
     return true;
  } // destroyP() 

  // set the FAS-modified system
  bool setModifiedSystem(bool flag, Epetra_Vector* fbar, Epetra_Vector* fxbar)
  { isFASmodfied_ = flag; fbar_ = fbar; fxbar_ = fxbar; return true; }
  
  // check the FAS flag
  bool isFAS() { return isFASmodfied_; }
  
  // restrict vector to next coarser level
  Epetra_Vector* restrict_to_next_coarser_level(Epetra_Vector* thisvec, int current, int next); 

  // prolong vector from next coarser level to current level
  Epetra_Vector* prolong_to_this_level(Epetra_Vector* coarsevec, int current, int next);                                  

  // apply constraints to the gradient
  void ApplyAllConstraints(Epetra_Vector& gradient);
  
  // return the level of this interface
  int Level() { return level_; }
  
  // return the rowmap of this level
  const Epetra_BlockMap* BlockMap() { return this_RowMap_; }
  
private:  
  // don't want copy-ctor and = operator
  Nox_CoarseProblem_Interface(ML_NOX::Nox_CoarseProblem_Interface& old);
  Nox_CoarseProblem_Interface operator = (const ML_NOX::Nox_CoarseProblem_Interface& old);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 
  ML_NOX::Ml_Nox_Fineinterface&    fineinterface_;
  int                              level_;         // the level this interface belongs to
  int                              nlevel_;        // total number of levels
  int                              ml_printlevel_; // ML-printlevel
  int                              nFcalls_;       // number of calls to computeF
  // NOTE: this class is NOT in charge of destroying P_!
  Epetra_CrsMatrix**               P_;             // ptr to prolongators extracted from the ml-hierarchy 
  const Epetra_BlockMap*           this_RowMap_;   // the RowMap of this levels operator
#if 1
  vector<Epetra_Vector*>           wvec_;          // working vectors for prolongation/restriction
#endif
  // NOTE: this class in not in charge of destroying fbar_, fxbar_ !
  bool                             isFASmodfied_;  // set the flag, that computeFshall return modfied answer 
  Epetra_Vector*                   fbar_;          // vector used to modify the answer when solving FAS-coarse grid
  Epetra_Vector*                   fxbar_;         // vector used to modify the answer when solving FAS-coarse grid
};
} // namespace ML_NOX


#endif // NOX_COARSEPROBLEM_INTERFACE_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
