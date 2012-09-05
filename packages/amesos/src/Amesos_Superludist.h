// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef AMESOS_SUPERLUDIST_H
#define AMESOS_SUPERLUDIST_H

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Amesos_Control.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_RCP.hpp"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif


// Amesos_Superludist_Pimpl contains a pointer to structures defined in 
// superlu_ddefs.h.  This prevents Amesos_Superludist.h 
// from having to include superludist.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Superlu_Pimpl ; 
#endif


//! Amesos_Superludist:  An object-oriented wrapper for Superludist.
/*!  Amesos_Superludist will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Superludist solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.
*/
class Amesos_Superludist: public Amesos_BaseSolver,
                          private Amesos_Time,
                          private Amesos_NoCopiable,
                          private Amesos_Utils,
                          private Amesos_Control,
                          private Amesos_Status 
{

public: 

  //@{ \name Constructor methods
  //! Amesos_Superludist Constructor.
  /*! Creates an Amesos_Superludist instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Superludist(const Epetra_LinearProblem& LinearProblem);

  //! Amesos_Superludist Destructor.
  /*! Completely deletes an Amesos_Superludist object.  
  */  ~Amesos_Superludist(void);
  //@}

  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();


  //@}
  //@{ \name Attribute set methods
  
  //!  Amesos_Superludist does not support transpose at this time.
  /*!  returns 0 if UseTranspose is set to false, else 1 (failure)
   */
  int SetUseTranspose(bool UseTranspose) { return( UseTranspose?1:0 );};

  //@}
  //@{ \name Attribute access functions

  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if SUPERLUDIST can handle this matrix shape 
  /*! Returns true if the matrix shape is one that SUPERLUDIST can
    handle. SUPERLUDIST only works with square matrices.  
  */
  bool MatrixShapeOK() const;

  //! Always returns false. Superludist doesn't support transpose solve
  bool UseTranspose() const {return(false);};
  //@}

  int SetParameters( Teuchos::ParameterList &ParameterList ) ;

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( Amesos_Status::NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( Amesos_Status::NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( Amesos_Status::NumSolve_ ); }

  //! Print various timig.
  void PrintTiming() const;
  
  //! Print various information about the parameters used by Superludist.
  void PrintStatus() const;

  //! Extracts timing information from the current solver and places it in the parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming(TimingParameterList); }
  
private:  
  inline const Epetra_Comm& Comm() const {return(GetProblem()->GetOperator()->Comm());};

  inline const Epetra_Import& Importer() const
  {
    return(*(Importer_.get()));
  }

  inline const Epetra_Map& UniformMap() const
  {
    return(*(UniformMap_.get()));
  }

  inline const Epetra_RowMatrix& UniformMatrix() const
  {
    return(*(UniformMatrix_.get()));
  }

  inline Epetra_CrsMatrix& CrsUniformMatrix()
  {
    return(*(CrsUniformMatrix_.get()));
  }

  //
  //  PrivateSuperluData_ contains pointers to data needed by klu whose
  //  data structures are defined by klu.h
  //
  Teuchos::RCP<Amesos_Superlu_Pimpl> PrivateSuperluData_; 

  int RedistributeA();

  int ReFactor();
  int Factor();
  
  const Epetra_LinearProblem* Problem_;
  Epetra_RowMatrix *RowMatrixA_ ;  // Problem_->GetOperator()

  RCP<Epetra_Map> UniformMap_;
  RCP<Epetra_CrsMatrix> CrsUniformMatrix_;  
  RCP<Epetra_RowMatrix> UniformMatrix_;  
  Teuchos::RCP<Epetra_Import> Importer_;

  //! Allows FactOption to be used on subsequent calls to pdgssvx from NumericFactorization
  bool ReuseSymbolic_; 
  //! redistribute the input matrix prior to calling Superludist
  bool Redistribute_ ; 

  //! \c true if the SuperLU_DIST's grid has been created (and has to be free'd)
  int GridCreated_ ; 
  int FactorizationDone_ ; 
  //! \c true if NumericFactorization() has been successfully called.
  bool FactorizationOK_ ;

  //! Global dimension of the matrix.
  int NumGlobalRows_; 

  // Ap, Ai, Aval form the compressed row storage used by SuperLU_DIST
  std::vector <int> Ap_;
  std::vector <int> Ai_;
  std::vector <double> Aval_;
  //! Contains the global ID of local columns.
  int* Global_Columns_;

  int nprow_;
  int npcol_;

  bool PrintNonzeros_;
  std::string ColPerm_;
  std::string RowPerm_;
  int* perm_c_;
  int* perm_r_;
  std::string IterRefine_;
  bool ReplaceTinyPivot_;
  bool Equil_;

  int MtxConvTime_, MtxRedistTime_, VecRedistTime_;
  int NumFactTime_, SolveTime_, OverheadTime_;
 
};  // End of  class Amesos_Superludist  
#endif /* AMESOS_SUPERLUDIST_H */
