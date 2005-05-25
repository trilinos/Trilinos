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
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"

#include "superlu_ddefs.h"
#include "supermatrix.h"
//  SuperLU defines Reduce to be a macro in util.h, this conflicts with Reduce() in Epetra_MultiVector.h
#undef Reduce

//! Amesos_Superludist:  An object-oriented wrapper for Superludist.
/*!  Amesos_Superludist will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Superludist solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.


  Superludist execution can be tuned through a variety of parameters.
  Amesos_Superludist.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.
    
*/
class Amesos_Superludist: public Amesos_BaseSolver,
                          private Amesos_Time,
                          private Amesos_NoCopiable,
                          private Amesos_Utils,
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
  */
  ~Amesos_Superludist(void);
  //@}

  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name Atribute set methods
  //!  Amesos_Superludist does not support transpose at this time.
  /*!  returns 0 if UseTranspose is set to false, else 1 (failure)
   */
  int SetUseTranspose(bool UseTranspose) {
    AMESOS_CHK_ERR(-1);
  }

  //@}
  
  //@{ \name Atribute access functions

  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if SUPERLUDIST can handle this matrix shape 
  /*! Returns true if the matrix shape is one that SUPERLUDIST can
    handle. SUPERLUDIST only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  bool UseTranspose() const {return(UseTranspose_);};

  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};
  //@}

  int SetParameters( Teuchos::ParameterList &ParameterList ) ;

  //! Print various timig.
  void PrintTiming() const;
  
  //! Print various information about the parameters used by Superludist.
  void PrintStatus() const;
  
private:  

  int RedistributeA();

  int ReFactor();
  int Factor();

  //
  //  Parameters set by the Parameter list
  //
  bool ReuseSymbolic_; // default is false ; Allows FactOption to be used on subsequent
                       // calls to pdgssvx from NumericFactorization
  bool AddZeroToDiag_; // default is false ; Adds zero to diagonal of redistributed matrix
                       // (in case Superludist chokes on a matrix with a partly empty diag)
  fact_t FactOption_;  // default is SamePattern_SameRow
  bool Redistribute_ ; //  default = true;  redistribute the input matrix 
                       //  prior to calling Superludist
  int MaxProcesses_;   // default is -1 ; If positive, distribute problem over
                       // MaxPricesses


  //
  //  These are used to determine what needs to be cleaned up.
  //
  const Epetra_LinearProblem * Problem_;
  int GridCreated_ ; 
  int FactorizationDone_ ; 

  int NumRows_; 
  int NumGlobalNonzeros_; 
  Epetra_Map *UniformMap_ ;    //  Uniform Map (SuperLUdist requires a linear map)
  Epetra_CrsMatrix *UniformMatrix_;  
  Epetra_Export *ExportToDist_; // Exporter from Input Matrix to UniformMatrix_

  Epetra_Import *ImportToDistributed_ ;
  Epetra_Import *ImportBackToOriginal_ ;

  //
  //  These variables are here just to keep the code a bit shorter
  //  They are set in NumericFactorization_
  //
  Epetra_RowMatrix *RowMatrixA_ ;  // Problem_->GetOperator()
  int iam_;                        // Comm_.MyPID() ;
  int NumProcs_;

  Epetra_RowMatrix *SuperluMat_ ;  // As passed to Superludist

  //
  //  Ap, Ai, Aval form the compressed row storage used by Klu
  //
  vector <int> Ap_;
  vector <int> Ai_;
  vector <double> Aval_;

  Epetra_MultiVector *vecBdistributed_; 
  Epetra_MultiVector *vecXdistributed_; 

  vector<int>ColIndicesV_;
  vector<double>RowValuesV_;
  vector <int>Global_Columns_; 


  //
  //  Here are the structures used by Superlu
  //
  SuperMatrix SuperluA_;
  ScalePermstruct_t ScalePermstruct_;
  LUstruct_t LUstruct_;
  SOLVEstruct_t SOLVEstruct_; 

  int nprow_;
  int npcol_;
  gridinfo_t grid_;                 // SuperLU's grid information
  superlu_options_t options_;

  bool PrintNonzeros_;
  string ColPerm_;
  string RowPerm_;
  int* perm_c_;
  int* perm_r_;
  string IterRefine_;
  bool ReplaceTinyPivot_;
  bool Equil_;
  
  bool FactorizationOK_ ;           // True if the matrix factorization has
                                    // been performed more recently than the
                                    // latest call to SymbolicFactorization()

  bool UseTranspose_;      // Set by SetUseTranpose() 
  
  //  int NumSymbolicFact_;  //  Amesos_Superludist_ does not separate Symbolic from Numeric Factorization
  int NumNumericFact_;
  int NumSolve_;
  
};  // End of  class Amesos_Superludist  
#endif /* AMESOS_SUPERLUDIST_H */
