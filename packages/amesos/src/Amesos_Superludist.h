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

#ifndef _AMESOS_SUPERLUDIST_H_
#define _AMESOS_SUPERLUDIST_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
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
class Amesos_Superludist: public Amesos_BaseSolver { 

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

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      </ul>

      postconditions:<ul>
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      </ul>
      
      bugs:<ul>
      <li>Construction and destruction of an Amesos_Superludist object leaks 24 bytes 
      (This happens in superlu_gridinit() but it could be that I am not calling the right
      destructor.) 
     </ul>


    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization (and symbolic
      factorization if necessary) on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      </ul>

      postconditions:<ul>
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Atribute set methods
  //!  Amesos_Superludist does not support transpose at this time.
  /*!  returns 0 if UseTranspose is set to false, else 1 (failure)
   */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return( UseTranspose_?1:0 );};

  //@}
  
  //@{ \name Atribute access functions

#if 0
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if SUPERLUDIST can handle this matrix shape 
  /*! Returns true if the matrix shape is one that SUPERLUDIST can
    handle. SUPERLUDIST only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};
  //@}

  //!  Updates internal variables. 
  /*!  
      <br \>Preconditions:<ul>
      <li>None.</li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>Internal variables controlling the factorization and solve will
      be updated and take effect on all subsequent calls to NumericFactorization() 
      and Solve().</li>
      <li>All parameters whose value are to differ from the default values must 
be included in ParameterList.  Parameters not specified in ParameterList 
revert to their default values.
      </ul>

    Amesos_Superludist accepts the following parameters:
    <ul>
    <li>"AddZeroToDiag" - boolean:false - Adds zero to the diagonal, only active if Redistribute is true 
    <li>"Redistribute" - boolean:true - Redistributes the matrix 
    </ul>



   */
  int SetParameters( const Teuchos::ParameterList &ParameterList ) = 0 ;


  //! Print various information about the parameters used by Superludist.
  int PrintStatistics();
  
 private:  

  int RedistributeA() ;

  int ReFactor() ;
  int Factor() ;

 protected:

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


  Epetra_RowMatrix *SuperluMat_ ;  // As passed to Superludist

  //
  //  Ap, Ai, Aval form the compressed row storage used by Klu
  //
  vector <int> Ap_;
  vector <int> Ai_;
  vector <double> Aval_;

  Epetra_MultiVector *vecXdistributed_; 
  Epetra_MultiVector *vecBdistributed_; 

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

  bool PrintStat_;
  string ColPerm_;
  int * perm_c_;
  string RowPerm_;
  int * perm_r_;
  string IterRefine_;
  bool ReplaceTinyPivot_;
  bool Equil_;
  
  bool FactorizationOK_ ;           // True if the matrix factorization has
                                    // been performed more recently than the
                                    // latest call to SymbolicFactorization()

  bool UseTranspose_;      // Set by SetUseTranpose() 
  const Epetra_LinearProblem * Problem_;
  bool PrintStatistics_;            // print some information in the destruction phase
                                    // (defaulted to false)
  
};  // End of  class Amesos_Superludist  
#endif /* _AMESOS_SUPERLUDIST_H_ */
