/*
  TODO:
  1)  Eliminate mb_ and grid_mb_  DONE 

  1)  Allow user level control over whether to use a 1D or 2D data 
      distribution.
  2)  Figure out how to redistribute the vectors.  (Allow no more than 
      nb vectors at a time for now).  Done - I think 
  3)  Move the code from software:MyExamples/TwodMap/TwoDMap.cpp in to
      this routine (for the 2D case) Done - I think
  4)  Create the ScaLAPACK 2D dense matrix.  Done - I think 
 */

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

#ifndef _AMESOS_SCALAPACK_H_
#define _AMESOS_SCALAPACK_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"

//! Amesos_Scalapack:  A serial and parallel dense solver.  For now, we implement only the unsymmetric ScaLAPACK solver.
/*!  Amesos_Scalapack, an object-oriented wrapper for LAPACK and ScaLAPACK, 
   will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the ScaLAPACK library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

<br /><br /><p>Amesos_Scalapack can be competitive for matrices 
that are not particularly sparse.  ScaLAPACK solves matrices
for which the fill-in is roughly 10% to 20% of the matrix size 
in time comparable to that achieve by other Amesos classes.  Amesos_Scalapack 
scales well and hence its performance advantage will be largest 
when large number of processes are involved.  

<br /><br /><p>Amesos_Scalapack uses the ScaLAPACK functions PDGETRF and PDGETRS 
if more than one process is used.  If only one process is used, Amesos_ScaLAPACK uses 
the LAPACK function PDGETRF and PDGETRS.

<br /><br /><p>AmesosScaLAPACK uses full partial pivoting and will
therefore provide answers that are at least as accurate as any 
direct sparse solver.

<br /><br /><p>AmesosScalapack makes sense under the following circumstances:
<ul>
<li>There is sufficient memory to store the entrie dense matrix.  8*n^2/p bytes 
will be required on each process.  -AND- one of the following
<ul>
<li>The matrix is relatively small and dense.  Amesos_Scalapack will solve 
matrices less than 100 by 100 faster than other Amesos classes unless the matrices are
very sparse.  

<li>The matrix is relatively dense and many processes are available.  If a thousand 
processes are available, Amesos_Scalapack should be competetive with other sparse direct 
solvers even for matrices whose L and U factors contain only 5% non-zeros.  

<li>The matrix is quite dense.  Amesos_Scalapack will be well on any
matrix whose L and U factors contain 20% or more non-zeros.

<li>Execution time is less important than robustness.  Amesos_Scalapack is among the most robust 
parallel direct solvers.

</ul> 

</ul> 

<h1>Common control parameters :</h1>
Amesos_Scalapack supports the following parameters which are common to accross multiple Amesos solvers:
<ul>
<li>ParamList.set("MaxProcs", int MaximumProcessesToUse );  <br>By default, this is set to -1, which causes Amesos_Scalapack to use a heuristic to determine how many processes to use.  If set to a postive value, MaximumProcessesToUse, Amesos_Scalapack will use MaximumProcessesToUse provided that there are that many processes available.  Testing should be performed with MaximumProcessesToUse set to some value larger than one to force parallel execution.   
<li><li>ParamList.set("PrintTiming", bool  );  <br>
<li><li>ParamList.set("PrintStatus", bool );  <br>
<li><li>ParamList.set("ComputeVectorNorms", bool ); <br> 
<li><li>ParamList.set("ComputeTrueResidual", bool ); <br> 
<li><li>ParamList.set("OutputLevel", int ); <br> 
<li><li>ParamList.set("DebugLevel", int ); <br> 
<li><li>ParamList.set("ComputeTrueResidual", bool ); <br> 
</ul> 
Amesos_Scalapack supports the following parameters specific to Amesos_Scalapack.
<br>
    Teuchos::ParameterList ScalapackParams = ParameterList.sublist("Scalapack") ;
<ul>
<li>ScalapackParams.set("2D distribution", bool );  <br>By default this is set "true".  In general, 
because a two dimensional data distribution generally produces faster results.  
However, in some cases, a one dimensional data distribution may provide faster 
execution time.  The code for the one dimensional data distribution uses a different data 
redistribution algorithm and uses the transpose of the matrix internally (all of which is transparent to the user).
<li>ScalapackParams.set("grid_nb", bool );  <br>By default this is set to 32.  On some 
machines, it may be possible to improve performance by up to 10% by changing the value 
of grid_nb.  (16,24,48,64 or 128) are reasonable values to try.  For testing on small 
matrices, small values of grid_nb will (if "MaxProcs" is set to a value greater than 1) 
force the code to execute in parallel.  
</ul>
<h1>Limitations:</h1>

<p>None of the following limitations would be particularly difficult to remove.

<br /><br /><p>The present implementation limits the number of right hand sides 
to the number of rows assigned to each process.  i.e. nrhs < n/p. 

<br /><br /><p>The present implementation does not take advantage of
symmetric or symmetric positive definite matrices, although ScaLAPACK has 
separate routines to  take advantages of such matrices.

*/

class Amesos_Scalapack: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Scalapack Constructor.
  /*! Creates an Amesos_Scalapack instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Scalapack( const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Scalapack Destructor.
  /*! Completely deletes an Amesos_Scalapack object.  
  */
  ~Amesos_Scalapack(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      There is no symbolic factorization phase in ScaLAPACK, as it operates
      only on dense matrices.  Hence, Amesos_Scalapack::SymbolicFactorization()
      takes no action.

    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization 
      on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)  NOT IMPLEMENTED
      <li>The non-zero structure of the matrix should not have changed
      since the last call to SymbolicFactorization().  Irrelevant for
      Amesos_Scalapack.
      <li>The distribution of the matrix should not have changed 
      since the last call to SymbolicFactorization(). Irrelevant for
      Amesos_Scalapack.
      </ul>

      postconditions:<ul>
      <li>nprow_, npcol_, DescA_ 
      <li>DenseA will be factored
      <li>Ipiv_ contains the pivots
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> X = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)   NOT IMPLEMENTED
      <li>X and B must have the same shape (NOT CHECKED)
      <li>X and B must have fewer than nb right hand sides.  EPETRA_CHK_ERR(-2)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      the Scalapack solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

#if 0
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if SCALAPACK can handle this matrix shape 
  /*! Returns true if the matrix shape is one that SCALAPACK can
    handle. SCALAPACK only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) is more efficient in Amesos_Scalapack
  /*! 
<ul>
  <li>If SetUseTranspose() is set to true, 
    <ul>
       <li><p class="code">A<sup>T</sup> X = B</p> is computed</li>
    </ul></li>
  <li>else
    <ul>
       <li><p class="code">A X = B</p> is computed</li>
    </ul></li>
</ul>
  */  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

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

    \return Integer error code, set to 0 if successful. 
   */
  int SetParameters( Teuchos::ParameterList &ParameterList ) ;

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();
  
  //@}

 private:  
  /*
  RedistributeA - Convert matrix to a dense ScaLAPACK matrix
    Preconditions:
      Problem_ must be set 
      SetParameters() 
      ScaLAPACK1DMap and ScaLAPACK1DMatrix must either be 0 or be pointers to 
        appropriatly allocate objects.  If they are non-zero, those objects
	will be deleted and recreated.  
	
    Postconditions:
      nprow_, npcol_, DescA_ 
      RowMatrixA_ 
      ScaLAPACK1DMap_ 
      ScaLAPACK1DMatrix_ 
      ImportToScaLAPACK1D_
      ImportBackToOriginal_

   */
  int RedistributeA();

  /*
    ConvertToScalapack - Convert matirx to form expected by Scalapack: Ai, Ap, Aval
    Preconditions:
      Problem_ 
    Postconditions:
      nprow_, npcol_, 
  */
  int ConvertToScalapack();     

  /*
    PerformNumericFactorization - Call Scalapack to perform numeric factorization
    Preconditions:

    Postconditions:
      DenseA_, DescA_
  */
  int PerformNumericFactorization(); 

 protected:

  int MaxProcesses_;                     // default is -1 ; If positive, distribute 
                                         // problem over MaxProcesses

  int iam_;                              //  Process number (i.e. Comm().MyPID() 
  
  int NumGlobalElements_;                //  Number of rows and columns in the Problem_->GetOperator()
  int NumGlobalNonzeros_; 

  //
  //  The following variables are required for the ScaLAPACK interface:
  //
  int nprow_ ;                           //  number of process rows: 1 for now
  int npcol_ ;                           //  number of process columns
  int ictxt_ ;                           //  BLACS context
  int m_per_p_;                          //  Number of columns per process
  int DescA_[10];                        //  ScaLAPACK array descriptor 

  Epetra_Map *ScaLAPACK1DMap_ ;          //  Points to a 1D Map which matches a ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  Epetra_CrsMatrix *ScaLAPACK1DMatrix_ ; //  Points to a  ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  Epetra_Map *VectorMap_ ;               //  Points to a Map for vectors X and B
  vector<double> DenseA_;                //  The data in a ScaLAPACK 1D blocked format
  vector<int> Ipiv_ ;                    //  ScaLAPACK pivot information
  int NumOurRows_ ;
  int NumOurColumns_ ;
 

  bool UseTranspose_;     
  const Epetra_LinearProblem * Problem_;
  
  bool PrintTiming_;
  bool PrintStatus_;
  bool ComputeVectorNorms_;
  bool ComputeTrueResidual_;
  
  int verbose_;
  int debug_;

  // some timing internal to MUMPS
  double ConTime_;                        // time to convert to MUMPS format
  double SymTime_;                        // time for symbolic factorization
  double NumTime_;                        // time for numeric factorization
  double SolTime_;                        // time for solution
  double VecTime_;                        // time to redistribute vectors
  double MatTime_;                        // time to redistribute matrix
  
  int NumSymbolicFact_;
  int NumNumericFact_;
  int NumSolve_;  



  //
  //  Control of the data distribution
  //
  bool TwoD_distribution_;  // True if 2D data distribution is used
  int grid_nb_;             // Row and Column blocking factor (only used in 2D distribution)  
  int mypcol_;              // Process column in the ScaLAPACK2D grid
  int myprow_;              // Process row in the ScaLAPACK2D grid
  Epetra_CrsMatrix* FatOut_;//

  //
  //  Blocking factors (For both 1D and 2D data distributions)
  //
  int nb_;
  int lda_;

  Epetra_Time * Time_;
  
};  // End of  class Amesos_Scalapack  
#endif /* _AMESOS_SCALAPACK_H_ */
