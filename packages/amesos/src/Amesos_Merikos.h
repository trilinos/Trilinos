/* 
Task list:
    Amesos_Merikos.h
    Amesos_Merikos.cpp
      Partition the matrix - store L as L^T?
      Build tree
      Initial data redistribution
      Change row and column ownership (pass them up to the parent)
    Amesos_Component_Solver.h
    Amesos_BTF.h 
    Amesos_BTF.cpp



  Communications issues/challenges:
  **  Redistributing the original matrix to the arrowhead form that we need, options:
      1)  Two matrices:  L^T and U
      2)  One matrix:  U | L^T 
      3)  Intermediate "fat" matrix - the way that I did Scalapack

  **  Adding the Schur complements SC_0 and SC_1 to the original 
          trailing matirx A_2, owned by the parent
      1)  Precompute the final size and shape of A_2 + SC_0 + SC_1
      2)  Perform A_2 + SC_0 + SC_1 in an empty matrix of size n by n 
          and extract the non-empty rows.
      CHALLENGES:
        A)  Only process 0/1 knows the size and map for SC_0/SC_1
        B)  It would be nice to allow SC_0 and SC_1 to be sent as soon as 
	    they are available
	C)  It would be nice to have just one copy of the matrix on the
            parent.  Hence, it would be nice to know the shape of 
	    A_2 + SC_0 + SC_1 in advance. 
        D)  An import would do the trick provided that we know both maps
            in advance.  But, neither map is available to us in advance. 
            The original map (which has SC_0 on process 0 and SC_1 on 
	    process 1) is not known 
      QUESTION:
        Should the maps be in some global address space or should they be 
	in a local address space?
	I'd like to keep them in the global address space as long as possible,
	but we can't do the import of the A_2 + SC_0 + SC_1 in a global 
	address space because that would require a map that changes at each 

  **  Redistributing the right hand side vector, b
      If we create a map that reflects the post-pivoting reality, assigning
      each row of U and each column of L to the process that owns the diagonal 
      entry, we can redistribute the right hand side vector, b, to the 
      processes where the values of b will first be used, in a single, efficient, 
      import operation.  


Observations:
1)  Although this algorithm is recursive, a non-recursive implementation 
    might be cleaner.  If it is done recursively, it should be done in place,
    i.e. any data movement of the matrix itself should have been done in 
    advance.
2)  There are two choices for the basic paradigm for parallelism.  Consider 
    a two level bisection of the matrix, yielding seven tasks or diaganol 
    blocks::  T0, T1, T01, T2, T3, T23 and T0123.  In both paradigms, 
    T0, T1, T2 and T3 would each 
    be handled by a single process.  Up until now, we have been assuming 
    that T01 would be handled by processes 0 and/or 1 while T23 would be 
    handled by processes 2 and/or 3.  The other option is to arrange the 
    tasks/diagonal blocks as follows:  T0, T1, T2, T3, T01, T23, T0123 and
    treat the last three blocks:  T01, T23 and T0123 as a single block to be
    handled by all four processes.  This second paradigm includes an
    additional synchronization, but may allow a better partitioning of 
    the remaining matrix because the resulting schur complement is now 
    known.   This improved partitioning will also improve the refactorization
    (i.e. pivotless factorization).  The second paradigm may also allow for 
    better load balancing.  For example, after using recursive minimum 
    degree bi-section (or any other scheme) to partition the matrix, one could 
    run a peephole optimization pass that would look for individuals blocks 
    that could be moved from the largest partition to a smaller one.  Finally,
    if it is clear that a given partition is going to be the slowest, it might
    make sense to shift some rows/columns off of it into the splitter just 
    for load balancing.  
3)  It seems possible that Merikos would be a cleaner code if rows
    which are shared by multiple processes are split such that each row
    resides entirely within a given process.

4)  Support for pivotless refactorization is important.
5)  There is no mention of the required row and column permutations.
6)  Amesos_Merikos only needs to support the Amesos_Component interface if
    it will call itself recursively on the subblocks.  
7)  Perhaps Amesos_Component.h should be an added interface.  Instead 
    of replacing Amesos_BaseSolver.h, maybe it should add functionality 
    to it.   
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

#ifndef _AMESOS_MERIKOS_H_
#define _AMESOS_MERIKOS_H_

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


//! Amesos_Merikos:  A parallel  divide and conquer solver.
/*!
<br /><br /><p>Merikos partitions the rows of a matrix into two or more 
disjoint submatrices.  i.e. if rows i and j are in different submatrices, 
A[i,j] == 0 == A[j,i].  Rows/columns not in any of the submatrices, i.e. the 
rows/columsn of the separator, are permuted to the bottom right.  

<br /><br /><p>Merikos factors each of the disjoint submatrices in parallel, 
(potentially by calling Amesos_Merikos() recursively), 
updating the rows and columns of the separator which belong to it and forming the schur complement of those rows and columns of the separator.  

<br /><br /><p>Merikos updates the trailing block of the matrix and then factors it.  

<br /><br /><p>Merikos is a Greek word for partial, reflecting the
fact that Amesos_Merikos uses a series of partial LU factorizations,
performed in parallel, to piece together the full LU decomposition.

*/

class Amesos_Merikos: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Merikos Constructor.
  /*! Creates an Amesos_Merikos instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

  */
  Amesos_Merikos(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Merikos Destructor.
  /*! Completely deletes an Amesos_Merikos object.  
  */
  ~Amesos_Merikos(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 

    SymbolicFactorization() takes no action in Amesos_Merikos().

    \return Integer error code, set to 0 if successful.
  */
    int RedistributeA() ;
    int ConvertToScalapack() ;
    int PerformNumericFactorization() ;
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!

<pre>
    Static pivoting (i.e. scale and permute the matrix to
      produce a zero-free diagonal and to minimize the need for
      pivoting later).
    Partition the matrix 
    Redistribute the matrix to match the partitioning
    Foreach subblock of the matrix do:
      Note:  this will happen in parallel
      Create an instance of an Amesos solver object (must  
        support the Amesos_Component interface)
      Call PartialFactorization 
      Add the Schur Complement into the trailing block of the matrix.
    Endfor
    Create an Amesos instance for the trailing block of the matrix.
    Call SymbolicFactorization on the trailing block 
    Call NumericFactorization on the trailing block 
</pre>



     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves L X = B 
    /*! 

<pre>
     | L11   0   0  |  X1      B1
     | L21 L22   0  |  X2   =  B2
     | L31 L32 L33  |  X3   =  B3

    Foreach subblock of the matrix do:
      Note:  this will happen in parallel
      Lsolve() 
        i.e. L11.Solve(X1, B1) and L22.Solve(X2, B2) 
      Update the elements of B corresponding to the seperator,
        i.e. B3 = B3 - L31 X1 - L32 X2 
    Endfor
    Perform a solve on the trailing matrix:
      i.e. L33.LSolve(X3,B3) 
</pre>    

     \return Integer error code, set to 0 if successful.
  */
    int LSolve();


    //! Solves U X = B 
    /*! 

<pre>
     | U11 U12 U13  |  X1      B1
     |   0 U22 U23  |  X2   =  B2
     |   0   0 U33  |  X3   =  B3

    Perform a solve on the trailing matrix:
      i.e. U33.USolve(X3,B3) 
    Foreach subblock of the matrix do:
      Note:  this will happen in parallel
      Update the elements of B corresponding to this block
        i.e. B2 = B2 - U23 X3 ; B1 = B1 - U13 X3 
      Usolve() 
        i.e. U11.Solve(X1, B1) and U22.Solve(X2, B2) 
    Endfor
</pre>    

     \return Integer error code, set to 0 if successful.
  */
    int USolve();

    //! Solves A X = B 
    /*! 

<pre>
     | L11\U11     U12     U13  |  X1      B1
     | L21     L22\U22     U23  |  X2   =  B2
     | L31     L32         A33  |  X3   =  B3

    Foreach subblock of the matrix do:
      Note:  this will happen in parallel
      Lsolve() 
        i.e. L11.Solve(X1, B1) and L22.Solve(X2, B2) 
      Update the elements of B corresponding to the seperator,
        i.e. B3 = B3 - L31 X1 - L32 X2 
    Endfor
    Perform a solve on the trailing matrix:
      i.e. A33.Solve(X3,B3) 

    B = X ;
    Foreach subblock of the matrix do:
      Note:  this will happen in parallel
      Update the elements of B corresponding to this block
        i.e. B2 = B2 - U23 X3 ; B1 = B1 - U13 X3 
      Usolve() 
        i.e. U11.Solve(X1, B1) and U22.Solve(X2, B2) 
    Endfor
</pre>    

     \return Integer error code, set to 0 if successful.
  */
    int Solve();



  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

    
  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if MERIKOS can handle this matrix shape 
  /*! Returns true if the matrix shape is one that MERIKOS can
    handle. MERIKOS only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose() controls whether to compute AX=B or A<sup>T</sup>X = B
  /*! 
  */  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //!  Updates internal variables. 
  /*!  
    \return Integer error code, set to 0 if successful. 
   */
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();
  
  //@}

 protected:

  bool UseTranspose_;
  const Epetra_LinearProblem * Problem_;

  Epetra_CrsMatrix *L;
  Epetra_CrsMatrix *U;

  bool PrintTiming_;
  bool PrintStatus_;
  bool ComputeVectorNorms_;
  bool ComputeTrueResidual_;
  
  int verbose_;
  int debug_;

  // some timing internal, copied from MUMPS
  double ConTime_;                        // time to convert to MERIKOS format
  double SymTime_;                        // time for symbolic factorization
  double NumTime_;                        // time for numeric factorization
  double SolTime_;                        // time for solution
  double VecTime_;                        // time to redistribute vectors
  double MatTime_;                        // time to redistribute matrix
  
  int NumSymbolicFact_;
  int NumNumericFact_;
  int NumSolve_;  

  Epetra_Time * Time_;



//
//  These allow us to use the Scalapack based Merikos code
//
  Epetra_Map *ScaLAPACK1DMap_ ;          //  Points to a 1D Map which matches a ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  Epetra_CrsMatrix *ScaLAPACK1DMatrix_ ; //  Points to a  ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  Epetra_Map *VectorMap_ ;               //  Points to a Map for vectors X and B
  vector<double> DenseA_;                //  The data in a ScaLAPACK 1D blocked format
  vector<int> Ipiv_ ;                    //  ScaLAPACK pivot information
  int NumOurRows_ ;
  int NumOurColumns_ ;


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

int iam_;
int nprow_;
int npcol_;
int NumGlobalElements_;
int m_per_p_;


  
};  // End of  class Amesos_Merikos  
#endif /* _AMESOS_MERIKOS_H_ */
