#ifndef _PETRA_RDP_DENSEMATRIX_H_
#define _PETRA_RDP_DENSEMATRIX_H_

//! Petra_RDP_DenseMatrix: A class for constructing and using general dense matrices.

/*! The Petra_RDP_DenseMatrix class enables the construction and use of real-valued, general, 
    double-precision dense matrices.  It is built on the BLAS and LAPACK and derives from the Petra_BLAS and 
    Petra_LAPACK classes. 

The Petra_RDP_DenseMatrix class is intended to provide full-featured support for solving linear and eigensystem
problems for general dense rectangular (or square) matrices.  It is written on top of BLAS and LAPACK and thus has excellent
performance and numerical capabilities.  Using this class, one can either perform simple factorizations and solves or
apply all the tricks available in LAPACK to get the best possible solution for very ill-conditioned problems.

<b>Petra_RDP_DenseMatrix vs. Petra_LAPACK</b>

The Petra_LAPACK class provides access to most of the same functionality as Petra_RDP_DenseMatrix.
The primary difference is that Petra_LAPACK is a "thin" layer on top of LAPACK and Petra_RDP_DenseMatrix
attempts to provide easy access to the more sophisticated aspects of solving dense linear and eigensystems.
<ul>
<li> When you should use Petra_LAPACK:  If you are simply looking for a convenient wrapper around the Fortran LAPACK
     routines and you have a well-conditioned problem, you should probably use Petra_LAPACK directly.
<li> When you should use Petra_RDP_DenseMatrix: If you want to (or potentially want to) solve ill-conditioned 
     problems or want to work with a more object-oriented interface, you should probably use Petra_RDP_DenseMatrix.
     
</ul>

<b>Constructing Petra_RDP_DenseMatrix Objects</b>

There are three Petra_RDP_DenseMatrix constructors.  The first constructs a zero-sized object which should be made
to appropriate length using the Shape() or Reshape() functions and then filled with the [] or () operators. 
The second is a constructor that accepts user
data as a 2D array, the third is a copy constructor. The second constructor has
two data access modes (specified by the Petra_DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

<b>Setting vectors used for linear solves</b>

Setting the X and B vectors (which are Petra_RDP_DenseMatrix objects) used for solving linear systems is 
done separately from the constructor.  
This allows a single matrix factor to be used for multiple solves.  Similar to the constructor, 
the vectors X and B can be copied or viewed using the Petra_DataAccess argument.

<b>Extracting Data from Petra_RDP_DenseMatrix Objects</b>

Once a Petra_RDP_DenseMatrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


<b>Vector and Utility Functions</b>

Once a Petra_RDP_DenseMatrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Factorizations.
  <li> Solves.
  <li> Condition estimates.
  <li> Equilibration.
  <li> Norms.
</ul>

The final useful function is Flops().  Each Petra_RDP_DenseMatrix object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Petra_Time class, one can get accurate parallel performance
numbers.

<b>Strategies for Solving Linear Systems</b>
In many cases, linear systems can be accurately solved by simply computing the LU factorization
of the matrix and then performing a forward back solve with a given set of right hand side vectors.  However,
in some instances, the factorization may be very poorly conditioned and this simple approach may not work.  In
these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
the factorization. 

Petra_RDP_DenseMatrix will use equilibration with the factorization if, once the object
is constructed and \e before it is factored, you call the function FactorWithEquilibration(true) to force 
equilibration to be used.  If you are uncertain if equilibration should be used, you may call the function
ShouldEquilibrate() which will return true if equilibration could possibly help.  ShouldEquilibrate() uses
guidelines specified in the LAPACK User Guide, namely if SCOND < 0.1 and AMAX < Underflow or AMAX > Overflow, to 
determine if equilibration \e might be useful. 
 
Petra_RDP_DenseMatrix will use iterative refinement after a forward/back solve if you call
SolveToRefinedSolution(true).  It will also compute forward and backward error estimates if you call
EstimateSolutionErrors(true).  Access to the forward (back) error estimates is available via FERR() (BERR()).

Examples using Petra_RDP_DenseMatrix can be found in the Petra test directories.

*/
#include "Petra_Petra.h" 
#include "Petra_Flops.h"
#include "Petra_BLAS.h"
#include "Petra_LAPACK.h"


//=========================================================================
class Petra_RDP_DenseMatrix : public Petra_Flops, public Petra_BLAS, public Petra_LAPACK{

  // Give ostream << function some access to private and protected data/functions.

  friend ostream& operator << (ostream& os, const Petra_RDP_DenseMatrix& A);
  public:
  
  //! Default constructor; defines a zero size object.
  /*!
    Petra_RDP_DenseMatrix objects defined by the default constructor should be sized with the 
    Shape() or Reshape functions.  
    Values should be defined by using the [] or () operators.
   */
  Petra_RDP_DenseMatrix(void);
  
  //! Set object values from two-dimensional array.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In
           A - Pointer to an array of double precision numbers.  The first vector starts at A.
	   The second vector starts at A+LDA, the third at A+2*LDA, and so on.
    \param In
           LDA - The "Leading Dimension", or stride between vectors in memory.
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   See Detailed Description section for further discussion.
  */
  Petra_RDP_DenseMatrix(Petra_DataAccess CV, double *A, int LDA, int NumRows, int NumCols);
  
  //! Petra_RDP_DenseMatrix copy constructor.
  
  Petra_RDP_DenseMatrix(const Petra_RDP_DenseMatrix& Source);
  
  //! Set dimensions of a Petra_RDP_DenseMatrix object; init values to zero.
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   Allows user to define the dimensions of a Petra_RDP_DenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Shape(int NumRows, int NumCols);
  
  //! Reshape a Petra_RDP_DenseMatrix object.
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   Allows user to define the dimensions of a Petra_RDP_DenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new shape.  If the new shape is smaller than the original, the upper left portion
	   of the original matrix (the principal submatrix) is copied to the new matrix.

    \return Integer error code, set to 0 if successful.
  */
  int Reshape(int NumRows, int NumCols);

  //! Petra_RDP_DenseMatrix destructor.  
  virtual ~Petra_RDP_DenseMatrix ();

  //! Sets the pointers for RHS and solution vectors 
  /* B should contain the RHS vectors.  X should be space allocated for the solution.  It is possible to have the
     solution returned in B by pasing in the same object for both B and X.
     \param B In 
            Right hand side as a Petra_RDP_DenseMatrix object.
     \param X In 
            Left hand side as a Petra_RDP_DenseMatrix object.

    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
   */
  int SetVectors(const Petra_RDP_DenseMatrix & B, Petra_RDP_DenseMatrix & X);
  
  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  int OneNorm();

  //! Causes equilibration to be called just before the matrix factorization as part of the call to Factor.
  /*! This function must be called before the factorization is performed. 
   */
  void FactorWithEquilibration(bool Flag) {Equilibrate_ = Flag; return;};

  //! If Flag is true, causes all subsequent function calls to work with the transpose of \e this matrix, otherwise not.
  void SolveWithTranspose(bool Flag) {Transpose_ = Flag; if (Flag) TRANS_ = 'T'; else TRANS_ = 'N'; return;};

  //! Causes all solves to compute solution to best ability using iterative refinement.
  void SolveToRefinedSolution(bool Flag) {RefineSolution_ = Flag; return;};

  //! Computes the in-place LU factorization of the matrix using the LAPACK routine \e DGETRF.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual int Factor(void);

  //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual int Solve(void);

  //! Causes all solves to estimate the forward and backward solution error. 
  /*! Error estimates will be in the arrays FERR and BERR, resp, after the solve step is complete.
      These arrays are accessible via the FERR() and BERR() access functions.
  */
  void EstimateSolutionErrors(bool Flag) {EstimateSolutionErrors_ = Flag; return;};

  //! Inverts the \e this matrix.
  /*! Note: This function works a little differently that DPOTRI in that it fills the entire
      matrix with the inverse, independent of the UPLO specification.

    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Invert(void);

  //! Computes the scaling vector S(i) = 1/sqrt(A(i,i) of the \e this matrix.
  /*! 
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int ComputeEquilibrateScaling(void);

  //! Equilibrates the \e this matrix.
  /*! 
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Equilibrate_A(void);

  //! Equilibrates the current RHS.
  /*! 
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Equilibrate_B(void);


  //! Apply Iterative Refinement.
  /*! 
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int ApplyRefinement(void);

  //! Unscales the solution vectors if equilibration was used to solve the system.
  /*! 
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Unequilibrate_X(void);

  //! Returns the reciprocal of the 1-norm condition number of the \e this matrix.
  /*! 
    \param Value Out
           On return contains the reciprocal of the 1-norm condition number of the \e this matrix.
    
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int ReciprocalConditionEstimate(double & Value);
  //! Matrix-Matrix multiplication, \e this = Scalar*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations.

  \param In
         TransA - Operate with the transpose of A if = 'T', else no transpose if = 'N'.
  \param In
         TransB - Operate with the transpose of B if = 'T', else no transpose if = 'N'.

  \param In
         ScalarAB - Scalar to multiply with A*B.
  \param In
         A - Dense Matrix.
  \param In
         B - Dense Matrix.
  \param In
         Scalar - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.
	 
  */
  int  Multiply (char TransA, char TransB, double ScalarAB, 
                 const Petra_RDP_DenseMatrix& A, 
                 const Petra_RDP_DenseMatrix& B,
                 double Scalar );

  //! Returns true if transpose of \e this matrix has and will be used.
  bool Transpose() {return(Transpose_);};

  //! Returns true if matrix is factored (factor available via AF() and LDAF()).
  bool Factored() {return(Factored_);};

  //! Returns true if factor is equilibrated (factor available via AF() and LDAF()).
  bool A_Equilibrated() {return(A_Equilibrated_);};

  //! Returns true if RHS is equilibrated (RHS available via B() and LDB()).
  bool B_Equilibrated() {return(B_Equilibrated_);};

  //! Returns true if the LAPACK general rules for equilibration suggest you should equilibrate the system.
  virtual bool ShouldEquilibrate() {ComputeEquilibrateScaling(); return(ShouldEquilibrate_);};

  //! Returns true if forward and backward error estimated have been computed (available via FERR() and BERR()).
  bool SolutionErrorsEstimated() {return(SolutionErrorsEstimated_);};

  //! Returns true if matrix inverse has been computed (inverse available via AF() and LDAF()).
  bool Inverted() {return(Inverted_);};

  //! Returns true if the condition number of the \e this matrix has been computed (value available via ReciprocalConditionEstimate()).
  bool ReciprocalConditionEstimated() {return(ReciprocalConditionEstimated_);};

  //! Returns true if the current set of vectors has been solved.
  bool Solved() {return(Solved_);};

  //! Returns true if the current set of vectors has been refined.
  bool SolutionRefined() {return(SolutionRefined_);};

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
    double& operator () (int RowIndex, int ColIndex);

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
    const double& operator () (int RowIndex, int ColIndex) const;

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
  */
    double* operator [] (int ColIndex);

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
  */
    const double* operator [] (int ColIndex) const;
    
  //! Returns row dimension of system.
  int M()  const {return(M_);};

  //! Returns column dimension of system.
  int N()  const {return(N_);};

  //! Returns pointer to the \e this matrix.
  double * A()  const {return(A_);};

  //! Returns the leading dimension of the \e this matrix.
  int LDA()  const {return(LDA_);};

  //! Returns pointer to current RHS.
  double * B()  const {return(B_);};

  //! Returns the leading dimension of the RHS.
  int LDB()  const {return(LDB_);};

  //! Returns the number of current right hand sides and solution vectors.
  int NRHS()  const {return(NRHS_);};

  //! Returns pointer to current solution.
  double * X()  const {return(X_);};

  //! Returns the leading dimension of the solution.
  int LDX()  const {return(LDX_);};

  //! Returns pointer to the factored matrix (may be the same as A() if factorization done in place).
  double * AF()  const {return(AF_);};

  //! Returns the leading dimension of the factored matrix.
  int LDAF()  const {return(LDAF_);};

  //! Returns pointer to pivot vector (if factorization has been computed), zero otherwise.
  int * IPIV()  const {return(IPIV_);};

  //! Returns the 1-Norm of the \e this matrix (returns -1 if not yet computed).
  double ANORM()  const {return(ANORM_);};

  //! Returns the reciprocal of the condition number of the \e this matrix (returns -1 if not yet computed).
  double RCOND()  const {return(RCOND_);};

  //! Ratio of smallest to largest row scale factors for the \e this matrix (returns -1 if not yet computed).
  /*! If ROWCND() is >= 0.1 and AMAX() is not close to overflow or underflow, then equilibration is not needed.
   */
  double ROWCND()  const {return(ROWCND_);};

  //! Ratio of smallest to largest column scale factors for the \e this matrix (returns -1 if not yet computed).
  /*! If COLCND() is >= 0.1 then equilibration is not needed.
   */
  double COLCND()  const {return(COLCND_);};

  //! Returns the absolute value of the largest entry of the \e this matrix (returns -1 if not yet computed).
  double AMAX()  const {return(AMAX_);};

  //! Returns a pointer to the forward error estimates computed by LAPACK.
  double * FERR()  const {return(FERR_);};

  //! Returns a pointer to the backward error estimates computed by LAPACK.
  double * BERR()  const {return(BERR_);};

  //! Returns a pointer to the row scaling vector used for equilibration.
  double * R()  const {return(R_);};

  //! Returns a pointer to the column scale vector used for equilibration.
  double * C()  const {return(C_);};
  

 protected:

  void CopyMat(double * A, int LDA, int NumRows, int NumCols, double * B, int LDB);
  void AllocateWORK(void) {if (WORK_==0) {LWORK_ = 4*N_; WORK_ = new double[LWORK_];} return;};
  void AllocateIWORK(void) {if (IWORK_==0) IWORK_ = new int[N_]; return;};
  void DeleteArrays(void);
  bool Equilibrate_;
  bool ShouldEquilibrate_;
  bool A_Equilibrated_;
  bool B_Equilibrated_;
  bool Transpose_;
  bool Factored_;
  bool EstimateSolutionErrors_;
  bool SolutionErrorsEstimated_;
  bool Solved_;
  bool Inverted_;
  bool ReciprocalConditionEstimated_;
  bool RefineSolution_;
  bool SolutionRefined_;
  bool A_Copied_;
  bool B_Copied_;
  bool X_Copied_;

  char TRANS_;

  int M_;
  int N_;
  int Min_MN_;
  int NRHS_;
  int LDA_;
  int LDAF_;
  int LDB_;
  int LDX_;
  int INFO_;
  int LWORK_;

  int * IPIV_;
  int * IWORK_;

  double ANORM_;
  double RCOND_;
  double ROWCND_;
  double COLCND_;
  double AMAX_;

  double * A_;
  double * FERR_;
  double * BERR_;
  double * AF_;
  double * WORK_;
  double * R_;
  double * C_;

  double * B_;
  double * X_;


};

const double Overflow_ = 1.79E308; // Used to test if equilibration should be done.
const double Underflow_ = 2.23E-308;
//! << operator will work for Petra_RDP_DenseMatrix objects.
ostream& operator << (ostream& os, const Petra_RDP_DenseMatrix& A);

#endif /* _PETRA_RDP_DENSEMATRIX_H_ */
