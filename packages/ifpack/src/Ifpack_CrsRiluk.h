#ifndef _IFPACK_CrsRiluk_H_
#define _IFPACK_CrsRiluk_H_

#include "Ifpack_ScalingType.h"
#include "Ifpack_IlukGraph.h"
#include "Epetra_CompObject.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_Vector;
class Epetra_MultiVector;

//! Ifpack_CrsRiluk: A class for constructing and using an incomplete lower/upper (ILU) factorization of a given Epetra_RowMatrix.

/*! The Ifpack_CrsRiluk class computes a "Relaxed" ILU factorization with level k fill 
    of a given Epetra_RowMatrix.  The factorization 
    that is produced is a function of several parameters:
<ol>
  <li> The pattern of the matrix - All fill is derived from the original matrix nonzero structure.  Level zero fill
       is defined as the original matrix pattern (nonzero structure), even if the matrix value at an entry is stored
       as a zero. (Thus it is possible to add entries to the ILU factors by adding zero entries the original matrix.)
  <li> Level of fill - Starting with the original matrix pattern as level fill of zero, the next level of fill is
       determined by analyzing the graph of the previous level and determining nonzero fill that is a result of combining
       entries that were from previous level only (not the current level).  This rule limits fill to entries that
       are direct decendents from the previous level graph.  Fill for level k is determined by applying this rule
       recursively.  For sufficiently large values of k, the fill would eventually be complete and an exact LU
       factorization would be computed.
  <li> Level of overlap - All Ifpack preconditioners work on parallel distributed memory computers by using
       the row partitioning the user input matrix to determine the partitioning for local ILU factors.  If the level of
       overlap is set to zero,
       the rows of the user matrix that are stored on a given processor are treated as a self-contained local matrix
       and all column entries that reach to off-processor entries are ignored.  Setting the level of overlap to one
       tells Ifpack to increase the size of the local matrix by adding rows that are reached to by rows owned by this
       processor.  Increasing levels of overlap are defined recursively in the same way.  For sufficiently large levels
       of overlap, the entire matrix would be part of each processor's local ILU factorization process.

       Once the factorization is computed, applying the factorization \(LUy = x\) 
       results in redundant approximations for any elements of y that correspond to 
       rows that are part of more than one local ILU factor.  The OverlapMode (changed by calling SetOverlapMode())
       defines how these redundancies are
       handled using the Epetra_CombineMode enum.  The default is to zero out all values of y for rows that
       were not part of the original matrix row distribution.

  <li> Degree of relaxation - Ifpack_CrsRiluk computes the ILU factorization row-by-row.  As entries at a given
       row are computed, some number of them will be dropped because they do match the prescribed sparsity pattern.
       The relaxation factor determines how these dropped values will be handled.  If the RelaxValue (changed by calling
       SetRelaxValue()) is zero, then these extra entries will by dropped.  This is a classical ILU approach.
       If the RelaxValue is 1, then the sum
       of the extra entries will be added to the diagonal.  This is a classical Modified ILU (MILU) approach.  If
       RelaxValue is between 0 and 1, then RelaxValue times the sum of extra entries will be added to the diagonal.

       For most situations, RelaxValue should be set to zero.  For certain kinds of problems, e.g., reservoir modeling,
       there is a conservation principle involved such that any operator should obey a zero row-sum property.  MILU 
       was designed for these cases and you should set the RelaxValue to 1.  For other situations, setting RelaxValue to
       some nonzero value may improve the stability of factorization, and can be used if the computed ILU factors
       are poorly conditioned.
  <li> Shift values - Prior to computing the factorization, it is possible to modify the diagonal entries of the matrix
       for which the factorization will be computing.  If the absolute and relative shift values are both zero, the
       factorization will be compute for the original user matrix A.  If these shift values are nonzero, the factorization
       will computed for a matrix that differs from the original user matrix in the diagonal values only as follows

       We often have difficulty computing usable incomplete
       factorizations for our problems.  The most common source of problems
       is that the factorization may encounter a small or zero pivot,
       in which case the factorization can fail, or even if the factorization
       succeeds, the factors may be so poorly conditioned that use of them in
       the iterative phase produces meaningless results.  Before we can fix
       this problem, we must be able to detect it.  To this end, we use a
       simple but effective condition number estimate for $(LU)^{-1}$.

       <b>Estimating Preconditioner Condition Numbers</b>

       The condition of a matrix $B$, called $cond_p(B)$, is defined as
       $cond_p(B)
       = \|B\|_p\|B^{-1}\|_p$ in some appropriate norm $p$.  $cond_p(B)$
       gives some indication of how many accurate floating point
       digits can be expected from operations involving the matrix and its
       inverse.  A condition number approaching the accuracy of a given
       floating point number system, about 15 decimal digits in IEEE double
       precision, means that any results involving $B$ or $B^{-1}$ may be
       meaningless.
       
       The $\infty$-norm of a vector $y$ is defined as the maximum of the
       absolute values of the vector entries, and the $\infty$-norm of a
       matrix C is defined as
       $\|C\|_\infty = \max_{\|y\|_\infty = 1} \|Cy\|_\infty$.
       A crude lower bound for the $cond_\infty(C)$ is
       $\|C^{-1}e\|_\infty$ where $e = (1, 1, \ldots, 1)^T$.  It is a
       lower bound because $cond_\infty(C) = \|C\|_\infty\|C^{-1}\|_\infty
       \ge \|C^{-1}\|_\infty \ge |C^{-1}e\|_\infty$.
       
       For our purposes, we want to estimate $cond_\infty(LU)$, where $L$ and
       $U$ are our incomplete factors.  Chow~\cite{Chow:97} demonstrates that
       $\|(LU)^{-1}e\|_\infty$ provides an effective estimate for
       $cond_\infty(LU)$.  Furthermore, since finding $z$ such that $LUz = y$
       is a basic kernel for applying the preconditioner, computing this
       estimate of $cond_\infty(LU)$ is performed by setting $y = e$, calling
       the solve kernel to compute $z$ and then
       computing $\|z\|_\infty$.  The condition number estimates reported
       in Section~\ref{s:test_problems} are obtained using this approach.
       

       <b>{\it A priori} Diagonal Perturbations</b>

       Given the above method to estimate the conditioning of the incomplete factors,
       if we detect that our factorization is too ill-conditioned
       we can improve the conditioning by perturbing the matrix diagonal and
       restarting the factorization using
       this more diagonally dominant matrix.  In order to apply perturbation,
       prior to starting
       the factorization, we compute a diagonal perturbation of our matrix
       $A$ in Eq.~\ref{e:axb} and perform the factorization on this perturbed
       matrix.  The overhead cost of perturbing the diagonal is minimal since
       the first step in computing the incomplete factors is to copy the
       matrix $A$ into the memory space for the incomplete factors.  We
       simply compute the perturbed diagonal at this point.  The actual
       perturbation values we use are discussed below.
       
       we replace the diagonal values $(d_1, d_2, \ldots, d_n)$
       with $d_i = \sign(d_i)\alpha + d_i\rho$, $i=1, 2, \ldots, n$, where
       $n$ is the matrix dimension and $\sign(d_i)$ returns
       the sign of the diagonal entry.  This has the effect of
       forcing the diagonal values to have minimal magnitude of $\alpha$ and
       to increase each by an amount proportional to $\rho$, and still keep
       the sign of the original diagonal entry.

<b>Constructing Ifpack_CrsRiluk objects</b>

Constructing Ifpack_CrsRiluk objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Ifpack_CrsRiluk instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
</ol>

Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
create new entries.

<b> Counting Floating Point Operations </b>

Each Ifpack_CrsRiluk object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Epetra_Map is required for the Ifpack_CrsRiluk constructor.

*/    


class Ifpack_CrsRiluk: public Epetra_CompObject {
      
  // Give ostream << function some access to private and protected data/functions.

  friend ostream& operator << (ostream& os, const Ifpack_CrsRiluk& A);

 public:
  //! Ifpack_CrsRiluk constuctor with variable number of indices per row.
  /*! Creates a Ifpack_CrsRiluk object and allocates storage.  
    
    \param In 
           A - User matrix to be factored.
    \param In
           Graph - Graph generated by Ifpack_IlukGraph.
  */
  Ifpack_CrsRiluk(const Epetra_CrsMatrix &A, const Ifpack_IlukGraph & Graph);
  
  //! Copy constructor.
  Ifpack_CrsRiluk(const Ifpack_CrsRiluk & Matrix);

  //! Ifpack_CrsRiluk Destructor
  virtual ~Ifpack_CrsRiluk();

  //! Initialize L and U with values from user matrix A.
  /*! Copies values from the user's matrix into the nonzero pattern of L and U.
   */
  int InitValues();

  //! If values have been initialized, this query returns true, otherwise it returns false.
  bool ValuesInitialized() const {return(ValuesInitialized_);};

  //! Set RILU(k) relaxation parameter
  void SetRelaxValue( double RelaxValue) {RelaxValue_ = RelaxValue; return;}

  //! Set diagonal shift parameter
  void SetShiftValue( double ShiftValue) {ShiftValue_ = ShiftValue; return;}

  //! Set overlap mode type
  void SetOverlapMode( Epetra_CombineMode OverlapMode) {OverlapMode_ = OverlapMode; return;}

  //! Compute ILU factors L and U using the specified graph, diagonal shift and relaxation parameters.
  /*! This function computes the RILU(k) factors L and U using the current:
    <ol>
    <li> Ifpack_IlukGraph specifying the structure of L and U.
    <li> Value for the RILU(k) relaxation parameter (default is 0,
         value can be set via SetRelaxValue()).
    <li> Value for the \e a \e priori diagonal shift value (default is 0, 
         value can be set via SetShiftValue()).
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  int Factor();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool Factored() const {return(Factored_);};
  

  // Mathematical functions.
  
  
  //! Returns the result of a Ifpack_CrsRiluk forward/back solve on a Epetra_Vector x in y.
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    x -A Epetra_Vector to solve for.
    \param Out
    y -A Epetra_Vector containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_Vector& x, Epetra_Vector& y);
  
  //! Returns the result of a Ifpack_CrsRiluk forward/back solve on a Epetra_MultiVector X in Y.
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    NumVectors -Number of vectors in X and Y.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y);
  
  // Atribute access functions
  
    
  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {return(NumGlobalRows_);};
  
  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {return(NumGlobalCols_);};
  
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {return(NumGlobalNonzeros_);};
  
  //! Returns the number of diagonal entries found in the global input graph.
  virtual int NumGlobalDiagonals() const {return(NumGlobalDiagonals_);};
  
  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(NumMyRows_);};
  
  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(NumMyCols_);};
  
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {return(NumMyNonzeros_);};
  
  //! Returns the number of diagonal entries found in the local input graph.
  virtual int NumMyDiagonals() const {return(NumMyDiagonals_);};
  
  //! Returns the index base for row and column indices for this graph.
  int IndexBase() const {return(IndexBase_);};
  
  //! Returns the address of the Ifpack_IlukGraph associated with this factored matrix.
  const Ifpack_IlukGraph & Graph() const {return(Graph_);};
  
  //! Returns the address of the L factor associated with this factored matrix.
  const Epetra_CrsMatrix & L() const {return(*L_);};
    
  //! Returns the address of the D factor associated with this factored matrix.
  const Epetra_Vector & D() const {return(*D_);};
    
  //! Returns the address of the L factor associated with this factored matrix.
  const Epetra_CrsMatrix & U() const {return(*U_);};

 protected:
  void SetFactored(bool Flag) {Factored_ = Flag;};
  void SetValuesInitialized(bool Flag) {ValuesInitialized_ = Flag;};
  bool Allocated() const {return(Allocated_);};
  int SetAllocated(bool Flag) {Allocated_ = Flag; return(0);};
  
 private:
  
  
  int Allocate();
    
  const Epetra_CrsMatrix &A_;
  const Ifpack_IlukGraph & Graph_;
  Epetra_CrsMatrix * L_;
  Epetra_CrsMatrix * U_;
  Epetra_Vector * D_;

  
  bool Allocated_;
  bool ValuesInitialized_;
  bool Factored_;
  double RelaxValue_;
  double ShiftValue_;
  
  int IndexBase_;

  int NumGlobalRows_;
  int NumGlobalCols_;
  int NumGlobalDiagonals_;
  int NumGlobalNonzeros_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;
  int NumMyNonzeros_;
  Epetra_MultiVector * OverlapX_;
  Epetra_MultiVector * OverlapY_;
  Epetra_CombineMode OverlapMode_;

#ifdef PETRA_LEVELSCHEDULING  
  
 public:
  //! Returns the result of a Ifpack_CrsRiluk forward/back solve on a Epetra_Vector x in y.
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    x -A Epetra_Vector to solve for.
    \param Out
    y -A Epetra_Vector containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int LevelSolve(bool Trans, const Epetra_Vector& x, Epetra_Vector& y);
  
  //! Returns the result of a Ifpack_CrsRiluk forward/back solve on a Epetra_MultiVector X in Y.
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    NumVectors -Number of vectors in X and Y.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int LevelSolve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y);

#endif

};

//! << operator will work for Ifpack_CrsRiluk.
ostream& operator << (ostream& os, const Ifpack_CrsRiluk& A);

#endif /* _IFPACK_CrsRiluk_H_ */
