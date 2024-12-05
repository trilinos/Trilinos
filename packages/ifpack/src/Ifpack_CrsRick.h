/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef _IFPACK_CRSRICK_H_
#define _IFPACK_CRSRICK_H_

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ScalingType.h"
#include "Ifpack_IlukGraph.h"
#include "Epetra_CompObject.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Object.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_Vector;
class Epetra_MultiVector;

namespace Teuchos {
  class ParameterList;
}

//! Ifpack_CrsRick: A class for constructing and using an incomplete Cholesky (IC) factorization of a given Epetra_CrsMatrix.

/*! The Ifpack_CrsRick class computes a "Relaxed" IC factorization with level k fill
    of a given Epetra_CrsMatrix.  The factorization
    that is produced is a function of several parameters:
<ol>
  <li> The pattern of the matrix - All fill is derived from the original matrix nonzero structure.  Level zero fill
       is defined as the original matrix pattern (nonzero structure), even if the matrix value at an entry is stored
       as a zero. (Thus it is possible to add entries to the IC factors by adding zero entries the original matrix.)

  <li> Level of fill - Starting with the original matrix pattern as level fill of zero, the next level of fill is
       determined by analyzing the graph of the previous level and determining nonzero fill that is a result of combining
       entries that were from previous level only (not the current level).  This rule limits fill to entries that
       are direct decendents from the previous level graph.  Fill for level k is determined by applying this rule
       recursively.  For sufficiently large values of k, the fill would eventually be complete and an exact Cholesky
       factorization would be computed.  Level of fill is defined during the construction of the Ifpack_IlukGraph object.

  <li> Level of overlap - All Ifpack preconditioners work on parallel distributed memory computers by using
       the row partitioning the user input matrix to determine the partitioning for local IC factors.  If the level of
       overlap is set to zero,
       the rows of the user matrix that are stored on a given processor are treated as a self-contained local matrix
       and all column entries that reach to off-processor entries are ignored.  Setting the level of overlap to one
       tells Ifpack to increase the size of the local matrix by adding rows that are reached to by rows owned by this
       processor.  Increasing levels of overlap are defined recursively in the same way.  For sufficiently large levels
       of overlap, the entire matrix would be part of each processor's local IC factorization process.
       Level of overlap is defined during the construction of the Ifpack_IlukGraph object.

       Once the factorization is computed, applying the factorization \(LL^{T}y = x\)
       results in redundant approximations for any elements of y that correspond to
       rows that are part of more than one local IC factor.  The OverlapMode (changed by calling SetOverlapMode())
       defines how these redundancies are
       handled using the Epetra_CombineMode enum.  The default is to zero out all values of y for rows that
       were not part of the original matrix row distribution.

  <li> Fraction of relaxation - Ifpack_CrsRick computes the IC factorization row-by-row.  As entries at a given
       row are computed, some number of them will be dropped because they do match the prescribed sparsity pattern.
       The relaxation factor determines how these dropped values will be handled.  If the RelaxValue (changed by calling
       SetRelaxValue()) is zero, then these extra entries will by dropped.  This is a classical IC approach.
       If the RelaxValue is 1, then the sum
       of the extra entries will be added to the diagonal.  This is a classical Modified IC (MIC) approach.  If
       RelaxValue is between 0 and 1, then RelaxValue times the sum of extra entries will be added to the diagonal.

       For most situations, RelaxValue should be set to zero.  For certain kinds of problems, e.g., reservoir modeling,
       there is a conservation principle involved such that any operator should obey a zero row-sum property.  MIC
       was designed for these cases and you should set the RelaxValue to 1.  For other situations, setting RelaxValue to
       some nonzero value may improve the stability of factorization, and can be used if the computed IC factors
       are poorly conditioned.

  <li> Diagonal perturbation - Prior to computing the factorization, it is possible to modify the diagonal entries of the matrix
       for which the factorization will be computing.  If the absolute and relative perturbation values are zero and one,
       respectively, the
       factorization will be compute for the original user matrix A.  Otherwise, the factorization
       will computed for a matrix that differs from the original user matrix in the diagonal values only.  Below we discuss
       the details of diagonal perturbations.
       The absolute and relative threshold values are set by calling SetAbsoluteThreshold() and SetRelativeThreshold(), respectively.
</ol>

<b> Estimating Preconditioner Condition Numbers </b>

For ill-conditioned matrices, we often have difficulty computing usable incomplete
factorizations.  The most common source of problems is that the factorization may encounter a small or zero pivot,
in which case the factorization can fail, or even if the factorization
succeeds, the factors may be so poorly conditioned that use of them in
the iterative phase produces meaningless results.  Before we can fix
this problem, we must be able to detect it.  To this end, we use a
simple but effective condition number estimate for \f$(LU)^{-1}\f$.

The condition of a matrix \f$B\f$, called \f$cond_p(B)\f$, is defined as
\f$cond_p(B) = \|B\|_p\|B^{-1}\|_p\f$ in some appropriate norm \f$p\f$.  \f$cond_p(B)\f$
gives some indication of how many accurate floating point
digits can be expected from operations involving the matrix and its
inverse.  A condition number approaching the accuracy of a given
floating point number system, about 15 decimal digits in IEEE double
precision, means that any results involving \f$B\f$ or \f$B^{-1}\f$ may be
meaningless.

The \f$\infty\f$-norm of a vector \f$y\f$ is defined as the maximum of the
absolute values of the vector entries, and the \f$\infty\f$-norm of a
matrix C is defined as
\f$\|C\|_\infty = \max_{\|y\|_\infty = 1} \|Cy\|_\infty\f$.
A crude lower bound for the \f$cond_\infty(C)\f$ is
\f$\|C^{-1}e\|_\infty\f$ where \f$e = (1, 1, \ldots, 1)^T\f$.  It is a
lower bound because \f$cond_\infty(C) = \|C\|_\infty\|C^{-1}\|_\infty
\ge \|C^{-1}\|_\infty \ge |C^{-1}e\|_\infty\f$.

For our purposes, we want to estimate \f$cond_\infty(LU)\f$, where \f$L\f$ and
\f$U\f$ are our incomplete factors.  Edmond in his Ph.D. thesis demonstrates that
\f$\|(LU)^{-1}e\|_\infty\f$ provides an effective estimate for
\f$cond_\infty(LU)\f$.  Furthermore, since finding \f$z\f$ such that \f$LUz = y\f$
is a basic kernel for applying the preconditioner, computing this
estimate of \f$cond_\infty(LU)\f$ is performed by setting \f$y = e\f$, calling
the solve kernel to compute \f$z\f$ and then
computing \f$\|z\|_\infty\f$.


<b>\e A \e priori Diagonal Perturbations</b>

Given the above method to estimate the conditioning of the incomplete factors,
if we detect that our factorization is too ill-conditioned
we can improve the conditioning by perturbing the matrix diagonal and
restarting the factorization using
this more diagonally dominant matrix.  In order to apply perturbation,
prior to starting
the factorization, we compute a diagonal perturbation of our matrix
\f$A\f$ and perform the factorization on this perturbed
matrix.  The overhead cost of perturbing the diagonal is minimal since
the first step in computing the incomplete factors is to copy the
matrix \f$A\f$ into the memory space for the incomplete factors.  We
simply compute the perturbed diagonal at this point.

The actual perturbation values we use are the diagonal values \f$(d_1, d_2, \ldots, d_n)\f$
with \f$d_i = sgn(d_i)\alpha + d_i\rho\f$, \f$i=1, 2, \ldots, n\f$, where
\f$n\f$ is the matrix dimension and \f$sgn(d_i)\f$ returns
the sign of the diagonal entry.  This has the effect of
forcing the diagonal values to have minimal magnitude of \f$\alpha\f$ and
to increase each by an amount proportional to \f$\rho\f$, and still keep
the sign of the original diagonal entry.

<b>Constructing Ifpack_CrsRick objects</b>

Constructing Ifpack_CrsRick objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Ifpack_CrsRick instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
</ol>

Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
create new entries.

<b> Counting Floating Point Operations </b>

Each Ifpack_CrsRick object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this
information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Epetra_Map is required for the Ifpack_CrsRick constructor.

*/


class Ifpack_CrsRick: public Epetra_Object, public Epetra_CompObject, public virtual Epetra_Operator {

  // Give ostream << function some access to private and protected data/functions.

  friend std::ostream& operator << (std::ostream& os, const Ifpack_CrsRick& A);

 public:
  //! Ifpack_CrsRick constuctor with variable number of indices per row.
  /*! Creates a Ifpack_CrsRick object and allocates storage.

    \param In
           A - User matrix to be factored.
    \param In
           Graph - Graph generated by Ifpack_IlukGraph.
  */
  Ifpack_CrsRick(const Epetra_CrsMatrix &A, const Ifpack_IlukGraph & Graph);

  //! Copy constructor.
  Ifpack_CrsRick(const Ifpack_CrsRick & Matrix);

  //! Ifpack_CrsRick Destructor
  virtual ~Ifpack_CrsRick();

  //! Initialize L and U with values from user matrix A.
  /*! Copies values from the user's matrix into the nonzero pattern of L and U.
   */
  int InitValues();

  //! If values have been initialized, this query returns true, otherwise it returns false.
  bool ValuesInitialized() const {return(ValuesInitialized_);};

  //! Set RIC(k) relaxation parameter
  void SetRelaxValue( double RelaxValue) {RelaxValue_ = RelaxValue; return;}

  //! Set absolute threshold value
  void SetAbsoluteThreshold( double Athresh) {Athresh_ = Athresh; return;}

  //! Set relative threshold value
  void SetRelativeThreshold( double Rthresh) {Rthresh_ = Rthresh; return;}

  //! Set overlap mode type
  void SetOverlapMode( Epetra_CombineMode OverlapMode) {OverlapMode_ = OverlapMode; return;}

  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes four parameter names: relax_value,
     absolute_threshold, relative_threshold and overlap_mode. These names
     are case insensitive. The ParameterEntry must have type double for all
     except overlap_mode, which must have type Epetra_CombineMode.
  */
  int SetParameters(const Teuchos::ParameterList& parameterlist,
                    bool cerr_warning_if_unused=false);

  //! Compute IC factor L using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the RIC(k) factor L using the current:
    <ol>
    <li> Ifpack_IlukGraph specifying the structure of L.
    <li> Value for the RIC(k) relaxation parameter.
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  int Factor();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool Factored() const {return(Factored_);};


  // Mathematical functions.


  //! Returns the result of a Ifpack_CrsRick forward/back solve on a Epetra_Vector x in y.
  /*!
    \param In
    Trans -If true, solve transpose problem.
    \param In
    x -A Epetra_Vector to solve for.
    \param Out
    y -A Epetra_Vector containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_Vector& x, Epetra_Vector& y) const;

  //! Returns the result of a Ifpack_CrsRick forward/back solve on a Epetra_MultiVector X in Y.
  /*!
    \param In
    Trans -If true, solve transpose problem.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of multiplying U, D and U^T in that order on an Epetra_MultiVector X in Y.
  /*!
    \param In
    Trans -If true, multiply by L^T, D and U^T in that order.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the maximum over all the condition number estimate for each local IC set of factors.
  /*! This functions computes a local condition number estimate on each processor and return the
      maximum over all processor of the estimate.
   \param In
    Trans -If true, solve transpose problem.
    \param Out
    ConditionNumberEstimate - The maximum across all processors of
    the infinity-norm estimate of the condition number of the inverse of LDU.
  */
  int Condest(bool Trans, double & ConditionNumberEstimate) const;
  // Attribute access functions


  //! Get RIC(k) relaxation parameter
  double GetRelaxValue() {return RelaxValue_;}

  //! Get absolute threshold value
  double GetAbsoluteThreshold() {return Athresh_;}

  //! Get relative threshold value
  double GetRelativeThreshold() {return Rthresh_;}

  //! Get overlap mode type
  Epetra_CombineMode GetOverlapMode() {return OverlapMode_;}

  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {return(Graph().NumGlobalRows());};

  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {return(Graph().NumGlobalCols());};

  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {return(Graph().NumGlobalNonzeros());};

  //! Returns the number of diagonal entries found in the global input graph.
  virtual int NumGlobalDiagonals() const {return(Graph().NumGlobalDiagonals());};

  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(Graph().NumMyRows());};

  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(Graph().NumMyCols());};

  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {return(Graph().NumMyNonzeros());};

  //! Returns the number of diagonal entries found in the local input graph.
  virtual int NumMyDiagonals() const {return(Graph().NumMyDiagonals());};

  //! Returns the index base for row and column indices for this graph.
  int IndexBase() const {return(Graph().IndexBase());};

  //! Returns the address of the Ifpack_IlukGraph associated with this factored matrix.
  const Ifpack_IlukGraph & Graph() const {return(Graph_);};

  //! Returns the address of the D factor associated with this factored matrix.
  const Epetra_Vector & D() const {return(*D_);};

  //! Returns the address of the U factor associated with this factored matrix.
  const Epetra_CrsMatrix & U() const {return(*U_);};

  //@{ \name Additional methods required to support the Epetra_Operator interface.

    //! Returns a character string describing the operator
    const char * Label() const {return(Epetra_Object::Label());};

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
        does not support transpose use, this method should return a value of -1.

    \param In
           UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! Note that this implementation of Apply does NOT perform a forward back solve with
        the LDU factorization.  Instead it applies these operators via multiplication with
        U, D and L respectively.  The ApplyInverse() method performs a solve.

    \param In
           X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
           Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Multiply(Ifpack_CrsRick::UseTranspose(), X, Y));};

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
        method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(),
        the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
        if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
        be accounted for when doing a triangular solve.

    \param In
           X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
           Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Solve(Ifpack_CrsRick::UseTranspose(), X, Y));};

    //! Returns 0.0 because this class cannot compute Inf-norm.
    double NormInf() const {return(0.0);};

    //! Returns false because this class cannot compute an Inf-norm.
    bool HasNormInf() const {return(false);};

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(UseTranspose_);};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(A_.DomainMap());};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const{return(A_.RangeMap());};
  //@}

 protected:
  void SetFactored(bool Flag) {Factored_ = Flag;};
  void SetValuesInitialized(bool Flag) {ValuesInitialized_ = Flag;};
  bool Allocated() const {return(Allocated_);};
  int SetAllocated(bool Flag) {Allocated_ = Flag; return(0);};

 private:


  int Allocate();

  const Epetra_CrsMatrix &A_;
  const Ifpack_IlukGraph & Graph_;
  Epetra_CrsMatrix * U_;
  Epetra_Vector * D_;
  bool UseTranspose_;


  bool Allocated_;
  bool ValuesInitialized_;
  bool Factored_;
  double RelaxValue_;
  double Athresh_;
  double Rthresh_;
  mutable double Condest_;

  mutable Epetra_MultiVector * OverlapX_;
  mutable Epetra_MultiVector * OverlapY_;
  Epetra_CombineMode OverlapMode_;


};

//! << operator will work for Ifpack_CrsRick.
std::ostream& operator << (std::ostream& os, const Ifpack_CrsRick& A);

#endif /* _IFPACK_CRSRICK_H_ */
