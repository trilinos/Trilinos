/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK2_CRSRILUK_DECL_HPP
#define IFPACK2_CRSRILUK_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_ScalingType.hpp"
#include "Ifpack2_IlukGraph.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"

#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace Ifpack2 {

/** \class RILUK
\brief ILU(k) (incomplete LU with fill level k) factorization of a given Tpetra::RowMatrix.
\tparam MatrixType A specialization of Tpetra::RowMatrix.

This class implements a "relaxed" incomplete ILU (ILU) factorization with level k fill.

\section Ifpack2_RILUK_Parameters Parameters

For a complete list of valid parameters, see the documentation of setParameters().

The computed factorization is a function of several parameters:
<ol>
  <li> The pattern of the matrix - All fill is derived from the original matrix nonzero structure.  Level zero fill
       is defined as the original matrix pattern (nonzero structure), even if the matrix value at an entry is stored
       as a zero. (Thus it is possible to add entries to the ILU factors by adding zero entries to the original matrix.)

  <li> Level of fill - Starting with the original matrix pattern as level fill of zero, the next level of fill is
       determined by analyzing the graph of the previous level and determining nonzero fill that is a result of combining
       entries that were from previous level only (not the current level).  This rule limits fill to entries that
       are direct decendents from the previous level graph.  Fill for level k is determined by applying this rule
       recursively.  For sufficiently large values of k, the fill would eventually be complete and an exact LU
       factorization would be computed.

  <li> Level of overlap - All Ifpack2 preconditioners work on parallel distributed-memory computers by using
       the row partitioning the user input matrix to determine the partitioning for local ILU factors.  If the level of
       overlap is set to zero,
       the rows of the user matrix that are stored on a given processor are treated as a self-contained local matrix
       and all column entries that reach to off-processor entries are ignored.  Setting the level of overlap to one
       tells Ifpack to increase the size of the local matrix by adding rows that are reached to by rows owned by this
       processor.  Increasing levels of overlap are defined recursively in the same way.  For sufficiently large levels
       of overlap, the entire matrix would be part of each processor's local ILU factorization process.
       Level of overlap is defined during the construction of the Ifpack2_IlukGraph object.

       Once the factorization is computed, applying the factorization \(LUy = x\)
       results in redundant approximations for any elements of y that correspond to
       rows that are part of more than one local ILU factor.  The OverlapMode (changed by calling SetOverlapMode())
       defines how these redundancies are
       handled using the Tpetra::CombineMode enum.  The default is to zero out all values of y for rows that
       were not part of the original matrix row distribution.

  <li> Fraction of relaxation - Ifpack2_RILUK computes the ILU factorization row-by-row.  As entries at a given
       row are computed, some number of them will be dropped because they do match the prescribed sparsity pattern.
       The relaxation factor determines how these dropped values will be handled.  If the RelaxValue (changed by calling
       setRelaxValue()) is zero, then these extra entries will by dropped.  This is a classical ILU approach.
       If the RelaxValue is 1, then the sum
       of the extra entries will be added to the diagonal.  This is a classical Modified ILU (MILU) approach.  If
       RelaxValue is between 0 and 1, then RelaxValue times the sum of extra entries will be added to the diagonal.

       For most situations, RelaxValue should be set to zero.  For certain kinds of problems, e.g., reservoir modeling,
       there is a conservation principle involved such that any operator should obey a zero row-sum property.  MILU
       was designed for these cases and you should set the RelaxValue to 1.  For other situations, setting RelaxValue to
       some nonzero value may improve the stability of factorization, and can be used if the computed ILU factors
       are poorly conditioned.

  <li> Diagonal perturbation - Prior to computing the factorization, it is possible to modify the diagonal entries of the matrix
       for which the factorization will be computing.  If the absolute and relative perturbation values are zero and one,
       respectively, the
       factorization will be compute for the original user matrix A.  Otherwise, the factorization
       will computed for a matrix that differs from the original user matrix in the diagonal values only.  Below we discuss
       the details of diagonal perturbations.
       The absolute and relative threshold values are set by calling SetAbsoluteThreshold() and SetRelativeThreshold(), respectively.
</ol>

\section Ifpack2_RILUK_CondEst Estimating preconditioner condition numbers

For ill-conditioned matrices, we often have difficulty computing
usable incomplete factorizations.  The most common source of problems
is that the factorization may encounter a small or zero pivot.  In
that case, the factorization may fail.  Even if the factorization
succeeds, the factors may be so poorly conditioned that use of them in
the iterative phase produces meaningless results.  Before we can fix
this problem, we must be able to detect it.  To this end, we use a
simple but effective condition number estimate for \f$(LU)^{-1}\f$.

The condition number of a matrix \f$B\f$, called \f$cond_p(B)\f$, is
defined as \f$cond_p(B) = \|B\|_p\|B^{-1}\|_p\f$ in some appropriate
norm \f$p\f$.  \f$cond_p(B)\f$ gives some indication of how many
accurate floating point digits can be expected from operations
involving the matrix and its inverse.  A condition number approaching
the accuracy of a given floating point number system, about 15 decimal
digits in IEEE double precision, means that any results involving
\f$B\f$ or \f$B^{-1}\f$ may be meaningless.

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

\section Ifpack2_RILUK_DiagPerturb A priori diagonal perturbations

If we detect using the above method that our factorization is too
ill-conditioned, we can improve the conditioning by perturbing the
matrix diagonal and restarting the factorization using this more
diagonally dominant matrix.  In order to apply perturbation, prior to
starting the factorization, we compute a diagonal perturbation of our
matrix \f$A\f$ and perform the factorization on this perturbed matrix.
The overhead cost of perturbing the diagonal is minimal since the
first step in computing the incomplete factors is to copy the matrix
\f$A\f$ into the memory space for the incomplete factors.  We simply
compute the perturbed diagonal at this point.

The actual perturbation values we use are the diagonal values \f$(d_1, d_2, \ldots, d_n)\f$
with \f$d_i = sgn(d_i)\alpha + d_i\rho\f$, \f$i=1, 2, \ldots, n\f$, where
\f$n\f$ is the matrix dimension and \f$sgn(d_i)\f$ returns
the sign of the diagonal entry.  This has the effect of
forcing the diagonal values to have minimal magnitude of \f$\alpha\f$ and
to increase each by an amount proportional to \f$\rho\f$, and still keep
the sign of the original diagonal entry.

\section Ifpack2_RILUK_Phases Phases of computation

Every Ifpack2 preconditioner has the following phases of computation:
1. initialize()
2. compute()
3. apply()

RILUK constructs the symbolic incomplete factorization (that is, the
structure of the incomplete factors) in the initialize() phase.  It
computes the numerical incomplete factorization (that is, it fills in
the factors' entries with their correct values) in the compute()
phase.  The apply() phase applies the incomplete factorization to a
given multivector using two triangular solves.

\section Ifpack2_RILUK_Measuring Measuring performance

Each RILUK object keeps track of both the time required for various
operations, and the number of times those operations have been applied
for that object.  The operations tracked include:
- initialize() (via getNumInitialize() and getInitializeTime())
- compute() (via getNumCompute() and getComputeTime())
- apply() (via getNumApply() and getApplyTime())

The <tt>getNum*</tt> methods return the number of times that operation
was called.  The <tt>get*Time</tt> methods return the number of
seconds spent in <i>all</i> invocations of that operation.  For
example, getApplyTime() returns the number of seconds spent in all
apply() calls.  For an average time per apply() call, divide by
getNumApply(), the total number of calls to apply().  
*/
template<class MatrixType>
class RILUK: public virtual Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {

 public:
  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //! RILUK constuctor with variable number of indices per row.
  /*! Creates a RILUK object and allocates storage.

    \param In
           Graph_in - Graph generated by IlukGraph.
  */
  RILUK(const Teuchos::RCP<const MatrixType>& A_in);

 private:
  //! Copy constructor.
  RILUK(const RILUK<MatrixType> & src);

 public:
  //! Ifpack2_RILUK Destructor
  virtual ~RILUK();

  //! Set RILU(k) relaxation parameter
  void SetRelaxValue( magnitude_type RelaxValue) {RelaxValue_ = RelaxValue;}

  //! Set absolute threshold value
  void SetAbsoluteThreshold( magnitude_type Athresh) {Athresh_ = Athresh;}

  //! Set relative threshold value
  void SetRelativeThreshold( magnitude_type Rthresh) {Rthresh_ = Rthresh;}

  //! Set overlap mode type
  void SetOverlapMode( Tpetra::CombineMode OverlapMode) {OverlapMode_ = OverlapMode;}

  /// Set parameters for the incomplete factorization.
  ///
  /// This preconditioner supports the following parameters:
  /// - "fact: iluk level-of-fill" (int)
  /// - "fact: absolute threshold" (magnitude_type)
  /// - "fact: relative threshold" (magnitude_type)
  /// - "fact: relax value" (magnitude_type)
  ///
  /// It will eventually also support the following parameter,
  /// although it currently does not:
  /// - "fact: iluk level-of-overlap" (int)
  void setParameters(const Teuchos::ParameterList& params);

  //! Initialize by computing the symbolic incomplete factorization.
  void initialize();

  //! Whether initialize() has been called.
  bool isInitialized() const {return isInitialized_;}

  //! How many times initialize() has been called for this object.
  int getNumInitialize() const {return numInitialize_;}

  /// \brief Compute the (numeric) incomplete factorization.
  ///
  /// This function computes the RILU(k) factors L and U using the current:
  /// - Ifpack2_IlukGraph specifying the structure of L and U.
  /// - Value for the RILU(k) relaxation parameter.
  /// - Value for the a priori diagonal threshold values.
  ///
  /// initialize() must be called first, before this method may be called.
  void compute();

  //! Whether compute() has been called.
  bool isComputed() const {return(Factored_);}

  //! How many times compute() has been called for this object.
  int getNumCompute() const {return numCompute_;}

  //! How many times apply() has been called for this object.
  int getNumApply() const {return numApply_;}

  double getInitializeTime() const {return -1;}
  double getComputeTime() const {return -1;}
  double getApplyTime() const {return -1;}

  // Mathematical functions.

  /// \brief Apply the (inverse of the) incomplete factorization to X, resulting in Y.
  ///
  /// In Matlab(tm) notation, if the incomplete factorization is \f$A \approx LDU\f$, 
  /// this method computes <tt>Y = beta*Y + alpha*(U \ (D \ (L \ X)))</tt> if mode=Teuchos::NO_TRANS, or 
  /// <tt>Y = beta*Y + alpha*(L^T \ (D^T \ (U^T \ X)))</tt> if mode=Teuchos::TRANS, or
  /// <tt>Y = beta*Y + alpha*(L^* \ (D^* \ (U^* \ X)))</tt> if mode=Teuchos::CONJ_TRANS.
  ///
  /// \param X [in] The input multivector.
  ///
  /// \param Y [in/out] The output multivector.
  ///
  /// \param mode [in] If Teuchos::TRANS resp. Teuchos::CONJ_TRANS,
  ///   apply the transpose resp. conjugate transpose of the incomplete
  ///   factorization.  Otherwise, don't apply the tranpose.
  ///
  /// \param alpha [in] Scaling factor for the result of applying the preconditioner.
  ///
  /// \param beta [in] Scaling factor for the initial value of Y.
  void apply(
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;


  /// \brief Apply the incomplete factorization (as a product) to X, resulting in Y.
  ///
  /// In Matlab(tm) notation, if the incomplete factorization is \f$A \approx LDU\f$, 
  /// this method computes <tt>Y = beta*Y + alpha*(L \ (D \ (U \ X)))</tt> mode=Teuchos::NO_TRANS, or 
  /// <tt>Y = beta*Y + alpha*(U^T \ (D^T \ (L^T \ X)))</tt> if mode=Teuchos::TRANS, or
  /// <tt>Y = beta*Y + alpha*(U^* \ (D^* \ (L^* \ X)))</tt> if mode=Teuchos::CONJ_TRANS.
  /// 
  /// \param X [in] The input multivector.
  ///
  /// \param Y [in/out] The output multivector.
  ///
  /// \param mode [in] If Teuchos::TRANS resp. Teuchos::CONJ_TRANS,
  ///   apply the transpose resp. conjugate transpose of the incomplete
  ///   factorization.  Otherwise, don't apply the tranpose.
  int Multiply(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                     Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Returns the maximum over all the condition number estimate for each local ILU set of factors.
  /*! This functions computes a local condition number estimate on each processor and return the
      maximum over all processors of the estimate.
   \param In
    Trans -If true, solve transpose problem.
    \param Out
    ConditionNumberEstimate - The maximum across all processors of
    the infinity-norm estimate of the condition number of the inverse of LDU.
  */
  magnitude_type computeCondEst(Teuchos::ETransp mode) const;
  magnitude_type computeCondEst(CondestType CT = Ifpack2::Cheap,
                               local_ordinal_type MaxIters = 1550,
                               magnitude_type Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &Matrix = Teuchos::null)
  {
    std::cerr << "Warning, Ifpack2::RILUK::computeCondEst currently does not use MaxIters/Tol/etc arguments..." << std::endl;
    return computeCondEst(Teuchos::NO_TRANS);
  }

  magnitude_type getCondEst() const {return Condest_;}

  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const
  {
    return A_;
  }

  // Attribute access functions

  //! Get RILU(k) relaxation parameter
  magnitude_type GetRelaxValue() const {return RelaxValue_;}

  //! Get absolute threshold value
  magnitude_type getAbsoluteThreshold() const {return Athresh_;}

  //! Get relative threshold value
  magnitude_type getRelativeThreshold() const {return Rthresh_;}

  int getLevelOfFill() const { return LevelOfFill_; }

  //! Get overlap mode type
  Tpetra::CombineMode getOverlapMode() {return OverlapMode_;}

  //! Returns the number of nonzero entries in the global graph.
  int getGlobalNumEntries() const {return(getL().getGlobalNumEntries()+getU().getGlobalNumEntries());}

  //! Returns the Ifpack2::IlukGraph associated with this factored matrix.
  const Teuchos::RCP<Ifpack2::IlukGraph<local_ordinal_type,global_ordinal_type,node_type> >& getGraph() const {return(Graph_);}

  //! Returns the L factor associated with this factored matrix.
  const MatrixType& getL() const {return(*L_);}

  //! Returns the D factor associated with this factored matrix.
  const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> & getD() const {return(*D_);}

  //! Returns the U factor associated with this factored matrix.
  const MatrixType& getU() const {return(*U_);}

  //@{ \name Additional methods required to support the Tpetra::Operator interface.

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& getDomainMap() const
  { return Graph_->getL_Graph()->getDomainMap(); }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& getRangeMap() const
  { return Graph_->getU_Graph()->getRangeMap(); }

  //@}

 protected:
  void setFactored(bool Flag) {Factored_ = Flag;}
  void setInitialized(bool Flag) {isInitialized_ = Flag;}
  bool isAllocated() const {return(isAllocated_);}
  void setAllocated(bool Flag) {isAllocated_ = Flag;}

 private:


  void allocate_L_and_U();
  void initAllValues(const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> & overlapA);
  void generateXY(Teuchos::ETransp mode,
                 const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Xin,
     const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Yin,
     Teuchos::RCP<const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Xout,
     Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Yout) const;
  bool isOverlapped_;
  Teuchos::RCP<Ifpack2::IlukGraph<local_ordinal_type,global_ordinal_type,node_type> > Graph_;
  const Teuchos::RCP<const MatrixType> A_;
  Teuchos::RCP<MatrixType> L_;
  Teuchos::RCP<MatrixType> U_;
  Teuchos::RCP<Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > D_;
  bool UseTranspose_;

  int LevelOfFill_;
  int LevelOfOverlap_;

  int NumMyDiagonals_;
  bool isAllocated_;
  bool isInitialized_;
  mutable int numInitialize_;
  mutable int numCompute_;
  mutable int numApply_;
  bool Factored_;
  magnitude_type RelaxValue_;
  magnitude_type Athresh_;
  magnitude_type Rthresh_;
  mutable magnitude_type Condest_;

  mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > OverlapX_;
  mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > OverlapY_;
  Tpetra::CombineMode OverlapMode_;
};

}//namespace Ifpack2

#endif /* IFPACK2_CRSRILUK_DECL_HPP */
