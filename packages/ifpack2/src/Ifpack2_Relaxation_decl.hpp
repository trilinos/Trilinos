/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_RELAXATION_DECL_HPP
#define IFPACK2_RELAXATION_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_CrsMatrix.hpp" // Don't need the definition here
#include "Tpetra_BlockCrsMatrix.hpp"
#include <type_traits>
#include <KokkosKernels_Handle.hpp>
#include "KokkosSparse_sor_sequential_impl.hpp"
#include "Tpetra_BlockView.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Ifpack2 {
namespace Details {

template<class TpetraOperatorType>
class ScaledDampedResidual; // forward declaration

} // namespace Details
} // namespace Ifpack2

namespace Teuchos {
  // forward declarations
  class ParameterList;
  class Time;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Ifpack2 {

/** \class Relaxation
\brief Relaxation preconditioners for Tpetra::RowMatrix and
  Tpetra::CrsMatrix sparse matrices.
\tparam MatrixType A specialization of Tpetra::RowMatrix.

\section Ifpack_Relaxation_Summary Summary

This class implements several different relaxation preconditioners for
Tpetra::RowMatrix or its subclass Tpetra::CrsMatrix.  Relaxation is
derived from Preconditioner, which is itself derived from
Tpetra::Operator.  Therefore this object may be used as a
preconditioner for Belos linear solvers, and for any linear solver
that treats preconditioners as instances of Tpetra::Operator.

This class implements the following relaxation methods:
- Richardson
- Jacobi
- Gauss-Seidel
- Symmetric Gauss-Seidel

All methods allow you to set an optional damping parameter.  The
"Gauss-Seidel" methods technically only perform Gauss-Seidel within an
MPI process, but Jacobi between processes.  To compensate, these
methods include an "L1" option, which can improve convergence by
weighting contributions near process boundaries differently.  For more
details, please refer to the following publication:

A. H. Baker, R. D. Falgout, T. V. Kolev, and U. M. Yang.
"Multigrid Smoothers for Ultraparallel Computing."
<i>SIAM J. Sci. Comput.</i>, Vol. 33, No. 5. (2011),
pp. 2864-2887.

\section Ifpack_Relaxation_Performance Performance

Richardson and Jacobi will always use your matrix's native sparse matrix-vector
multiply kernel.  This should give good performance, since we have
spent a lot of effort tuning Tpetra's kernels.  Depending on the node_type
type of your Tpetra matrix, it may also exploit threads for additional
parallelism within each MPI process.  In contrast, Gauss-Seidel and
symmetric Gauss-Seidel are intrinsically sequential methods within an
MPI process.  This prevents us from exposing more parallelism via
threads.  The difference should become more apparent as your code
moves away from a "one MPI process per core" model, to a "one MPI
process per socket or node" model, assuming that you are using a
thread-parallel node_type type.

Relaxation works with any Tpetra::RowMatrix input.  If your
Tpetra::RowMatrix happens to be a Tpetra::CrsMatrix, the Gauss-Seidel
and symmetric Gauss-Seidel relaxations may be able to exploit this for
better performance.  You don't have to do anything to figure this out
(we test via \c dynamic_cast).

\section Ifpack_Relaxation_Create Creating a Relaxation preconditioner

The following code snippet shows how to create a Relaxation
preconditioner.

\code
#include "Ifpack2_Relaxation.hpp"

...
using Teuchos::ParameterList;
using Teuchos::RCP;
typedef Tpetra::CrsMatrix<double> crs_matrix_type;
typedef Ifpack2::Preconditioner<double> precond_type;
...

// Create the sparse matrix A somehow.  It must be fill complete
// before you may create an Ifpack2 preconditioner from it.
RCP<crs_matrix_type> A = ...;

// Create the relaxation.  You could also do this using
// Ifpack2::Factory (the preconditioner factory) if you like.
precond_type prec (A);

// Make the list of relaxation parameters.
Teuchos::ParameterList params;
// Do symmetric SOR / Gauss-Seidel.
params.set ("relaxation: type", "Symmetric Gauss-Seidel");
// Two sweeps (of symmetric SOR / Gauss-Seidel) per apply() call.
params.set ("relaxation: sweeps", 2);
// ... Set any other parameters you want to set ...

// Set parameters.
prec.setParameters (params);

// Prepare the relaxation instance for use.
prec.initialize ();
prec.compute ();

// Now prec may be used as a preconditioner or smoother,
// by calling its apply() method, just like any Tpetra::Operator.
\endcode

\section Ifpack_Relaxation_Algorithms Algorithms

We now briefly describe the relaxation algorithms this class
implements.  Consider the linear system \f$Ax=b\f$, where \f$A\f$ is a
square matrix, and \f$x\f$ and \f$b\f$ are two vectors of compatible
dimensions.  Suppose that \f$x^{(0)}\f$ is the starting vector and
\f$x^{(k)}\f$ is the approximate solution for \f$x\f$ computed by
iteration $k+1$ of whatever relaxation method we are using.  Here,
\f$x^{(k)}_i\f$ is the $i$-th element of vector \f$x^{(k)}\f$.

The Richardson method computes
\f[
x^{(k+1)}_i = x_^{(k)}_i + alpha ( b_i - \sum_{j} A_{ij} x^{(k)}_j ).
\f]

The Jacobi method computes
\f[
x^{(k+1)}_i = A_{ii}^{-1} ( b_i - \sum_{j \neq i} A_{ij} x^{(k)}_j ).
\f]
The "damped" Jacobi method generalizes Jacobi.  It introduces a
damping parameter \f$\omega \f$, and computes
\f[
x^{(k+1)}_i = (1 - \omega) x^{(k)}_i + \omega A_{ii}^{-1} ( b_i - \sum_{j \neq i} A_{ij} x^{(k)}_j ).
\f]

The "damped Gauss-Seidel method" is actually successive over-relaxation
(SOR), with Gauss-Seidel as a special case when the damping parameter
\f$\omega = 1\f$.  We implement has two different sweep directions: Forward and
Backward.  The Forward sweep direction computes
\f[
x^{(k+1)}_i = (1 - \omega) x^{(k)}_i + \omega A_{ii}^{-1} ( b_i - \sum_{j < i} A_{ij} x^{(k+1)}_j - \sum_{j > i} A_{ij} x^{(k)}_j ),
\f]
and the Backward sweep direction computes
\f[
x^{(k+1)}_i = (1 - \omega) x^{(k)}_i + \omega A_{ii}^{-1} ( b_i - \sum_{j > i} A_{ij} x^{(k+1)}_j - \sum_{j < i} A_{ij} x^{(k)}_j ),
\f]
Users may set the sweep direction via the "relaxation: backward mode"
option.  See the documentation of setParameters() for details.

Gauss-Seidel / SOR also comes in a symmetric version.  This method
first does a Forward sweep, then a Backward sweep.  Only the symmetric
version of this preconditioner is guaranteed to be symmetric (or Hermitian,
if the matrix's data are complex).

Users may set the relaxation method via the "relaxation: type"
parameter.  For all relaxation methods, users may specify the number
of sweeps per call to apply() and the damping parameter \f$\omega \f$.
For a list of all supported parameters, please refer to the
documentation of the setParameters() method.  For advice on picking
\f$\omega \f$ for a preconditioner, please refer to the following
book: "Templates for the Solution of Linear Systems: Building Blocks
for Iterative Methods, 2nd Edition," R. Barrett et al., SIAM, 1994.

\note This class does not actually use the formulae above to apply
Jacobi or SOR.  For example, the computational kernels for the above SOR
sweeps actually do not require branches in the inner loop to distinguish
between the lower triangle, diagonal, and upper triangle of A.  One can
see this by multiplying through the forward sweep expression by \f$A_{ii}\f$
and combining terms, then dividing through again by \f$A_{ii}\f$.  This
results in the expression
\f[
x^{(k+1)}_i = x^{(k)}_i + \omega b_i - \frac{\omega}{A_{ii}} ( \sum_{j \geq i} A_{ij} x^{(k)}_j + \sum_{j < i} x^{(k+1)}_j ).
\f]
Executing this expression in a forward sweep does not require
distinguishing between the lower and upper triangle of A.  The
same thing holds for the backward sweep.
*/
template<class MatrixType>
class Relaxation :
  virtual public Ifpack2::Preconditioner<
    typename MatrixType::scalar_type,
    typename MatrixType::local_ordinal_type,
    typename MatrixType::global_ordinal_type,
    typename MatrixType::node_type>,
  virtual public Ifpack2::Details::CanChangeMatrix<
    Tpetra::RowMatrix<typename MatrixType::scalar_type,
                      typename MatrixType::local_ordinal_type,
                      typename MatrixType::global_ordinal_type,
                      typename MatrixType::node_type> >
{
public:
  //! @name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The node_type type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The Kokkos device type used by the input MatrixType.
  typedef typename MatrixType::node_type::device_type device_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Tpetra::RowMatrix specialization used by this class.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::Relaxation: Please use MatrixType = Tpetra::RowMatrix.  This saves build times, library sizes, and executable sizes.  Don't worry, this class still works with CrsMatrix and BlockCrsMatrix; those are both subclasses of RowMatrix.");

  //@}
  //! @name Constructors and destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param A [in] The matrix for which to make the constructor.
  ///   Tpetra::RowMatrix is the base class of Tpetra::CrsMatrix, so
  ///   you may give either a Tpetra::RowMatrix or a Tpetra::CrsMatrix
  ///   here.
  ///
  /// The results of apply() are undefined if you change the diagonal
  /// entries of the sparse matrix after invoking this constructor,
  /// without first calling compute().  In particular, the compute()
  /// method may extract the diagonal entries and precompute their
  /// inverses, in order to speed up any of the relaxation methods that
  /// this class implements.
  ///
  /// The "explicit" keyword just means that you must invoke the
  /// Relaxation constructor explicitly; you aren't allowed to use it
  /// as an implicit conversion ("cast").  For example, you may do
  /// this (namespaces and Tpetra template parameters omitted for
  /// brevity):
  /// \code
  /// RCP<const CrsMatrix<...> > A = ...;
  /// Relaxation<CrsMatrix<...> > R (A);
  /// \endcode
  /// but you may not do this:
  /// \code
  /// void foo (const Relaxation<CrsMatrix<...> >& R);
  ///
  /// RCP<const CrsMatrix<...> > A = ...;
  /// foo (A);
  /// \endcode
  explicit Relaxation (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor.
  virtual ~Relaxation () = default;

  //@}
  //! @name Preconditioner computation methods
  //@{

  /// \brief Set the relaxation / preconditioner parameters.
  ///
  /// \warning All parameters are case sensitive.  We make no attempt
  ///   to check the spelling of parameter names in your input list.
  ///
  /// The "relaxation: type" (string) parameter sets the relaxation /
  /// preconditioner method you want to use.  It currently accepts the
  /// following values (the default is "Jacobi"):
  /// - "Richardson"
  /// - "Jacobi"
  /// - "Gauss-Seidel"
  /// - "Symmetric Gauss-Seidel"
  ///
  /// The "relaxation: sweeps" (int) parameter sets the number of
  /// sweeps, that is, the number of times to apply the relaxation on
  /// each invocation of apply().  The default number of sweeps is 1.
  ///
  /// The "relaxation: damping factor" (scalar_type -- the type of the
  /// entries of the matrix) parameter is the value of the damping
  /// factor \f$\omega \f$.  The main documentation of this class
  /// explains how we use this value.  The default value is 1.0.
  ///
  /// The "relaxation: zero starting solution" (bool) parameter
  /// governs whether or not we use the existing values in the output
  /// multivector Y when applying the relaxation.  Its default value
  /// is true, meaning that we fill Y with zeros before applying
  /// relaxation sweeps.  If false, we use the existing values in Y.
  ///
  /// If the "relaxation: backward mode" (bool) parameter is true, we
  /// perform Gauss-Seidel in reverse mode.  The default value is
  /// false, meaning that we do forward-mode Gauss-Seidel.  This only
  /// affects standard Gauss-Seidel, not symmetric Gauss-Seidel.
  ///
  /// The "relaxation: fix tiny diagonal entries" (bool) parameter
  /// defaults to false.  If true, the compute() method will do extra
  /// work (computation only, no MPI communication) to "fix" diagonal
  /// entries that are less than or equal to the threshold given by
  /// the (magnitude of the) "relaxation: min diagonal value"
  /// parameter.  The default behavior imitates that of Aztec, which
  /// does not do any special modification of the diagonal.
  ///
  /// The "relaxation: min diagonal value" (scalar_type) parameter
  /// only matters if "relaxation: fix tiny diagonal entries" (see
  /// above) is true.  This parameter limits how close to zero the
  /// diagonal elements of the matrix are allowed to be.  If the
  /// magnitude of a diagonal element of the matrix is less than the
  /// magnitude of this value, then we set that diagonal element to
  /// this value.  (We don't actually modify the matrix; we just
  /// remember the diagonal values.)  The use of magnitude rather than
  /// the value itself makes this well defined if scalar_type is
  /// complex.  The default value of this parameter is zero, in which
  /// case we will replace diagonal entries that are exactly equal to
  /// zero with a small nonzero value (machine precision for the given
  /// \c Scalar type) before inverting them.  Note that if
  /// "relaxation: fix tiny diagonal entries" is false, the default
  /// value, this parameter does nothing.)
  ///
  /// The "relaxation: check diagonal entries" (bool) parameter
  /// defaults to false.  If true, the compute() method will do extra
  /// work (both computation and communication) to count diagonal
  /// entries that are zero, have negative real part, or are small in
  /// magnitude.  The describe() method will then print this
  /// information for you.  You may find this useful for checking
  /// whether your input matrix has issues that make Jacobi or
  /// Gauss-Seidel a poor choice of preconditioner.
  ///
  /// The last two parameters govern the L1 variant of Gauss-Seidel.
  /// The "relaxation: use l1" (bool) parameter, if true, turns on the
  /// L1 variant.  (In "l1", the first character is a lower-case L,
  /// and the second character is the numeral 1 (one).)  This
  /// parameter's value is false by default.  The "relaxation: l1 eta"
  /// (magnitude_type) parameter is the \f$\eta \f$ parameter
  /// associated with that method; its default value is 1.5.  Recall
  /// that "magnitude_type" is the type of the absolute value of a
  /// scalar_type value.  This is the same as scalar_type for
  /// real-valued floating-point types (like \c float and \c double).
  /// If scalar_type is <tt>std::complex<T></tt> for some type \c T,
  /// then magnitude_type is \c T.
  void setParameters (const Teuchos::ParameterList& params);

  //! Return a list of all the parameters that this class accepts.
  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters () const;

  //! Initialize the preconditioner ("symbolic setup").
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return isInitialized_;
  }

  //! Compute the preconditioner ("numeric setup");
  void compute ();

  //! Return true if compute() has been called.
  inline bool isComputed() const {
    return IsComputed_;
  }

  //@}
  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param A [in] The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method with a matrix different than the current
  /// matrix resets the preconditioner's state.  After calling this
  /// method with a nonnull input, you must first call initialize()
  /// and compute() (in that order) before you may call apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  virtual void
  setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //@}
  //! @name Implementation of the Tpetra::Operator interface
  //@{

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method computes Y = beta*Y + alpha*M*X, where M*X
  /// represents the action of the preconditioner on the input
  /// multivector X.
  ///
  /// \param X [in] The multivector input of the preconditioner.
  /// \param Y [in/out] The multivector output of the preconditioner.
  /// \param mode [in] Whether to apply the transpose or conjugate
  ///   transpose of the preconditioner.  Not all preconditioners
  ///   support options other than the default (no transpose); please
  ///   call hasTransposeApply() to determine whether nondefault
  ///   options are supported.
  /// \param alpha [in] Scaling factor for the preconditioned input.
  /// \param beta [in] Scaling factor for the output.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getDomainMap () const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getRangeMap () const;

  /// \brief Whether apply() and applyMat() let you apply the
  ///   transpose or conjugate transpose.
  bool hasTransposeApply () const;

  /// \brief Apply the input matrix to X, returning the result in Y.
  ///
  /// \param X [in] The multivector input of the preconditioner.
  /// \param Y [in/out] The multivector output of the preconditioner.
  /// \param mode [in] Whether to apply the transpose or conjugate
  ///   transpose of the matrix.
  void
  applyMat (const Tpetra::MultiVector<
              scalar_type,
              local_ordinal_type,
              global_ordinal_type,
              node_type>& X,
            Tpetra::MultiVector<
              scalar_type,
              local_ordinal_type,
              global_ordinal_type,
              node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! @name Attribute accessor methods
  //@{

  //! The communicator over which the matrix and vectors are distributed.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! Total number of floating-point operations over all calls to compute().
  double getComputeFlops() const;

  //! Total number of floating-point operations over all calls to apply().
  double getApplyFlops() const;

  //! Total number of calls to initialize().
  int getNumInitialize() const;

  //! Total number of calls to compute().
  int getNumCompute() const;

  //! Total number of calls to apply().
  int getNumApply() const;

  //! Total time in seconds spent in all calls to initialize().
  double getInitializeTime() const;

  //! Total time in seconds spent in all calls to compute().
  double getComputeTime() const;

  //! Total time in seconds spent in all calls to apply().
  double getApplyTime() const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}
  //! @name Implementation of Teuchos::Describable interface
  //@{

  /// \brief A simple one-line description of this object.
  ///
  /// Be aware that this will print a very long line, because some
  /// users really like to see all the attributes of the object in a
  /// single line.  If you prefer multiple lines of output, you should
  /// call describe() instead.
  std::string description () const;

  /// \brief Print the object's attributes to the given output stream.
  ///
  /// This method will print a constant amount of information (not
  /// proportional to the matrix's dimensions or number of entries) on
  /// Process 0 of the communicator over which this object is
  /// distributed.
  ///
  /// You may wrap an std::ostream in a Teuchos::FancyOStream by
  /// including "Teuchos_FancyOStream.hpp" and calling
  /// Teuchos::getFancyOStream().  For example:
  /// \code
  /// using Teuchos::RCP;
  /// using Teuchos::rcpFromRef;
  /// using Teuchos::FancyOStream;
  ///
  /// // Wrap std::cout in a FancyOStream.
  /// RCP<FancyOStream> wrappedCout = getFancyOStream (rcpFromRef (std::cout));
  ///
  /// // Wrap an output file in a FancyOStream.
  /// RCP<std::ofstream> outFile (new std::ofstream ("myFile.txt"));
  /// RCP<FancyOStream> wrappedFile = getFancyOStream (outFile);
  /// \endcode
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  //! \name Internal typedefs
  //@{

  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

    typedef typename Kokkos::ArithTraits<scalar_type>::val_type impl_scalar_type;

  /// \brief Tpetra::CrsMatrix specialization used by this class.
  ///
  /// We use this for dynamic casts to dispatch to the most efficient
  /// implementation of various relaxation kernels.
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> crs_matrix_type;
  typedef Tpetra::CrsGraph<local_ordinal_type,
                            global_ordinal_type, node_type> crs_graph_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> multivector_type;
  typedef Tpetra::BlockCrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> block_crs_matrix_type;
  typedef Tpetra::BlockMultiVector<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> block_multivector_type;

  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;


  //@}
  //! \name Implementation of multithreaded Gauss-Seidel.
  //@{

  typedef typename crs_matrix_type::local_matrix_type local_matrix_type;
  typedef typename local_matrix_type::StaticCrsGraphType::row_map_type lno_row_view_t;
  typedef typename local_matrix_type::StaticCrsGraphType::entries_type lno_nonzero_view_t;
  typedef typename local_matrix_type::values_type scalar_nonzero_view_t;
  typedef typename local_matrix_type::StaticCrsGraphType::device_type TemporaryWorkSpace;
  typedef typename local_matrix_type::StaticCrsGraphType::device_type PersistentWorkSpace;
  typedef typename local_matrix_type::StaticCrsGraphType::execution_space MyExecSpace;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle
      <typename lno_row_view_t::const_value_type, local_ordinal_type,typename scalar_nonzero_view_t::value_type,
      MyExecSpace, TemporaryWorkSpace,PersistentWorkSpace > mt_kernel_handle_type;
  Teuchos::RCP<mt_kernel_handle_type> mtKernelHandle_;

  //@}
  //! \name Unimplemented methods that you are syntactically forbidden to call.
  //@{

  //! Copy constructor (not implemented; you are not allowed to call this).
  Relaxation (const Relaxation<MatrixType>& RHS);

  //! Assignment operator (not implemented; you are not allowed to call this).
  Relaxation<MatrixType>& operator= (const Relaxation<MatrixType>& RHS);

  //@}
  //! @name Internal methods
  //@{

  /// \brief Variant of setParameters() that takes a nonconst Teuchos::ParameterList.
  ///
  /// This variant fills in default values for any valid parameters
  /// that are not in the input list.
  void setParametersImpl (Teuchos::ParameterList& params);

 //! Apply Richardson to X, returning the result in Y.
  void ApplyInverseRichardson(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Jacobi to X, returning the result in Y.
  void ApplyInverseJacobi(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Jacobi to X, returning the result in Y.
  void ApplyInverseJacobi_BlockCrsMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply multi-threaded Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseMTGS_CrsMatrix(
          const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;


  //! Apply Gauss-Seidel for a Tpetra::RowMatrix specialization.
  void ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Gauss-Seidel for a Tpetra::CrsMatrix specialization.
  void
  ApplyInverseGS_CrsMatrix (const crs_matrix_type& A,
                            const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Gauss-Seidel for a Tpetra::BlockCrsMatrix specialization.
  void
  ApplyInverseGS_BlockCrsMatrix (const block_crs_matrix_type& A,
                            const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y);

  //! Apply symmetric Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseSGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply symmetric multi-threaded Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseMTSGS_CrsMatrix(
          const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  void MTGaussSeidel (
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& B,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
      const Tpetra::ESweepDirection direction) const;

  //! Apply symmetric Gauss-Seidel for a Tpetra::RowMatrix specialization.
  void ApplyInverseSGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply symmetric Gauss-Seidel for a Tpetra::CrsMatrix specialization.
  void
  ApplyInverseSGS_CrsMatrix (const crs_matrix_type& A,
                             const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                             Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply symmetric Gauss-Seidel for a Tpetra::BlockCrsMatrix specialization.
  void
  ApplyInverseSGS_BlockCrsMatrix (const block_crs_matrix_type& A,
                             const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                             Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y);

  void computeBlockCrs ();

  //! A service routine for updating the cached MultiVector
  void updateCachedMultiVector(const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& map, size_t numVecs) const;


  //@}
  //! @name Internal data and parameters
  //@{

  /// \brief List of valid parameters.
  ///
  /// This is created on demand by getValidParameters().  That method
  /// is const to help comply with the Teuchos::ParameterListAcceptor
  /// interface (which we might like to use later), which is why we
  /// have to declare this field \c mutable.
  mutable Teuchos::RCP<const Teuchos::ParameterList> validParams_;

  //! The matrix for which to construct the preconditioner or smoother.
  Teuchos::RCP<const row_matrix_type> A_;
  //! Time object to track timing (setup).
  //! Importer for parallel Gauss-Seidel and symmetric Gauss-Seidel.
  Teuchos::RCP<const import_type> Importer_;
  //! Importer for block multivector versions of GS and SGS
  Teuchos::RCP<const import_type> pointImporter_;
  //! Contains the diagonal elements of \c A_.
  Teuchos::RCP<Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Diagonal_;
  //! MultiVector for caching purposes (so apply doesn't need to allocate one on each call)
  mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > cachedMV_;

  typedef Kokkos::View<typename block_crs_matrix_type::impl_scalar_type***,
                       typename block_crs_matrix_type::device_type> block_diag_type;
  typedef Kokkos::View<typename block_crs_matrix_type::impl_scalar_type***,
                       typename block_crs_matrix_type::device_type,
                       Kokkos::MemoryUnmanaged> unmanaged_block_diag_type;

  /// \brief Storage of the BlockCrsMatrix's block diagonal.
  ///
  /// This is only allocated and used if the input matrix is a
  /// Tpetra::BlockCrsMatrix.  In that case, Ifpack2::Relaxation does
  /// block relaxation, using the (small dense) blocks in the
  /// BlockCrsMatrix.
  ///
  /// "Block diagonal" means the blocks in the BlockCrsMatrix
  /// corresponding to the diagonal entries of the BlockCrsMatrix's
  /// graph.  To get the block corresponding to local (graph
  /// a.k.a. "mesh") row index i, do the following:
  /// \code
  /// auto D_ii = Kokkos::subview (blockDiag_, i, Kokkos::ALL (), Kokkos::ALL ());
  /// \endcode
  block_diag_type blockDiag_;

  Teuchos::RCP<block_multivector_type> yBlockColumnPointMap_;

  //! How many times to apply the relaxation per apply() call.
  int NumSweeps_ = 1;
  //! Number of inner-sweeps for the two-stage Gauss Seidel
  int NumInnerSweeps_ = 1;
  //! Which relaxation method to use.
  Details::RelaxationType PrecType_ = Ifpack2::Details::JACOBI;
  //! Damping factor
  scalar_type DampingFactor_ = STS::one();
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_ = true;
  //! If true, do backward-mode Gauss-Seidel.
  bool DoBackwardGS_ = false;
  //! If true, do the L1 version of Jacobi, Gauss-Seidel, or symmetric Gauss-Seidel.
  bool DoL1Method_ = false;
  //! Eta parameter for modified L1 method
  magnitude_type L1Eta_ = Teuchos::as<magnitude_type>(1.5);
  //! Minimum diagonal value
  scalar_type MinDiagonalValue_ = STS::zero();
  //! Whether to fix up zero or tiny diagonal entries.
  bool fixTinyDiagEntries_ = false;
  //! Whether to spend extra effort and all-reduces checking diagonal entries.
  bool checkDiagEntries_ = false;
  //! Whether to use sparse-triangular solve instead of inner-iterations
  bool InnerSpTrsv_ = false;
  //! For MTSGS, the cluster size (use point coloring if equal to 1)
  int clusterSize_ = 1;

  //!Wheter the provided matrix is structurally symmetric or not.
  bool is_matrix_structurally_symmetric_ = false;

  //!Whether to write the given input file
  bool ifpack2_dump_matrix_ = false;


  //! If \c true, the preconditioner has been initialized successfully.
  bool isInitialized_ = false;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_ = false;
  //! The number of successful calls to initialize().
  int NumInitialize_ = 0;
  //! the number of successful calls to compute().
  int NumCompute_ = 0;
  //! The number of successful calls to apply().
  mutable int NumApply_ = 0;
  //! Total time in seconds for all successful calls to initialize().
  double InitializeTime_ = 0.0;
  //! Total time in seconds for all successful calls to compute().
  double ComputeTime_ = 0.0;
  //! Total time in seconds for all successful calls to apply().
  mutable double ApplyTime_ = 0.0;
  //! The total number of floating-point operations for all successful calls to compute().
  double ComputeFlops_ = 0.0;
  //! The total number of floating-point operations for all successful calls to apply().
  mutable double ApplyFlops_ = 0.0;

  //! Global magnitude of the diagonal entry with the minimum magnitude.
  magnitude_type globalMinMagDiagEntryMag_ = STM::zero();
  //! Global magnitude of the diagonal entry with the maximum magnitude.
  magnitude_type globalMaxMagDiagEntryMag_ = STM::zero();
  //! Global number of small (in magnitude) diagonal entries detected by compute().
  size_t globalNumSmallDiagEntries_ = 0;
  //! Global number of zero diagonal entries detected by compute().
  size_t globalNumZeroDiagEntries_ = 0;
  //! Global number of negative (real part) diagonal entries detected by compute().
  size_t globalNumNegDiagEntries_ = 0;
  /// \brief Absolute two-norm difference between computed and actual inverse diagonal.
  ///
  /// "Actual inverse diagonal" means the result of 1/diagonal,
  /// without any protection against zero or small diagonal entries.
  magnitude_type globalDiagNormDiff_ = STM::zero();

  /// \brief Precomputed offsets of local diagonal entries of the matrix.
  ///
  /// These are only used if the matrix has a const ("static") graph.
  /// In that case, the offsets of the diagonal entries will never
  /// change, even if the values of the diagonal entries change.
  Kokkos::View<size_t*, typename node_type::device_type> diagOffsets_;

  /// \brief Whether we have precomputed offsets of diagonal entries.
  ///
  /// We need this flag because it is not enough just to test if
  /// diagOffsets_ has size zero.  It is perfectly legitimate for the
  /// matrix to have zero rows on the calling process.
  bool savedDiagOffsets_ = false;

  bool hasBlockCrsMatrix_ = false;

  /// \brief In case of local/reordered smoothing, the unknowns to use
  Teuchos::ArrayRCP<local_ordinal_type> localSmoothingIndices_;

  //@}
  
  static Teuchos::RCP<multivector_type> getRowMapMultiVector(const crs_matrix_type* A, const multivector_type& U_rangeMap, bool force)
  {
    const size_t numVecs = U_rangeMap.getNumVectors ();
    auto exporter = A->getGraph ()->getExporter ();
    Teuchos::RCP<const map_type> rowMap = A->getRowMap ();
    Teuchos::RCP<multivector_type> U_rowMap; // null by default
    if (! exporter.is_null () || force) {
      U_rowMap = rcp (new multivector_type(rowMap, numVecs));
    }
    return U_rowMap;
  }

  static Teuchos::RCP<multivector_type> getColumnMapMultiVector(const crs_matrix_type* A, const multivector_type& U_domainMap, bool force)
  {
    const size_t numVecs = U_domainMap.getNumVectors ();
    auto importer = A->getGraph ()->getImporter();
    Teuchos::RCP<const map_type> colMap = A->getColMap();
    Teuchos::RCP<multivector_type> U_colMap; // null by default
    if (! importer.is_null () || force) {
      U_colMap = rcp (new multivector_type(colMap, numVecs));
    }
    return U_colMap;
  }

  void CRS_localGaussSeidel (
      const crs_matrix_type* A,
      const multivector_type &B,
      multivector_type &X,
      const multivector_type &D,
      const scalar_type& dampingFactor,
      const Tpetra::ESweepDirection direction) const
  {
    typedef typename node_type::device_type::memory_space dev_mem_space;
    typedef typename multivector_type::dual_view_type::t_host::device_type host_mem_space;
    typedef typename crs_matrix_type::crs_graph_type::local_graph_type k_local_graph_type;
    typedef typename k_local_graph_type::size_type offset_type;
    const char prefix[] = "Tpetra::CrsMatrix::localGaussSeidel: ";

    TEUCHOS_TEST_FOR_EXCEPTION
      (! A->isFillComplete (), std::runtime_error,
       prefix << "The matrix is not fill complete.");
    const size_t lclNumRows = A->getNodeNumRows ();
    const size_t numVecs = B.getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getNumVectors () != numVecs, std::invalid_argument,
       prefix << "B.getNumVectors() = " << numVecs << " != "
       "X.getNumVectors() = " << X.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.getLocalLength () != lclNumRows, std::invalid_argument,
       prefix << "B.getLocalLength() = " << B.getLocalLength ()
       << " != this->getNodeNumRows() = " << lclNumRows << ".");

    // mfh 28 Aug 2017: The current local Gauss-Seidel kernel only
    // runs on host.  (See comments below.)  Thus, we need to access
    // the host versions of these data.
    const_cast<multivector_type&> (B).sync_host ();
    X.sync_host ();
    X.modify_host ();
    const_cast<multivector_type&> (D).sync_host ();

    auto B_lcl = B.template getLocalView<host_mem_space> ();
    auto X_lcl = X.template getLocalView<host_mem_space> ();
    auto D_lcl = D.template getLocalView<host_mem_space> ();

    offset_type B_stride[8], X_stride[8], D_stride[8];
    B_lcl.stride (B_stride);
    X_lcl.stride (X_stride);
    D_lcl.stride (D_stride);

    local_matrix_type lclMatrix = A->getLocalMatrix ();
    k_local_graph_type lclGraph = lclMatrix.graph;
    typename local_matrix_type::row_map_type ptr = lclGraph.row_map;
    typename local_matrix_type::index_type ind = lclGraph.entries;
    typename local_matrix_type::values_type val = lclMatrix.values;
    const offset_type* const ptrRaw = ptr.data ();
    const local_ordinal_type* const indRaw = ind.data ();
    const impl_scalar_type* const valRaw = val.data ();

    const std::string dir ((direction == Tpetra::Forward) ? "F" : "B");
    // NOTE (mfh 28 Aug 2017) This assumes UVM.  We can't get around
    // that on GPUs without using a GPU-based sparse triangular
    // solve to implement Gauss-Seidel.  This exists in cuSPARSE,
    // but we would need to implement a wrapper with a fall-back
    // algorithm for unsupported Scalar and LO types.
    KokkosSparse::Impl::Sequential::gaussSeidel<local_ordinal_type, offset_type, impl_scalar_type, impl_scalar_type, impl_scalar_type>(
        lclNumRows,
        numVecs,
        ptrRaw, indRaw, valRaw,
        B_lcl.data (), B_stride[1],
        X_lcl.data (), X_stride[1],
        D_lcl.data (),
        static_cast<impl_scalar_type> (dampingFactor),
        dir.c_str ());
    const_cast<multivector_type&> (B).template sync<dev_mem_space> ();
    X.template sync<dev_mem_space> ();
    const_cast<multivector_type&> (D).template sync<dev_mem_space> ();
  }

  void CRS_reorderedLocalGaussSeidel (
      const crs_matrix_type* A,
      const multivector_type& B,
      multivector_type& X,
      const multivector_type& D,
      const Teuchos::ArrayView<local_ordinal_type>& rowIndices,
      const scalar_type& dampingFactor,
      const Tpetra::ESweepDirection direction) const
  {
    typedef typename node_type::device_type::memory_space dev_mem_space;
    typedef typename multivector_type::dual_view_type::t_host::device_type host_mem_space;
    typedef typename crs_matrix_type::crs_graph_type::local_graph_type k_local_graph_type;
    typedef typename k_local_graph_type::size_type offset_type;
    const char prefix[] = "Tpetra::CrsMatrix::reorderedLocalGaussSeidel: ";

    TEUCHOS_TEST_FOR_EXCEPTION
      (! A->isFillComplete (), std::runtime_error,
       prefix << "The matrix is not fill complete.");
    const size_t lclNumRows = A->getNodeNumRows ();
    const size_t numVecs = B.getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getNumVectors () != numVecs, std::invalid_argument,
       prefix << "B.getNumVectors() = " << numVecs << " != "
       "X.getNumVectors() = " << X.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.getLocalLength () != lclNumRows, std::invalid_argument,
       prefix << "B.getLocalLength() = " << B.getLocalLength ()
       << " != this->getNodeNumRows() = " << lclNumRows << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<size_t> (rowIndices.size ()) < lclNumRows,
       std::invalid_argument, prefix << "rowIndices.size() = "
       << rowIndices.size () << " < this->getNodeNumRows() = "
       << lclNumRows << ".");

    // mfh 28 Aug 2017: The current local Gauss-Seidel kernel only
    // runs on host.  (See comments below.)  Thus, we need to access
    // the host versions of these data.
    const_cast<multivector_type&> (B).sync_host ();
    X.sync_host ();
    X.modify_host ();
    const_cast<multivector_type&> (D).sync_host ();

    auto B_lcl = B.template getLocalView<host_mem_space> ();
    auto X_lcl = X.template getLocalView<host_mem_space> ();
    auto D_lcl = D.template getLocalView<host_mem_space> ();

    offset_type B_stride[8], X_stride[8], D_stride[8];
    B_lcl.stride (B_stride);
    X_lcl.stride (X_stride);
    D_lcl.stride (D_stride);

    local_matrix_type lclMatrix = A->getLocalMatrix ();
    auto lclGraph = lclMatrix.graph;
    typename local_matrix_type::index_type ind = lclGraph.entries;
    typename local_matrix_type::row_map_type ptr = lclGraph.row_map;
    typename local_matrix_type::values_type val = lclMatrix.values;
    const offset_type* const ptrRaw = ptr.data ();
    const local_ordinal_type* const indRaw = ind.data ();
    const impl_scalar_type* const valRaw = val.data ();

    const std::string dir = (direction == Tpetra::Forward) ? "F" : "B";
    // NOTE (mfh 28 Aug 2017) This assumes UVM.  We can't get around
    // that on GPUs without using a GPU-based sparse triangular
    // solve to implement Gauss-Seidel, and also handling the
    // permutations correctly.
    KokkosSparse::Impl::Sequential::reorderedGaussSeidel<local_ordinal_type, offset_type, impl_scalar_type, impl_scalar_type, impl_scalar_type>(
      lclNumRows,
      numVecs,
      ptrRaw, indRaw, valRaw,
      B_lcl.data (),
      B_stride[1],
      X_lcl.data (),
      X_stride[1],
      D_lcl.data (),
      rowIndices.getRawPtr (),
      lclNumRows,
      static_cast<impl_scalar_type> (dampingFactor),
      dir.c_str ());
    const_cast<multivector_type&> (B).template sync<dev_mem_space> ();
    X.template sync<dev_mem_space> ();
    const_cast<multivector_type&> (D).template sync<dev_mem_space> ();
  }

  void CRS_gaussSeidel(
      const crs_matrix_type* A,
      const multivector_type& B,
      multivector_type& X,
      const multivector_type& D,
      const scalar_type& dampingFactor,
      const Tpetra::ESweepDirection direction,
      const int numSweeps) const
  {
    CRS_reorderedGaussSeidel (A, B, X, D, Teuchos::null, dampingFactor, direction, numSweeps);
  }

  void CRS_reorderedGaussSeidel (
      const crs_matrix_type* A,
      const multivector_type& B,
      multivector_type& X,
      const multivector_type& D,
      const Teuchos::ArrayView<local_ordinal_type>& rowIndices,
      const scalar_type& dampingFactor,
      const Tpetra::ESweepDirection direction,
      const int numSweeps) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;

    TEUCHOS_TEST_FOR_EXCEPTION(
      A->isFillComplete() == false, std::runtime_error,
      "Tpetra::CrsMatrix::gaussSeidel: cannot call this method until "
      "fillComplete() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numSweeps < 0,
      std::invalid_argument,
      "Tpetra::CrsMatrix::gaussSeidel: The number of sweeps must be , "
      "nonnegative but you provided numSweeps = " << numSweeps << " < 0.");

    // Translate from global to local sweep direction.
    // While doing this, validate the input.
    Tpetra::ESweepDirection localDirection;
    if (direction == Tpetra::Forward) {
      localDirection = Tpetra::Forward;
    }
    else if (direction == Tpetra::Backward) {
      localDirection = Tpetra::Backward;
    }
    else if (direction == Tpetra::Symmetric) {
      // We'll control local sweep direction manually.
      localDirection = Tpetra::Forward;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Tpetra::CrsMatrix::gaussSeidel: The 'direction' enum does not have "
        "any of its valid values: Forward, Backward, or Symmetric.");
    }

    if (numSweeps == 0) {
      return; // Nothing to do.
    }

    // We don't need the Export object because this method assumes
    // that the row, domain, and range Maps are the same.  We do need
    // the Import object, if there is one, though.
    RCP<const import_type> importer = A->getGraph()->getImporter();

    RCP<const map_type> domainMap = A->getDomainMap ();
    RCP<const map_type> rangeMap = A->getRangeMap ();
    RCP<const map_type> rowMap = A->getGraph ()->getRowMap ();
    RCP<const map_type> colMap = A->getGraph ()->getColMap ();

    {
      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "B be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "D be in the row Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rowMap->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the domain Map and "
        "the range Map of the matrix be the same.");
    }

    // If B is not constant stride, copy it into a constant stride
    // multivector.  We'l handle the right-hand side B first and deal
    // with X right before the sweeps, to improve locality of the
    // first sweep.  (If the problem is small enough, then that will
    // hopefully keep more of the entries of X in cache.  This
    // optimizes for the typical case of a small number of sweeps.)
    RCP<const multivector_type> B_in;
    if (B.isConstantStride()) {
      B_in = Teuchos::rcpFromRef (B);
    }
    else {
      // The range Map and row Map are the same in this case, so we
      // can use the (possibly cached) row Map multivector to store a
      // constant stride copy of B.  We don't have to copy back, since
      // Gauss-Seidel won't modify B.
      RCP<multivector_type> B_in_nonconst = getRowMapMultiVector (A, B, true);
      deep_copy (*B_in_nonconst, B); // Copy from B into B_in(_nonconst).
      B_in = rcp_const_cast<const multivector_type> (B_in_nonconst);

      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (),
        std::runtime_error,
        "gaussSeidel: The current implementation of the Gauss-Seidel kernel "
        "requires that X and B both have constant stride.  Since B does not "
        "have constant stride, we had to make a copy.  This is a limitation of "
        "the current implementation and not your fault, but we still report it "
        "as an efficiency warning for your information.");
    }

    // If X is not constant stride, copy it into a constant stride
    // multivector.  Also, make the column Map multivector X_colMap,
    // and its domain Map view X_domainMap.  (X actually must be a
    // domain Map view of a column Map multivector; exploit this, if X
    // has constant stride.)

    RCP<multivector_type> X_domainMap;
    RCP<multivector_type> X_colMap;
    bool copiedInput = false;

    if (importer.is_null ()) { // Domain and column Maps are the same.
      if (X.isConstantStride ()) {
        X_domainMap = Teuchos::rcpFromRef (X);
        X_colMap = X_domainMap;
        copiedInput = false;
      }
      else {
        // Get a temporary column Map multivector, make a domain Map
        // view of it, and copy X into the domain Map view.  We have
        // to copy here because we won't be doing Import operations.
        X_colMap = getColumnMapMultiVector (A, X, true);
        X_domainMap = X_colMap; // Domain and column Maps are the same.
        deep_copy (*X_domainMap, X); // Copy X into the domain Map view.
        copiedInput = true;
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (), std::runtime_error,
          "Tpetra::CrsMatrix::gaussSeidel: The current implementation of the "
          "Gauss-Seidel kernel requires that X and B both have constant "
          "stride.  Since X does not have constant stride, we had to make a "
          "copy.  This is a limitation of the current implementation and not "
          "your fault, but we still report it as an efficiency warning for "
          "your information.");
      }
    }
    else { // We will be doing Import operations in the sweeps.
      if (X.isConstantStride ()) {
        X_domainMap = Teuchos::rcpFromRef (X);
        // This kernel assumes that X is a domain Map view of a column
        // Map multivector.  We will only check if this is valid if
        // the CMake configure Teuchos_ENABLE_DEBUG is ON.
        X_colMap = X_domainMap->offsetViewNonConst (colMap, 0);

        // FIXME (mfh 19 Mar 2013) Do we need to fill the remote
        // entries of X_colMap with zeros?  Do we need to fill all of
        // X_domainMap initially with zeros?  Ifpack
        // (Ifpack_PointRelaxation.cpp, line 906) creates an entirely
        // new MultiVector each time.

        // Do the first Import for the first sweep.  This simplifies
        // the logic in the sweeps.
        X_colMap->doImport (X, *importer, Tpetra::INSERT);
        copiedInput = false;
      }
      else {
        // Get a temporary column Map multivector X_colMap, and make a
        // domain Map view X_domainMap of it.  Instead of copying, we
        // do an Import from X into X_domainMap.  This saves us a
        // copy, since the Import has to copy the data anyway.
        X_colMap = getColumnMapMultiVector (A, X, true);
        X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);
        X_colMap->doImport (X, *importer, Tpetra::INSERT);
        copiedInput = true;
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (), std::runtime_error,
          "Tpetra::CrsMatrix::gaussSeidel: The current implementation of the "
          "Gauss-Seidel kernel requires that X and B both have constant stride.  "
          "Since X does not have constant stride, we had to make a copy.  "
          "This is a limitation of the current implementation and not your fault, "
          "but we still report it as an efficiency warning for your information.");
      }
    }

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
      if (! importer.is_null () && sweep > 0) {
        // We already did the first Import for the zeroth sweep.
        X_colMap->doImport (*X_domainMap, *importer, Tpetra::INSERT);
      }

      // Do local Gauss-Seidel.
      if (direction != Tpetra::Symmetric) {
        if (rowIndices.is_null ()) {
          CRS_localGaussSeidel(A, *B_in, *X_colMap, D, dampingFactor, localDirection);
        }
        else {
          CRS_reorderedLocalGaussSeidel(A, *B_in, *X_colMap, D, rowIndices, dampingFactor, localDirection);
        }
      }
      else { // direction == Symmetric
        const bool doImportBetweenDirections = false;
        if (rowIndices.is_null ()) {
          CRS_localGaussSeidel(A, *B_in, *X_colMap, D, dampingFactor, Tpetra::Forward);
          // mfh 18 Mar 2013: Aztec's implementation of "symmetric
          // Gauss-Seidel" does _not_ do an Import between the forward
          // and backward sweeps.  This makes sense, because Aztec
          // considers "symmetric Gauss-Seidel" a subdomain solver.
          if (doImportBetweenDirections) {
            // Communicate again before the Backward sweep.
            if (! importer.is_null ()) {
              X_colMap->doImport (*X_domainMap, *importer, Tpetra::INSERT);
            }
          }
          CRS_localGaussSeidel(A, *B_in, *X_colMap, D, dampingFactor, Tpetra::Backward);
        }
        else {
          CRS_reorderedLocalGaussSeidel(A, *B_in, *X_colMap, D, rowIndices, dampingFactor, Tpetra::Forward);
          if (doImportBetweenDirections) {
            // Communicate again before the Backward sweep.
            if (! importer.is_null ()) {
              X_colMap->doImport (*X_domainMap, *importer, Tpetra::INSERT);
            }
          }
          CRS_reorderedLocalGaussSeidel(A, *B_in, *X_colMap, D, rowIndices, dampingFactor, Tpetra::Backward);
        }
      }
    }

    if (copiedInput) {
      deep_copy (X, *X_domainMap); // Copy back from X_domainMap to X.
    }
  }

  void CRS_gaussSeidelCopy (
      const crs_matrix_type* A,
      multivector_type& X,
      const multivector_type& B,
      const multivector_type& D,
      const scalar_type& dampingFactor,
      const Tpetra::ESweepDirection direction,
      const int numSweeps,
      const bool zeroInitialGuess) const
  {
    CRS_reorderedGaussSeidelCopy (A, X, B, D, Teuchos::null, dampingFactor, direction,
                              numSweeps, zeroInitialGuess);
  }

  void CRS_reorderedGaussSeidelCopy (
      const crs_matrix_type* A,
      multivector_type& X,
      const multivector_type& B,
      const multivector_type& D,
      const Teuchos::ArrayView<local_ordinal_type>& rowIndices,
      const scalar_type& dampingFactor,
      const Tpetra::ESweepDirection direction,
      const int numSweeps,
      const bool zeroInitialGuess) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using Teuchos::rcp_const_cast;
    typedef scalar_type Scalar;
    const char prefix[] = "Tpetra::CrsMatrix::(reordered)gaussSeidelCopy: ";
    const scalar_type ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! A->isFillComplete (), std::runtime_error,
      prefix << "The matrix is not fill complete.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numSweeps < 0, std::invalid_argument,
      prefix << "The number of sweeps must be nonnegative, "
      "but you provided numSweeps = " << numSweeps << " < 0.");

    // Translate from global to local sweep direction.
    // While doing this, validate the input.
    Tpetra::ESweepDirection localDirection;
    if (direction == Tpetra::Forward) {
      localDirection = Tpetra::Forward;
    }
    else if (direction == Tpetra::Backward) {
      localDirection = Tpetra::Backward;
    }
    else if (direction == Tpetra::Symmetric) {
      // We'll control local sweep direction manually.
      localDirection = Tpetra::Forward;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        prefix << "The 'direction' enum does not have any of its valid "
        "values: Forward, Backward, or Symmetric.");
    }

    if (numSweeps == 0) {
      return;
    }

    RCP<const import_type> importer = A->getGraph ()->getImporter ();

    RCP<const map_type> domainMap = A->getDomainMap ();
    RCP<const map_type> rangeMap = A->getRangeMap ();
    RCP<const map_type> rowMap = A->getGraph ()->getRowMap ();
    RCP<const map_type> colMap = A->getGraph ()->getColMap ();

    {
      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "B be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "D be in the row Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rowMap->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the domain Map and "
        "the range Map of the matrix be the same.");
    }

    // Fetch a (possibly cached) temporary column Map multivector
    // X_colMap, and a domain Map view X_domainMap of it.  Both have
    // constant stride by construction.  We know that the domain Map
    // must include the column Map, because our Gauss-Seidel kernel
    // requires that the row Map, domain Map, and range Map are all
    // the same, and that each process owns all of its own diagonal
    // entries of the matrix.

    RCP<multivector_type> X_colMap;
    RCP<multivector_type> X_domainMap;
    bool copyBackOutput = false;
    if (importer.is_null ()) {
      if (X.isConstantStride ()) {
        X_colMap = Teuchos::rcpFromRef (X);
        X_domainMap = Teuchos::rcpFromRef (X);
        // Column Map and domain Map are the same, so there are no
        // remote entries.  Thus, if we are not setting the initial
        // guess to zero, we don't have to worry about setting remote
        // entries to zero, even though we are not doing an Import in
        // this case.
        if (zeroInitialGuess) {
          X_colMap->putScalar (ZERO);
        }
        // No need to copy back to X at end.
      }
      else { // We must copy X into a constant stride multivector.
        // Just use the cached column Map multivector for that.
        // force=true means fill with zeros, so no need to fill
        // remote entries (not in domain Map) with zeros.
        X_colMap = getColumnMapMultiVector (A, X, true);
        // X_domainMap is always a domain Map view of the column Map
        // multivector.  In this case, the domain and column Maps are
        // the same, so X_domainMap _is_ X_colMap.
        X_domainMap = X_colMap;
        if (! zeroInitialGuess) { // Don't copy if zero initial guess
          try {
            deep_copy (*X_domainMap , X); // Copy X into constant stride MV
          } catch (std::exception& e) {
            std::ostringstream os;
            os << "Tpetra::CrsMatrix::reorderedGaussSeidelCopy: "
              "deep_copy(*X_domainMap, X) threw an exception: "
               << e.what () << ".";
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what ());
          }
        }
        copyBackOutput = true; // Don't forget to copy back at end.
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (),
          std::runtime_error,
          "gaussSeidelCopy: The current implementation of the Gauss-Seidel "
          "kernel requires that X and B both have constant stride.  Since X "
          "does not have constant stride, we had to make a copy.  This is a "
          "limitation of the current implementation and not your fault, but we "
          "still report it as an efficiency warning for your information.");
      }
    }
    else { // Column Map and domain Map are _not_ the same.
      X_colMap = getColumnMapMultiVector (A, X, true);
      X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);

#ifdef HAVE_TPETRA_DEBUG
      auto X_colMap_host_view = X_colMap->getLocalViewHost ();
      auto X_domainMap_host_view = X_domainMap->getLocalViewHost ();

      if (X_colMap->getLocalLength () != 0 && X_domainMap->getLocalLength ()) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (X_colMap_host_view.data () != X_domainMap_host_view.data (),
           std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: Pointer to "
           "start of column Map view of X is not equal to pointer to start of "
           "(domain Map view of) X.  This may mean that Tpetra::MultiVector::"
           "offsetViewNonConst is broken.  "
           "Please report this bug to the Tpetra developers.");
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_colMap_host_view.extent (0) < X_domainMap_host_view.extent (0) ||
        X_colMap->getLocalLength () < X_domainMap->getLocalLength (),
        std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
        "X_colMap has fewer local rows than X_domainMap.  "
        "X_colMap_host_view.extent(0) = " << X_colMap_host_view.extent (0)
        << ", X_domainMap_host_view.extent(0) = "
        << X_domainMap_host_view.extent (0)
        << ", X_colMap->getLocalLength() = " << X_colMap->getLocalLength ()
        << ", and X_domainMap->getLocalLength() = "
        << X_domainMap->getLocalLength ()
        << ".  This means that Tpetra::MultiVector::offsetViewNonConst "
        "is broken.  Please report this bug to the Tpetra developers.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_colMap->getNumVectors () != X_domainMap->getNumVectors (),
        std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
        "X_colMap has a different number of columns than X_domainMap.  "
        "X_colMap->getNumVectors() = " << X_colMap->getNumVectors ()
        << " != X_domainMap->getNumVectors() = "
        << X_domainMap->getNumVectors ()
        << ".  This means that Tpetra::MultiVector::offsetViewNonConst "
        "is broken.  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      if (zeroInitialGuess) {
        // No need for an Import, since we're filling with zeros.
        X_colMap->putScalar (ZERO);
      } else {
        // We could just copy X into X_domainMap.  However, that
        // wastes a copy, because the Import also does a copy (plus
        // communication).  Since the typical use case for
        // Gauss-Seidel is a small number of sweeps (2 is typical), we
        // don't want to waste that copy.  Thus, we do the Import
        // here, and skip the first Import in the first sweep.
        // Importing directly from X effects the copy into X_domainMap
        // (which is a view of X_colMap).
        X_colMap->doImport (X, *importer, Tpetra::INSERT);
      }
      copyBackOutput = true; // Don't forget to copy back at end.
    } // if column and domain Maps are (not) the same

    // The Gauss-Seidel / SOR kernel expects multivectors of constant
    // stride.  X_colMap is by construction, but B might not be.  If
    // it's not, we have to make a copy.
    RCP<const multivector_type> B_in;
    if (B.isConstantStride ()) {
      B_in = Teuchos::rcpFromRef (B);
    }
    else {
      // Range Map and row Map are the same in this case, so we can
      // use the cached row Map multivector to store a constant stride
      // copy of B.
      RCP<multivector_type> B_in_nonconst = getRowMapMultiVector (A, B, true);
      try {
        deep_copy (*B_in_nonconst, B);
      } catch (std::exception& e) {
        std::ostringstream os;
        os << "Tpetra::CrsMatrix::reorderedGaussSeidelCopy: "
          "deep_copy(*B_in_nonconst, B) threw an exception: "
           << e.what () << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what ());
      }
      B_in = Teuchos::rcp_const_cast<const multivector_type> (B_in_nonconst);

      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (),
        std::runtime_error,
        "gaussSeidelCopy: The current implementation requires that B have "
        "constant stride.  Since B does not have constant stride, we had to "
        "copy it into a separate constant-stride multivector.  This is a "
        "limitation of the current implementation and not your fault, but we "
        "still report it as an efficiency warning for your information.");
    }

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
      if (! importer.is_null () && sweep > 0) {
        // We already did the first Import for the zeroth sweep above,
        // if it was necessary.
        X_colMap->doImport (*X_domainMap, *importer, Tpetra::INSERT);
      }

      // Do local Gauss-Seidel.
      if (direction != Tpetra::Symmetric) {
        if (rowIndices.is_null ()) {
          CRS_localGaussSeidel(A, *B_in, *X_colMap, D, dampingFactor, localDirection);
        }
        else {
          CRS_reorderedLocalGaussSeidel(A, *B_in, *X_colMap, D, rowIndices, dampingFactor, localDirection);
        }
      }
      else { // direction == Tpetra::Symmetric
        if (rowIndices.is_null ()) {
          CRS_localGaussSeidel(A, *B_in, *X_colMap, D, dampingFactor, Tpetra::Forward);
          CRS_localGaussSeidel(A, *B_in, *X_colMap, D, dampingFactor, Tpetra::Backward);
        }
        else {
          CRS_reorderedLocalGaussSeidel(A, *B_in, *X_colMap, D, rowIndices, dampingFactor, Tpetra::Forward);
          CRS_reorderedLocalGaussSeidel(A, *B_in, *X_colMap, D, rowIndices, dampingFactor, Tpetra::Backward);
        }
      }
    }

    if (copyBackOutput) {
      try {
        deep_copy (X , *X_domainMap); // Copy result back into X.
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, prefix << "deep_copy(X, *X_domainMap) "
          "threw an exception: " << e.what ());
      }
    }
  }

  void BlockCRS_localGaussSeidel (
      const block_crs_matrix_type* A,
      const block_multivector_type& B,
      block_multivector_type& X,
      const Kokkos::View<impl_scalar_type***, device_type, Kokkos::MemoryUnmanaged>& D_inv,
      const scalar_type& omega,
      const Tpetra::ESweepDirection direction) const
  {
    using Kokkos::ALL;
    using little_host_vec_type = typename block_crs_matrix_type::little_host_vec_type;
    using const_little_block_type = typename block_crs_matrix_type::const_little_block_type;
    const impl_scalar_type zero =
      Kokkos::Details::ArithTraits<impl_scalar_type>::zero ();
    const impl_scalar_type one =
      Kokkos::Details::ArithTraits<impl_scalar_type>::one ();
    const local_ordinal_type numLocalMeshRows = A->getNodeNumRows();
    const local_ordinal_type numVecs = X.getNumVectors ();

    const local_ordinal_type blockSize = A->getBlockSize ();
    const size_t bs2 = blockSize * blockSize;
    Teuchos::Array<impl_scalar_type> localMem (blockSize);
    //X_lcl is scratch for computing the row-vec product in GS
    little_host_vec_type X_lcl (localMem.getRawPtr (), blockSize);

    if (direction == Tpetra::Symmetric) {
      BlockCRS_localGaussSeidel (A, B, X, D_inv, omega, Tpetra::Forward);
      BlockCRS_localGaussSeidel (A, B, X, D_inv, omega, Tpetra::Backward);
      return;
    }

    bool forward = direction == Tpetra::Forward;

    const scalar_type one_minus_omega = Teuchos::ScalarTraits<scalar_type>::one()-omega;
    const scalar_type minus_omega = -omega;

    Kokkos::fence();

    const crs_graph_type& graph = A->getCrsGraph();
    auto localGraph = graph.getLocalGraph();
    const_cast<block_crs_matrix_type*>(A)->sync_host();
    auto rowmap = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localGraph.row_map);
    auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localGraph.entries);
    auto values = A->getValuesHost();

    if (numVecs == 1) {
      for (local_ordinal_type i = 0; i < numLocalMeshRows; i++) {
        const local_ordinal_type actlRow = forward ? i : numLocalMeshRows - i - 1;

        little_host_vec_type B_cur = B.getLocalBlock (actlRow, 0);
        Tpetra::COPY (B_cur, X_lcl);
        Tpetra::SCAL (static_cast<impl_scalar_type> (omega), X_lcl);

        const size_t meshBeg = rowmap[actlRow];
        const size_t meshEnd = rowmap[actlRow+1];
        for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
          const local_ordinal_type meshCol = entries[absBlkOff];
          const_little_block_type A_cur(values.data() + bs2 * absBlkOff, blockSize, blockSize);
          little_host_vec_type X_cur = X.getLocalBlock (meshCol, 0);

          // X_lcl += alpha*A_cur*X_cur
          const scalar_type alpha = meshCol == actlRow ? one_minus_omega : minus_omega;
          //X_lcl.matvecUpdate (alpha, A_cur, X_cur);
          Tpetra::GEMV (static_cast<impl_scalar_type> (alpha), A_cur, X_cur, X_lcl);
        } // for each entry in the current local row of the matrix

        // NOTE (mfh 20 Jan 2016) The two input Views here are
        // unmanaged already, so we don't have to take unmanaged
        // subviews first.
        auto D_lcl = Kokkos::subview (blockDiag_, actlRow, ALL (), ALL ());
        little_host_vec_type X_update = X.getLocalBlock (actlRow, 0);
        Tpetra::FILL (X_update, zero);
        Tpetra::GEMV (one, D_lcl, X_lcl, X_update); // overwrite X_update
      } // for each local row of the matrix
    }
    else {
      for (local_ordinal_type i = 0; i < numLocalMeshRows; i++) {
        const local_ordinal_type actlRow = forward ? i : numLocalMeshRows - i - 1;
        for (local_ordinal_type j = 0; j < numVecs; ++j) {

          little_host_vec_type B_cur = B.getLocalBlock (actlRow, j);
          Tpetra::COPY (B_cur, X_lcl);
          Tpetra::SCAL (static_cast<impl_scalar_type> (omega), X_lcl);

          const size_t meshBeg = rowmap[actlRow];
          const size_t meshEnd = rowmap[actlRow+1];
          for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
            const local_ordinal_type meshCol = entries[absBlkOff];
            const_little_block_type A_cur(values.data() + bs2 * absBlkOff, blockSize, blockSize);
            little_host_vec_type X_cur = X.getLocalBlock (meshCol, j);

            // X_lcl += alpha*A_cur*X_cur
            const scalar_type alpha = meshCol == actlRow ? one_minus_omega : minus_omega;
            Tpetra::GEMV (static_cast<impl_scalar_type> (alpha), A_cur, X_cur, X_lcl);
          } // for each entry in the current local row of the matrx

          auto D_lcl = Kokkos::subview (blockDiag_, actlRow, ALL (), ALL ());
          auto X_update = X.getLocalBlock (actlRow, j);
          Tpetra::FILL (X_update, zero);
          Tpetra::GEMV (one, D_lcl, X_lcl, X_update); // overwrite X_update
        } // for each entry in the current local row of the matrix
      } // for each local row of the matrix
    }
  }


}; //class Relaxation

}//namespace Ifpack2

#endif // IFPACK2_RELAXATION_DECL_HPP
