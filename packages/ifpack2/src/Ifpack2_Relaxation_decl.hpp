// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "Ifpack2_Details_GaussSeidel.hpp"
#include "Ifpack2_Details_InverseDiagonalKernel.hpp"

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
spent a lot of effort tuning Tpetra's kernels.  Depending on the Node
type of your Tpetra matrix, it may also exploit threads for additional
parallelism within each MPI process.  In contrast, Gauss-Seidel and
symmetric Gauss-Seidel are intrinsically sequential methods within an
MPI process.  This prevents us from exposing more parallelism via
threads.  The difference should become more apparent as your code
moves away from a "one MPI process per core" model, to a "one MPI
process per socket or node" model, assuming that you are using a
thread-parallel Node type.

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
if the matrix data are complex).

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

  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The Kokkos device type used by the input MatrixType.
  typedef typename MatrixType::node_type::device_type device_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Tpetra::RowMatrix specialization used by this class.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;

  //! Tpetra::Operator specialization used by this class.
  typedef Tpetra::Operator<scalar_type,
                           local_ordinal_type,
                           global_ordinal_type,
                           node_type> op_type;

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

  bool supportsZeroStartingSolution() { return true; }

  void setZeroStartingSolution (bool zeroStartingSolution) { ZeroStartingSolution_ = zeroStartingSolution; };

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

  typedef typename crs_matrix_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename crs_matrix_type::nonconst_values_host_view_type nonconst_values_host_view_type;

  Teuchos::RCP<Ifpack2::Details::InverseDiagonalKernel<op_type> > invDiagKernel_;

  //@}
  //! \name Implementation of multithreaded Gauss-Seidel.
  //@{

  typedef typename crs_matrix_type::local_matrix_device_type local_matrix_device_type;
  typedef typename local_matrix_device_type::StaticCrsGraphType::row_map_type lno_row_view_t;
  typedef typename local_matrix_device_type::StaticCrsGraphType::entries_type lno_nonzero_view_t;
  typedef typename local_matrix_device_type::values_type scalar_nonzero_view_t;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type TemporaryWorkSpace;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type PersistentWorkSpace;
  typedef typename local_matrix_device_type::StaticCrsGraphType::execution_space MyExecSpace;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle
      <typename lno_row_view_t::const_value_type, local_ordinal_type,typename scalar_nonzero_view_t::value_type,
      MyExecSpace, TemporaryWorkSpace,PersistentWorkSpace > mt_kernel_handle_type;
  Teuchos::RCP<mt_kernel_handle_type> mtKernelHandle_;

  //@}
  //! \name Implementation of serial Gauss-Seidel.
  //@{

  using SerialGaussSeidel = Ifpack2::Details::GaussSeidel<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  Teuchos::RCP<SerialGaussSeidel> serialGaussSeidel_;

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

  //! Apply multi-threaded Gauss-Seidel for solving AX = B.
  void ApplyInverseMTGS_CrsMatrix(
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& B,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
      Tpetra::ESweepDirection direction) const;

  //! Apply Gauss-Seidel to X, returning the result in Y. Calls one of the Crs/Row/BlockCrs specializations below.
  void ApplyInverseSerialGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
              Tpetra::ESweepDirection direction) const;

  //! Apply Gauss-Seidel for a Tpetra::CrsMatrix specialization.
  void ApplyInverseSerialGS_CrsMatrix (const crs_matrix_type& A,
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& B,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
      Tpetra::ESweepDirection direction) const;

  //! Apply Gauss-Seidel for a Tpetra::RowMatrix specialization.
  void ApplyInverseSerialGS_RowMatrix(
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
      Tpetra::ESweepDirection direction) const;

  //! Apply Gauss-Seidel for a Tpetra::BlockCrsMatrix specialization.
  void ApplyInverseSerialGS_BlockCrsMatrix(
      const block_crs_matrix_type& A,
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
      Tpetra::ESweepDirection direction);

  //! Apply symmetric multi-threaded Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseMTSGS_CrsMatrix(
          const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  void MTGaussSeidel (
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& B,
      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
      const Tpetra::ESweepDirection direction) const;

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
  //! For MTSGS, the cluster size (use point coloring if equal to 1)
  int clusterSize_ = 1;
  //! For MTSGS, the threshold for long/bulk rows (rows with at least this many nonzeros)
  int longRowThreshold_ = 0;
  //! For MTSGS, the coloring algorithm to use (as a string)
  KokkosGraph::ColoringAlgorithm mtColoringAlgorithm_ = KokkosGraph::COLORING_DEFAULT;

  //! Number of outer-sweeps for the two-stage Gauss Seidel
  int NumOuterSweeps_ = 1;
  //! Number of inner-sweeps for the two-stage Gauss Seidel
  int NumInnerSweeps_ = 1;
  //! Whether to use sparse-triangular solve instead of inner-iterations
  bool InnerSpTrsv_ = false;
  //! Damping factor for inner-sweeps
  scalar_type InnerDampingFactor_ = STS::one();
  //! Whether to use compact form of recurrence for the two-stage Gauss Seidel
  bool CompactForm_ = false;

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
  //! Whether apply should register and use a timer
  bool TimerForApply_ = true;
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
}; //class Relaxation

}//namespace Ifpack2

#endif // IFPACK2_RELAXATION_DECL_HPP
