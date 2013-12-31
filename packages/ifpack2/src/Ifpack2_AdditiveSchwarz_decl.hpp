/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_ADDITIVESCHWARZ_DECL_HPP
#define IFPACK2_ADDITIVESCHWARZ_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"

#include "Ifpack2_LocalFilter_decl.hpp"
#include "Ifpack2_OverlappingRowMatrix_decl.hpp"
#include "Ifpack2_ReorderFilter_decl.hpp"
#include "Ifpack2_SingletonFilter_decl.hpp"

#include "Ifpack2_Details_NestedPreconditioner.hpp"


namespace Ifpack2 {

/** \class AdditiveSchwarz
\brief Additive Schwarz domain decomposition for Tpetra sparse matrices
\tparam MatrixType A specialization of Tpetra::RowMatrix
\tparam LocalInverseType The type of the solver for the local
  (subdomain) problem.  This must be a specialization of a concrete
  subclass of Ifpack2::Preconditioner.  (This template parameter will
  be DEPRECATED soon; refer to discussion below.)

\section Ifpack2_AdditiveSchwarz_Summary Summary

This class implements Additive Schwarz domain decomposition, with
optional overlap.  It operates on a given Tpetra::RowMatrix.  Each
subdomain corresponds to exactly one MPI process in the given matrix's
MPI communicator.

This class implements Tpetra::Operator, like all other subclasses of
Preconditioner.  Thus, the apply() method applies the preconditioner
to a multivector.

\section Ifpack2_AdditiveSchwarz_Alg Algorithm

One-level overlapping domain decomposition preconditioners use local
solvers of Dirichlet type. This means that the inverse of the local
matrix (possibly with overlap) is applied to the residual to be
preconditioned.

The preconditioner can be written as:
\f[
P_{AS}^{-1} = \sum_{i=1}^M P_i A_i^{-1} R_i,
\f]
where \f$M\f$ is the number of subdomains (in this case, the number of
(MPI) processes in the computation), \f$R_i\f$ is an operator that
restricts the global vector to the vector lying on subdomain \f$i\f$,
\f$P_i\f$ is the prolongator operator, and
\f[
A_i = R_i A P_i.
\f]

Constructing a Schwarz preconditioner takes two steps:
<ol>
<li> Definition of the restriction and prolongation operators </li>
<li> Definition of a solver for linear systems involving \f$A_i\f$ </li>
</ol>

The definition of the restriction and prolongation operators \f$R_i\f$
and \f$R_i^T\f$ depends on the level of overlap.  If the overlap level
is zero (no overlap), their implementation is trivial; \f$R_i\f$ will
return all the local components.  For nonzero overlap, Tpetra's data
redistribution facilities (Tpetra::Import) will be used to bring in
the required data.  Users may control how these data are combined with
existing data, by setting the combine mode parameter.  Valid combine
modes include "ADD", "INSERT", "REPLACE", "ABSMAX", and "ZERO".  These
correspond to the valid values of Tpetra::CombineMode.

To solve linear systems involving \f$A_i\f$ on each subdomain, the
user can adopt any subclass of Preconditioner.  Currently, users
control this at compile time by setting the second template parameter
\c LocalInverseType.  Soon, this option will be removed in favor of
run-time control of the subdomain solver.

The local matrix \f$A_i\f$ can be filtered, to eliminate singletons,
and reordered. At the present time, RCM and METIS can be used to
reorder the local matrix.

\section Ifpack2_AdditiveSchwarz_CombineMode Combine modes

AdditiveSchwarz supports different <i>combine modes</i>.  These are
rules for combining overlapping entries from the subdomain solves'
results, to form the final additive Schwarz result.  The two modes
that users will likely find most useful are "Add" and "Zero".  "Add"
sums overlapping results, and "Zero" ignores them.

The Zero combine mode can be a bit confusing.  It helps to compare the
Add combine mode without overlap, to the Zero combine mode with overlap.
Consider the following 4 x 4 linear system \f$Ax = b\f$, where
\f[
A =
\left(
\begin{array}{rrrr}
  2 & -1 &  0 &  0 \\
 -1 &  2 & -1 &  0 \\
  0 & -1 &  2 & -1 \\
  0 &  0 & -1 &  2
\end{array} \right)
\f]
and
\f[
b =
\left(
\begin{array}{rrrr}
  1 \\
  4 \\
  9 \\
 16
\end{array} \right)
Suppose that we give the first two rows of A and b to Process 0, and
the last two rows of A and b to Process 1.

If we use additive Schwarz without overlap, and use the (default) Add
combine mode, then each process must solve a linear system with the
following 2 x 2 matrix:
\f[
A_{local} =
\begin{array}{rr}
  2 & -1 \\
 -1 &  2
\end{array} \right).
\f]
The inverse of this matrix is
\f[
A_{local}^{-1} =
\frac{1}{3} \cdot
\begin{array}{rr}
  2 &  1 \\
  1 &  2
\end{array} \right).
\f]

Process 0 gets the right-hand side \f$(1, 4)^T\f$ and thus the local
solution \f$(2, 3)^T\f$.  Process 1 gets the right-hand side \f$(9,
16)^T\f$ and thus the local solution \f$(34/3, 41/3)^T\f$.  The Add
combine mode sums entries corresponding to overlap, but in this case
there is no overlap, so we just concatenate the local results.  Thus,
the result of applying additive Schwarz once is
\f[
x_{Add} =
\begin{array}{r}
   2 \\
   3 \\
34/3 \\
41/3
\end{array}.
\f]

If we introduce one level of overlap on each of these processes, and
use the Zero combine mode with additive Schwarz, then each process has
to solve a linear system with the following 3 x 3 matrix:
\f[
A_{local} =
\begin{array}{rrr}
  2 & -1 &  0 \\
 -1 &  2 & -1 \\
  0 & -1 &  2
\end{array} \right).
The inverse of this matrix is
\f[
A_{local}^{-1} =
\frac{1}{4} \cdot
\begin{array}{rrr}
  3 &  2 &  1 \\
  2 &  4 &  2 \\
  1 &  2 &  3
\end{array} \right).
\f]
Process 0 gets the right-hand side \f$(1, 4, 9)^T\f$ and thus the
local solution \f$(5, 9, 9)^T\f$.  Process 1 gets the right-hand side
\f$(4, 9, 16)^T\f$ and thus the local solution \f$(23/2, 19,
35/2)^T\f$.  The Zero combine mode ignores "remote entries," so that
the result of applying additive Schwarz once is
\f[
x_{Zero} =
\begin{array}{r}
   5 \\
   9 \\
  19 \\
35/2
\end{array}.
\f]

Even though both of the above examples combine the local results in
the same way, Zero produces a different final result than Add, because
the Zero example uses overlap 1, whereas the Add example uses overlap 0.

\section Ifpack2_AdditiveSchwarz_DevNotes Notes to Ifpack2 developers

The second template parameter (\c LocalInverseType) is being
DEPRECATED.  It causes a lot of trouble for explicit template
instantiation, and we can perfectly well support an arbitrary
subdomain solver type by run-time polymorphism.  For example, a
subclass of Preconditioner could expose a scalar_type (of the input
and output vectors of apply()) different than its internal storage
type, via a mechanism like that of Tpetra::ApplyOp.  The only issue is
that the preconditioner would need to be creatable by Factory, but we
could get around that by adding a method to AdditiveSchwarz that lets
the user supply an arbitrary Preconditioner subclass instance for the
subdomain solver.
*/
template<class MatrixType, class LocalInverseType>
class AdditiveSchwarz :
    virtual public Preconditioner<typename MatrixType::scalar_type,
                                  typename MatrixType::local_ordinal_type,
                                  typename MatrixType::global_ordinal_type,
                                  typename MatrixType::node_type>,
    virtual public Details::NestedPreconditioner<Preconditioner<typename MatrixType::scalar_type,
                                                                typename MatrixType::local_ordinal_type,
                                                                typename MatrixType::global_ordinal_type,
                                                                typename MatrixType::node_type> >
{
public:
  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  /// \brief The Tpetra::RowMatrix specialization matching MatrixType.
  ///
  /// MatrixType may be either a Tpetra::CrsMatrix specialization or a
  /// Tpetra::RowMatrix specialization.  This typedef will always be a
  /// Tpetra::RowMatrix specialization which is either the same as
  /// MatrixType, or the parent class of MatrixType.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  //@}
  // \name Deprecated typedefs
  //@{

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type         Scalar;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type  LocalOrdinal;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type           Node;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //! Preserved only for backwards compatibility.  Please use "row_matrix_type".
  TEUCHOS_DEPRECATED typedef typename Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type>  LocalMatrixType;

  //@}
  // \name Constructors and destructor
  //@{

  /// \brief Constructor that takes a matrix.
  ///
  /// \param A [in] The matrix to be preconditioned.
  AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A);

  /// \brief Constructor that takes a matrix and the level of overlap.
  ///
  /// \warning This version of the constructor is DEPRECATED, because
  ///   the single-argument version suffices; users may specify the
  ///   overlap level via the "schwarz: overlap level" parameter.
  ///
  /// \param A [in] The matrix to be preconditioned.
  /// \param overlapLevel [in] The level of overlap.  Must be
  ///   nonnegative.  Zero means no overlap.
  AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                   const int overlapLevel);

  //! Destructor
  virtual ~AdditiveSchwarz();

  //@}
  //! \name Implementation of Tpetra::Operator
  //@{

  //! The domain Map of this operator.
  virtual Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getDomainMap() const;

  //! The range Map of this operator.
  virtual Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getRangeMap() const;

  //! Apply the preconditioner to X, putting the result in Y.
  virtual void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //@}

  /// \brief Set the inner preconditioner.
  ///
  /// Most users do not need to call this method.  Factory will handle
  /// calling this method if necessary.  One case in which users might
  /// need to call this method, is if they implement their own
  /// subdomain solver as a subclass of Preconditioner, and want to
  /// use that subdomain solver in this class.
  ///
  /// \param innerPrec [in/out] The inner preconditioner.  Its matrix
  ///   (if it has one) may be replaced by a matrix specified by the
  ///   outer (this) preconditioner.
  ///
  /// \warning CIRCULAR DEPENDENCIES ARE FORBIDDEN.  You may NOT give
  ///   this object (<tt>*this</tt>) to itself as an inner solver.
  ///   You MAY use an inner solver of the same TYPE as this object
  ///   (as long as this makes sense mathematically), but it must be a
  ///   different instance.
  ///
  /// \pre <tt>&*innerPrec != this</tt>.
  virtual void
  setInnerPreconditioner (const Teuchos::RCP<Preconditioner<scalar_type,
                                                            local_ordinal_type,
                                                            global_ordinal_type,
                                                            node_type> >& innerPrec);

  //! The input matrix.
  virtual Teuchos::RCP<const row_matrix_type> getMatrix() const;

  /// \brief Set the preconditioner's parameters.
  ///
  /// This version of the method takes a const list, as required by
  /// the Preconditioner interface.  setParameterList() takes a
  /// nonconst pointer to a list, in order to match the
  /// Teuchos::ParameterListAcceptor interface.
  ///
  /// Accepted parameters include the following:
  ///   - "schwarz: compute condest" (\c bool): If true, estimate the
  ///     condition number each time compute() is called.
  ///   - "schwarz: combine mode" (\c std::string): The rule for
  ///     combining incoming data with existing data in overlap
  ///     regions.  Valid values include "ADD", "INSERT", "REPLACE",
  ///     "ABSMAX", and "ZERO".  These correspond to the valid values
  ///     of Tpetra::CombineMode.
  ///   - "schwarz: overlap level" (\c int): The level of overlap.
  ///   - "schwarz: use reordering" (\c bool): Whether to use Zoltan2
  ///     to do reordering.  If true, then Trilinos must have been
  ///     built with Zoltan2 and Xpetra enabled.
  ///   - "schwarz: subdomain id" (\c int): This option does not
  ///     currently work.
  ///   - "schwarz: filter singletons" (\c bool): If true, exclude
  ///     rows with just a single entry on the calling process.
  virtual void setParameters (const Teuchos::ParameterList& List);

  /// \brief Set the preconditioner's parameters.
  ///
  /// This version of the method takes a nonconst pointer to a list,
  /// in order to match the Teuchos::ParameterListAcceptor interface.
  /// setParameters() takes a const list, as required by the
  /// Preconditioner interface.
  ///
  /// Accepted parameters include the following:
  ///   - "schwarz: compute condest" (\c bool): If true, estimate the
  ///     condition number each time compute() is called.
  ///   - "schwarz: combine mode" (\c std::string): The rule for
  ///     combining incoming data with existing data in overlap
  ///     regions.  Valid values include "ADD", "INSERT", "REPLACE",
  ///     "ABSMAX", and "ZERO".  These correspond to the valid values
  ///     of Tpetra::CombineMode.
  ///   - "schwarz: overlap level" (\c int): The level of overlap.
  ///   - "schwarz: use reordering" (\c bool): Whether to use Zoltan2
  ///     to do reordering.  If true, then Trilinos must have been
  ///     built with Zoltan2 and Xpetra enabled.
  ///   - "schwarz: subdomain id" (\c int): This option does not
  ///     currently work.
  ///   - "schwarz: filter singletons" (\c bool): If true, exclude
  ///     rows with just a single entry on the calling process.
  void
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist);

  /// \brief Get a list of the preconditioner's default parameters.
  ///
  /// See the documentation of setParameters() for a list of the
  /// parameters that AdditiveSchwarz accepts.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const;

  //! Computes all (graph-related) data necessary to initialize the preconditioner.
  virtual void initialize();

  //! Returns true if the  preconditioner has been successfully initialized, false otherwise.
  virtual bool isInitialized() const;

  //! Computes all (coefficient) data necessary to apply the preconditioner.
  virtual void compute();

  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  virtual bool isComputed() const;

  /// \brief Compute the condition number estimate and return its value.
  ///
  /// \warning This method is DEPRECATED.  It was inherited from
  ///   Ifpack, and Ifpack never clearly stated what this method
  ///   computes.  Furthermore, Ifpack's method just estimates the
  ///   condition number of the matrix A, and ignores the
  ///   preconditioner -- which is probably not what users thought it
  ///   did.  If there is sufficient interest, we might reintroduce
  ///   this method with a different meaning and a better algorithm.
  virtual magnitude_type TEUCHOS_DEPRECATED
  computeCondEst (CondestType CT = Ifpack2::Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const row_matrix_type> &Matrix = Teuchos::null);

  /// \brief Return the computed condition number estimate, or -1 if not computed.
  ///
  /// \warning This method is DEPRECATED.  See warning for computeCondEst().
  virtual magnitude_type TEUCHOS_DEPRECATED getCondEst() const;

  //! Returns the number of calls to initialize().
  virtual int getNumInitialize() const;

  //! Returns the number of calls to compute().
  virtual int getNumCompute() const;

  //! Returns the number of calls to apply().
  virtual int getNumApply() const;

  //! Returns the time spent in initialize().
  virtual double getInitializeTime() const;

  //! Returns the time spent in compute().
  virtual double getComputeTime() const;

  //! Returns the time spent in apply().
  virtual double getApplyTime() const;

  //! \name Implementation of Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& print(std::ostream& os) const;

  //! Returns the level of overlap.
  virtual int getOverlapLevel() const;


private:
  //! Specialization of Tpetra::Map.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  //! Specialization of Tpetra::Import.
  typedef Tpetra::Import<local_ordinal_type,
                         global_ordinal_type,
                         node_type> import_type;
  //! Specialization of Tpetra::MultiVector.
  typedef Tpetra::MultiVector<scalar_type,
                              local_ordinal_type,
                              global_ordinal_type,
                              node_type> MV;
  //! Specialization of Preconditioner.
  typedef Preconditioner<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> prec_type;

  //! Copy constructor (unimplemented; do not use)
  AdditiveSchwarz (const AdditiveSchwarz& RHS);

  //! Set up the localized matrix and the singleton filter.
  void setup ();

  //! The current inner preconditioner name.
  std::string innerPrecName () const;

  //! The default inner preconditioner name.
  static std::string defaultInnerPrecName ();

  /// \brief The matrix to be preconditioned.
  ///
  /// This is the same matrix as the input argument to this class' constructor.
  Teuchos::RCP<const row_matrix_type> Matrix_;

  //! The overlapping matrix.
  Teuchos::RCP<OverlappingRowMatrix<row_matrix_type> > OverlappingMatrix_;

  //! The reordered matrix.
  Teuchos::RCP<ReorderFilter<row_matrix_type> > ReorderedLocalizedMatrix_;

  //! The "inner matrix" given to the inner (subdomain) solver.
  Teuchos::RCP<row_matrix_type> innerMatrix_;

  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! If true, overlapping is used
  bool IsOverlapping_;
  //! Level of overlap among the processors.
  int OverlapLevel_;

  //! Store a copy of the list given in setParameters()
  Teuchos::ParameterList List_;
  //! Valid (default) parameters; computed and cached in getValidParameters().
  mutable Teuchos::RCP<const Teuchos::ParameterList> validParams_;

  //! Combine mode for off-process elements (only if overlap is used)
  Tpetra::CombineMode CombineMode_;
  //! Contains the estimated condition number.
  magnitude_type Condest_;
  //! If \c true, compute the condition number estimate each time Compute() is called.
  bool ComputeCondest_;
  //! If \c true, reorder the local matrix.
  bool UseReordering_;
  //! Record reordering for output purposes.
  std::string ReorderingAlgorithm_;
  //! If true, subdomain filtering is used
  bool UseSubdomain_;
  //! Whether to filter singleton rows.
  bool FilterSingletons_;
  //! Matrix from which singleton rows have been filtered.
  Teuchos::RCP<SingletonFilter<row_matrix_type> > SingletonMatrix_;

  //! The total number of successful calls to initialize().
  int NumInitialize_;
  //! The total number of successful calls to compute().
  int NumCompute_;
  //! The total number of successful calls to apply().
  mutable int NumApply_;
  //! The total time in seconds of all successful calls to initialize().
  double InitializeTime_;
  //! The total time in seconds of all successful calls to compute().
  double ComputeTime_;
  //! The total time in seconds of all successful calls to apply().
  mutable double ApplyTime_;
  //! The total number of floating-point operations over all initialize() calls.
  double InitializeFlops_;
  //! The total number of floating-point operations over all compute() calls.
  double ComputeFlops_;
  //! The total number of floating-point operations over all apply() calls.
  mutable double ApplyFlops_;
  //! The inner (that is, per subdomain local) solver.
  Teuchos::RCP<prec_type> Inverse_;
  //! SerialMap for filtering multivector with no overlap.
  Teuchos::RCP<const map_type> SerialMap_;
  //! Distributed map for filtering multivector with no overlap.
  Teuchos::RCP<const map_type> DistributedMap_;
  //! Local distributed map for filtering multivector with no overlap.
  Teuchos::RCP<const map_type> LocalDistributedMap_;

  /// \brief Import object used in apply().
  ///
  /// If the domain decomposition is <i>not</i> overlapping, then
  /// apply() needs to redistribute the X input vector from the domain
  /// Map of this operator to SerialMap_ (see above).  Computing the
  /// Import is expensive, and the Import does not change as long as
  /// the overlap level does not change, so it makes sense to keep the
  /// Import object around.  apply() creates this on demand if
  /// necessary, which explains why this is marked \c mutable.
  mutable Teuchos::RCP<const import_type> SerialImporter_;

  /// \brief Import object used in apply().
  ///
  /// Import object from the domain Map of this operator to
  /// DistributedMap_ (see above).  This is used in apply(), where it
  /// is created on demand if necessary.  This explains why this is
  /// marked \c mutable.
  mutable Teuchos::RCP<const import_type> DistributedImporter_;
}; // class AdditiveSchwarz

}// end namespace

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
