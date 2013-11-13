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
#include "Tpetra_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Ifpack2_ReorderFilter.hpp"
#include "Ifpack2_SingletonFilter.hpp"
#include "Ifpack2_OverlappingRowMatrix.hpp"

#ifdef HAVE_IFPACK2_ZOLTAN2
#include "Zoltan2_OrderingSolution.hpp"
#endif


namespace Ifpack2 {

/** \class AdditiveSchwarz
\brief Additive Schwarz domain decomposition for Tpetra sparse matrices.

\section Ifpack2_AdditiveSchwarz_Summary Summary

This class implements an Additive Schwarz (one-level overlapping
domain decomposition) preconditioner.  It operates on a given
Tpetra::RowMatrix.  This class implements Tpetra::Operator, like all
other subclasses of Preconditioner.  Thus, the apply() method applies
the preconditioner to a multivector.

\section Ifpack2_AdditiveSchwarz_Alg Algorithm

One-level overlapping domain decomposition preconditioners use local
solvers of Dirichlet type. This means that the inverse of the local
matrix (with minimal or wider overlap) is applied to the residual to
be preconditioned.

The preconditioner can be written as:
\f[
P_{AS}^{-1} = \sum_{i=1}^M P_i A_i^{-1} R_i,
\f]
where \f$M\f$ is the number of subdomains (in this case, the number of
processors in the computation), \f$R_i\f$ is an operator that
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
and \f$R_i^T\f$ depends on the level of overlap. If minimal overlap is
chosen, their implementation is trivial; \f$R_i\f$ will return all the
local components.  For wider overlap, Tpetra::Import and
Tpetra::Export will be used to import resp. export data.  The user
must provide both the matrix to be preconditioned (whose which must
have minimal overlap) and the matrix with wider overlap.

To solve linear systems involving \f$A_i\f$ on each subdomain, the
user can adopt any subclass of Preconditioner. This can be easily
accomplished, as AdditiveSchwarz is templated with the solver for each
subdomain.

The local matrix \f$A_i\f$ can be filtered, to eliminate singletons,
and reordered. At the present time, RCM and METIS can be used to
reorder the local matrix.

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
template<class MatrixType,class LocalInverseType>
class AdditiveSchwarz :
    virtual public Preconditioner<typename MatrixType::scalar_type,
                                  typename MatrixType::local_ordinal_type,
                                  typename MatrixType::global_ordinal_type,
                                  typename MatrixType::node_type> {
public:
  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type         scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type  local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type           node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  /// \brief The Tpetra::RowMatrix specialization matching MatrixType.
  ///
  /// MatrixType may be either a Tpetra::CrsMatrix specialization or a
  /// Tpetra::RowMatrix specialization.  This typedef will always be a
  /// Tpetra::RowMatrix specialization which is either the same as
  /// MatrixType, or the parent class of MatrixType.
  typedef typename Tpetra::RowMatrix<scalar_type,
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

  //! The input matrix.
  virtual Teuchos::RCP<const row_matrix_type> getMatrix() const;

  /// \brief Set the preconditioner's parameters.
  ///
  /// Accepted parameters include the following:
  ///   - "schwarz: compute condest" (\c bool): If true, estimate the
  ///     condition number each time compute() is called.
  ///   - "schwarz: combine mode" (\c std::string): The (Tpetra)
  ///     CombineMode used for combining incoming data with existing
  ///     data in overlap regions.  Valid values include "Add",
  ///     "Insert", "Replace", and "AbsMax".
  ///   - "schwarz: overlap level" (\c int): The level of overlap.
  ///   - "schwarz: use reordering" (\c bool): Whether to use Zoltan2
  ///     to do reordering.  If true, then Trilinos must have been
  ///     built with Zoltan2 and Xpetra enabled.
  ///   - "schwarz: subdomain id" (\c int): This option does not
  ///     currently work.
  ///   - "schwarz: filter singletons" (\c bool): If true, exclude
  ///     rows with just a single entry on the calling process.
  virtual void setParameters (const Teuchos::ParameterList& List);

  //! Computes all (graph-related) data necessary to initialize the preconditioner.
  virtual void initialize();

  //! Returns true if the  preconditioner has been successfully initialized, false otherwise.
  virtual bool isInitialized() const;

  //! Computes all (coefficient) data necessary to apply the preconditioner.
  virtual void compute();

  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  virtual bool isComputed() const;

  //! Computes the condition number estimate and returns its value.
  virtual magnitude_type
  computeCondEst (CondestType CT = Ifpack2::Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const row_matrix_type> &Matrix = Teuchos::null);

  //! Returns the computed condition number estimate, or -1.0 if not computed.
  virtual magnitude_type getCondEst() const;

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
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> map_type;

protected:
  //! Copy constructor (unimplemented; do not use)
  AdditiveSchwarz (const AdditiveSchwarz& RHS);

  //! Set up the localized matrix and the singleton filter.
  void setup ();

  /// \brief The matrix to be preconditioned.
  ///
  /// This is the same matrix as the input argument to this class' constructor.
  Teuchos::RCP<const row_matrix_type> Matrix_;

  //! The overlapping matrix.
  Teuchos::RCP<OverlappingRowMatrix<row_matrix_type> > OverlappingMatrix_;

  //! Localized version of Matrix_ or OverlappingMatrix_.
  // CMS: Probably not a great idea, but this will remove the Local/Subdomain conflict here
  Teuchos::RCP<row_matrix_type> LocalizedMatrix_;
  //! The reordered matrix.
  Teuchos::RCP<ReorderFilter<row_matrix_type> > ReorderedLocalizedMatrix_;

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
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to apply().
  mutable int NumApply_;
  //! Contains the time for all successful calls to initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to apply().
  mutable double ApplyTime_;
  //! Contains the number of flops for Initialize().
  double InitializeFlops_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyFlops_;
  //! Object used for timing purposes.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! Pointer to the local solver.
  Teuchos::RCP<LocalInverseType> Inverse_;
  //! SerialMap for filtering multivector with no overlap.
  Teuchos::RCP<const map_type> SerialMap_;
  //! Distributed map for filtering multivector with no overlap.
  Teuchos::RCP<const map_type> DistributedMap_;
  //! Local distributed map for filtering multivector with no overlap.
  Teuchos::RCP<const map_type> LocalDistributedMap_;
}; // class AdditiveSchwarz

}// end namespace

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
