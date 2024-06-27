// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_Hiptmair_decl.hpp
/// \brief Declaration of Ifpack2::Hiptmair class.
/// \author Paul Tsuji

#ifndef IFPACK2_HIPTMAIR_DECL_HPP
#define IFPACK2_HIPTMAIR_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Tpetra_Map_fwd.hpp"
#include <type_traits>

namespace Teuchos {
  class ParameterList; // forward declaration
}

namespace Ifpack2 {

  /// \class Hiptmair
  /// \brief Wrapper for Hiptmair smoothers.
  /// \tparam A specialization of Tpetra::RowMatrix.
  ///
  /// Ifpack2::Hiptmair does smoothing on two spaces; a primary space
  /// and an auxiliary space. This situation arises when preconditioning
  /// Maxwell's equations discretized by edge elements.
  ///
  /// For a list of all run-time parameters that this class accepts,
  /// see the documentation of setParameters().
  template<class MatrixType>
  class Hiptmair :
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type>
  {
  public:
    // \name Public typedefs
    //@{

    //! The type of the entries of the input MatrixType.
    typedef typename MatrixType::scalar_type scalar_type;

    //! The type of local indices in the input MatrixType.
    typedef typename MatrixType::local_ordinal_type local_ordinal_type;

    //! The type of global indices in the input MatrixType.
    typedef typename MatrixType::global_ordinal_type global_ordinal_type;

    //! The Node type used by the input MatrixType.
    typedef typename MatrixType::node_type node_type;

    //! The type of the magnitude (absolute value) of a matrix entry.
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

    //! Type of the Tpetra::RowMatrix specialization that this class uses.
    typedef Tpetra::RowMatrix<scalar_type,
                              local_ordinal_type,
                              global_ordinal_type,
                              node_type> row_matrix_type;

    static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::Hiptmair: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

    //! Type of the Ifpack2::Preconditioner specialization from which this class inherits.
    typedef Ifpack2::Preconditioner<scalar_type,
                                    local_ordinal_type,
                                    global_ordinal_type,
                                    node_type> prec_type;

    //@}
    // \name Constructors and Destructors
    //@{

    //! Constructor that takes 1 Tpetra matrix (assumes we'll get the rest off the parameter list)
    explicit Hiptmair (const Teuchos::RCP<const row_matrix_type>& A);


    //! Constructor that takes 3 Tpetra matrices.
    explicit Hiptmair (const Teuchos::RCP<const row_matrix_type>& A,
                       const Teuchos::RCP<const row_matrix_type>& PtAP,
                       const Teuchos::RCP<const row_matrix_type>& P,
                       const Teuchos::RCP<const row_matrix_type>& Pt=Teuchos::null);

    //! Destructor
    virtual ~Hiptmair ();

    //@}
    //! \name Methods for setting parameters and initialization
    //@{

    /// \brief Set the preconditioner's parameters.
    ///
    /// This preconditioner accepts the following parameters:
    ///   - "hiptmair: smoother type 1" (\c std::string)
    ///   - "hiptmair: smoother type 2" (\c std::string)
    ///   - "hiptmair: smoother list 1" (\c Teuchos::ParameterList)
    ///   - "hiptmair: smoother list 2" (\c Teuchos::ParameterList)
    ///   - "hiptmair: pre or post" (\c std::string)
    ///   - "hiptmair: zero starting solution" (\c bool)
    void setParameters (const Teuchos::ParameterList& params);

    bool supportsZeroStartingSolution() { return true; }

    void setZeroStartingSolution (bool zeroStartingSolution) { ZeroStartingSolution_ = zeroStartingSolution; };

    //! Do any initialization that depends on the input matrix's structure.
    void initialize ();

    //! Return \c true if initialize() completed successfully, else \c false.
    inline bool isInitialized () const {
      return IsInitialized_;
    }

    //! Do any initialization that depends on the input matrix's values.
    void compute ();

    //! Return \c true if compute() completed successfully, else \c false.
    inline bool isComputed() const {
      return IsComputed_;
    }

    //@}
    //! @name Implementation of Tpetra::Operator
    //@{

    //! Apply the preconditioner to X, putting the result in Y.
    void
    apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
           scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

    void
    applyHiptmairSmoother(const Tpetra::MultiVector<typename MatrixType::scalar_type,
                          typename MatrixType::local_ordinal_type,
                          typename MatrixType::global_ordinal_type,
                          typename MatrixType::node_type>& X,
                          Tpetra::MultiVector<typename MatrixType::scalar_type,
                          typename MatrixType::local_ordinal_type,
                          typename MatrixType::global_ordinal_type,
                          typename MatrixType::node_type>& Y) const;

    //! A service routine for updating the cached MultiVectors
    void updateCachedMultiVectors(const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& map1,
                                  const Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>>& map2,
                                  size_t numVecs) const;

    //! Tpetra::Map representing the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getDomainMap() const;

    //! Tpetra::Map representing the range of this operator.
    Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getRangeMap() const;

    //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
    bool hasTransposeApply() const;

    //@}
    //! \name Mathematical functions.
    //@{

    //! Returns the operator's communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

    //! Returns a reference to the matrix to be preconditioned.
    Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const;

    //! Returns the number of calls to Initialize().
    int getNumInitialize() const;

    //! Returns the number of calls to Compute().
    int getNumCompute() const;

    //! Returns the number of calls to apply().
    int getNumApply() const;

    //! Returns the time spent in Initialize().
    double getInitializeTime() const;

    //! Returns the time spent in Compute().
    double getComputeTime() const;

    //! Returns the time spent in apply().
    double getApplyTime() const;

    //! Returns prec 1
    Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> > getPrec1();

    //! Returns prec 2
    Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> > getPrec2();

    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

  private:
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    //! Copy constructor (should never be used)
    Hiptmair (const Hiptmair<MatrixType>& RHS);

    //! operator= (should never be used)
    Hiptmair<MatrixType>& operator= (const Hiptmair<MatrixType>& RHS);

    //! The 3 matrices necessary:
    //  A - matrix in primary space
    //  PtAP - matrix in auxiliary space
    //  P  - prolongator matrix
    Teuchos::RCP<const row_matrix_type> A_, PtAP_, P_, Pt_;

    //! Preconditioner types
    std::string precType1_, precType2_, preOrPost_;

    //! If true, the starting solution is always the zero vector.
    bool ZeroStartingSolution_;

    //! If false, explicitely apply Pt
    bool ImplicitTranspose_;

    //! Preconditioner parameters.
    Teuchos::ParameterList precList1_, precList2_;

    //! \c true if \c this object has been initialized
    bool IsInitialized_;
    //! \c true if \c this object has been computed
    bool IsComputed_;
    //! Contains the number of successful calls to Initialize().
    int NumInitialize_;
    //! Contains the number of successful call to Compute().
    int NumCompute_;
    //! Contains the number of successful call to apply().
    mutable int NumApply_;
    //! Contains the time for all successful calls to Initialize().
    double InitializeTime_;
    //! Contains the time for all successful calls to Compute().
    double ComputeTime_;
    //! Contains the time for all successful calls to apply().
    mutable double ApplyTime_;
    //! preconditioners for the two spaces
    Teuchos::RCP<prec_type> ifpack2_prec1_, ifpack2_prec2_;
    //! MultiVectors for caching purposes (so apply doesn't need to allocate them on each call)
    mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > cachedResidual1_;
    mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > cachedSolution1_;
    mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > cachedResidual2_;
    mutable Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > cachedSolution2_;
  };

} // namespace Ifpack2

#endif // IFPACK2_HIPTMAIR_DECL_HPP
