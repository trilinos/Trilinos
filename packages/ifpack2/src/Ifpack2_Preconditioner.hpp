// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// @file Ifpack2_Preconditioner.hpp

#ifndef IFPACK2_PRECONDITIONER_HPP
#define IFPACK2_PRECONDITIONER_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <iostream>

namespace Ifpack2 {

/// \class Preconditioner
/// \brief Interface for all Ifpack2 preconditioners
/// \tparam Scalar Type of the matrix's entries; same as the first
///   template parameter of Tpetra::RowMatrix
/// \tparam LocalOrdinal Type of the matrix's local indices; same as the
///   second template parameter of Tpetra::RowMatrix
/// \tparam GlobalOrdinal Type of the matrix's global indices; same as the
///   third template parameter of Tpetra::RowMatrix
/// \tparam Node The matrix's Node type; same as the fourth template
///   parameter of Tpetra::RowMatrix
///
/// The Preconditioner class defines the interface that all Ifpack2
/// preconditioners must implement.  Preconditioner inherits from
/// Tpetra::Operator.  Its apply() method applies the preconditioner.
///
/// If you are familiar with the IFPACK package, please be aware that
/// this is different from IFPACK.  In IFPACK, the ApplyInverse()
/// method applies or "solves with" the preconditioner \f$M^{-1}\f$,
/// and the Apply() method "applies" the preconditioner \f$M\f$.  In
/// Ifpack2, the apply() method applies or "solves with" the
/// preconditioner \f$M^{-1}\f$.  Ifpack2 has no method comparable to
/// IFPACK's Apply().
///
/// Preconditioner provides the following methods
///   - initialize() performs all operations based on the graph of the
///     matrix (without considering the numerical values)
///   - isInitialized() returns true if the preconditioner has been
///     successfully initialized
///   - compute() computes everything required to apply the
///     preconditioner, using the matrix's values (and assuming that the
///     graph structure of the matrix has not changed)
///   - isComputed() returns true if the preconditioner has been
///     successfully computed, false otherwise.
///   - getMatrix() returns a reference to the matrix to be preconditioned
///
/// Implementations of compute() must internally call initialize() if
/// isInitialized() returns false. The preconditioner is applied by
/// apply() (which returns if isComputed() is false). Every time that
/// initialize() is called, the object destroys all the previously
/// allocated information, and reinitializes the preconditioner. Every
/// time compute() is called, the object recomputes the actual values of
/// the preconditioner.
template<class Scalar =
           Tpetra::Operator<>::scalar_type,
         class LocalOrdinal =
           typename Tpetra::Operator<Scalar>::local_ordinal_type,
         class GlobalOrdinal =
           typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
         class Node =
           typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class Preconditioner :
  virtual public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:
  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  //! Destructor.
  virtual ~Preconditioner(){}

  /// @name Methods implementing Tpetra::Operator.
  //@{

  /// @brief The domain Map of this operator.
  ///
  /// The domain Map describes the distribution of valid input vectors
  /// X to the apply() method.
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getDomainMap () const = 0;

  /// @brief The range Map of this operator.
  ///
  /// The range Map describes the distribution of valid output vectors
  /// Y to the apply() method.
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getRangeMap () const = 0;

  /// @brief Apply the preconditioner to X, putting the result in Y.
  ///
  /// If the result of applying this preconditioner to a vector X is
  /// \f$F \cdot X$\f$, then this method computes \f$\beta Y + \alpha F \cdot X\f$.
  /// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  virtual void
  apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
         Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const = 0;
  //@}

  //! Set this preconditioner's parameters.
  virtual void setParameters (const Teuchos::ParameterList& List) = 0;

  virtual bool supportsZeroStartingSolution () { return false; };

  //! Set this preconditioner's parameters.
  virtual void setZeroStartingSolution (bool zeroStartingSolution) { TEUCHOS_ASSERT(false); };

  /// @brief Set up the graph structure of this preconditioner.
  ///
  /// If the graph structure of the constructor's input matrix has
  /// changed, or if you have not yet called initialize(), you must
  /// call initialize() before you may call compute() or apply().
  ///
  /// Thus, initialize() corresponds to the "symbolic factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  virtual void initialize() = 0;

  //! True if the preconditioner has been successfully initialized, else false.
  virtual bool isInitialized() const = 0;

  /// @brief Set up the numerical values in this preconditioner.
  ///
  /// If the values of the constructor's input matrix have changed, or
  /// if you have not yet called compute(), you must call compute()
  /// before you may call apply().
  ///
  /// Thus, compute() corresponds to the "numeric factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  virtual void compute() = 0;

  //! True if the preconditioner has been successfully computed, else false.
  virtual bool isComputed() const = 0;

  //! The input matrix given to the constructor.
  virtual Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const = 0;

  //! The number of calls to initialize().
  virtual int getNumInitialize() const = 0;

  //! The number of calls to compute().
  virtual int getNumCompute() const = 0;

  //! The number of calls to apply().
  virtual int getNumApply() const = 0;

  //! The time (in seconds) spent in initialize().
  virtual double getInitializeTime() const = 0;

  //! The time (in seconds) spent in compute().
  virtual double getComputeTime() const = 0;

  //! The time (in seconds) spent in apply().
  virtual double getApplyTime() const = 0;
};

}//namespace Ifpack2

#endif // IFPACK2_PRECONDITIONER_HPP
