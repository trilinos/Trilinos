// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PITCHFORK_MOORESPENCE_EXTENDEDMULTIVECTOR_H
#define LOCA_PITCHFORK_MOORESPENCE_EXTENDEDMULTIVECTOR_H

#include "LOCA_Extended_MultiVector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Pitchfork {
    namespace MooreSpence {
      class ExtendedVector;
    }
  }
}

namespace LOCA {

  namespace Pitchfork {

    namespace MooreSpence {

      /*!
       * \brief %MultiVector class to hold solution vectors, %Newton vectors,
       * etc.for the Moore-Spence pitchfork formulation.
       */
      /*!
       * This class uses the LOCA::Extended::MultiVector implementation to
       * store the solution, null, parameter, and slack components of the
       * pitchfork multivector and merely provides an interface for naming
       * which components of the multivector these quantities correspond to.
      */
      class ExtendedMultiVector : public LOCA::Extended::MultiVector {

    /*!
     * \brief Declare LOCA::Pitchfork::MooreSpece::ExtendedVector as
     * a friend class so it can call protected methods.
     */
    friend class ExtendedVector;

      public:

    //! Constructor.
    /*!
     * Generates a multivector with \em nColumns columns from \em cloneVec
     */
    ExtendedMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& cloneVec,
            int nColumns);

    //! Constructor.
    /*!
     * Construct the multivector from \em xVec, \em nullVec, and
     * \em bifParams
     */
    ExtendedMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::MultiVector& xVec,
            const NOX::Abstract::MultiVector& nullVec,
            const NOX::Abstract::MultiVector::DenseMatrix& slacks,
            const NOX::Abstract::MultiVector::DenseMatrix& bifParams);

    //! Copy constructor.
    ExtendedMultiVector(const ExtendedMultiVector& source,
                NOX::CopyType type = NOX::DeepCopy);

    /*!
     * \brief Copy constructor that creates a new multivector with
     * \em nColumns columns.
     */
    ExtendedMultiVector(const ExtendedMultiVector& source, int nColumns);

    /*!
     * \brief Copy constructor that creates a sub copy or view of the
     * given multivector
     */
    ExtendedMultiVector(const ExtendedMultiVector& source,
                const std::vector<int>& index, bool view);

    //! Destructor.
    virtual ~ExtendedMultiVector();

    //! Assignment operator
    virtual ExtendedMultiVector&
    operator=(const ExtendedMultiVector& y);

    //! Assignment operator
    virtual LOCA::Extended::MultiVector&
    operator=(const LOCA::Extended::MultiVector& y);

    //! Assignment operator
    virtual NOX::Abstract::MultiVector&
    operator=(const NOX::Abstract::MultiVector& y);

    /*!
     * \brief Create a new multi-vector of the same underlying type by
     * cloning "this", and return a pointer to the new vector.
     */
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

    //! Creates a new multi-vector with \c numvecs columns
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    clone(int numvecs) const;

    /*!
     * \brief Creates a new multi-vector with \c index.size() columns whose
     * columns are copies of the columns of \c *this given by \c index.
     */
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    subCopy(const std::vector<int>& index) const;

    /*!
     * \brief Creates a new multi-vector with \c index.size() columns that
     * shares the columns of \c *this given by \c index.
     */
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    subView(const std::vector<int>& index) const;

    //! Returns the solution vector component of extended multivector
    virtual Teuchos::RCP<const NOX::Abstract::MultiVector>
    getXMultiVec() const;

    //! Returns the solution vector component of extended multivector
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    getXMultiVec();

    //! Returns the null vector component of extended multivector
    virtual Teuchos::RCP<const NOX::Abstract::MultiVector>
    getNullMultiVec() const;

    //! Returns the null vector component of extended multivector
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    getNullMultiVec();

    //! Returns slack component of the extended multivector
    virtual
    Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
    getSlacks() const;

    //! Returns slack component of the extended multivector
    virtual
    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
    getSlacks();

    //! Returns bifurcation parameter component of the extended multivector
    virtual
    Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
    getBifParams() const;

    //! Returns bifurcation parameter component of the extended multivector
    virtual
    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
    getBifParams();

    //! Returns ith column as an extended vector
    virtual
    Teuchos::RCP<LOCA::Pitchfork::MooreSpence::ExtendedVector>
    getColumn(int i);

    //! Returns ith column as an extended vector
    virtual
    Teuchos::RCP<const LOCA::Pitchfork::MooreSpence::ExtendedVector>
    getColumn(int i) const;

      protected:

    //! Constructor.
    /*!
     * Creates an empty multivector with \em nColumns columns
     */
    ExtendedMultiVector(
             const Teuchos::RCP<LOCA::GlobalData>& global_data,
             int nColumns);

    //! Generate a derived extended vector
    /*!
     * Returns a vector of type
     * LOCA::Pitchfork::MooreSpence::ExtendedVector
     */
    virtual Teuchos::RCP<LOCA::Extended::Vector>
    generateVector(int nVecs, int nScalarRows) const;

      }; // class ExtendedVector
    } // namespace MooreSpence
  } // namespace Pitchfork
} // namespace LOCA

#endif
