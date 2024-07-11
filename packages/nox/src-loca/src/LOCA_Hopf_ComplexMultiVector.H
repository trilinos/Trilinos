// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_COMPLEXMULTIVECTOR_H
#define LOCA_HOPF_COMPLEXMULTIVECTOR_H

#include "LOCA_Extended_MultiVector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Hopf {
    class ComplexVector;
  }
}

namespace LOCA {

  namespace Hopf {

    /*!
     * \brief Multi-vector class to hold two multi-vectors to represent a
     * complex multi-vector.
     */
    /*!
     * This is not a true complex multi-vector.  Operations like dot()
     * and multiply() are not correct for complex vectors.  This class
     * exists to make some aspects of the real-equivalent formulation of
     * complex linear algebra simpler to implement.
     */
    class ComplexMultiVector : public LOCA::Extended::MultiVector {

      /*!
       * \brief Declare LOCA::Hopf::ComplexVector as
       * a friend class so it can call protected methods.
       */
      friend class ComplexVector;

    public:

      //! Constructor.
      /*!
       * Generates a multivector with \em nColumns columns from \em cloneVec
       */
      ComplexMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& cloneVec,
            int nColumns);

      //! Constructor.
      ComplexMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::MultiVector& realVec,
            const NOX::Abstract::MultiVector& imagVec);

      //! Copy constructor.
      ComplexMultiVector(const ComplexMultiVector& source,
             NOX::CopyType type = NOX::DeepCopy);

      /*!
       * \brief Copy constructor that creates a new multivector with
       * \em nColumns columns.
       */
      ComplexMultiVector(const ComplexMultiVector& source, int nColumns);

      /*!
       * \brief Copy constructor that creates a sub copy or view of the
       * given multivector
       */
      ComplexMultiVector(const ComplexMultiVector& source,
             const std::vector<int>& index, bool view);

      //! Destructor.
      virtual ~ComplexMultiVector();

      //! Assignment operator
      virtual ComplexMultiVector&
      operator=(const ComplexMultiVector& y);

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
       * \brief Creates a new multi-vector with \c index.size() columns
       * whose columns are copies of the columns of \c *this given by
       * \c index.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      subCopy(const std::vector<int>& index) const;

      /*!
       * \brief Creates a new multi-vector with \c index.size() columns that
       * shares the columns of \c *this given by \c index.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      subView(const std::vector<int>& index) const;

      //! Returns the real component of extended multivector
      virtual Teuchos::RCP<const NOX::Abstract::MultiVector>
      getRealMultiVec() const;

      //! Returns the real component of extended multivector
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      getRealMultiVec();

      //! Returns the imaginary component of extended multivector
      virtual Teuchos::RCP<const NOX::Abstract::MultiVector>
      getImagMultiVec() const;

      //! Returns the imaginary component of extended multivector
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      getImagMultiVec();

      //! Returns ith column as an extended vector
      virtual
      Teuchos::RCP<LOCA::Hopf::ComplexVector>
      getColumn(int i);

      //! Returns ith column as an extended vector
      virtual
      Teuchos::RCP<const LOCA::Hopf::ComplexVector>
      getColumn(int i) const;

    protected:

      //! Constructor.
      /*!
       * Creates an empty multivector with \em nColumns columns
       */
      ComplexMultiVector(
             const Teuchos::RCP<LOCA::GlobalData>& global_data,
             int nColumns);

      //! Generate a derived extended vector
      /*!
       * Returns a vector of type
       * LOCA::Hopf::ComplexVector
       */
      virtual Teuchos::RCP<LOCA::Extended::Vector>
      generateVector(int nVecs, int nScalarRows) const;

    }; // class ComplexMultiVector
  } // namespace Hopf
} // namespace LOCA

#endif
