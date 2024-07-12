// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_EXTENDEDVECTOR_H
#define LOCA_MULTICONTINUATION_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    /*!
      \brief %Vector class to hold solution vectors, %Newton vectors, etc.
      for continuation equations.
    */
    /*!
       This class uses the LOCA::Extended::Vector implementation to store the
       solution and parameter components of the continuation vector
       and merely provides an interface for naming which components of the
       multivector these quantities correspond to.
    */
    class ExtendedVector : public LOCA::Extended::Vector {

      /*!
       * \brief Declare LOCA::MultiContinuation::ExtendedMultiVector as a
       * friend class so it can call protected methods.
       */
      friend class ExtendedMultiVector;

    public:

      //! Constructor.
      ExtendedVector(const Teuchos::RCP<LOCA::GlobalData>& global_data,
             const NOX::Abstract::Vector& xVec, int nScalars);

      //! Copy constructor.
      ExtendedVector(const ExtendedVector& source,
             NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ExtendedVector();

      //! Assignment operator
      virtual ExtendedVector&
      operator=(const ExtendedVector& y);

      //! Assignment operator
      virtual LOCA::Extended::Vector&
      operator=(const LOCA::Extended::Vector& y);

      //! Assignment operator
      virtual NOX::Abstract::Vector&
      operator=(const NOX::Abstract::Vector& y);

      //! Assignment operator
      virtual Teuchos::RCP<NOX::Abstract::Vector>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Returns the solution vector component of extended vector
      virtual Teuchos::RCP<const NOX::Abstract::Vector>
      getXVec() const;

      //! Returns the solution vector component of extended vector
      virtual Teuchos::RCP<NOX::Abstract::Vector> getXVec();

    protected:

      //! Constructor
      ExtendedVector(const Teuchos::RCP<LOCA::GlobalData>& global_data,
             int nScalars);

      //! generate an extended multi-vector
      /*!
       * Returns an empty multivector of type
       * LOCA::MultiContinuation::ExtendedMultiVector
       */
      virtual Teuchos::RCP<LOCA::Extended::MultiVector>
      generateMultiVector(int nColumns, int nVectorRows,
              int nScalarRows) const;

    }; // class ExtendedVector
  } // namespace Continuation
} // namespace LOCA

#endif
