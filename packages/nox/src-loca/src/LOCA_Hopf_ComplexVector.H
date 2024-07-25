// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_COMPLEXVECTOR_H
#define LOCA_HOPF_COMPLEXVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Hopf {
    class ComplexMultiVector;
  }
}

namespace LOCA {

  namespace Hopf {

    /*!
     * \brief %Vector class to hold two vectors to represent a complex vector.
     */
    /*!
     * This is not a true complex vector.  Operations like innerProduct()
     * are not correct for complex vectors.  This class
     * exists to make some aspects of the real-equivalent formulation of
     * complex linear algebra simpler to implement.
     */
    class ComplexVector : public LOCA::Extended::Vector {

      /*!
       * \brief Declare LOCA::Hopf::ComplexMultiVector
       * as a friend class so it can call protected methods.
       */
      friend class ComplexMultiVector;

    public:

      //! Constructor.
      ComplexVector(const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& realVec,
            const NOX::Abstract::Vector& imagVec);

      //! Copy constructor.
      ComplexVector(const ComplexVector& source,
            NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ComplexVector();

      //! Assignment operator
      virtual NOX::Abstract::Vector&
      operator=(const NOX::Abstract::Vector& y);

      //! Assignment operator
      virtual LOCA::Extended::Vector&
      operator=(const LOCA::Extended::Vector& y);

      //! Assignment operator
      virtual ComplexVector&
      operator=(const ComplexVector& y);

      //! Cloning function
      virtual Teuchos::RCP<NOX::Abstract::Vector>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Sets the vector by setting its components
      virtual void setVec(const NOX::Abstract::Vector& realVec,
              const NOX::Abstract::Vector& imagVec);

      //! Returns the real component of extended vector
      virtual Teuchos::RCP<const NOX::Abstract::Vector>
      getRealVec() const;

      //! Returns the imaginary component of extended vector
      virtual Teuchos::RCP<const NOX::Abstract::Vector>
      getImagVec() const;

      //! Returns the real component of extended vector
      virtual Teuchos::RCP<NOX::Abstract::Vector>
      getRealVec();

      //! Returns the imaginary component of extended vector
      virtual Teuchos::RCP<NOX::Abstract::Vector>
      getImagVec();

    protected:

      //! Default Constructor
      ComplexVector(
           const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! generate an extended multi-vector
      /*!
       * Returns an empty multivector of type
       * LOCA::Hopf::ComplexMultiVector
       */
      virtual Teuchos::RCP<LOCA::Extended::MultiVector>
      generateMultiVector(int nColumns, int nVectorRows,
              int nScalarRows) const;

    }; // class ComplexVector
  } // namespace Hopf
} // namespace LOCA

#endif
