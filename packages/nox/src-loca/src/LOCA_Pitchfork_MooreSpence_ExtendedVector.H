// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PITCHFORK_MOORESPENCE_EXTENDEDVECTOR_H
#define LOCA_PITCHFORK_MOORESPENCE_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Pitchfork {
    namespace MooreSpence {
      class ExtendedMultiVector;
    }
  }
}

namespace LOCA {

  namespace Pitchfork {

    namespace MooreSpence {

      /*!
       * \brief %Vector class to hold solution vectors, %Newton vectors, etc.
       * for the Moore-Spence turning point formulation.
       */
      /*!
       * This class uses the LOCA::Extended::Vector implementation to store the
       * solution, null, and parameter components of the turning point vector
       * and merely provides an interface for naming which components of the
       * vector these quantities correspond to.
      */
      class ExtendedVector : public LOCA::Extended::Vector {

    /*!
     * \brief Declare LOCA::Pitchfork::MooreSpence::ExtendedMultiVector
     * as a friend class so it can call protected methods.
     */
    friend class ExtendedMultiVector;

      public:

    //! %Constructor.
    ExtendedVector(
           const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const NOX::Abstract::Vector& xVec,
           const NOX::Abstract::Vector& nullVec,
           double slack,
           double bifParam);

    //! Copy constructor.
    ExtendedVector(const ExtendedVector& source,
               NOX::CopyType type = NOX::DeepCopy);

    //! Destructor.
    virtual ~ExtendedVector();

    //! Assignment operator
    virtual NOX::Abstract::Vector&
    operator=(const NOX::Abstract::Vector& y);

    //! Assignment operator
    virtual LOCA::Extended::Vector&
    operator=(const LOCA::Extended::Vector& y);

    //! Assignment operator
    virtual ExtendedVector& operator=(const ExtendedVector& y);

    //! Cloning function
    virtual Teuchos::RCP<NOX::Abstract::Vector>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

    //! Sets the %Vector by setting its three components
    virtual void setVec(const NOX::Abstract::Vector& xVec,
                const NOX::Abstract::Vector& nullVec,
                double slack,
                double bifPar);

    //! Returns the solution vector component of extended vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getXVec() const;

    //! Returns the null vector component of extended vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getNullVec() const;

    //! Get slack component
    virtual double getSlack() const;

    //! Get %Bifurcation parameter
    virtual double getBifParam() const;

    //! Returns the solution vector component of extended vector
    virtual Teuchos::RCP<NOX::Abstract::Vector> getXVec();

    //! Returns the null vector component of extended vector
    virtual Teuchos::RCP<NOX::Abstract::Vector> getNullVec();

    //! Get slack component
    virtual double& getSlack();

    //! Get %Bifurcation parameter
    virtual double& getBifParam();

      protected:

    //! Default Constructor
    ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data);

    //! generate an extended multi-vector
    /*!
     * Returns an empty multivector of type
     * LOCA::Pitchfork::MooreSpence::ExtendedMultiVector
     */
    virtual Teuchos::RCP<LOCA::Extended::MultiVector>
    generateMultiVector(int nColumns, int nVectorRows,
                int nScalarRows) const;

      }; // class ExtendedVector
    } // namepsace MooreSpence
  } // namespace Pitchfork
} // namespace LOCA

#endif
