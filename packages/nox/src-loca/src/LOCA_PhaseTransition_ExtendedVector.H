// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PHASE_TRANSITION_EXTENDEDVECTOR_H
#define LOCA_PHASE_TRANSITION_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace PhaseTransition {
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  namespace PhaseTransition {

      /*!
       * \brief %Vector class to hold solution vectors, %Newton vectors, etc.
       * for the Phase Transition tracking formulation.
       */
      /*!
       * This class uses the LOCA::Extended::Vector implementation to store the
       * solution1, solution2, and parameter components of the phase transition vector
       * and merely provides an interface for naming which components of the
       * vector these quantities correspond to.
      */
      class ExtendedVector : public LOCA::Extended::Vector {

    /*!
     * \brief Declare LOCA::PhaseTransition::MooreSpence::ExtendedMultiVector
     * as a friend class so it can call protected methods.
     */
    friend class ExtendedMultiVector;

      public:

    //! %Constructor.
    ExtendedVector(
           const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const NOX::Abstract::Vector& x1Vec,
           const NOX::Abstract::Vector& x2Vec,
           double ptp);

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
                double bifPar);

    //! Returns the solution1 vector component of extended vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    X1() const;

    //! Returns the solution2 vector component of extended vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    X2() const;

    //! Get %Bifurcation parameter
    virtual double PTP() const;

    //! Returns the solution vector component of extended vector
    virtual Teuchos::RCP<NOX::Abstract::Vector> X1();

    //! Returns the null vector component of extended vector
    virtual Teuchos::RCP<NOX::Abstract::Vector> X2();

    //! Get %Bifurcation parameter
    virtual double& PTP();

      protected:

    //! Default Constructor
    ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data);

    //! generate an extended multi-vector
    /*!
     * Returns an empty multivector of type
     * LOCA::PhaseTransition::ExtendedMultiVector
     */
    virtual Teuchos::RCP<LOCA::Extended::MultiVector>
    generateMultiVector(int nColumns, int nVectorRows,
                int nScalarRows) const;

      }; // class ExtendedVector
  } // namespace PhaseTransition
} // namespace LOCA

#endif
