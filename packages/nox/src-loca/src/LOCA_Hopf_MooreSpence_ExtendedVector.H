// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MOORESPENCE_EXTENDEDVECTOR_H
#define LOCA_HOPF_MOORESPENCE_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Hopf {
    namespace MooreSpence {
      class ExtendedMultiVector;
    }
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MooreSpence {

      /*!
       * \brief %Vector class to hold solution vectors, %Newton vectors, etc.
       * for Moore-Spence Hopf equations.
       */
      /*!
       * This class uses the LOCA::Extended::Vector implementation to store the
       * solution, real and imaginary eigenvector, frequency and parameter
       * components of the Hopf vector and merely provides an interface for
       * naming which components of the multivector these quantities
       * correspond to.
       */
      class ExtendedVector : public LOCA::Extended::Vector {

    /*!
     * \brief Declare LOCA::Hopf::MooreSpence::ExtendedMultiVector
     * as a friend class so it can call protected methods.
     */
    friend class ExtendedMultiVector;

      public:

    //! Constructor.
    ExtendedVector(
           const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const NOX::Abstract::Vector& xVec,
           const NOX::Abstract::Vector& realEigenVec,
           const NOX::Abstract::Vector& imagEigenVec,
           double frequency, double bifParam);

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
    virtual ExtendedVector&
    operator=(const ExtendedVector& y);

    //! Cloning function
    virtual Teuchos::RCP<NOX::Abstract::Vector>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

    //! Sets the Hopf vector by setting its five components
    virtual void setVec(const NOX::Abstract::Vector& xVec,
                const NOX::Abstract::Vector& realEigenVec,
                const NOX::Abstract::Vector& imagEigenVec,
                double frequency,
                double bifPar);

    //! Returns the solution vector component of extended vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getXVec() const;

    /*!
     * \brief Returns the real part of the eigenvector component of
     * extended vector
     */
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getRealEigenVec() const;

    /*!
     * \brief Returns the imaginary part of the eigenvector component of
     * extended vector
     */
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getImagEigenVec() const;

    //! Returns the frequency component of the extended vector
    virtual double getFrequency() const;

    //! Get Bifurcation parameter
    virtual double getBifParam() const;

    //! Returns the solution vector component of extended vector
    virtual Teuchos::RCP<NOX::Abstract::Vector> getXVec();

    /*!
     * \brief Returns the real part of the eigenvector component of
     * extended vector
     */
    virtual Teuchos::RCP<NOX::Abstract::Vector> getRealEigenVec();

    /*!
     * \brief Returns the imaginary part of the eigenvector component of
     * extended vector
     */
    virtual Teuchos::RCP<NOX::Abstract::Vector> getImagEigenVec();

    //! Returns the frequency component of the extended vector
    virtual double& getFrequency();

    //! Get Bifurcation parameter
    virtual double& getBifParam();

      protected:

    //! Default Constructor
    ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data);

    //! generate an extended multi-vector
    /*!
     * Returns an empty multivector of type
     * LOCA::Hopf::MooreSpence::ExtendedMultiVector
     */
    virtual Teuchos::RCP<LOCA::Extended::MultiVector>
    generateMultiVector(int nColumns, int nVectorRows,
                int nScalarRows) const;

      }; // class ExtendedVector
    } // namespace MooreSpence
  } // namespace Hopf
} // namespace LOCA

#endif
