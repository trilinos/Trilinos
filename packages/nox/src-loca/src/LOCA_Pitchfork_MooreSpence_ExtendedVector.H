// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
