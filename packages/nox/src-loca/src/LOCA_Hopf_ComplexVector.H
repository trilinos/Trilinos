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
