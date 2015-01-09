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
