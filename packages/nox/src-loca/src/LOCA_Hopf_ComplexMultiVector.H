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
