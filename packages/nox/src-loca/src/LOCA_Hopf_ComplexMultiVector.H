// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns);
  
      //! Constructor.  
      ComplexMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
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
			 const vector<int>& index, bool view);

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
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
      clone(NOX::CopyType type = NOX::DeepCopy) const;
      
      //! Creates a new multi-vector with \c numvecs columns
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
      clone(int numvecs) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns 
       * whose columns are copies of the columns of \c *this given by 
       * \c index.
       */
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
      subCopy(const vector<int>& index) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns that 
       * shares the columns of \c *this given by \c index.
       */
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
      subView(const vector<int>& index) const;

      //! Returns the real component of extended multivector
      virtual Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> 
      getRealMultiVec() const;

      //! Returns the real component of extended multivector
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector> 
      getRealMultiVec();

      //! Returns the imaginary component of extended multivector
      virtual Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> 
      getImagMultiVec() const;

      //! Returns the imaginary component of extended multivector
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector> 
      getImagMultiVec();

      //! Returns ith column as an extended vector
      virtual 
      Teuchos::RefCountPtr<LOCA::Hopf::ComplexVector>
      getColumn(int i);

      //! Returns ith column as an extended vector
      virtual 
      Teuchos::RefCountPtr<const LOCA::Hopf::ComplexVector>
      getColumn(int i) const;

    protected:

      //! Constructor.  
      /*!
       * Creates an empty multivector with \em nColumns columns
       */
      ComplexMultiVector(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		     int nColumns);

      //! Generate a derived extended vector
      /*!
       * Returns a vector of type 
       * LOCA::Hopf::ComplexVector
       */
      virtual Teuchos::RefCountPtr<LOCA::Extended::Vector>
      generateVector(int nVecs, int nScalarRows) const;

    }; // class ComplexMultiVector
  } // namespace Hopf
} // namespace LOCA

#endif
