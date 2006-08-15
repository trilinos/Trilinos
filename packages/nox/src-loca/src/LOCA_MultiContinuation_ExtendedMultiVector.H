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

#ifndef LOCA_MULTICONTINUATION_EXTENDEDMULTIVECTOR_H
#define LOCA_MULTICONTINUATION_EXTENDEDMULTIVECTOR_H

#include "LOCA_Extended_MultiVector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedVector;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    /*! 
      \brief %MultiVector class to hold solution vectors, %Newton vectors, etc.
      for continuation equations.
    */
    /*! 
     * This class uses the LOCA::Extended::MultiVector implementation to 
     * store the solution and parameter components of the continuation vector
     * and merely provides an interface for naming which components of the
     * multivector these quantities correspond to.
    */
    class ExtendedMultiVector : public LOCA::Extended::MultiVector {

      //! Declare LOCA::Continuation::ExtendedVector as a friend class
      /*!
       * So it can call protected methods.
       */
      friend class ExtendedVector;

    public:

      //! Constructor.  
      /*!
       * Generates a multivector with nColumns from \em xVec amd
       * \em nScalarRows of zeros.
       */
      ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& xVec, int nColumns,
		    int nScalarRows, NOX::CopyType type = NOX::DeepCopy);
  
      //! Constructor.  
      /*!
       * Initializes the scalar matrix to \em nScalarRows rows and 
       * \em xVec.numVectors() columns of zeros
       */
      ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::MultiVector& xVec, 
		    int nScalarRows);

      //! Constructor.  
      /*!
       * Sets the scalar matrix explicitly
       */
      ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::MultiVector& xVec,
		    const NOX::Abstract::MultiVector::DenseMatrix& params);

      //! Copy constructor.
      ExtendedMultiVector(const ExtendedMultiVector& source, 
		     NOX::CopyType type = NOX::DeepCopy);

      /*!
       * \brief Copy constructor that creates a new multivector with
       * \em nColumns columns.
       */
      ExtendedMultiVector(const ExtendedMultiVector& source, int nColumns);

      /*!
       * \brief Copy constructor that creates a sub copy or view of the 
       * given multivector
       */
      ExtendedMultiVector(const ExtendedMultiVector& source, 
			  const vector<int>& index, bool view);

      //! Destructor.
      virtual ~ExtendedMultiVector();

      //! Assignment operator
      virtual ExtendedMultiVector& 
      operator=(const ExtendedMultiVector& y);

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
       * \brief Creates a new multi-vector with \c index.size() columns whose 
       * columns are copies of the columns of \c *this given by \c index.
       */
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
      subCopy(const vector<int>& index) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns that 
       * shares the columns of \c *this given by \c index.
       */
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
      subView(const vector<int>& index) const;

      //! Returns the solution vector component of extended multivector
      virtual Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> 
      getXMultiVec() const;

      //! Returns the solution vector component of extended multivector
      virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector> getXMultiVec();

    protected:

      //! Constructor.  
      /*!
       * Creates an empty multivector with \em nColumns columns
       * and \em nScalarRows scalar rows 
       */
      ExtendedMultiVector(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		     int nColumns, int nScalarRows);

      //! Generate a derived extended vector
      /*!
       * Returns a vector of type LOCA::MultiContinuation::ExtendedVector
       */
      virtual Teuchos::RefCountPtr<LOCA::Extended::Vector>
      generateVector(int nVecs, int nScalarRows) const;

    }; // class ExtendedVector
  } // namespace Continuation
} // namespace LOCA

#endif
