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

#ifndef LOCA_HOPF_MOORESPENCE_EXTENDEDMULTIVECTOR_H
#define LOCA_HOPF_MOORESPENCE_EXTENDEDMULTIVECTOR_H

#include "LOCA_Extended_MultiVector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Hopf {
    namespace MooreSpence {
      class ExtendedVector;
    }
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MooreSpence {

      /*! 
       * \brief Multi-vector class to hold solution vectors, %Newton vectors, 
       * etc.for the Moore-Spence Hopf eqautions.
       */
      /*! 
       * This class uses the LOCA::Extended::MultiVector implementation to 
       * store the solution, real and imaginary eigenvector, frequency and 
       * parameter components of the Hopf multi vector and merely provides an 
       * interface for naming which components of the multivector these 
       * quantities correspond to.
      */
      class ExtendedMultiVector : public LOCA::Extended::MultiVector {

	/*!
	 * \brief Declare LOCA::Hopf::MooreSpence::ExtendedVector as 
	 * a friend class so it can call protected methods.
	 */
	friend class ExtendedVector;

      public:

	//! Constructor.  
	/*!
	 * Generates a multivector with \em nColumns columns from \em cloneVec
	 */
	ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns);
  
	//! Constructor.  
	ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::MultiVector& xVec,
		    const NOX::Abstract::MultiVector& realEigenVec,
		    const NOX::Abstract::MultiVector& imagEigenVec,
		    const NOX::Abstract::MultiVector::DenseMatrix& freqs,
		    const NOX::Abstract::MultiVector::DenseMatrix& bifParams);

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

	//! Returns the solution vector component of extended multivector
	virtual Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> 
	getXMultiVec() const;

	//! Returns the solution vector component of extended multivector
	virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector> 
	getXMultiVec();

	/*! 
	 * \brief Returns the real part of the eigenvector component of 
	 * extended multivector
	 */
	virtual Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> 
	getRealEigenMultiVec() const;

	/*! 
	 * \brief Returns the real part of the eigenvector component of 
	 * extended multivector
	 */
	virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector> 
	getRealEigenMultiVec();

	/*! 
	 * \brief Returns the imaginary part of the eigenvector component of 
	 * extended multivector
	 */
	virtual Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> 
	getImagEigenMultiVec() const;

	/*! 
	 * \brief Returns the imaginary part of the eigenvector component of 
	 * extended multivector
	 */
	virtual Teuchos::RefCountPtr<NOX::Abstract::MultiVector> 
	getImagEigenMultiVec();

	//! Returns frequency component of extended multi vector
	virtual 
	Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>
	getFrequencies() const;

	//! Returns frequency component of extended multi vector
	virtual 
	Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix>
	getFrequencies();

	//! Returns bifurcation parameter component of extended multi vector
	virtual 
	Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>
	getBifParams() const;

	//! Returns bifurcation parameter component of extended multi vector
	virtual 
	Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix>
	getBifParams();

	//! Returns ith column as an extended vector
	virtual 
	Teuchos::RefCountPtr<LOCA::Hopf::MooreSpence::ExtendedVector>
	getColumn(int i);

	//! Returns ith column as an extended vector
	virtual 
	Teuchos::RefCountPtr<const LOCA::Hopf::MooreSpence::ExtendedVector>
	getColumn(int i) const;

      protected:

	//! Constructor.  
	/*!
	 * Creates an empty multivector with \em nColumns columns
	 */
	ExtendedMultiVector(
		     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		     int nColumns);

	//! Generate a derived extended vector
	/*!
	 * Returns a vector of type 
	 * LOCA::Hopf::MooreSpence::ExtendedVector
	 */
	virtual Teuchos::RefCountPtr<LOCA::Extended::Vector>
	generateVector(int nVecs, int nScalarRows) const;

      }; // class ExtendedVector
    } // namespace MooreSpence
  } // namespace Hopf
} // namespace LOCA

#endif
