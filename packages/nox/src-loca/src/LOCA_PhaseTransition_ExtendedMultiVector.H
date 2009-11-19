// $Id: LOCA_PhaseTransition_MooreSpence_ExtendedMultiVector.H,v 1.8 2007/06/21 16:22:52 rhoope Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_MooreSpence_ExtendedMultiVector.H,v $ 

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
//  $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_MooreSpence_ExtendedMultiVector.H,v $
//  $Author: rhoope $
//  $Date: 2007/06/21 16:22:52 $
//  $Revision: 1.8 $
// ************************************************************************
//@HEADER

#ifndef LOCA_PHASETRANSITION_EXTENDEDMULTIVECTOR_H
#define LOCA_PHASETRANSITION_EXTENDEDMULTIVECTOR_H

#include "LOCA_Extended_MultiVector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace PhaseTransition {
      class ExtendedVector;
  }
}

namespace LOCA {

  namespace PhaseTransition {

      /*! 
       * \brief %MultiVector class to hold solution vectors, %Newton vectors, 
       * etc. for the phase transition tracking algorithm.
       */
      /*! 
       * This class uses the LOCA::Extended::MultiVector implementation to 
       * store the solution1, solution2, and parameter components of the phase 
       * transition multivector and merely provides an interface for naming which 
       * components of the multivector these quantities correspond to.
      */
      class ExtendedMultiVector : public LOCA::Extended::MultiVector {

	/*!
	 * \brief Declare LOCA::PhaseTransition::MooreSpence::ExtendedVector as 
	 * a friend class so it can call protected methods.
	 */
	friend class ExtendedVector;

      public:

	//! Constructor.  
	/*!
	 * Generates a multivector with \em nColumns columns from \em cloneVec
	 */
	ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns);
  
	//! Constructor.  
	/*!
	 * Construct the multivector from \em xVec, \em nullVec, and 
	 * \em bifParams
	 */
	ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::MultiVector& xVec,
		    const NOX::Abstract::MultiVector& nullVec,
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
	virtual Teuchos::RCP<NOX::Abstract::MultiVector>
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Creates a new multi-vector with \c numvecs columns
	virtual Teuchos::RCP<NOX::Abstract::MultiVector>
	clone(int numvecs) const;

	/*! 
	 * \brief Creates a new multi-vector with \c index.size() columns whose 
	 * columns are copies of the columns of \c *this given by \c index.
	 */
	virtual Teuchos::RCP<NOX::Abstract::MultiVector>
	subCopy(const vector<int>& index) const;

	/*! 
	 * \brief Creates a new multi-vector with \c index.size() columns that 
	 * shares the columns of \c *this given by \c index.
	 */
	virtual Teuchos::RCP<NOX::Abstract::MultiVector>
	subView(const vector<int>& index) const;

	//! Returns the solution vector component of extended multivector
	virtual Teuchos::RCP<const NOX::Abstract::MultiVector> 
	getXMultiVec() const;

	//! Returns the solution vector component of extended multivector
	virtual Teuchos::RCP<NOX::Abstract::MultiVector> 
	getXMultiVec();

	//! Returns the null vector component of extended multivector
	virtual Teuchos::RCP<const NOX::Abstract::MultiVector> 
	getNullMultiVec() const;

	//! Returns the null vector component of extended multivector
	virtual Teuchos::RCP<NOX::Abstract::MultiVector> 
	getNullMultiVec();

	//! Returns ith column as an extended vector
	virtual 
	Teuchos::RCP<LOCA::PhaseTransition::ExtendedVector>
	getColumn(int i);

	//! Returns ith column as an extended vector
	virtual 
	Teuchos::RCP<const LOCA::PhaseTransition::ExtendedVector>
	getColumn(int i) const;

      protected:

	//! Constructor.  
	/*!
	 * Creates an empty multivector with \em nColumns columns
	 */
	ExtendedMultiVector(
		     const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     int nColumns);

	//! Generate a derived extended vector
	/*!
	 * Returns a vector of type 
	 * LOCA::PhaseTransition::MooreSpence::ExtendedVector
	 */
	virtual Teuchos::RCP<LOCA::Extended::Vector>
	generateVector(int nVecs, int nScalarRows) const;

      }; // class ExtendedVector
  } // namespace PhaseTransition
} // namespace LOCA

#endif
