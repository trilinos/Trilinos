// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSYSTEM_MANAGER_H
#define LOCA_BORDEREDSYSTEM_MANAGER_H

#include "LOCA_BorderedSystem_Generic.H"  // base class

namespace LOCA {

  namespace BorderedSystem {
    
    /*! 
     * \brief %Manager for all bordered system solver methods.  See
     * LOCA::BorderedSystem::Generic for a description of what these
     * classes compute.
     */
    /*! 
     * The parameters passed to the constructor or reset should specify
     * the "Bordered Solver Method", as described below, as well as any 
     * additional parameters for that particular method.
     *
     * <ul>
     * <li> "Bordered Solver Method" - Name of the method.  Valid choices are:
     * <ul> 
     * <li> "%Bordering" (LOCA::BorderedSystem::Bordering) [Default]
     * </ul>
     * </ul>
     */  
    class Manager : public LOCA::BorderedSystem::Generic {

    public:

      //! Constructor. 
      Manager(NOX::Parameter::List& params);

      //! Copy constructor
      Manager(const Manager& source);

      //! Destructor
      virtual ~Manager();

      //! Assignment operator
      virtual Manager& 
      operator = (const Manager& source);

      //! Clone function
      virtual LOCA::BorderedSystem::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::BorderedSystem::Generic& 
      operator = (const LOCA::BorderedSystem::Generic& source);

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);

      //! Set flags indicating if any blocks are zero
      virtual void setIsZero(bool flagA, bool flagB, bool flagC);

      //! Set blocks
      virtual void setBlocks(
			const NOX::Abstract::Group* group,
			const NOX::Abstract::MultiVector* blockA,
			const NOX::Abstract::MultiVector* blockB,
			const NOX::Abstract::MultiVector::DenseMatrix* blockC);
  
      /*! 
       * \brief Computed extended matrix-multivector product
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::MultiVector& X,
	    const NOX::Abstract::MultiVector::DenseMatrix& Y,
	    NOX::Abstract::MultiVector& U,
	    NOX::Abstract::MultiVector::DenseMatrix& V) const;

      /*! 
       * \brief Computed extended matrix transpose-multivector product
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyTranspose(const NOX::Abstract::MultiVector& X,
		     const NOX::Abstract::MultiVector::DenseMatrix& Y,
		     NOX::Abstract::MultiVector& U,
		     NOX::Abstract::MultiVector::DenseMatrix& V) const;

      /*! 
       * \brief Solves the extended system as defined above using bordering
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(NOX::Parameter::List& params,
		   bool isZeroF,
		   bool isZeroG,
		   bool contiguousRHS,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      //! stores bordere system method
      string method;

      //! stores pointer to current bordered system object
      LOCA::BorderedSystem::Generic *borderedPtr;

    };
  } // namespace BorderedSystem
} // namespace LOCA

#endif
