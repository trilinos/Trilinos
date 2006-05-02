// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_OPERATOR_HPP
#define TPETRA_OPERATOR_HPP

#include "Tpetra_Object.hpp"

namespace Tpetra {

  template<typename OrdinalType, typename ScalarType> class Vector;
  template<typename OrdinalType, typename ScalarType> class VectorSpace;

  /** \brief Abstract interface for linear operators that accept Tpetra
   * vectors.
   */
  template<typename OrdinalType, typename ScalarType>
	class Operator : virtual public Object {
	public:
  
		/** \name Pure virtual functions to be overridden by subclasses. */
    //@{

		//! Returns the VectorSpace associated with the domain of this linear operator.
		virtual VectorSpace<OrdinalType,ScalarType> const& getDomainDist() const = 0;
  
		//! Returns the VectorSpace associated with the range of this linear operator.
		virtual VectorSpace<OrdinalType,ScalarType> const& getRangeDist() const = 0;

    //! Computes the matrix-vector multiplication y = Ax.
		virtual void apply(Vector<OrdinalType,ScalarType> const& x, Vector<OrdinalType, ScalarType> & y, bool transpose=false) const = 0;
    
    //@}
    
	};

} // Tpetra namespace

#endif // TPETRA_OPERATOR_HPP
