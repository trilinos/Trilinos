// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_VECTORSPACE_DECL_HPP
#define THYRA_VECTORSPACE_DECL_HPP

#include "Teuchos_Handle.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_VectorSpaceBaseDecl.hpp"
#include "Thyra_VectorBaseDecl.hpp"
#include "Thyra_VectorDecl.hpp"

namespace Thyra
{

  /** 
   *
   */
  template <class Scalar>
  class VectorSpace : public Teuchos::ConstHandle<VectorSpaceBase<Scalar> >
  {
  public:
    TEUCHOS_CONST_HANDLE_CTORS(VectorSpace<Scalar>, VectorSpaceBase<Scalar>);

    /** Create a new element of this vector space */
    Vector<Scalar>  createMember() const ;

    /** Return the dimension of the space */
    int dim() const {return this->ptr()->dim();}

    /** Return the lowest global index accessible on this processor */
    int lowestLocallyOwnedIndex() const ;

    /** Return the number of elements owned by this processor */
    int numLocalElements() const ;

    /** Check compatibility with another space. Implementation note: 
     * we don't know if the argument vec space is a handle to another
     * vector space or the contents of a handle, and we want the operation
     * to work the same in either case. We can make this work as
     * follows: have the argument check compatibility with the contents
     * of this handle. If the argument is a handle, the process 
     * will be repeated, interchanging places again so that both handles
     * are dereferenced. If the argument is not a handle, then it
     * ends up comparing to the concrete contents of this handle, giving the
     * same results. */
    bool isCompatible(const VectorSpace<Scalar>& vecSpc) const; 


    /** Tell if vectors of this space are in core  */
    bool isInCore() const {return this->ptr()->isInCore();}

   

    /** test equality between two spaces */
    bool operator==(const VectorSpace<Scalar>& other) const ;


    /** test inequality of two spaces */
    bool operator!=(const VectorSpace<Scalar>& other) const ;


    /** test whether the space contains a given vector */
    bool contains(const Vector<Scalar>& vec) const ;


    /** return the number of subblocks. */
    int numBlocks() const ;

    /** get the i-th subblock */
    VectorSpace<Scalar> getBlock(int i) const ;


    /** set the i-th subblock */
    void setBlock(int i, const VectorSpace<Scalar>& space);


  };

  /** */
  template <class Scalar>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > 
  productSpace(const Teuchos::Array<VectorSpace<Scalar> >& spaces);

  /** */
  template <class Scalar>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1);

  /** */
  template <class Scalar>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1, 
               VectorSpace<Scalar>& s2);

  /** */
  template <class Scalar>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1,VectorSpace<Scalar>& s2,
               VectorSpace<Scalar>& s3);
}

#endif
