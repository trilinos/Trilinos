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
#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_VectorBaseDecl.hpp"
#include "Thyra_VectorDecl.hpp"

namespace Thyra
{

  /** \brief Handle class for <tt>VectorSpaceBase</tt>.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
  template <class Scalar>
  class VectorSpace : public Teuchos::ConstHandle<VectorSpaceBase<Scalar> >
  {
  public:
    TEUCHOS_CONST_HANDLE_CTORS(VectorSpace<Scalar>, VectorSpaceBase<Scalar>);

    /** \brief Create a new element of this std::vector space. */
    Vector<Scalar>  createMember() const ;

    /** \brief Create a multivector in this std::vector space.
     *
     * Note: There is not yet a handle class for multi-vectors yet and
     * therefore this function returns <tt>Teuchos::RCP</tt> object to
     * the created raw <tt>MultiVectorBase</tt> object.  In the future, this
     * will be replaced to return a true <tt>Thyra::MultiVector</tt> handle
     * object onces this class is created.
     */
    Teuchos::RCP<MultiVectorBase<Scalar> > createMembers(int n) const ;

    /** \brief Return the dimension of the space. */
    int dim() const {return this->constPtr()->dim();}

    /** \brief Check compatibility with another space.
     *
     * Implementation note: we don't know if the argument vec space is a
     * handle to another std::vector space or the contents of a handle, and we want
     * the operation to work the same in either case. We can make this work as
     * follows: have the argument check compatibility with the contents of
     * this handle. If the argument is a handle, the process will be repeated,
     * interchanging places again so that both handles are dereferenced. If
     * the argument is not a handle, then it ends up comparing to the concrete
     * contents of this handle, giving the same results. */
    bool isCompatible(const VectorSpace<Scalar>& vecSpc) const; 

    /** \brief Tell if vectors of this space are in core.  */
    bool isInCore() const {return this->constPtr()->isInCore();}

    /** \brief Test equality between two spaces */
    bool operator==(const VectorSpace<Scalar>& other) const ;

    /** \brief Test inequality of two spaces */
    bool operator!=(const VectorSpace<Scalar>& other) const ;

    /** \brief Test whether the space contains a given std::vector */
    bool contains(const Vector<Scalar>& vec) const ;

    /** \brief Return the number of std::vector space subblocks. */
    int numBlocks() const ;

    /** \brief Get the k-th subblock where (<tt>0 <= k < numBlocks()</tt>) */
    VectorSpace<Scalar> getBlock(int k) const ;

    /** \brief Set the k-th subblock where (<tt>0 <= k < numBlocks()</tt>) */
    void setBlock(int k, const VectorSpace<Scalar>& space);

  };

  /** \brief Create a member.
   *
   * \relates VectorSpace
   */
  template <class Scalar>
  Vector<Scalar> createMember( const VectorSpace<Scalar> &space )
  { return space.createMember(); }

  /** \brief Create a product space given an array of std::vector spaces blocks.
   *
   * \relates VectorSpace
   */
  template <class Scalar>
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(const Teuchos::Array<VectorSpace<Scalar> >& spaces);

  /** \brief Create a product space given a single std::vector space block.
   *
   * \relates VectorSpace
   */
  template <class Scalar>
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1);

  /** \brief Create a product space given two std::vector space blocks.
   *
   * \relates VectorSpace
   */
  template <class Scalar>
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1, 
               VectorSpace<Scalar>& s2);

  /** \brief Create a product space given three std::vector space blocks.
   *
   * \relates VectorSpace
   */
  template <class Scalar>
  Teuchos::RCP<const VectorSpaceBase<Scalar> > 
  productSpace(VectorSpace<Scalar>& s1,VectorSpace<Scalar>& s2,
               VectorSpace<Scalar>& s3);

  /** \brief Returns the lowest global index accessible on this processor.
   * This function will throw an std::exception if the space is not a
   * SpmdVectorSpaceBase.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>isSPMD(s)==true</tt>
   * </ul>
   *
   * \relates VectorSpace 
   */
  template <class Scalar>
  int lowestLocallyOwnedIndex(const VectorSpace<Scalar>& s) ;

  /** \brief Return the number of elements owned by this processor.  This
   * function will throw an std::exception if the space is not a
   * SpmdVectorSpaceBase.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>isSPMD(s)==true</tt>
   * </ul>
   *
   * \relates VectorSpace 
   */
  template <class Scalar>
  int numLocalElements(const VectorSpace<Scalar>& s) ;

  /** \brief Test whether a specified index is local to this processor.
   *
   * \relates VectorSpace 
   */
  template <class Scalar>
  bool indexIsLocal(const VectorSpace<Scalar>& s, Index i);

  /** \brief Test whether a std::vector space is an atomic SPMD object.
   *
   * \relates VectorSpace 
   */
  template <class Scalar>
  bool isSPMD(const VectorSpace<Scalar>& s);

}

#endif
