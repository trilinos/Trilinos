// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SPMD_VECTOR_BASE_DECL_HPP
#define THYRA_SPMD_VECTOR_BASE_DECL_HPP


#include "Thyra_VectorBase.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"


namespace Thyra {


/** \brief Base class for SPMD vectors that can provide views of contiguous
 * elements in a process.
 *
 * By inheriting from this base class, vector implementations allow their
 * vector objects to be seamlessly combined with other SPMD vector objects (of
 * potentially different concrete types) in <tt>applyOp()</tt>.  A big part of
 * this protocol is that every vector object can expose an
 * <tt>SpmdVectorSpaceBase</tt> object through the function
 * <tt>spmdSpace()</tt>.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class SpmdVectorBase :
    virtual public VectorBase<Scalar>,
    virtual public SpmdMultiVectorBase<Scalar>
{
public:

  /** \brief . */
  using SpmdMultiVectorBase<Scalar>::getNonconstLocalData;
  /** \brief . */
  using SpmdMultiVectorBase<Scalar>::getLocalData;

  /** @name Public non-virtual interface functions */
  //@{

  /** \brief Get a non-const generalized view of local vector data.
   */
  RTOpPack::SubVectorView<Scalar> getNonconstLocalSubVector()
    { return getNonconstLocalSubVectorImpl(); }

  /** \brief Get a const generalized view of local vector data.
   */
  RTOpPack::ConstSubVectorView<Scalar> getLocalSubVector() const
    { return getLocalSubVectorImpl(); }

  /** \brief Returns a non-<tt>const</tt> pointer to the beginning of the
   * local vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to an
   * array of the local values.
   *
   * Preconditions:<ul>
   * <li> <tt>nonnull(localValues)==true</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>nonnull(*localValues)==true</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be freed by removing
   * all of the <tt>ArrayRCP</tt> objects (or setting them to null).
   */
  void getNonconstLocalData(const Ptr<ArrayRCP<Scalar> > &localValues)
    { this->getNonconstLocalVectorDataImpl(localValues); }

  /** \brief Returns a <tt>const</tt> pointer to the beginning of the local
   * vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to an
   * array of the local values.
   *
   * Preconditions:<ul>
   * <li> <tt>nonnull(localValues)==true</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>nonnull(*localValues)==true</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be freed by removing
   * all of the <tt>ArrayRCP</tt> objects (or setting them to null).
   */
  void getLocalData(const Ptr<ArrayRCP<const Scalar> > &localValues) const
    { this->getLocalVectorDataImpl(localValues); }

  //@}

protected:

  /** \name Protected pure virtual functions to be overridden */
  //@{

  /** \brief Virtual implementation for getNonconstLocalSubVector(). */
  virtual RTOpPack::SubVectorView<Scalar>
  getNonconstLocalSubVectorImpl() = 0;

  /** \brief Virtual implementation for getLocalSubVector(). */
  virtual RTOpPack::ConstSubVectorView<Scalar>
  getLocalSubVectorImpl() const = 0;

  /** \brief Implementation of getNonconstLocalData() */
  virtual void getNonconstLocalVectorDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues) = 0;

  /** \brief Implementation of getLocalData() */
  virtual void getLocalVectorDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues) const = 0;

  //@}


}; // end class SpmdVectorBase


} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_BASE_DECL_HPP
