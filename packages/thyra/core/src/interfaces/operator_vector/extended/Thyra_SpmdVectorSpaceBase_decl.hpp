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

#ifndef THYRA_SPMD_VECTOR_SPACE_BASE_BASE_DECL_HPP
#define THYRA_SPMD_VECTOR_SPACE_BASE_BASE_DECL_HPP

#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Teuchos_Comm.hpp"


namespace Thyra {


/** \brief Base abstract <tt>%VectorSpaceBase</tt> class for all SPMD-based
 * vector spaces.
 *
 * This interface defines a basic mechanism for the interoperability of SPMD
 * <tt>VectorBase</tt> and <tt>MultiVectorBase</tt> objects.  See the base
 * classes <tt>SPMDVetorBase</tt> and <tt>SpmdMultiVectorBase</tt> to
 * demonstrate how this interface is used to achieve universal
 * interoperability of vector and multi-vector objects.
 *
 * Specifically, these classes are designed to handle three different
 * use case:
 * <ul>
 * <li> Serial vectors in one process with <tt>localSubDim()==dim()</tt>.
 *      In this case, the <tt>SPMD_Comm</tt> returned from <tt>getComm()</tt> can
 *      be <tt>SPMD_COMM_NULL</tt>.
 * <li> Distributed parallel vectors in two or more processes with
 *      <tt>localSubDim() < dim()</tt>.  This implementation assumes
 *      that all of the elements are stored contiguously in each
 *      process and that there are no ghost elements as described
 *      below.
 * <li> Locally replicated vectors in one or more process.  This
 *      case is similar to serial vectors except that we have
 *      <tt>localSubDim()==dim()</tt> even if there is more than
 *      one process.  This case is checked for so that the reduction
 *      operations are performed correctly.
 * </ul>
 *
 * This interface provides all the information necessary to implement the
 * <tt>SpmdVectorBase::applyOp()</tt> and
 * <tt>SpmdMultiVectorBase::applyOp()</tt> functions in all of the above
 * described use cases.  This interface returns an SPMD communicator (of which
 * all compatible vector spaces must have the same communicator obviously)
 * through the method <tt>getComm()</tt>.
 *
 * <A NAME="SpmdVectorSpaceBase_Vector_layout"></A>
 * <b>%VectorBase data layout:</b>
 *
 * For the case of a distributed parallel vector, this interface base
 * class assumes that vector data is partitioned to processes in
 * contiguous chunks of dense subvectors.  To spell this out, let
 * <tt>v</tt> be the local vector that is sorted on this process and
 * let <tt>g</tt> be the global vector.  Then these two vectors are
 * related (using zero-based indexing) as:

 \verbatim
    v(k) == g(this->localOffset() + k), for k = 0...this->localSubDim()-1
 \endverbatim

 * Any type of mapping of vector data to processes that can not be
 * interpreted in this way can not rely on this base class for
 * interoperability.  Note that as long as the elements in a process
 * are partitioned to unique processes and no ghost elements are
 * present, the actual indexes used by the application with these
 * vectors is immaterial.  The indexes associated with this set of
 * interfaces, however, are only meaningful to abstract numerical
 * algorithms and provide an arbitrary label for certain types of
 * coordinate-dependent operations (like required in an active-set
 * method for optimization).  Therefore, as long as the underlying
 * vector represents a unique partitioning of elements, these classes
 * can be used.  There is a default implementation of
 * <tt>localOffset()</tt> that automatically assumes this contiguous
 * mapping of elements to processes and in general this should not be
 * changed.
 *
 * <b>Notes to subclass developers:</b>
 *
 * The pure virtual methods <tt>getComm()</tt> and <tt>localSubDim()</tt>
 * declared in this interface along with the pure virtual methods
 * <tt>dim()</tt> and <tt>createMember()</tt> are the only methods that must
 * be overridden.
 *
 * If <tt>this</tt> this is in an uninitialized state then
 * <tt>localSubDim()</tt> should return <tt>0</tt>.
 *
 * If optimized implementations of multi-vectors can be supported,
 * then the <tt>createMembers()</tt> method should also be overridden.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class SpmdVectorSpaceBase : virtual public VectorSpaceBase<Scalar> {
public:
  
  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Returns the SPMD communicator.
   */
  virtual Teuchos::RCP<const Teuchos::Comm<Ordinal> > getComm() const = 0;

  /** \brief Returns the number of local elements stored on this process.
   *
   * If <tt>this</tt> this is uninitialized then <tt>localSubDim()</tt>
   * returns <tt>0</tt>.
   */
   virtual Ordinal localSubDim() const = 0;

  /** \brief Returns the offset for the local sub-vector stored on this
   * process.
   */
  virtual Ordinal localOffset() const = 0;

  /** \brief Returns the code for the mapping of elements to processes.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->localSubDim() > 0</tt>] <tt>this->mapCode() > 0</tt>.
   * <li> [<tt>this->localSubDim() <= 0</tt>] <tt>this->mapCode() <= 0</tt>.
   * </ul>
   */
  virtual Ordinal mapCode() const = 0;

  //@}

private:
  
  // Not defined and not to be called
  SpmdVectorSpaceBase<Scalar>&
  operator=(const SpmdVectorSpaceBase<Scalar>&);
  
};


} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_SPACE_BASE_BASE_DECL_HPP
