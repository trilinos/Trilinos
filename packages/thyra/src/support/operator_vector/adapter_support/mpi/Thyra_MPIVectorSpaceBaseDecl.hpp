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

#ifndef THYRA_MPI_VECTOR_SPACE_BASE_DECL_HPP
#define THYRA_MPI_VECTOR_SPACE_BASE_DECL_HPP

#include "Thyra_ScalarProdVectorSpaceBaseDecl.hpp"
#include "RTOp_MPI_config.h"

namespace Thyra {

/** \brief Base <tt>%VectorSpaceBase</tt> class for all MPI-based vector
 * spaces with contiguous local storage.
 *
 * This interface defines a basic mechanism for the interoperability
 * of MPI SPMD <tt>VectorBase</tt> and <tt>MultiVectorBase</tt> objects.  See
 * the base classes <tt>MPIVetorBase</tt> and
 * <tt>MPIMultiVectorBase</tt> to demonstrate how this interface is
 * used to achieve universal interoperability of vector and
 * multi-vector objects.
 *
 * Specifically, these classes are designed to handle three different
 * use case:
 * <ul>
 * <li> Serial vectors on one processor with <tt>localSubDim()==dim()</tt>.
 *      In this case, the <tt>MPI_Comm</tt> returned from <tt>mpiComm()</tt> can
 *      be <tt>MPI_COMM_NULL</tt>.
 * <li> Distributed parallel vectors on two or more processors with
 *      <tt>localSubDim() < dim()</tt>.  This implementation assumes
 *      that all of the elements are stored contiguously on each
 *      processor and that there are no ghost elements as described
 *      below.
 * <li> Locally replicated vectors on one or more processor.  This
 *      case is similar to serial vectors except that we have
 *      <tt>localSubDim()==dim()</tt> even if there is more than
 *      one processor.  This case is checked for so that the reduction
 *      operations are performed correctly.
 * </ul>
 *
 * This interface provides all the information necessary to implement
 * <tt>MPIVectorBase::applyOp()</tt> in all of the above described use
 * cases.  This interface returns an MPI communicator (of which all
 * compatible vector spaces must have the same communicator obviously)
 * through the method <tt>mpiComm()</tt>.
 *
 * <A NAME="MPIVectorSpaceBase_Vector_layout"></A>
 * <b>%VectorBase data layout:</b>
 *
 * For the case of a distributed parallel vector, this interface base
 * class assumes that vector data is partitioned to processors in
 * contiguous chunks of dense subvectors.  To spell this out, let
 * <tt>v</tt> be the local vector that is sorted on this processor and
 * let <tt>g</tt> be the global vector.  Then these two vectors are
 * related (using zero-based indexing) as:

 \verbatim
    v(k) == g(this->localOffset() + k), for k = 0...this->localSubDim()-1

 \endverbatim

 * Any type of mapping of vector data to processors that can not be
 * interpreted in this way can not rely on this base class for
 * interoperability.  Note that as long as the elements in a processor
 * are partitioned to unique processors and no ghost elements are
 * present, the actual indexes used by the application with these
 * vectors is immaterial.  The indexes associated with this set of
 * interfaces, however, are only meaningful to abstract numerical
 * algorithms and provide an arbitrary label for certain types of
 * coordinate-dependent operations (like required in an active-set
 * method for optimization).  Therefore, as long as the underlying
 * vector represents a unique partitioning of elements, these classes
 * can be used.  There is a default implementation of
 * <tt>localOffset()</tt> that automatically assumes this contiguous
 * mapping of elements to processors and in general this should not be
 * changed.
 *
 * <b>Notes to subclass developers:</b>
 *
 * The pure virtual methods <tt>mpiComm()</tt> and
 * <tt>localSubDim()</tt> defined in this interface along with the
 * pure virtual methods <tt>dim()</tt> and <tt>createMember()</tt> are
 * the only methods that must be overridden.
 *
 * If <tt>this</tt> this is in an uninitialized state then
 * <tt>localSubDim()</tt> should return <tt>0</tt>.
 *
 * It should never be necessary to override the virtual functions
 * <tt>mapCode()</tt> and <tt>isCompatible()</tt> as these functions
 * have very good and very general implementations.  do.
 *
 * If optimized implementations of multi-vectors can be supported,
 * then the <tt>createMembers()</tt> method should also be overridden.
 *
 * This class defines a very general default implementation for
 * <tt>smallVecSpcFcty()</tt> that returns a
 * <tt>MPIVectorSpaceFactoryStd</tt> object.  This returned object
 * creates <tt>MPIVectorSpaceStd</tt>
 * objects. <tt>MPIVectorSpaceStd</tt> creates <tt>MPIVectorStd</tt>
 * and <tt>MPIMultiVectorStd</tt>.  This implementation is very
 * general should be very appropriate for many different concrete
 * implementations.
 *
 * <b>Note:</b> It is very important that subclasses call the
 * <tt>updateState()</tt> function whenever the state of
 * <tt>*this</tt> changes in a way that might affect the behavior of
 * any of the public member functions.  For example, if a different
 * value of <tt>localSubDim()</tt> will be returned the next time it
 * is called by a client, then <tt>%updateState()</tt> needs to be
 * called by the subclass.  External clients should never need to
 * worry about this function and that is why <tt>%updateState()</tt>
 * is declared protected.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_support_grp
 */
template<class Scalar>
class MPIVectorSpaceBase : public ScalarProdVectorSpaceBase<Scalar> {
public:

  /** \brief . */
  MPIVectorSpaceBase();
  
  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Returns the MPI communicator.
   */
  virtual MPI_Comm mpiComm() const = 0;
  /** \brief Returns the number of local elements stored on this processor.
   *
   * If <tt>this</tt> this is uninitialized then <tt>localSubDim()</tt>
   * returns <tt>0</tt>.
   */
   virtual Index localSubDim() const = 0;

  //@}

  /** @name Virtual methods with default implementations */
  //@{

  /** \brief Returns the offset for the local sub-vector stored on this
   * processor.
   *
   * This method has a default implementation which just assigns this offset
   * based on counting up <tt>localSubDim()</tt> on each processor and then
   * setting <tt>localOffset()</tt> by the rank of the processor.  For
   * example, if there are 5 elements in process 0 and 4 elements in process
   * rank, then <tt>localOffset</tt> on each of these processors will be set
   * as: <tt>localOffset=0</tt> on process 0, <tt>localOffset=5</tt> on
   * process 1, <tt>localOffset=9</tt> on process 2 and so on.
   */
  virtual Index localOffset() const;
  /** \brief Returns the code for the mapping of elements to processors.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->localSubDim() > 0</tt>] <tt>this->mapCode() > 0</tt>.
   * <li> [<tt>this->localSubDim() <= 0</tt>] <tt>this->mapCode() <= 0</tt>.
   * </ul>
   *
   * This method takes the data <tt>mpiComm()</tt>, <tt>numProc</tt>
   * (where <tt>numProc</tt> is returned from
   * <tt>MPI_Comm_size(this->mpiComm(),&numProc)</tt>,
   * <tt>localOffset()</tt> or <tt>localSubDim()</tt> on each
   * processor and then uses it to compute a value for
   * <tt>mapCode</tt> (using a single global reduction if
   * <tt>numProc > 1</tt>) which is returned from this function.
   *
   * The value returned from this default implementation of this
   * method must not be changed or this approach breaks down.  The
   * only reason for overriding this method is for the subclass to
   * be alerted of <em>when</em> this method is called but not
   * <em>what</em> is returned from this method.  If a subclass
   * developer does not understand what this means then <b>don't</b>
   * override this method!
   *
   * The default implementation will always return <tt>return >
   * 0</tt> (unless <tt>this</tt> is uninitialized) so that if this
   * method is overridden to return <tt>return <= </tt> then this is
   * a flag that the underlying vector map does not satisfy the
   * assumptions of this vector space interface and vectors that are
   * in <tt>*this</tt> vector space can not collaborate with other
   * MPI-based vector implementations.
   */
  virtual Index mapCode() const;

  //@}

  /** @name Overridden from VectorSpaceBase */
  //@{

  /** \brief Returns the sum of the local number of elements on every process.
   */
  Index dim() const;

  /** \brief Returns true if all of the elements are stored on one processor.
   *
   * Postconditions:<ul>
   * <li> <tt>return == (dim()==localSubDim())</tt>.
   * </ul>
   */
   bool isInCore() const;

  /** \brief Returns a <tt>MPIVectorSpaceFactoryStd</tt> object that has been given <tt>mpiComm()</tt>.
   */
  Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const;

  /** \brief Checks the general compatibility of parallel (or serial on one
   * processor) MPI-based vector spaces.
   *
   * @return Returns true if <tt>*this</tt> and <tt>vecSpace</tt>
   * are both serial in-core vectors or if <tt>vecSpc</tt> is of
   * type <tt>MPIVectorSpaceBase<Scalar></tt> and both <tt>*this</tt>
   * and <tt>vecSpc</tt> have the same MPI communicators and the same
   * mapping of elements to processors.
   *
   * Postconditions:<ul>
   * <li> <tt>return = ( this->isInCore() &&
   *      vecSpc->isInCore() ) || ( (mpiVecSpc = dynamic_cast<const
   *      MPIVectorSpaceBase<Scalar>*>(&vecSpc)) && this->mpiComm() ==
   *      mpiVecSpc->mpiComm() && this->mapCode() ==
   *      mpiVecSpc->mpiComm())</tt>.
   *
   * </ul>
   *
   * If the mapping of vector elements to processors is not as
   * described
   * <A HREF="classThyra_1_1MPIVectorSpaceBase.html#MPIVectorSpaceBase_Vector_layout">above</A>
   * then this method should be overridden in a way that is specific
   * to the vector implementation.
   */
   bool isCompatible(const VectorSpaceBase<Scalar>& vecSpc) const;
  
  //@}

protected:

  /** \brief This function must be called whenever the state of
   * <tt>this</tt> changes and some internal state must be updated.
   *
   * @param  globalDim  [in] If <tt>globalDim > 0</tt> then this determines
   *                    the global dimension of the vector space.  If <tt>globalDim==this->localSubDim()</tt>
   *                    then this is a locally replicated vector space.  If <tt>globalDim < 0</tt> then
   *                    the global dimension is computed using a global reduction.
   *                    If <tt>MPI_Comm_size(this->mpiComm(),&numProc)</tt> returns <tt>numProc==1</tt>
   *                    then this argument is ignored.
   *
   * Note that calling this function will involve one or more global
   * reductions being called if this is parallel vector space so it
   * should only be called when needed by subclasses.
   *
   * Usually, this operation only needs to be called once for every
   * *new* parallel vector space constructed and very few parallel
   * vector spaces will be created per application usually.
   */
  virtual void updateState( const Index globalDim );

private:

  // //////////////////////////////////////
  // Private data members

  Index     mapCode_;    // < 0 is a flag that everything needs initialized
  bool      isInCore_;
  Index     defaultLocalOffset_;
  Index     defaultGlobalDim_;

  Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> >  smallVecSpcFcty_;
  
}; // end class MPIVectorSpaceBase

} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_SPACE_BASE_DECL_HPP
