// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP
#define TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "mpi.h"

namespace Teuchos {

/** \brief Base class for an MPI-compatible reduction operator.
 *
 * The base class allows clients to create a reduction callback as an object
 * instead of just as a global function.  The single extra level of
 * indirection should not be a performance problem in must cases.
 *
 * Note, <tt>HAVE_MPI</tt> must be defined to use this class!.
 */
class MpiReductionOpBase : Describable {
public:
  /** \brief . */
  virtual void reduce(
    void              *invec
    ,void             *inoutvec
    ,int              *len
    ,MPI_Datatype     *datatype
    ) const = 0;
};

/** \brief Standard subclass implementation for <tt>MpiReductionOpBase</tt>
 * in terms of a templated <tt>ReductionOp<Ordinal,char></tt> object.
 */
template<typename Ordinal>
class MpiReductionOp : public MpiReductionOpBase {
public:
  /** \brief . */
  MpiReductionOp( const RefCountPtr<const ValueTypeReductionOp<Ordinal,char> > &reductOp );
  /** \brief . */
  void reduce(
    void              *invec
    ,void             *inoutvec
    ,int              *len
    ,MPI_Datatype     *datatype
    ) const;
private:
  RefCountPtr<const ValueTypeReductionOp<Ordinal,char> > reductOp_;
  // Not defined and not to be called
  MpiReductionOp();
  MpiReductionOp(const MpiReductionOp&);
  MpiReductionOp& operator=(const MpiReductionOp&);
};

/** \brief Create an <tt>MpiReductionOp</tt> object given an
 * <tt>ReductionOp</tt> object.
 */
template<typename Ordinal>
RefCountPtr<const MpiReductionOp<Ordinal> >
mpiReductionOp( const RefCountPtr<const ValueTypeReductionOp<Ordinal,char> > &reductOp );

/** \brief Utility class for setting an MPI-compatible reduction object and
 * using it to create an <tt>MPI_Op</tt> object.
 *
 * The destructor to this object will remove the set MPI-compatible reduction
 * operation.
 *
 * Note, this object can only be allocated on the stack and should be used
 * directly before a call to any MPI function that takes an <tt>MPI_Op</tt>
 * object.  For example:
 
 \code

  ???

 \endcode

 *
 * Note that this class can only be used in a program where MPI is called
 * from only one thread.
 *
 * Note, <tt>HAVE_MPI</tt> must be defined to use this class!.
 */
class MpiReductionOpSetter {
public:

  /** \brief Construct a new <tt>MPI_Op</tt> object which uses the passed in
   * <tt>reduct_op</tt> objects.
   *
   * Preconditions:<ul>
   * <li><tt>reduct_op.get()!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->mpi_op()</tt> gives access to the created <tt>MPI_Op</tt>
   *     object that can be used with MPI.
   * </ul>
   */
  MpiReductionOpSetter( const Teuchos::RefCountPtr<const MpiReductionOpBase>& reduct_op );

  /** \brief . */
  ~MpiReductionOpSetter();

  /** \breif Return the created <tt>MPI_Op</tt> reduction object that can be used
   * by MPI.
   *
   * Note, this reduction function object will only be valid while <tt>*this</tt> is still
   * in scope.  Therefore, it is recommended that clients only directly call this
   * function to pass in the returned <tt>MPI_Op</tt> object.
   */
  MPI_Op mpi_op() const;
  
private:
  // Not defined and not to be called!
  MpiReductionOpSetter();
  MpiReductionOpSetter(const MpiReductionOpSetter&);
  MpiReductionOpSetter& operator=(const MpiReductionOpSetter&);
};

// ///////////////////////////////////////
// Template implemenations

template<typename Ordinal>
MpiReductionOp<Ordinal>::MpiReductionOp(
  const RefCountPtr<const ValueTypeReductionOp<Ordinal,char> > &reductOp
  )
  :reductOp_(reductOp)
{}

template<typename Ordinal>
void MpiReductionOp<Ordinal>::reduce(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,MPI_Datatype     *datatype
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!invec);
  TEST_FOR_EXCEPT(!inoutvec);
  TEST_FOR_EXCEPT(!len);
  TEST_FOR_EXCEPT(!(*len>0));
  TEST_FOR_EXCEPT(!datatype);
  //TEST_FOR_EXCEPT(!(*datatype==MPI_CHAR));
  // We also allow datatypes that are blocks of chars!
#endif
  reductOp_->reduce(
    *len,reinterpret_cast<char*>(invec),reinterpret_cast<char*>(inoutvec)
    );
}

template<typename Ordinal>
RefCountPtr<const MpiReductionOp<Ordinal> >
mpiReductionOp( const RefCountPtr<const ValueTypeReductionOp<Ordinal,char> > &reductOp )
{
  return rcp(new MpiReductionOp<Ordinal>(reductOp));
}

} // namespace Teuchos

#endif // TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP
