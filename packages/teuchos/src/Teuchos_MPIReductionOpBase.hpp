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

#ifndef _TEUCHOS_MPI_REDUCTION_OP_BASE_HPP_
#define _TEUCHOS_MPI_REDUCTION_OP_BASE_HPP_

#include "Teuchos_RefCountPtr.hpp"

#ifdef HAVE_MPI

#include "mpi.h"

extern "C" {

/** \brief Function that is used to create a <tt>MPI_Op</tt> object to be used
 * in MPI calls.
 *
 * This global static function calls whatever reduction object that is set by
 * <tt>set_reduct_op()</tt> and accessed using <tt>get_reduct_op()</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>get_reduct_op()!=null<tt>
 * </ul>
 */
void Teuchos_MPI_reduction_op(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,MPI_Datatype    *datatype
  );

} // extern "C"

namespace Teuchos {

/** \brief Base class for an MPI-compatible reduction operator.
 *
 * Note, <tt>HAVE_MPI</tt> must be defined to use this class!.
 */
class MPIReductionOpBase {
public:
  /** \brief . */
  virtual ~MPIReductionOpBase() {}
  /** \brief . */
  virtual void reduce(
    void              *invec
    ,void             *inoutvec
    ,int              *len
    ,MPI_Datatype     *datatype
    ) const = 0;
};

/** \brief Sets reduction object that is called by <tt>Teuchos_MPI_reduction_op()</tt>.
 *
 * Preconditions:<ul>
 * <li>[<tt>get_reduct_op()!=null</tt>] <tt>reduct_ob==NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>get_reduct_op()==reduct_op</tt>
 * </ul>
 */
void set_reduct_op( const RefCountPtr<const MPIReductionOpBase>& reduct_op );

/** \brief Return the reduction object that is set by <tt>set_reduct_op()</tt>.
 *
 */
RefCountPtr<const MPIReductionOpBase> get_reduct_op();

/** \brief Utility class for setting an MPI-compatible reduction object and
 * using it to create an MPI_Op object.
 *
 * The destructor to this object will remove the set MPI-compatible reduction
 * operation.
 *
 * Note, this object can only be allocated on the stack and should be used
 * directly before a call to <tt>MPI_allreduce(...)</tt> or
 * <tt>MPI_reduce(...)</tt>.  For example:
 
 \code

  ???

 \endcode

 *
 * Note that this class can only be used in an program where MPI is called
 * from only one thread.
 *
 * Note, <tt>HAVE_MPI</tt> must be defined to use this class!.
 */
class MPIReductionOpCreator {
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
  MPIReductionOpCreator( const Teuchos::RefCountPtr<const MPIReductionOpBase>& reduct_op );
  /** \brief Destroy the created <tt>MPI_Op</tt> object and unset it. */
  ~MPIReductionOpCreator();
  /** \breif Return the created <tt>MPI_Op</tt> reduction object that can be used
   * by MPI.
   *
   * Note, this reduction function object will only be valid while <tt>*this</tt> is still
   * in scope.  Therefore, it is recommended that clients only directly call this
   * function to pass in the returned <tt>MPI_Op</tt> object.
   */
  const MPI_Op& mpi_op() const;
private:
  MPI_Op mpi_op_;
  MPIReductionOpCreator(); // Not defined and not to be called!
};

// /////////////////////////
// Inline defintions

inline
const MPI_Op& MPIReductionOpCreator::mpi_op() const
{
  return mpi_op_;
}

} // namespace Teuchos

#endif // HAVE_MPI

#endif // _TEUCHOS_MPI_REDUCTION_OP_BASE_HPP_
