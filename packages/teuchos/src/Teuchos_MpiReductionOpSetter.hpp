// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP
#define TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP

#include "Teuchos_RCP.hpp"
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
class TEUCHOS_LIB_DLL_EXPORT MpiReductionOpBase : virtual public Describable {
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
  MpiReductionOp( const RCP<const ValueTypeReductionOp<Ordinal,char> > &reductOp );
  /** \brief . */
  void reduce(
    void              *invec
    ,void             *inoutvec
    ,int              *len
    ,MPI_Datatype     *datatype
    ) const;
private:
  RCP<const ValueTypeReductionOp<Ordinal,char> > reductOp_;
  // Not defined and not to be called
  MpiReductionOp();
  MpiReductionOp(const MpiReductionOp&);
  MpiReductionOp& operator=(const MpiReductionOp&);
};

/** \brief Create an <tt>MpiReductionOp</tt> object given an
 * <tt>ReductionOp</tt> object.
 */
template<typename Ordinal>
RCP<const MpiReductionOp<Ordinal> >
mpiReductionOp( const RCP<const ValueTypeReductionOp<Ordinal,char> > &reductOp );

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
class TEUCHOS_LIB_DLL_EXPORT MpiReductionOpSetter {
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
  MpiReductionOpSetter( const Teuchos::RCP<const MpiReductionOpBase>& reduct_op );

  /** \brief . */
  ~MpiReductionOpSetter();

  /** \brief Return the created <tt>MPI_Op</tt> reduction object that can be used
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
  const RCP<const ValueTypeReductionOp<Ordinal,char> > &reductOp
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
  (void)datatype;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!invec);
  TEUCHOS_TEST_FOR_EXCEPT(!inoutvec);
  TEUCHOS_TEST_FOR_EXCEPT(!len);
  TEUCHOS_TEST_FOR_EXCEPT(!(*len>0));
  TEUCHOS_TEST_FOR_EXCEPT(!datatype);
  //TEUCHOS_TEST_FOR_EXCEPT(!(*datatype==MPI_CHAR));
  // We also allow datatypes that are blocks of chars!
#endif
  int sz;
  MPI_Type_size(*datatype, &sz);
  reductOp_->reduce(
    *len*sz,reinterpret_cast<char*>(invec),reinterpret_cast<char*>(inoutvec)
    );
}

template<typename Ordinal>
RCP<const MpiReductionOp<Ordinal> >
mpiReductionOp( const RCP<const ValueTypeReductionOp<Ordinal,char> > &reductOp )
{
  return rcp(new MpiReductionOp<Ordinal>(reductOp));
}

} // namespace Teuchos

#endif // TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP
