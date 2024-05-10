//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_Partitioner2D_hpp_
#define _Isorropia_Partitioner2D_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Isorropia_Partitioner.hpp>

namespace Isorropia {

/** Interface (abstract base class) for computing a new 2D partitioning and
  describing the layout of elements in the new partitions.

  If the methods which describe the new partitioning (e.g., 
  newPartitionNumber(), etc.) are called before compute_partitioning()
  has been called, behavior is not well defined. Implementations will
  either return empty/erroneous data, or throw an exception. In most
  cases, implementations will probably call compute_partitioning()
  internally in a constructor or factory method, so this won't usually
  be an issue.
*/
class Partitioner2D : public Isorropia::Partitioner
{
public:

  /** Destructor */
  virtual ~Partitioner2D() {}


  /** Method which does the work of computing a new partitioning.
     Implementations of this interface will typically be constructed
     with an object or information describing the existing ('old')
     partitioning. This method computes a 'new' rebalanced
     partitioning for that input data.

     \param force_repartitioning Optional argument defaults to false.
        Depending on the implementation, compute_partitioning() should
        only perform a repartitioning the first time it is called, and
        subsequent repeated calls are no-ops. If the user's intent is
        to re-compute the partitioning (e.g., if parameters or other
        inputs have been changed), then setting this flag to true
        will force a new partitioning to be computed.
   */
  virtual void partition(bool force_repartitioning=false) = 0;



  /** Return the number of LOCAL elements in a given part.  

      \param[in] part the part ID we want to know the number of local
      elements.

      \return number of local elements that belongs to the
      given part.

      \sa Isorropia::Operator::numElemsWithProperty()
  */
  virtual  int numElemsInPart(int part) const = 0;


  /** Fill user-allocated list (of length len) with the
      local element ids to be located in the given part

      \param[in] part the part ID we consider

      \param[out] elementList array of elements that belongs to this
      part ID, must be allocated by user with size at least @c len

      \param[in] len maximum number of elements we can put in the
      array. Usually, may be the result of
      Isorropia::Partitioner::numElemsInPart().  .

      \sa Isorropia::Operator::elemsWithProperty()
  */
  virtual  void elemsInPart(int part,
                            int* elementList,
                            int len) const = 0;



//   /** Return the new partition ID for a given element that
//      resided locally in the old partitioning.
//   */

//   virtual int newPartitionNumber(int myElem) const = 0;





};//class Partitioner2D

}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

