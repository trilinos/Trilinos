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

#ifndef _Isorropia_Orderer_hpp_
#define _Isorropia_Orderer_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Isorropia_Operator.hpp>

namespace Isorropia {

/** Interface (abstract base class) for computing a new ordering and
  describing the layout of elements in the new order.

  If the methods which describe the new ordering (e.g., operator[],
  etc.) are called before order() has been called, behavior is not
  well defined. Implementations will either return empty/erroneous
  data, or throw an exception. In most cases, implementations will
  probably call order() internally in a constructor or factory method,
  so this won't usually be an issue.
*/
class Orderer : virtual public Operator {
public:

  /** Destructor */
  virtual ~Orderer() {}

  /** Method which does the work of computing a new ordering.

     \param forceOrdering Optional argument defaults to false.
        Depending on the implementation, compute_partitioning() should
        only perform a repartitioning the first time it is called, and
        subsequent repeated calls are no-ops. If the user's intent is
        to re-compute the partitioning (e.g., if parameters or other
        inputs have been changed), then setting this flag to true
        will force a new partitioning to be computed.
   */
  virtual void order(bool forceOrdering=false) = 0;

   /** Give access of the "direct" permutation vector that is owned by the current
      processor.

      \param size [out] Number of elements in the array.

      \param array [out] Pointer to the the part assignements array inside
                        the object.

      \remark This pointer is only significant if the object still exists.
      Otherwise, you must use \see Isorropia::Operator::extractPartsCopy()

      \sa Isorropia::Operator::extractPropertiesView()
   */
  virtual int extractPermutationView(int& size,
			       const int*& array) const {
    return extractPropertiesView(size, array);
  }


  /** Copy a part of the "direct" permutation vector.

      \param len [in] of the array given by the user.

      \param size [out] Number of elements in the array.

      \param array [out] Direct permutation vector. Allocated by the user with
                        a size of at least @c len elements.

      \remark Memory space which is not useful in the array is not
      initialized or used in this method.

      \sa Isorropia::Operator::extractPropertiesCopy()
   */
  virtual int extractPermutationCopy(int len,
			       int& size,
			       int* array) const {
    return extractPropertiesCopy(len, size, array);
  }


};//class Orderer

}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

