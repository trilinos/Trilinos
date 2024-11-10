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

#ifndef _Isorropia_LevelScheduler_hpp_
#define _Isorropia_LevelScheduler_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Operator.hpp>

namespace Isorropia {

/** Interface (abstract base class) for an operator that computes
    a partitioning of local elements into levels.  On each process the
    levels begin at 1 and increase by 1.  All elements in the same 
    level can be used in calculation concurrently.  All elements in 
    level i are dependent upon the completion of the calculations involving 
    the elements in level i-1.
*/
class LevelScheduler : virtual public Operator {
public:

  /** Destructor */
  virtual ~LevelScheduler() {}


  /** Method which does the work of computing a new level schedule.

     \param forceScheduling Optional argument defaults to false.
       Depending on the implementation, schedule() should only perform a
       scheduling the first time it is called, and subsequent repeated
       calls are no-ops. If the user's intent is to re-compute the
       scheduling (e.g., if parameters or other inputs have been
       changed), then setting this flag to true will force a new
       scheduling to be computed.
   */
  virtual void schedule(bool forceScheduling=false) = 0;


  /** Method which returns the number of levels. 

      \return The number of levels on the local process. 
              Levels begin at 1 and increase by 1.

      \sa Isorropia::Operator::numProperties()
   */
  virtual int numLevels() const {
      return numLocalProperties(); }

  /** Return the number of \b elements in a given level.

      \param level The wanted level.

      \return The number of \b elements in the level.

      \sa Isorropia::Operator::numElemsWithProperty()
   */
  virtual int numElemsWithLevel(int level) const
  { return numElemsWithProperty(level); }


  /** Fill user-allocated list (of length len) with the local ID for each
    * element in the given level.

      \param level the wanted level

      \param elementList an array to receive local elements of the given level

      \param len the number of elements wanted

      \sa Isorropia::Operator::elemsWithProperty()
   */
  virtual void elemsWithLevel(int level,
			      int* elementList,
			      int len) const {
    return elemsWithProperty(level, elementList, len);}

};//class LevelScheduler

}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

