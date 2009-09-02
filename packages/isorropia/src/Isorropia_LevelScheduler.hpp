//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
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

