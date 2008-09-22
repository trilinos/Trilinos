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

#ifndef _Isorropia_Colorer_hpp_
#define _Isorropia_Colorer_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Isorropia_Operator.hpp>

namespace Isorropia {

/** Interface (abstract base class) for computing a new coloring and
    describing the result.

    The colors returned have values between 1 and C, where C is the number of colors used.

*/
class Colorer : virtual public Operator {
public:

  /** Destructor */
  virtual ~Colorer() {}


  /** Method which does the work of computing a new coloring.

     \param forceColoring Optional argument defaults to false.
       Depending on the implementation, color() should only perform a
       coloring the first time it is called, and subsequent repeated
       calls are no-ops. If the user's intent is to re-compute the
       coloring (e.g., if parameters or other inputs have been
       changed), then setting this flag to true will force a new
       coloring to be computed.
   */
  virtual void color(bool forceColoring=false) = 0;


  /** Method which returns the number (global) of colors used.

      \return The overall number of colors used. All colors used for
      all vertices are between 1 and this value (included).

      \sa Isorropia::Operator::numProperties()
   */
  virtual int numColors() const {
      return numProperties(); }


  /** Return the number of \b local elements of a given color.

      \param color The wanted color.

      \return The number of \b local of the asked color.

      \sa Isorropia::Operator::numElemsWithProperty()
   */
  virtual int numElemsWithColor(int color) const
  { return numElemsWithProperty(color); }


  /** Fill user-allocated list (of length len) with the
   *  local element ids for LOCAL elements of the given color.

      \param color the wanted color

      \param elementList an array to receive local elements of the given color

      \param len the number of elements wanted

      \sa Isorropia::Operator::elemsWithProperty()
   */
  virtual void elemsWithColor(int color,
			      int* elementList,
			      int len) const {
    return elemsWithProperty(color, elementList, len);}

};//class Colorer

}//namespace Isorropia

#endif

