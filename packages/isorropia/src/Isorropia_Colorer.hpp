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

      \param theColor The wanted color.

      \return The number of \b local of the asked color.

      \sa Isorropia::Operator::numElemsWithProperty()
   */
  virtual int numElemsWithColor(int theColor) const
  { return numElemsWithProperty(theColor); }


  /** Fill user-allocated list (of length len) with the
   *  local element ids for LOCAL elements of the given color.

      \param theColor the wanted color

      \param elementList an array to receive local elements of the given color

      \param len the number of elements wanted

      \sa Isorropia::Operator::elemsWithProperty()
   */
  virtual void elemsWithColor(int theColor,
                              int* elementList,
                              int len) const {
    return elemsWithProperty(theColor, elementList, len);}

  /** Give access of the color assignments array that is owned by the current
      processor.

      \param size [out] Number of elements in the array.

      \param array [out] Pointer to the color assignements array inside
                        the object.

      \remark This pointer is only significant if the object still exists.
      Otherwise, you must use \see Isorropia::Operator::extractPartsCopy()

      \sa Isorropia::Operator::extractPropertiesView()
   */
  virtual int extractColorsView(int& size,
                               const int*& array) const {
    return extractPropertiesView(size, array);
  }


  /** Copy a part of the color assignments array.

      \param len [in] of the array given by the user.

      \param size [out] Number of elements in the array.

      \param array [out] Array of color assignments. Allocated by the user with
                        a size of at least @c len elements.

      \remark Memory space which is not useful in the array is not
      initialized or used in this method.

      \sa Isorropia::Operator::extractPropertiesCopy()
   */
  virtual int extractColorsCopy(int len,
                               int& size,
                               int* array) const {
    return extractPropertiesCopy(len, size, array);
  }


};//class Colorer

}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

