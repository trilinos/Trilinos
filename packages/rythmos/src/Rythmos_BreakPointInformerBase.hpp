//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_BREAK_POINT_INFORMER_H
#define Rythmos_BREAK_POINT_INFORMER_H

namespace Rythmos {

/** \brief Interface for using breakpoints.   */
template<class Scalar>
class BreakPointInformerBase 
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<BreakPointInformerBase<Scalar> >
{
  public:
    /** \brief Basic breakpoint interface. 
     *
     *  Returns true if there is a breakpoint between time0 and time0+dt
     *
     */ 
    virtual bool testForBreakPoint(Scalar& time0, Scalar& dt) const =0;

    /** \brief Get next break point. 
     *
     *  Returns the next breakpoint after time0.
     *
     */
    virtual Scalar getNextBreakPoint(Scalar& time0) const =0;

    /** \brief Remove breakpoint from list. 
     *
     *  Removes the next breakpoint after time0.
     */
    virtual void removeNextBreakPoint(Scalar& time0) =0;

};

} // namespace Rythmos

#endif //Rythmos_BREAK_POINT_INFORMER_H
