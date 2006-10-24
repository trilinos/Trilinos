//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_INTERPOLATOR_BASE_H
#define Rythmos_INTERPOLATOR_BASE_H

#include "Rythmos_DataStore.hpp"

namespace Rythmos {

/** \brief Base strategy class for interpolation functionality. */
template<class Scalar> 
class InterpolatorBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolatorBase<Scalar> >
{
  public:

    /// Destructor
    virtual ~InterpolatorBase() {};

    /// Interpolation:
    /// This function must support passing node values back out directly,
    /// handling when only one node value is passed in,
    /// and dealing with xdot == Teuchos::null. 
    /// There is no guarantee at this time that data_out will be in the same order as t_values
    virtual bool interpolate(
        const typename DataStore<Scalar>::DataStoreVector_t &data_in
        ,const std::vector<Scalar> &t_values
        ,typename DataStore<Scalar>::DataStoreVector_t *data_out
        ) const =0;

    /// Order of interpolation:
    virtual int order() const =0;

};

} // namespace Rythmos

#endif //Rythmos_INTERPOLATOR_BASE_H
