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

#ifndef Rythmos_INTERPOLATION_BUFFER_BASE_H
#define Rythmos_INTERPOLATION_BUFFER_BASE_H

#include "Thyra_VectorBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Rythmos {

/** \brief Base class for defining interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBufferBase 
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolationBufferBase<Scalar> >
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /// Destructor
    virtual ~InterpolationBufferBase() {};

    /// Add points to buffer
    virtual bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec
      ) =0; 

    /// Get values from buffer
    virtual bool GetPoints(
      const std::vector<Scalar>& time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
      ,std::vector<ScalarMag>* accuracy_vec
      ) const =0;

    /// Fill data in from another interpolation buffer
    /// The time points in the interval [time_lower,time_upper] will be
    /// inserted into the interpolation buffer.
    virtual bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar>& IB
      ) =0;

    /// Get interpolation nodes
    virtual bool GetNodes(std::vector<Scalar>* time_vec) const =0;

    /// Remove interpolation nodes
    virtual bool RemoveNodes(std::vector<Scalar>& time_vec) =0;

    /// Get order of interpolation
    virtual int GetOrder() const =0;

};

} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_BASE_H
