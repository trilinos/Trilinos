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

#ifndef Rythmos_LINEAR_INTERPOLATOR_DECL_H
#define Rythmos_LINEAR_INTERPOLATOR_DECL_H

#include "Rythmos_InterpolatorBase.hpp"
#include "Rythmos_Types.hpp"

namespace Rythmos {


/** \brief Concrete implemenation of <tt>InterpolatorBase</tt> just just does
 * simple linear interploation.
 */
template<class Scalar>
class LinearInterpolator : virtual public InterpolatorBase<Scalar>
{
public:

  /** \brief . */
  ~LinearInterpolator() {};

  /** \brief . */
  LinearInterpolator();

  /** \brief . */
  bool supportsCloning() const;

  /** \brief . */
  RCP<InterpolatorBase<Scalar> > cloneInterpolator() const;

  /** \brief . */
  void setNodes(
    const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes
    );

  /** \brief . */
  void interpolate(
    const Array<Scalar> &t_values,
    typename DataStore<Scalar>::DataStoreVector_t *data_out
    ) const;

  /** \brief . */
  int order() const; 

  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<ParameterList> unsetParameterList();

  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

private:

  RCP<const typename DataStore<Scalar>::DataStoreVector_t> nodes_;

  RCP<ParameterList> parameterList_;

};

// non-member constructor
template<class Scalar>
RCP<LinearInterpolator<Scalar> > linearInterpolator();


} // namespace Rythmos


#endif // Rythmos_LINEAR_INTERPOLATOR_DECL_H
