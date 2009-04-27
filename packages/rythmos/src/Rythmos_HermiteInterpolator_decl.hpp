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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_HERMITE_INTERPOLATOR_DECL_H
#define Rythmos_HERMITE_INTERPOLATOR_DECL_H

#include "Rythmos_InterpolatorBase.hpp"
#include "Rythmos_Types.hpp"

namespace Rythmos {
/** 
  This class implements piecewise Hermite interpolation on each interval where the data is:
  (t0,x(t0)), (t1,x(t1)), (t0,x'(t0)), (t1,x'(t1))
  The Hermite Interpolation polynomial is:
  H_3(t) = x[z0] + x[z0,z1](t-t0) + x[z0,z1,z2](t-t0)^2 + x[z0,z1,z2,z3](t-t0)^2(t-t1)
  where z0 = z1 = t0 and z2 = z3 = t1 and
  x[z0,z1] = x'(t0) and x[z2,z3] = x'(t1)
  This reduces to:
  H_3(t) = x(t0) + x'(t0)(t-t0) + ((x(t1)-x(t0))/(t1-t0) - x'(t0))(t-t0)^2/(t1-t0)
           +(x'(t1) - 2(x(t1)-x(t0))/(t1-t0) + x'(t0))(t-t0)^2(t-t1)/(t1-t0)^2
  With derivative:
  H_3'(t) =        x'(t0) + 2*((x(t1)-x(t0))/(t1-t0) - x'(t0))(t-t0)/(t1-t0)
           +(x'(t1) - 2(x(t1)-x(t0))/(t1-t0) + x'(t0))[2*(t-t0)(t-t1) + (t-t0)^2]/(t1-t0)^2
  With the error expression:
  x(t) - H_3(t) = (x^{(3)}(\xi(t))/(4!))(t-t0)^2(t-t1)^2
  The Hermite Interpolant will match 3rd degree polynomials exactly with both function values and derivatives.
  **/
template<class Scalar>
class HermiteInterpolator : virtual public InterpolatorBase<Scalar>
{
  public:

    /// Destructor
    ~HermiteInterpolator() {};

    /// Constructor
    HermiteInterpolator();

    void setNodes(
      const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes
      );

    /// Interpolation:
    /** \brief Hermite interpolation function.
     *
     * <b>Preconditions:</b><ul>
     * <li>Preconditions of InterpolatorBase<Scalar> apply 
     * <li><tt>(*nodes_)[i].xdot != Teuchos::null</tt> for all <tt>i=0..(*nodes_).size()-1</tt>
     * </ul>
     */
    void interpolate(
      const Array<Scalar> &t_values,
      typename DataStore<Scalar>::DataStoreVector_t *data_out
      ) const;

    /// Order of interpolation:
    int order() const; 

    /// Inherited from Teuchos::Describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream &out,
      const Teuchos::EVerbosityLevel verbLevel
      ) const;

    /// Redefined from ParameterListAcceptor
    /** \brief . */
    void setParameterList(RCP<ParameterList> const& paramList);

    /** \brief . */
    RCP<ParameterList> getNonconstParameterList();

    /** \brief . */
    RCP<ParameterList> unsetParameterList();

    /** \brief. */
    RCP<const Teuchos::ParameterList> getValidParameters() const;

    void assertInterpolatePreconditions(
        const typename DataStore<Scalar>::DataStoreVector_t &data_in
        ,const Array<Scalar> &t_values
        ,typename DataStore<Scalar>::DataStoreVector_t *data_out
        ) const;

  private:

    RCP<const typename DataStore<Scalar>::DataStoreVector_t> nodes_;

    RCP<ParameterList> parameterList_;

};

// non-member constructor
template<class Scalar>
RCP<HermiteInterpolator<Scalar> > hermiteInterpolator();

} // namespace Rythmos

#endif // Rythmos_HERMITE_INTERPOLATOR_DECL_H
