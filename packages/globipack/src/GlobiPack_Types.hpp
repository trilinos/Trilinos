/*
// @HEADER
// ***********************************************************************
// 
//    GlobiPack: Collection of Scalar 1D globalizaton utilities
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef GLOBIPACK_TYPES_HPP
#define GLOBIPACK_TYPES_HPP


#include "GlobiPack_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace Teuchos {

/** \brief . */
class ParameterList;

} // namespace Teuchos


namespace GlobiPack {


/** \brief . */
using Teuchos::RCP;
/** \brief . */
using Teuchos::Ptr;
/** \brief . */
using Teuchos::Array;
/** \brief . */
using Teuchos::ArrayView;
/** \brief . */
using Teuchos::ScalarTraits;
/** \brief . */
using Teuchos::ParameterList;


/** \brief Represents the evaluation point of the merit function
 * <tt>phi(alpha)</tt> and/or is derivative <tt>Dphi(alpha)</tt>.
 *
 * If a value has not been set it will be equal to <tt>valNotGiven()</tt>.
 */
template<class Scalar>
struct PointEval1D {
  /** \brief . */
  static Scalar valNotGiven() { return std::numeric_limits<Scalar>::max(); }
  /** \brief . */
  PointEval1D()
    : alpha(valNotGiven()), phi(valNotGiven()), Dphi(valNotGiven())
    {}
  /** \brief . */
  PointEval1D( const Scalar &alpha_in, const Scalar &phi_in,
    const Scalar &Dphi_in = valNotGiven())
    : alpha(alpha_in), phi(phi_in), Dphi(Dphi_in)
    {}
  /** \brief The value of the unknown <tt>alpha</tt>. */
  Scalar alpha;
  /** \brief The value of the merit function <tt>phi(alpha)</tt>. */
  Scalar phi;
  /** \brief The value of the derivative of the merit function
   * <tt>Dphi(alpha)</tt>. */
  Scalar Dphi;
};


namespace Exceptions {


/** \brief Thrown if search direction not a descent direction for the merit
 * function.
 */
class NotDescentDirection : public std::logic_error
{public: NotDescentDirection(const std::string& what_arg) : std::logic_error(what_arg) {}};


} // namespace Exceptions


} // namespace GlobiPack


#endif // GLOBIPACK_TYPES_HPP
