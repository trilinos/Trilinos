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
