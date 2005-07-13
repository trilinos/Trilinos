// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_TESTING_TOOLS_HPP
#define THYRA_TESTING_TOOLS_HPP

#include "Thyra_TestingToolsDecl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"

// Utilities

template <bool hasMachineParameters, class Scalar>
class relErrSmallNumber {
public:
  static Scalar smallNumber()
    {
      return Teuchos::ScalarTraits<Scalar>::ThisShouldNotCompile();
    }
};

template <class Scalar>
class relErrSmallNumber<false,Scalar> {
public:
  static Scalar smallNumber()
    {
      return Scalar(1e-8);
    }
};

template <class Scalar>
class relErrSmallNumber<true,Scalar> {
public:
  static Scalar smallNumber()
    {
      return Teuchos::ScalarTraits<Scalar>::eps();
    }
};

//

template <class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::relErr( const Scalar &s1, const Scalar &s2 )
{
  return
    Teuchos::ScalarTraits<Scalar>::magnitude( s1 - s2 )
    / (
      Teuchos::ScalarTraits<Scalar>::magnitude(
        relErrSmallNumber<Teuchos::ScalarTraits<Scalar>::hasMachineParameters,Scalar>::smallNumber()
        )
      + std::max( Teuchos::ScalarTraits<Scalar>::magnitude(s1), Teuchos::ScalarTraits<Scalar>::magnitude(s1) )
      );
}

template<class Scalar>
bool Thyra::testMaxErr(
  const std::string                                             &error_name
  ,const Scalar                                                 &error
  ,const std::string                                            &max_error_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &max_error
  ,const std::string                                            &max_warning_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &max_warning
  ,std::ostream                                                 *out
  ,const std::string                                            &li
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const bool success = ( error <= max_error );
  if(out) {
    *out
      << std::endl
      << li << "Check: " << error_name << " = " << error
      << " <= " << max_error_name << " = " << max_error << " : " << passfail(success) << std::endl;
    if( success && error >= max_warning ) {
      *out
        << li << "Warning! " << error_name << " = " << error
        << " >= " << max_warning_name << " = " << max_warning << "!\n";
    }
  }
  return success;
}

template<class Scalar>
bool Thyra::testRelErr(
  const std::string                                             &v1_name
  ,const Scalar                                                 &v1
  ,const std::string                                            &v2_name
  ,const Scalar                                                 &v2
  ,const std::string                                            &maxRelErr_error_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr_error
  ,const std::string                                            &maxRelErr_warning_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr_warning
  ,std::ostream                                                 *out
  ,const std::string                                            &li
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const ScalarMag rel_err = relErr( v1, v2 );
  const bool success = ( rel_err <= maxRelErr_error );
  if(out) {
    *out
      << std::endl
      << li << "Check: rel_err(" << v1_name << "," << v2_name << ")\n"
      << li << "       = rel_err(" << v1 << "," << v2 << ") "
      << "= " << rel_err << std::endl
      << li << "         <= " << maxRelErr_error_name << " = " << maxRelErr_error << " : " << passfail(success) << std::endl;
    if( success && rel_err >= maxRelErr_warning ) {
      *out
        << li << "Warning! rel_err(" << v1_name << "," << v2_name << ")\n"
        << li << "       = rel_err(" << v1 << "," << v2 << ") " << "= " << rel_err << std::endl
        << li << "         >= " << maxRelErr_warning_name << " = " << maxRelErr_warning << "!\n";
    }
  }
  return success;
}

template<class Scalar>
std::ostream& Thyra::operator<<( std::ostream& o, const VectorBase<Scalar>& v )
{
  return o << describe(v,Teuchos::VERB_HIGH);
}

template<class Scalar>
std::ostream& Thyra::operator<<( std::ostream& o, const LinearOpBase<Scalar>& M )
{
  return o << describe(M,Teuchos::VERB_EXTREME);
}

#endif // THYRA_TESTING_TOOLS_HPP
