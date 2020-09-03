// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
//
// @HEADER

#ifndef SACADO_FAD_EXP_GENERALFADTESTINGHELPERS_HPP

#ifndef TEUCHOS_TESTING_HELPERS_HPP
#include "Teuchos_TestingHelpers.hpp"
#endif

#define SACADO_FAD_EXP_GENERALFADTESTINGHELPERS_HPP

namespace Sacado {
  namespace Fad {
  namespace Exp {
    template <typename Storage>
    typename GeneralFad<Storage>::scalar_type relErr( const GeneralFad<Storage> &s1, const GeneralFad<Storage> &s2 );
  }
  }
}

namespace Teuchos {
  using Sacado::Fad::Exp::relErr;
}

namespace Sacado {
  namespace Fad {
  namespace Exp {

    template <typename Scalar1, typename Scalar2>
    typename Scalar1::scalar_type
    relErrFadImpl( const Scalar1 &s1, const Scalar2 &s2 )
    {
      typename Scalar1::scalar_type maxRelErr = Teuchos::relErr(s1.val(), s2.val());
      for (int i=0; i<s1.size(); ++i) {
        typename Scalar1::scalar_type tmpRelErr = Teuchos::relErr(s1.dx(i), s2.dx(i));
        if (tmpRelErr >= maxRelErr)
          maxRelErr = tmpRelErr;
      }
      return maxRelErr;
    }

    template <typename Storage>
    typename GeneralFad<Storage>::scalar_type relErr( const GeneralFad<Storage> &s1, const GeneralFad<Storage> &s2 ) { return relErrFadImpl(s1.derived(),s2.derived()); }
  }
  }
}

namespace Teuchos {
  template <typename T1, typename T2>
  struct TestRelErr<T1,T2,
                    typename std::enable_if<Sacado::IsExpr<T1>::value ||
                                            Sacado::IsExpr<T2>::value >::type >
  {
    typedef typename Sacado::ScalarType<T1>::type scalarType1;
    typedef typename Sacado::ScalarType<T2>::type scalarType2;
    typedef typename std::common_type<scalarType1,scalarType2>::type magnitudeType;
    static bool eval(
      const std::string &v1_name,
      const T1 &v1,
      const std::string &v2_name,
      const T2 &v2,
      const std::string &maxRelErr_error_name,
      const magnitudeType &maxRelErr_error,
      const std::string &maxRelErr_warning_name,
      const magnitudeType &maxRelErr_warning,
      const Ptr<std::ostream> &out
    )
    {
      using std::endl;
      typedef Teuchos::ScalarTraits<magnitudeType> SMT;
      const auto rel_err = Teuchos::relErr( v1.derived(), v2.derived() );
      const bool success = ( !SMT::isnaninf(rel_err) && !SMT::isnaninf(maxRelErr_error)
        && rel_err <= maxRelErr_error );
      if (!is_null(out)) {
        *out
          << endl
          << "Check: rel_err(" << v1_name << ", " << v2_name << ")\n"
          << "       = rel_err(" << v1 << ", " << v2 << ") "
          << "= " << rel_err << endl
          << "         <= " << maxRelErr_error_name
          << " = " << maxRelErr_error << " : " << Teuchos::passfail(success) << endl;
        if( success && rel_err >= maxRelErr_warning ) {
          *out
            << "Warning! rel_err(" << v1_name << ", " << v2_name << ")\n"
            << "       = rel_err(" << v1 << ", " << v2 << ") "
            << "= " << rel_err << endl
            << "         >= " << maxRelErr_warning_name
            << " = " << maxRelErr_warning << "!\n";
        }
      }
      return success;
    }
  };
}

#define TEUCHOS_TEST_FLOATING_NOT_EQUALITY( v1, v2, tol, out, success ) \
  { \
    const bool l_result = Teuchos::testRelErr( \
      #v1, v1, #v2, v2, "tol", tol, "tol", tol, Teuchos::outArg(out) ); \
    if (l_result) (success) = false; \
  }

#endif
