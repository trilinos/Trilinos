// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_GENERALFADTESTINGHELPERS_HPP

#include "Sacado_ConfigDefs.h"
#ifdef HAVE_SACADO_TEUCHOSCORE

#define SACADO_FAD_GENERALFADTESTINGHELPERS_HPP

namespace Sacado {
 namespace Fad {
   // relErrFadImpl version for non-Fad types
   template <typename T1, typename T2>
   typename std::enable_if<
     !IsExpr<T1>::value && !IsExpr<T2>::value,
     typename Teuchos::ScalarTraits< typename std::common_type<T1,T2>::type >::magnitudeType
   >::type
   relErrFadImpl(const T1& v1, const T2& v2) { return Teuchos::relErr(v1,v2); }

   // relErrFadImpl version for Fad types
   template <typename T1, typename T2>
   typename std::enable_if<
     IsExpr<T1>::value || IsExpr<T2>::value,
     typename Teuchos::ScalarTraits<
       typename std::common_type< typename ScalarType<T1>::type, typename ScalarType<T2>::type >::type
     >::magnitudeType
   >::type
   relErrFadImpl(const T1& v1, const T2& v2)
    {
      typedef typename Teuchos::ScalarTraits<typename std::common_type< typename ScalarType<T1>::type, typename ScalarType<T2>::type >::type>::magnitudeType magnitudeType;
      magnitudeType maxRelErr = relErrFadImpl(v1.val(), v2.val());
      for (int i=0; i<v1.size(); ++i) {
        magnitudeType tmpRelErr = relErrFadImpl(v1.dx(i), v2.dx(i));
        if (tmpRelErr >= maxRelErr)
          maxRelErr = tmpRelErr;
      }
      return maxRelErr;
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
    typedef typename std::common_type<scalarType1,scalarType2>::type scalarType;
    typedef typename Teuchos::ScalarTraits<scalarType>::magnitudeType magnitudeType;
    static bool eval(
      const std::string &v1_name,
      const T1 &vv1,
      const std::string &v2_name,
      const T2 &vv2,
      const std::string &maxRelErr_error_name,
      const magnitudeType &maxRelErr_error,
      const std::string &maxRelErr_warning_name,
      const magnitudeType &maxRelErr_warning,
      const Ptr<std::ostream> &out
    )
    {
      using std::endl;
      typedef Teuchos::ScalarTraits<magnitudeType> SMT;
      typename Sacado::BaseExprType<T1>::type v1(vv1);
      typename Sacado::BaseExprType<T2>::type v2(vv2);
      const magnitudeType rel_err = Sacado::Fad::relErrFadImpl( v1, v2);
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

#endif