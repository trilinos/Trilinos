// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_SCALARTRAITSIMP_HPP
#define SACADO_FAD_SCALARTRAITSIMP_HPP

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_SACADO_TEUCHOSCORE

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {

  namespace Fad {

    //! Implementation for Teuchos::ScalarTraits for all Fad types
    template <typename FadType>
    struct ScalarTraitsImp {
      typedef typename Sacado::ValueType<FadType>::type ValueT;

      typedef typename mpl::apply<FadType,typename Teuchos::ScalarTraits<ValueT>::magnitudeType>::type magnitudeType;
      typedef typename mpl::apply<FadType,typename Teuchos::ScalarTraits<ValueT>::halfPrecision>::type halfPrecision;
      typedef typename mpl::apply<FadType,typename Teuchos::ScalarTraits<ValueT>::doublePrecision>::type doublePrecision;

      static const bool isComplex = Teuchos::ScalarTraits<ValueT>::isComplex;
      static const bool isOrdinal = Teuchos::ScalarTraits<ValueT>::isOrdinal;
      static const bool isComparable =
        Teuchos::ScalarTraits<ValueT>::isComparable;
      static const bool hasMachineParameters =
        Teuchos::ScalarTraits<ValueT>::hasMachineParameters;
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType eps() {
        return Teuchos::ScalarTraits<ValueT>::eps();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType sfmin() {
        return Teuchos::ScalarTraits<ValueT>::sfmin();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType base()  {
        return Teuchos::ScalarTraits<ValueT>::base();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType prec()  {
        return Teuchos::ScalarTraits<ValueT>::prec();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType t()     {
        return Teuchos::ScalarTraits<ValueT>::t();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType rnd()   {
        return Teuchos::ScalarTraits<ValueT>::rnd();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType emin()  {
        return Teuchos::ScalarTraits<ValueT>::emin();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType rmin()  {
        return Teuchos::ScalarTraits<ValueT>::rmin();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType emax()  {
        return Teuchos::ScalarTraits<ValueT>::emax();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType rmax()  {
        return Teuchos::ScalarTraits<ValueT>::rmax();
      }
      static magnitudeType magnitude(const FadType& a) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
          a, "Error, the input value to magnitude(...) a = " << a <<
          " can not be NaN!" );
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(a) == false, std::runtime_error,
                           "Complex magnitude is not a differentiable "
                           "function of complex inputs.");
#endif
        //return std::fabs(a);
        magnitudeType b(a.size(),
                        Teuchos::ScalarTraits<ValueT>::magnitude(a.val()));
        if (Teuchos::ScalarTraits<ValueT>::real(a.val()) >= 0)
          for (int i=0; i<a.size(); i++)
            b.fastAccessDx(i) =
              Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessDx(i));
        else
          for (int i=0; i<a.size(); i++)
            b.fastAccessDx(i) =
              -Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessDx(i));
        return b;
      }
      static ValueT zero()  {
        return ValueT(0.0);
      }
      static ValueT one()   {
        return ValueT(1.0);
      }

      // Conjugate is only defined for real derivative components
      static FadType conjugate(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(x) == false, std::runtime_error,
                           "Complex conjugate is not a differentiable "
                           "function of complex inputs.");
#endif
        FadType y = x;
        y.val() = Teuchos::ScalarTraits<ValueT>::conjugate(x.val());
        return y;
      }

      // Real part is only defined for real derivative components
      static FadType real(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(x) == false, std::runtime_error,
                           "Real component is not a differentiable "
                           "function of complex inputs.");
#endif
        FadType y = x;
        y.val() = Teuchos::ScalarTraits<ValueT>::real(x.val());
        return y;
      }

      // Imaginary part is only defined for real derivative components
      static FadType imag(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(x) == false, std::runtime_error,
                           "Imaginary component is not a differentiable "
                           "function of complex inputs.");
#endif
        return FadType(Teuchos::ScalarTraits<ValueT>::imag(x.val()));
      }

      static ValueT nan() {
        return Teuchos::ScalarTraits<ValueT>::nan();
      }
      static bool isnaninf(const FadType& x) {
        if (Teuchos::ScalarTraits<ValueT>::isnaninf(x.val()))
          return true;
        for (int i=0; i<x.size(); i++)
          if (Teuchos::ScalarTraits<ValueT>::isnaninf(x.dx(i)))
            return true;
        return false;
      }
      static void seedrandom(unsigned int s) {
        Teuchos::ScalarTraits<ValueT>::seedrandom(s);
      }
      static ValueT random() {
        return Teuchos::ScalarTraits<ValueT>::random();
      }
      static std::string name() {
        return Sacado::StringName<FadType>::eval();
      }
      static FadType squareroot(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
          x, "Error, the input value to squareroot(...) a = " << x <<
          " can not be NaN!" );
#endif
        return std::sqrt(x);
      }
      static FadType pow(const FadType& x, const FadType& y) {
        return std::pow(x,y);
      }

      // Helper function to determine whether a complex value is real
      static bool is_complex_real(const ValueT& x) {
        return
          Teuchos::ScalarTraits<ValueT>::magnitude(x-Teuchos::ScalarTraits<ValueT>::real(x)) == 0;
      }

      // Helper function to determine whether a Fad type is real
      static bool is_fad_real(const FadType& x) {
        if (x.size() == 0)
          return true;
        if (Teuchos::ScalarTraits<ValueT>::isComplex) {
          if (!is_complex_real(x.val()))
            return false;
          for (int i=0; i<x.size(); i++)
            if (!is_complex_real(x.fastAccessDx(i)))
              return false;
        }
        return true;
      }

    }; // class ScalarTraitsImp

  } // namespace Fad

} // namespace Sacado

#endif // HAVE_SACADO_TEUCHOSCORE

#endif // SACADO_FAD_SCALARTRAITSIMP_HPP
