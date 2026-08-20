// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORY_BELOS_HPP
#define BELOSSOLVERFACTORY_BELOS_HPP

#include "Belos_Details_registerSolverFactory.hpp"

// Note that this file is currently included by BelosSolverFactory.hpp
// to maintain backwards compatibility. We don't include it here because
// gcc won't resolve the circular includes for namespacing
// #include "BelosSolverFactory.hpp"

#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"

#ifdef HAVE_TEUCHOSCORE_KOKKOS
#include "Kokkos_Core.hpp"
#endif

namespace Belos {

/** \example epetra/example/SolverFactory/SolverFactoryEpetraGaleriEx.cpp 
    This is an example of how to use the Belos::SolverFactory with Epetra.
*/
/** \example tpetra/example/SolverFactory/SolverFactoryTpetraGaleriEx.cpp 
    This is an example of how to use the Belos::SolverFactory with Tpetra.
*/

  class BelosSolverFactory : public Impl::SolverFactoryParent<double,MultiVec<double, DefaultDenseMatrix<int,double>>,Operator<double, DefaultDenseMatrix<int,double>>,
                                                          DefaultDenseMatrix<int,double>>
{
  public:
    BelosSolverFactory() {
      Details::registerSolverFactory();
    };
};

  class BelosFloatSolverFactory : public Impl::SolverFactoryParent<float,MultiVec<float, DefaultDenseMatrix<int,float>>,Operator<float, DefaultDenseMatrix<int,float>>,
                                                          DefaultDenseMatrix<int,float>>
{
  public:
    BelosFloatSolverFactory() {
      Details::registerSolverFactory();
    };
};

#if defined(HAVE_TEUCHOS_HALF) && \
    defined(HAVE_TEUCHOSCORE_KOKKOS) && defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
class BelosHalfSolverFactory : public Impl::SolverFactoryParent<Kokkos::Experimental::half_t,MultiVec<Kokkos::Experimental::half_t>,Operator<Kokkos::Experimental::half_t>>
{
  public:
    BelosHalfSolverFactory() {
      Details::registerSolverFactory();
    };
};
#endif

namespace Impl {

template<>
class SolverFactorySelector<double,MultiVec<double, DefaultDenseMatrix<int,double>>,Operator<double, DefaultDenseMatrix<int,double>>,DefaultDenseMatrix<int,double>> {
  public:
    typedef BelosSolverFactory type;
};

template<>
class SolverFactorySelector<float,MultiVec<float, DefaultDenseMatrix<int,float>>,Operator<float, DefaultDenseMatrix<int,float>>,DefaultDenseMatrix<int,float>> {
  public:
    typedef BelosFloatSolverFactory type;
};

#if defined(HAVE_TEUCHOS_HALF) && \
    defined(HAVE_TEUCHOSCORE_KOKKOS) && defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
template<>
class SolverFactorySelector<Kokkos::Experimental::half_t,MultiVec<Kokkos::Experimental::half_t>,Operator<Kokkos::Experimental::half_t>> {
  public:
    typedef BelosHalfSolverFactory type;
};
#endif

#ifdef HAVE_TEUCHOS_COMPLEX
class BelosComplexSolverFactory : public Impl::SolverFactoryParent<std::complex<double>,MultiVec<std::complex<double>>,
  Operator<std::complex<double>, DefaultDenseMatrix<int,std::complex<double>>>, DefaultDenseMatrix<int,std::complex<double>>>
{
  public:
    BelosComplexSolverFactory() {
      Details::registerSolverFactory();
    };
};

template<>
class SolverFactorySelector<std::complex<double>,MultiVec<std::complex<double>>,Operator<std::complex<double>, DefaultDenseMatrix<int,std::complex<double>>>,
                                    DefaultDenseMatrix<int,std::complex<double>>>
{
  public:
    typedef BelosComplexSolverFactory type;
};

class BelosFloatComplexSolverFactory : public Impl::SolverFactoryParent<std::complex<float>,MultiVec<std::complex<float>>,
                                                                        Operator<std::complex<float>, DefaultDenseMatrix<int,std::complex<float>>>, DefaultDenseMatrix<int,std::complex<float>>>
{
  public:
    BelosFloatComplexSolverFactory() {
      Details::registerSolverFactory();
    };
};

template<>
class SolverFactorySelector<std::complex<float>,MultiVec<std::complex<float>>,Operator<std::complex<float>, DefaultDenseMatrix<int,std::complex<float>>>,
                                    DefaultDenseMatrix<int,std::complex<float>>>
{
  public:
    typedef BelosFloatComplexSolverFactory type;
};
#endif

} // namespace Impl
} // namespace Belos

#endif // BELOSSOLVERFACTORY_BELOS_HPP
