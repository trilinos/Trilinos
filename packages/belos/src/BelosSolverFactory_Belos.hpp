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

namespace Belos {

/** \example epetra/example/SolverFactory/SolverFactoryEpetraGaleriEx.cpp 
    This is an example of how to use the Belos::SolverFactory with Epetra.
*/
/** \example tpetra/example/SolverFactory/SolverFactoryTpetraGaleriEx.cpp 
    This is an example of how to use the Belos::SolverFactory with Tpetra.
*/

class BelosSolverFactory : public Impl::SolverFactoryParent<double,MultiVec<double>,Operator<double>>
{
  public:
    BelosSolverFactory() {
      Details::registerSolverFactory();
    };
};

class BelosFloatSolverFactory : public Impl::SolverFactoryParent<float,MultiVec<float>,Operator<float>>
{
  public:
    BelosFloatSolverFactory() {
      Details::registerSolverFactory();
    };
};

namespace Impl {

template<>
class SolverFactorySelector<double,MultiVec<double>,Operator<double>> {
  public:
    typedef BelosSolverFactory type;
};

template<>
class SolverFactorySelector<float,MultiVec<float>,Operator<float>> {
  public:
    typedef BelosFloatSolverFactory type;
};

#ifdef HAVE_TEUCHOS_COMPLEX
class BelosComplexSolverFactory : public Impl::SolverFactoryParent<std::complex<double>,MultiVec<std::complex<double>>,Operator<std::complex<double>>>
{
  public:
    BelosComplexSolverFactory() {
      Details::registerSolverFactory();
    };
};

template<>
class SolverFactorySelector<std::complex<double>,MultiVec<std::complex<double>>,Operator<std::complex<double>>> {
  public:
    typedef BelosComplexSolverFactory type;
};

class BelosFloatComplexSolverFactory : public Impl::SolverFactoryParent<std::complex<float>,MultiVec<std::complex<float>>,Operator<std::complex<float>>>
{
  public:
    BelosFloatComplexSolverFactory() {
      Details::registerSolverFactory();
    };
};

template<>
class SolverFactorySelector<std::complex<float>,MultiVec<std::complex<float>>,Operator<std::complex<float>>> {
  public:
    typedef BelosFloatComplexSolverFactory type;
};
#endif

} // namespace Impl
} // namespace Belos

#endif // BELOSSOLVERFACTORY_BELOS_HPP
