// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <BelosConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <complex>

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "MyOperator.hpp"

//
// mfh 20 Jan 2014: This test ensures that Belos::SolverFactory can
// compile whether its ScalarType (first) template parameter is real
// or complex.
//
// This test requires that Trilinos was compiled with complex
// arithmetic support enabled.
//

TEUCHOS_UNIT_TEST( Factory, Real )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef double ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::SolverFactory<ST, MV, OP> factory_type;

  factory_type factory;
}

TEUCHOS_UNIT_TEST( Factory, Complex )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef std::complex<double> ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::SolverFactory<ST, MV, OP> factory_type;

  factory_type factory;
}
