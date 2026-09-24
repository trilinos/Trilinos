// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosSolverManager.hpp"

#include "BelosBiCGStabIter.hpp"
#include "BelosBlockCGIter.hpp"
#include "BelosBlockFGmresIter.hpp"
#include "BelosCGIter.hpp"
#include "BelosCGSingleRedIter.hpp"
#include "BelosCurrentSolutionProvider.hpp"
#include "BelosFGCRODRIter.hpp"
#include "BelosFixedPointIter.hpp"
#include "BelosLSQRIter.hpp"
#include "BelosMinresIter.hpp"
#include "BelosPCPGIter.hpp"
#include "BelosPseudoBlockCGIter.hpp"
#include "BelosPseudoBlockGmresIter.hpp"
#include "BelosPseudoBlockStochasticCGIter.hpp"
#include "BelosPseudoBlockTFQMRIter.hpp"
#include "BelosRCGIter.hpp"
#include "BelosTFQMRIter.hpp"
#include "BelosTpetraAdapter.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"

#include <cstdlib>
#include <iostream>
#include <type_traits>

namespace {

template<class Iter, class Provider>
void assertProvider()
{
  static_assert(std::is_base_of<Provider, Iter>::value,
                "Instrumented iteration must implement CurrentSolutionProvider");
}

template<class SC>
void checkProviderTypes()
{
  using MV = Tpetra::MultiVector<SC>;
  using OP = Tpetra::Operator<SC>;
  using DM = Teuchos::SerialDenseMatrix<int, SC>;
  using Provider = Belos::CurrentSolutionProvider<SC, MV, OP, DM>;

  assertProvider<Belos::BiCGStabIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::BlockCGIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::BlockFGmresIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::CGIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::CGSingleRedIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::FGCRODRIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::FixedPointIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::LSQRIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::MinresIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::PCPGIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::PseudoBlockCGIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::PseudoBlockGmresIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::PseudoBlockStochasticCGIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::PseudoBlockTFQMRIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::RCGIter<SC, MV, OP, DM>, Provider>();
  assertProvider<Belos::TFQMRIter<SC, MV, OP, DM>, Provider>();
}

} // namespace

int main(int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  checkProviderTypes<double>();
  if (Tpetra::getDefaultComm()->getRank() == 0)
    std::cout << "End Result: TEST PASSED" << std::endl;
  return EXIT_SUCCESS;
}
