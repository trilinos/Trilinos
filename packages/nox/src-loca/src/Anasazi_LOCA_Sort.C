// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Anasazi_LOCA_Sort.H"
#include "LOCA_EigenvalueSort_Strategies.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

Anasazi::LOCASort::LOCASort(
 const Teuchos::RCP<LOCA::GlobalData>& global_data,
 const Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>& strategy_)
  : globalData(global_data),
    strategy(strategy_)
{
}

Anasazi::LOCASort::~LOCASort()
{
}

void
Anasazi::LOCASort::sort(std::vector<double>& evals,
                              Teuchos::RCP<std::vector<int> > perm,
                        int n) const
{
  if (n == -1) {
    n = evals.size();
  }
  NOX::Abstract::Group::ReturnType res = strategy->sort(n, &evals[0], perm.get());
  globalData->locaErrorCheck->checkReturnType(res, "Anasazi::LOCASort::sort()");
}

void
Anasazi::LOCASort::sort(std::vector<double>& r_evals,
                              std::vector<double>& i_evals,
                              Teuchos::RCP<std::vector<int> > perm,
                        int n) const
{
  if (n == -1) {
    n = r_evals.size();
  }
  NOX::Abstract::Group::ReturnType res =
    strategy->sort(n, &r_evals[0], &i_evals[0], perm.get());
  globalData->locaErrorCheck->checkReturnType(res,
                          "Anasazi::LOCASort::sort()");
}
