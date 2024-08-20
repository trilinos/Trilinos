// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_DataLayout.hpp"

//**********************************************************************
std::ostream& PHX::operator<<(std::ostream& os, const PHX::DataLayout& t)
{
  t.print(os);
  return os;
}

//**********************************************************************
