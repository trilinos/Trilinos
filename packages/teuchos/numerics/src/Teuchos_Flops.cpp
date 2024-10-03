// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Flops.hpp"

namespace Teuchos
{

Flops::Flops(void) : flops_(0.0)
{
}

// 2007/11/26: rabartl: Below, is it correct that flops_in does not have its
// flops copied into the data member flops_?
Flops::Flops(const Flops& /* flops_in */) : flops_(0.0)
{
}

Flops::~Flops(void)
{
  flops_ = 0.0;
}

}  // namespace Teuchos
