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

// Constructor

#include "Teuchos_CompObject.hpp"

namespace Teuchos
{

CompObject::CompObject() : flopCounter_(0)
{
}

// Copy Constructor

CompObject::CompObject(const CompObject& source) : flopCounter_(source.flopCounter_)
{
}

// Destructor

CompObject::~CompObject()
{
  flopCounter_ = 0;
}

} // namespace Teuchos
