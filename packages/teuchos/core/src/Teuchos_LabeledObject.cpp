// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_LabeledObject.hpp"


namespace Teuchos {


LabeledObject::LabeledObject()
  : objectLabel_("")
{}


LabeledObject::~LabeledObject()
{}


void LabeledObject::setObjectLabel( const std::string &objectLabel )
{
  objectLabel_ = objectLabel;
}


std::string LabeledObject::getObjectLabel() const
{
  return objectLabel_;
}


} // namespace Teuchos
