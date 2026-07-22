// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Describable.hpp"
#include "Teuchos_TypeNameTraits.hpp"


namespace Teuchos {


const EVerbosityLevel Describable::verbLevel_default = VERB_DEFAULT;


std::string Describable::description () const
{
  const std::string objectLabel = this->getObjectLabel ();
  std::ostringstream oss;
  if (objectLabel.length ()) {
    oss << "\"" << objectLabel << "\": ";
  }
  oss << typeName (*this);
  return oss.str ();
}

void
Describable::describe (FancyOStream& out_arg,
                       const EVerbosityLevel /* verbLevel */) const
{
  RCP<FancyOStream> out = rcpFromRef (out_arg);
  OSTab tab (out);
  *out << this->description () << std::endl;
}

void
Describable::describe (std::ostream& out,
                       const EVerbosityLevel verbLevel) const
{
  RCP<FancyOStream> fancyOut = getFancyOStream (rcpFromRef (out));
  this->describe (*fancyOut, verbLevel);
}

Describable::~Describable () {}

} // namespace Teuchos
