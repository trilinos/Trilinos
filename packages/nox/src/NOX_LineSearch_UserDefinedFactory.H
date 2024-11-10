// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_USERDEFINED_FACTORY_H
#define NOX_LINESEARCH_USERDEFINED_FACTORY_H

#include "Teuchos_RCP.hpp"
#include "NOX_Common.H"

namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class GlobalData;
}

namespace NOX {
namespace LineSearch {

class Generic;

/*!  \brief Pure virtual interface for users to supply their own line search objects.

*/
class UserDefinedFactory {

public:

  //! Constructor
  UserDefinedFactory() {}

  //! Destructor
  virtual ~UserDefinedFactory() {}

  /*! \brief Builds a user defined line search object.

    @param gd A global data pointer that contains the top level
    parameter list.  Without storing this inside the line
    searchobject, there is no guarantee that the second parameter \c
    params will still exist.  It can be deleted by the top level RCP.

    @param params Sublist with line search construction parameters.

  */
  virtual Teuchos::RCP<NOX::LineSearch::Generic>
  buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params) const = 0;

};
} // namespace LineSearch
} // namespace NOX
#endif
