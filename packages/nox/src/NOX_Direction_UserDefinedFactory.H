// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_DIRECTION_USERDEFINED_FACTORY_H
#define NOX_DIRECTION_USERDEFINED_FACTORY_H

#include "Teuchos_RCP.hpp"
#include "NOX_Common.H"

namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class GlobalData;
}

namespace NOX {
namespace Direction {

class Generic;

/*!  \brief Pure virtual interface for users to supply their own direction objects.

*/
class UserDefinedFactory {

public:

  //! Constructor
  UserDefinedFactory() {}

  //! Destructor
  virtual ~UserDefinedFactory() {}

  /*! \brief Builds a user defined direction object.

    @param gd A global data pointer that contains the top level
    parameter list.  Without storing this inside the direction object,
    there is no guarantee that the second parameter \c params will
    still exist.  It can be deleted by the top level RCP.

    @param params Sublist with direction construction parameters.

  */
  virtual Teuchos::RCP<NOX::Direction::Generic>
  buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params) const = 0;

};
} // namespace Direction
} // namespace NOX
#endif
