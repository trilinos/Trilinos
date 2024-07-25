// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_DIRECTION_USERDEFINED_FACTORY_TEMPLATE_H
#define NOX_DIRECTION_USERDEFINED_FACTORY_TEMPLATE_H

#include "Teuchos_RCP.hpp"
#include "NOX_Direction_UserDefinedFactory.H"
#include "NOX_Common.H"

namespace NOX {
namespace Direction {

/*!  \brief Concrete instantiation of a NOX::Direction::UserDefinedFactory object that uses the base objects only for constuction.

   If the user writes their own direction and that object has the same
   constructor arguments as the nox directions (the \c gd and \c
   params as in the buildDirection method), then users can use this
   object instead of having to write their own factory.

   For example, if a user writes their own direction object:
   \code
   class MyDirection : public NOX::Direction::Generic {

     // Ctor that takes the standard direction arguments.
     MyDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
                 Teuchos::ParameterList& params);

     .
     .
     .

   };
   \endcode

   They can build that object using this factory and do not
   have to write their own factory

   \code
   ParameterList dl& = nox_params.sublist("Direction");
   RCP<NOX::Direction::UserDefinedFactory> uddf =
     rcp(new NOX::Direction::UserDefinedFactoryT<MyDirection>);
   dl.set("User Defined Direction Factory", uddf);
   \endcode

   It is critical that the user defined factory be set in the
   parameter list as a base class type object:
   NOX::Direction::UserDefinedFactory.

 */
template <typename T>
class UserDefinedFactoryT : public NOX::Direction::UserDefinedFactory {

public:

  //! Constructor
  UserDefinedFactoryT() {}

  //! Destructor
  ~UserDefinedFactoryT() {}

  Teuchos::RCP<NOX::Direction::Generic>
  buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params) const
  {
    using namespace Teuchos;
    RCP<NOX::Direction::Generic> ls = rcp(new T(gd, params));
    return ls;
  }

};
} // namespace Direction
} // namespace NOX
#endif
