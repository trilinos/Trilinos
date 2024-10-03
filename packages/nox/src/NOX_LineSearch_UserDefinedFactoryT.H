// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_USERDEFINED_FACTORY_TEMPLATE_H
#define NOX_LINESEARCH_USERDEFINED_FACTORY_TEMPLATE_H

#include "Teuchos_RCP.hpp"
#include "NOX_LineSearch_UserDefinedFactory.H"
#include "NOX_Common.H"

namespace NOX {
namespace LineSearch {

/*!  \brief Concrete instantiation of a NOX::LineSearch::UserDefinedFactory object that uses the base objects only for constuction.

   If the user writes their own line search and that object has the same
   constructor arguments as the nox line searches (the \c gd and \c
   params as in the buildDirection method), then users can use this
   object instead of having to write their own factory.

   For example, if a user writes their own line search object:
   \code
   class MyLineSearch : public NOX::LineSearch::Generic {

     // Ctor that takes the standard line search arguments.
     MyLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
                  Teuchos::ParameterList& params);

     .
     .
     .

   }
   \endcode

   They can build that object using this factory and do not
   have to write their own factory

   \code
   ParameterList dl& = nox_params.sublist("Line Search");
   RCP<NOX::LineSearch::UserDefinedFactory> udlsf =
     rcp(new NOX::LineSearch::UserDefinedFactoryT<MyLineSearch>);
   dl.set("User Defined LineSearch Factory", udlsf);
   \endcode

   It is critical that the user defined factory be set in the
   parameter list as a base class type object:
   NOX::LineSearch::UserDefinedFactory.

 */
template <typename T>
class UserDefinedFactoryT : public NOX::LineSearch::UserDefinedFactory {

public:

  //! Constructor
  UserDefinedFactoryT() {}

  //! Destructor
  ~UserDefinedFactoryT() {}

  Teuchos::RCP<NOX::LineSearch::Generic>
  buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params) const
  {
    using namespace Teuchos;
    RCP<NOX::LineSearch::Generic> ls = rcp(new T(gd, params));
    return ls;
  }

};
} // namespace LineSearch
} // namespace NOX
#endif
