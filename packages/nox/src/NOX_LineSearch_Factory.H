// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_LINESEARCH_FACTORY_H
#define NOX_LINESEARCH_FACTORY_H

#include "Teuchos_RCP.hpp"
#include "NOX_Common.H"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class Utils;
}

namespace NOX {

class GlobalData;

namespace LineSearch {

class Generic;

/*!

  \brief Factory to build line search objects derived from NOX::LineSearch::Generic.

  <b>Parameters</b>

   <ul>
   <li> "Method" <std::string> Name of the line search. Valid choices are:
     <ul>
     <li> "Full Step" (NOX::LineSearch::FullStep)
     <li> "Backtrack" (NOX::LineSearch::Backtrack)
     <li> "Polynomial" (NOX::LineSearch::Polynomial)
     <li> "More'-Thuente" (NOX::LineSearch::MoreThuente)
     <li> "User Defined" - see below
     </ul>

   <li> "User Defined Constructor" - see below

   <li> "Full Step" <sublist> Parameters to build a NOX::LineSearch::FullStep sublist.
   <li> "Backtrack" <sublist> Parameters to build a NOX::LineSearch::Backtrack sublist.
   <li> "Polynomial" <sublist> Parameters to build a NOX::LineSearch::Polynomial sublist.
   <li> "More'-Thuente" <sublist> Parameters to build a NOX::LineSearch::MoreThuente sublist.
   <li> "User Defined Line Search Factory" < RCP<NOX::LineSearch::UserDefinedFactory> > RCP to a NOX::LineSearch::UserDefinedFactory derived object.  This factory object is used to build user defined line search objects.
   </ul>

<b>Using a User-Defined Line Search</b>

The user has the option of passing in a user-defined line search.  First, they must implement their own line search, deriving from the base class interface NOX::LineSearch::Generic:

\code
   class MyLineSearch : public NOX::LineSearch::Generic {

     // Ctor that takes the standard line search arguments.
     MyLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
                  Teuchos::ParameterList& params);

     .
     .
     .

   };
\endcode

Next they must write a factory to build their object, deriving from the NOX::LineSearch::UserDefinedFactory base class interface:

\code
class MyFactory {
  MyDirFactory();
  ~MyDirFactory();
   Teuchos::RCP<NOX::LineSearch::Generic> buildLineSearch(const Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)
   .
   .
   .
};
\endcode

Then under the "Line Search" parameter sublist, they need to set the
method to "User Defined" and register the factory:

\code
using namespace Teuchos; // for RCP and ParameterList
ParameterList& ls_params = p.sublist("Line Search");
RCP<NOX::LineSearch::UserDefinedFactory> ls_facotry = rcp(new MyLineSearchFactory);
ls_params.set("Method", "User Defined");
ls_params.set("User Defined Line Search Factory", ls_factory);
\endcode

It is critical that the user defined factory be set in the parameter
list as a base class type object: NOX::LineSearch::UserDefinedFactory.

 */
class Factory {

public:

  //! Constructor
  Factory();

  //! Destructor
  ~Factory();

  /*! \brief Factory to build a line search object.

      @param gd A global data pointer that contains the top level parameter list.  Without storing this inside the line searchobject, there is no guarantee that the second parameter \c params will still exist.  It can be deleted by the top level RCP.
      @param params Sublist with line search construction parameters.

  */
  Teuchos::RCP<NOX::LineSearch::Generic>
  buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);

};

  /*! Nonmember function to build a line search object.

  \relates NOX::LineSearch::Factory

  */
  Teuchos::RCP<NOX::LineSearch::Generic>
  buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);


} // namespace LineSearch
} // namespace NOX
#endif
