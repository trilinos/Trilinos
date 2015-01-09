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

#ifndef NOX_GLOBALDATA_H
#define NOX_GLOBALDATA_H

#include "Teuchos_RCP.hpp"

// Global data objects forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class Utils;
  namespace MeritFunction {
    class Generic;
  }
}

namespace NOX {

  //! Container class to hold "global" %NOX objects
  /*!%GlobalData is a container class that holds ref-count pointers to
   * "global" objects, i.e., objects that nearly every NOX object will
   * need access to.  By putting them all in one container class, the
   * container class can be stored in each NOX object, and if a new
   * global object is needed, it can be added here without modifying
   * the rest of the code.  This is an alternative to true global or
   * static objects which are note safe in many contexts (threading).
   * In particular, this approach allows multiple %NOX "invocations"
   * to be in memory at the same time.
   */
  class GlobalData {

  public:

    //! Consturctor using the top level NOX parameter list.
    GlobalData(const Teuchos::RCP<Teuchos::ParameterList>& noxParams);

    //! Constructor taking a ref-count pointer to each global object.
    GlobalData(const Teuchos::RCP<NOX::Utils>& utils,
           const Teuchos::RCP<NOX::MeritFunction::Generic>& mf);

    //! Destructor.
    virtual ~GlobalData();

    //! Returns the print utils object.
    Teuchos::RCP<NOX::Utils> getUtils() const;

    //! Returns the merit function object.
    Teuchos::RCP<NOX::MeritFunction::Generic> getMeritFunction() const;

    /*! \brief Returns the top-level nox parameter list input by the user.

       This list is kept in global data so that any sublists of the
       main parameters list that objects may keep a refernece to is
       still valid.  The line searches and directions usually store
       data in an output sublist for the users to query.  These
       sublists are NOT wrapped in reference counted smart pointers,
       so if the base list is deleted, the references to the sublist
       will no longer be valid.  To remedy this, any object that
       stores a reference to a sublist should also store the global
       data object.
    */
    Teuchos::RCP<Teuchos::ParameterList> getNoxParameterList() const;

  private:

    //! Copy constructor is private to preclude copying.
    GlobalData(const GlobalData&);

    //! Assignment operator is private to preclude copying.
    GlobalData& operator = (const GlobalData&);

  private:

    //! Ref-count pointer to NOX::Utils object.
    /* All output is controlled by the NOX::Utils object, so anything
       that has the potential to output (especially for warnings and
       errors) needs this object.
    */
    Teuchos::RCP<NOX::Utils> utilsPtr;

    //! Ref-count pointer to a NOX::MeritFunction::Generic object.
    /* This is global because both line searches and directions
       *might* need a merit function for their tasks.  We don't want
       to require a merit function for all line search and direction
       ctors since only a few line searches and the steepest descent
       direction needs a merit function.  But we need all objects to
       have a consistent merit function.
    */
    Teuchos::RCP<NOX::MeritFunction::Generic> meritFunctionPtr;

    /*! \brief The top-level nox parameter list input by the user.

       This list is kept in global data so that any sublists of the
       main parameters list that objects may keep a refernece to is
       still valid.  The line searches and directions usually store
       data in an output sublist for the users to query.  These
       sublists are NOT wrapped in reference counted smart pointers,
       so if the base list is deleted, the references to the sublist
       will no longer be valid.  To remedy this, any object that
       stores a reference to a sublist should also store the global
       data object.
    */
    Teuchos::RCP<Teuchos::ParameterList> paramListPtr;

  }; // Class GlobalData

} // Namespace NOX

#endif
