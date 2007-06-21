// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
