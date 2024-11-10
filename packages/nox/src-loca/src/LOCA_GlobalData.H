// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_GLOBALDATA_H
#define LOCA_GLOBALDATA_H

#include "Teuchos_RCP.hpp"

// Global data objects forward declarations
namespace Teuchos {
    class ParameterList;
}
namespace NOX {
  class Utils;
}
namespace LOCA {
  class ErrorCheck;
  class Factory;
  namespace Abstract {
    class Factory;
  }
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  //! Container class to hold "global" %LOCA objects
  /*!
   * %GlobalData is a container class that holds ref-count pointers to
   * "global" objects, i.e., objects that nearly every LOCA object will need
   * access to.  By putting them all in one container class, the container
   * class can be stored in each LOCA object, and if a new global object
   * is needed, it can be added here without modifying the rest of the
   * code.  This is an alternative to true global or static objects which are
   * note safe in many contexts.  In particular, this approach allows multiple
   * %LOCA "invocations" to be in memory at the same time. Note that all data
   * members are declared public.
   */
  class GlobalData {

  public:

    //! Constructor taking a ref-count pointer to each global object
    GlobalData(const Teuchos::RCP<NOX::Utils>& loca_utils,
           const Teuchos::RCP<LOCA::ErrorCheck>& loca_error_check,
           const Teuchos::RCP<LOCA::Factory>& loca_factory);

    //! Destructor
    virtual ~GlobalData();

    //! Ref-count pointer to NOX::Utils object
    Teuchos::RCP<NOX::Utils> locaUtils;

    //! Ref-count pointer to LOCA::ErrorCheck object
    Teuchos::RCP<LOCA::ErrorCheck> locaErrorCheck;

    //! Ref-count pointer to LOCA::Factory object
    Teuchos::RCP<LOCA::Factory> locaFactory;

    //! Ref-count pointer to parsed parameter list
    Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

  private:

    //! Copy constructor is private to preclude copying
    GlobalData(const GlobalData&);

    //! Assignment operator is private to preclude copying
    GlobalData& operator = (const GlobalData&);

  }; // Class GlobalData

  //! Creates and initializes a LOCA::GlobalData object
  Teuchos::RCP<LOCA::GlobalData>
  createGlobalData(const Teuchos::RCP<Teuchos::ParameterList>& paramList,
           const Teuchos::RCP<LOCA::Abstract::Factory>& userFactory = Teuchos::null);

  //! De-initializes a LOCA::GlobalData object for destruction
  /*!
   * Sets the data members to Teuchos::null to remove circular references
   */
  void
  destroyGlobalData(const Teuchos::RCP<LOCA::GlobalData>& globalData);

} // Namespace LOCA

#endif
