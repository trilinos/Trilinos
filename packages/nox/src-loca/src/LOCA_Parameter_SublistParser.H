// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PARAMETER_SUBLISTPARSER_H
#define LOCA_PARAMETER_SUBLISTPARSER_H

#include "Teuchos_RCP.hpp"

#include "NOX_Common.H"
#include <string>
#include <map>

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Parameter {

    //! Class to parse a parameter list for sublists
    /*!
     * This class parses a supplied parameter list and locates various
     * sublists.  This saves the code from having to traverse the parameter
     * list to find sublists itself, and puts in one location the
     * hard-coded structure of the parameter list.
     */
    class SublistParser {

      //! Map type to hold sublists
      typedef std::map<std::string,Teuchos::RCP<Teuchos::ParameterList> > SublistMap;

      //! Iterator type to traverse SublistMap
      typedef SublistMap::iterator SublistMapIterator;

    public:

      //! Constructor
      SublistParser(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      ~SublistParser();

      //! Parse parameter list to find sublists
      void parseSublists(
         const Teuchos::RCP<Teuchos::ParameterList>& topLevelParams);

      //! Return sublist of name \c name
      Teuchos::RCP<Teuchos::ParameterList>
      getSublist(const std::string& name);

    private:

      //! Private to prohibit copying
      SublistParser(const SublistParser&);

      //! Private to prohibit copying
      SublistParser& operator = (const SublistParser&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Map to hold sublists
      SublistMap sublistMap;

    }; // Class SublistParser

  } // Namespace Parameter

} // Namespace LOCA

#endif
