// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_IOSSDatabaseTypeManager_hpp__
#define __Panzer_IOSSDatabaseTypeManager_hpp__

#include <string>
#include <set>

#include <Teuchos_TestForException.hpp>
#include <Ioss_IOFactory.h>

namespace panzer_ioss {

class IossDatabaseTypeManager {

public:

	/** Get valid IOSS database types.
	  *
      * \returns a vector of valid database types
      */
	static std::vector<std::string> getValidTypes() {
	  std::vector<std::string> validTypes;
      Ioss::IOFactory::describe(&validTypes);
      return validTypes;
	}

	/** Determine whether a sting names a valid IOSS database type.
	  *
	  * \param[in] type The IOSS database type name
	  *
	  * \returns true if the string names a valid IOSS database type
	  */
	static bool isValidType(const std::string inputType) {
		bool validType = false;
		std::vector<std::string> validTypes = getValidTypes();
		for (std::string type : validTypes) {
		  if (type == inputType) {
		    validType = true;
		    break;
		  }
		}
		return validType;
	}

	/** List valid IOSS database types.
	  *
      * \returns a string listing valid database types
      */
	static std::string validTypeList() {
	  std::string validTypeList = "";
	  std::vector<std::string> validTypes = getValidTypes();
	  for (std::string type : validTypes) {
	    validTypeList += type;
	    validTypeList += ", ";
	  }
	  return validTypeList;
	}

	/** Determine whether the IOSS database type supports multiple open databases.
		  *
		  * \param[in] type The IOSS database type name
		  *
		  * \returns true if this IOSS database type supports multiple open databases
		  */
	static bool supportsMultipleOpenDatabases(const std::string inputType) {
		TEUCHOS_TEST_FOR_EXCEPTION(!isValidType(inputType), std::logic_error,
				 "Error, " << inputType
				 << " is an invalid IOSS database type." << std::endl
				 << "Valid types are " << validTypeList());
		std::set<std::string> multipleNotSupported;
		multipleNotSupported.insert("pamgen");
		return multipleNotSupported.find(inputType) == multipleNotSupported.end();
	}


};

} // end namespace panzer_ioss

#endif
