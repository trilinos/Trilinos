// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.hpp
    \brief Example of how to supply ROL with parameters from a JSON file.
           Requires that 
           <a href="https://github.com/open-source-parsers/jsoncpp">json-cpp</a> be installed.

           To build this example, add the following to the cmake call in your trilinos build script
           -D TPL_ENABLE_JSONCPP=ON \
           -D JSONCPP_INCLUDE_DIRS:PATH=/usr/include/jsoncpp \
           -D JSONCPP_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu \
           -D JSONCPP_LIBRARY_NAMES:STRING="jsoncpp" \

           These example paths above are default for Ubuntu 64 bit if jsoncpp is installed using 
           sudo apt-get install libjsoncpp-dev  

           Possible build failures:
           Unable to find json.h - JSONCPP_INCLUDE_DIRS variable is incorrect
           Undefined reference to json::value - JSONCPP_LIBRARY variables are incorrect
    
    \author Created by Greg von Winckel

*/


#include "json/json.h"                    // For JSON definitions
#include "Teuchos_ParameterList.hpp"      // For Parameter List
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include <fstream>
#include <string>


namespace ROL {

void addJSONBlockToPL(const Json::Value& block,Teuchos::ParameterList& parlist);
void addJSONPairToPL(const Json::Value& block, const std::string &key,Teuchos::ParameterList& parlist);


/** \brief Given a JSON block and a key, get the value and insert the key-value pair into a
           Teuchos::ParameterList. If the value is itself a block, recursively iterate.

    @param[in]     block     is a block from a JSON object
    @param[in]     key       is a string key 
    @param[in/out] parlist   is a Teuchos::ParameterList
*/
void addJSONPairToPL(const Json::Value& block,
                     const std::string &key,
                     Teuchos::ParameterList& parlist) {

    Json::Value val = block[key];
 
    if(val.isString()) {
        parlist.set(key,val.asString());
    } else if(val.isBool()) {
        parlist.set(key,val.asBool());
    } else if(val.isInt()) {
        parlist.set(key,val.asInt());
    } else if(val.isUInt()) {
        parlist.set(key,val.asUInt());
    } else if(val.isDouble()) {
        parlist.set(key,val.asDouble());  
    } else if(val.isObject()) { // This is a block. Iterate over its pairs
        addJSONBlockToPL(val,parlist);
    }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, ">>> ERROR (addJSONPairToPL, "
                                                             "json value has unsupported type.");
    }
}


/** \brief Iterate over a block and insert key-value pairs into the Teuchos::ParameterList

    @param[in]     block     is a block from a JSON object
    @param[in/out] parlist   is a Teuchos::ParameterList
*/
void addJSONBlockToPL(const Json::Value& block,
                      Teuchos::ParameterList& parlist) {
    if(block.size()>0) {
        for(Json::ValueIterator itr = block.begin(); itr != block.end(); ++itr) {
            addJSONPairToPL(block,itr.key().asString(),parlist);
        }        
    }
}



/** \brief  Read a JSON file and store all parameters in a Teuchos::ParameterList.
            Checks for a key called "Algorithm" which has a string value which can
            specify a Step Type (Linesearch or Trust-Region) and either a
            Descent Type or a Trust-Region Subproblem Solver Type.
    
    @param[in]     block     is a block from a JSON object
    @param[in/out] parlist   is a Teuchos::ParameterList
*/
void JSON_Parameters(const std::string& jsonFileName, 
		     Teuchos::ParameterList& parlist) {

Json::Value json;   
std::ifstream stream(jsonFileName, std::ifstream::binary);
stream >> json;

if(json.isMember("ROL")) {
    Json::Value rolBlock = json["ROL"];

    // Make a flat parameter list from the ROL JSON block
    addJSONBlockToPL(rolBlock,parlist);    

    // Check for an algorithm  
    if(rolBlock.isMember("Algorithm")) {
	std::string rolAlgorithm = rolBlock["Algorithm"].asString();

	if(rolAlgorithm.find("Trust-Region") != std::string::npos) {

	    parlist.set("Step Type","Trust-Region");

	    // Set subproblem solver
	    if(rolAlgorithm.find("Cauchy Point") != std::string::npos) {
		parlist.set("Trust-Region Subproblem Solver Type","Cauchy Point");
	    } 
	    else if(rolAlgorithm.find("Double Dogleg") != std::string::npos) {
		parlist.set("Trust-Region Subproblem Solver Type","Double Dogleg");
	    }
	    else if(rolAlgorithm.find("Dogleg") != std::string::npos) {
		parlist.set("Trust-Region Subproblem Solver Type","Dogleg");
	    }
	    else if(rolAlgorithm.find("Truncated CG") != std::string::npos) {
		parlist.set("Trust-Region Subproblem Solver Type","Truncated CG");
	    }

	}         
	else { // Use Linesearch
	    parlist.set("Step Type","Linesearch");

	    // Set descent type
	    if(rolAlgorithm.find("Steepest Descent") != std::string::npos) {
		parlist.set("Descent Type","Steepest Descent");
	    } 
	    else if(rolAlgorithm.find("Quasi-Newton") != std::string::npos) {
		parlist.set("Descent Type","Quasi-Newton Method");
	    }
	    else if(rolAlgorithm.find("Newton-Krylov") != std::string::npos) {
		parlist.set("Descent Type","Newton-Krylov");
	    }
	    else if(rolAlgorithm.find("Nonlinear CG") != std::string::npos) { 
		    parlist.set("Descent Type","Nonlinear CG"); 
		}

	    }
	    
	} 
	else { // No algorithm block found - use defaults
	    parlist.set("Step Type","Linesearch");
	    parlist.set("Descent Type","Nonlinear CG"); 
	}

    }
}


/** \brief  A minimalist step factory which specializes the Step Type depending on 
            whether a Trust-Region or Linesearch has been selected.     
    @param[in]     parlist   is a Teuchos::ParameterList
    @param[in/out] step      is a ref count pointer to a ROL::Step 
*/
template <class Real>
void stepFactory(Teuchos::ParameterList &parlist,Teuchos::RCP<ROL::Step<Real> > &step) {
     
    if(parlist.get("Step Type","Linesearch")=="Trust-Region") {
	step = Teuchos::rcp(new ROL::TrustRegionStep<Real>(parlist));
    }
    else {
	step = Teuchos::rcp(new ROL::LineSearchStep<Real>(parlist));
	 
    }
}

}
