// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
