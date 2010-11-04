//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef __Belos_OrthoManagerFactory_hpp
#define __Belos_OrthoManagerFactory_hpp

#include <BelosConfigDefs.hpp>
#include <BelosTsqrOrthoManager.hpp>
#include <BelosICGSOrthoManager.hpp>
#include <BelosIMGSOrthoManager.hpp>
#include <BelosDGKSOrthoManager.hpp>

#include <Teuchos_StandardCatchMacros.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept> // #include <string>
#include <vector>

//////////////////////////////////////////////////////////////////////

namespace Belos {

  /// \class OrthoManagerFactory
  /// \brief Enumeration of all valid Belos (Mat)OrthoManager classes
  ///
  template< class Scalar, class MV, class OP >
  class OrthoManagerFactory {
  private:
    /// List of valid OrthoManager names
    ///
    std::vector< std::string > theList_;

  public:
    /// Number of valid command-line parameter values for the OrthoManager
    /// subclass to test.  (Must be at least one.)
    static int numOrthoManagers () { return 4; }

    /// Constructor
    ///
    OrthoManagerFactory () : theList_ (numOrthoManagers())
    {
      int index = 0;
      theList_[index++] = "TSQR";
      theList_[index++] = "ICGS";
      theList_[index++] = "IMGS";
      theList_[index++] = "DGKS";
    }

    /// Valid names of (Mat)OrthoManagers.  Useful as a list of valid
    /// command-line parameter values for choosing a (Mat)OrthoManager
    /// subclass to test.  
    ///
    /// \note To implementers: Anasazi and Belos currently implement
    ///   different sets of (Mat)OrthoManagers.  This method returns a
    ///   valid list of Belos (Mat)OrthoManagers.
    const std::vector< std::string >& 
    validNames () const { return theList_; }

    /// Return whether name names a valid OrthoManager
    ///
    bool
    isValidName (const std::string& name) const 
    {
      return (std::find (theList_.begin(), theList_.end(), name) != theList_.end());
    }

    /// Print the list of valid OrthoManager names to the given ostream
    ///
    void
    printValidNames (std::ostream& out) const
    {
      const int numValid = numOrthoManagers();
      TEST_FOR_EXCEPTION( numValid <= 0,
			  std::logic_error,
			  "Invalid number " 
			  << numValid 
			  << " of valid OrthoManager names" );
      if (numValid > 1)
	{
	  for (int k = 0; k < numValid - 1; ++k)
	    out << "\"" << theList_[k] << "\", ";
	  out << "or ";
	}
      out << "\"" << theList_[numValid-1] << "\"";
    }

    /// Return a list (as a string) of valid command-line parameter values
    /// for the OrthoManager subclass to test.
    std::string
    validNamesString () const
    {
      std::ostringstream os;
      printValidNames (os);
      return os.str();
    }

    /// Name of the "default" OrthoManager (e.g., for tests).
    ///
    const std::string& defaultName () const { return theList_[0]; }

    /// Instantiate and return an RCP to the specified MatOrthoManager
    /// subclass.
    ///
    /// \param ortho [in] Name of MatOrthoManager to instantiate
    /// \param M [in] Inner product operator.  If Teuchos::null,
    ///   orthogonalize with respect to the standard Euclidean 
    ///   inner product.
    /// \param params [in] Optional list of parameters for 
    ///   setting up the specific MatOrthoManager subclass
    ///
    /// \return (Smart pointer to an) MatOrthoManager instance
    ///   
    Teuchos::RCP< Belos::OrthoManager< Scalar, MV > >
    makeOrthoManager (const std::string& ortho, 
		      const Teuchos::RCP< const OP >& M,
		      const Teuchos::RCP< Teuchos::ParameterList >& params)
    {
      using Belos::TsqrMatOrthoManager;
      using Belos::ICGSOrthoManager;
      using Belos::IMGSOrthoManager;
      using Belos::DGKSOrthoManager;
      using Teuchos::RCP;
      using Teuchos::rcp;

      Teuchos::ParameterList theParams;
      if (! Teuchos::is_null (params))
	theParams = *params;

      const std::string label("Belos");
      if (ortho == "TSQR") {
	return rcp (new TsqrMatOrthoManager< Scalar, MV, OP > (theParams, label, M));
      }
      else if (ortho == "ICGS") {
	return rcp (new ICGSOrthoManager< Scalar, MV, OP >(label, M));
      }
      else if (ortho == "IMGS") {
	return rcp (new IMGSOrthoManager< Scalar, MV, OP >(label, M));
      }
      else if (ortho == "DGKS") {
	return rcp (new DGKSOrthoManager< Scalar, MV, OP >(label, M));
      }
      else {
	TEST_FOR_EXCEPTION( true, std::invalid_argument, 
			    "Invalid value for command-line parameter \"ortho\":"
			    " valid values are " << validNamesString() 
			    << "." );
      }
    }
  };

} // namespace Belos

#endif // __Belos_OrthoManagerFactory_hpp

