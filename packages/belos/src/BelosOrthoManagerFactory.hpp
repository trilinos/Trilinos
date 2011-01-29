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
#ifdef HAVE_BELOS_TSQR
#  include <BelosTsqrOrthoManager.hpp>
#endif // HAVE_BELOS_TSQR
#include <BelosICGSOrthoManager.hpp>
#include <BelosIMGSOrthoManager.hpp>
#include <BelosDGKSOrthoManager.hpp>
#include <BelosSimpleOrthoManager.hpp>

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
    //! List of valid OrthoManager names
    std::vector< std::string > theList_;

  public:
    /// Number of valid command-line parameter values for the OrthoManager
    /// subclass to test.  (Must be at least one.)
    static int numOrthoManagers () {
#ifdef HAVE_BELOS_TSQR
      return 5; 
#else
      return 4;
#endif // HAVE_BELOS_TSQR
    }

    /// Return true if and only if the given OrthoManager name is that
    /// of an OrthoManager subclass with rank-revealing capability.
    static bool isRankRevealing (const std::string& name) {
#ifdef HAVE_BELOS_TSQR
      // Currently only TSQR has a full rank-revealing capability.
      return (name == "TSQR");
#else
      return false;
#endif // HAVE_BELOS_TSQR
    }

    //! Constructor
    OrthoManagerFactory () : theList_ (numOrthoManagers())
    {
      int index = 0;
#ifdef HAVE_BELOS_TSQR
      theList_[index++] = "TSQR";
#endif // HAVE_BELOS_TSQR
      theList_[index++] = "ICGS";
      theList_[index++] = "IMGS";
      theList_[index++] = "DGKS";
      theList_[index++] = "Simple";
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

    //! Whether 'name' names a valid OrthoManager
    bool
    isValidName (const std::string& name) const 
    {
      return (std::find (theList_.begin(), theList_.end(), name) != theList_.end());
    }

    //! Print the list of valid OrthoManager names to the given ostream
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

    //! Name of the "default" OrthoManager (e.g., for tests).
    const std::string& defaultName () const { return theList_[0]; }

    /// Default parameters for the given orthogonalization manager
    ///
    /// \param name [in] Orthogonalization manager name.  Must be
    ///   supported by this factory (i.e., validName(name)==true).
    ///
    /// \return Default parameters for the given orthogonalization
    ///   manager
    /// 
    /// \warning This method may not be reentrant.
    Teuchos::RCP<const Teuchos::ParameterList> 
    getDefaultParameters (const std::string& name)
    {
      if (name == "DGKS") {
	return getDefaultDgksParameters<Scalar>();
      }
#ifdef HAVE_BELOS_TSQR
      else if (name == "TSQR") {
	return TsqrMatOrthoManager<Scalar, MV, OP>::getDefaultParameters();
      }
#endif // HAVE_BELOS_TSQR
      else if (name == "ICGS") {
	return getDefaultIcgsParameters<Scalar>();
      }
      else if (name == "IMGS") {
	return getDefaultImgsParameters<Scalar>();
      }
      else if (name == "Simple") {
	return SimpleOrthoManager<Scalar, MV>::getDefaultParameters();
      }
      else {
	TEST_FOR_EXCEPTION(true, std::invalid_argument, 
			   "Invalid orthogonalization manager name \"" << name 
			   << "\": Valid names are " << validNamesString() 
			   << ".  For many of the test executables, the "
			   "orthogonalization manager name often corresponds "
			   "to the \"ortho\" command-line argument.");
	// Placate the compiler if necessary; we should never reach
	// this point.
	return Teuchos::null; 
      }
    }


    /// \brief "Fast" parameters for the given (Mat)OrthoManager
    ///
    /// "Fast" parameters for the given orthogonalization manager.
    /// "Fast" usually means that accuracy and/or robustness with
    /// respect to rank deficiency might be compromised for the sake
    /// of better performance.
    ///
    /// \param name [in] Orthogonalization manager name.  Must be
    ///   supported by this factory (i.e., validName(name)==true).
    ///
    /// \return "Fast" parameters for the given orthogonalization
    ///   manager
    /// 
    /// \warning This method may not be reentrant.
    Teuchos::RCP<const Teuchos::ParameterList> 
    getFastParameters (const std::string& name)
    {
      if (name == "DGKS") {
	return getFastDgksParameters<Scalar>();
      }
#ifdef HAVE_BELOS_TSQR
      else if (name == "TSQR") {
	return TsqrMatOrthoManager<Scalar, MV, OP>::getFastParameters();
      }
#endif // HAVE_BELOS_TSQR
      else if (name == "ICGS") {
	return getFastIcgsParameters<Scalar>();
      }
      else if (name == "IMGS") {
	return getFastImgsParameters<Scalar>();
      }
      else if (name == "Simple") {
	return SimpleOrthoManager<Scalar, MV>::getFastParameters();
      }
      else {
	TEST_FOR_EXCEPTION(true, std::invalid_argument, 
			   "Invalid orthogonalization manager name \"" << name 
			   << "\": Valid names are " << validNamesString() 
			   << ".  For many of the test executables, the "
			   "orthogonalization manager name often corresponds "
			   "to the \"ortho\" command-line argument.");
	// Placate the compiler if necessary; we should never reach
	// this point.
	return Teuchos::null; 
      }
    }


    /// Create and return the specified MatOrthoManager subclass
    ///
    /// \param ortho [in] Name of MatOrthoManager subclass.  The \c
    ///   validNames() method returns a list of the supported names.
    /// \param M [in] Inner product operator.  If Teuchos::null,
    ///   orthogonalize with respect to the standard Euclidean 
    ///   inner product.
    /// \param label [in] Label for timers
    /// \param params [in] Optional (null) list of parameters for
    ///   setting up the specific MatOrthoManager subclass.  A default
    ///   parameter list with embedded documentation is available for
    ///   each MatOrthoManager subclass that this factory knows how to
    ///   make.
    ///
    /// \return (Smart pointer to a) MatOrthoManager instance
    ///   
    Teuchos::RCP< Belos::MatOrthoManager<Scalar, MV, OP> >
    makeMatOrthoManager (const std::string& ortho, 
			 const Teuchos::RCP<const OP>& M,
			 const std::string& label,
			 const Teuchos::RCP<const Teuchos::ParameterList>& params)
    {
#ifdef HAVE_BELOS_TSQR
      using Belos::TsqrMatOrthoManager;
#endif // HAVE_BELOS_TSQR
      using Belos::ICGSOrthoManager;
      using Belos::IMGSOrthoManager;
      using Belos::DGKSOrthoManager;
      using Belos::SimpleOrthoManager;
      using Teuchos::rcp;
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

      TEST_FOR_EXCEPTION(ortho == "Simple", std::logic_error,
			 "SimpleOrthoManager does not yet support the "
			 "MatOrthoManager interface");

      if (ortho == "DGKS") {
	int maxNumOrthogPasses;
	magnitude_type blkTol, depTol, singTol;
	readDgksParameters<Scalar> (params, maxNumOrthogPasses, blkTol, depTol, singTol);
	return rcp (new DGKSOrthoManager<Scalar, MV, OP> (label, M, maxNumOrthogPasses, 
							  blkTol, depTol, singTol));
      }
#ifdef HAVE_BELOS_TSQR
      else if (ortho == "TSQR") {
	// mfh 12 Jan 2011: TSQR knows how to read its own parameters.
	// I didn't want to change the other OrthoManager subclasses'
	// public interfaces to accept a parameter list input.
	return rcp (new TsqrMatOrthoManager<Scalar, MV, OP> (params, label, M));
      }
#endif // HAVE_BELOS_TSQR
      else if (ortho == "ICGS") {
	int maxNumOrthogPasses;
	magnitude_type blkTol, singTol;
	readIcgsParameters<Scalar> (params, maxNumOrthogPasses, blkTol, singTol);
	return rcp (new ICGSOrthoManager<Scalar, MV, OP>(label, M, maxNumOrthogPasses, 
							 blkTol, singTol));
      }
      else if (ortho == "IMGS") {
	int maxNumOrthogPasses;
	magnitude_type blkTol, singTol;
	readImgsParameters<Scalar> (params, maxNumOrthogPasses, blkTol, singTol);
	return rcp (new IMGSOrthoManager<Scalar, MV, OP>(label, M, maxNumOrthogPasses,
							 blkTol, singTol));
      }
      else {
	TEST_FOR_EXCEPTION(true, std::invalid_argument, 
			   "Invalid orthogonalization manager name: Valid names"
			   " are " << validNamesString() << ".  For many of "
			   "the test executables, the orthogonalization manager"
			   " name often corresponds to the \"ortho\" command-"
			   "line argument.");
      }
    }

    /// Create and return the specified OrthoManager subclass
    ///
    /// \param ortho [in] Name of OrthoManager subclass.  The \c
    ///   validNames() method returns a list of the supported names.
    /// \param M [in] Inner product operator.  If Teuchos::null,
    ///   orthogonalize with respect to the standard Euclidean 
    ///   inner product.
    /// \param label [in] Label for timers
    /// \param params [in] Optional list of parameters for 
    ///   setting up the specific OrthoManager subclass
    ///
    /// \return (Smart pointer to a) OrthoManager instance
    ///   
    Teuchos::RCP<Belos::OrthoManager<Scalar, MV> >
    makeOrthoManager (const std::string& ortho, 
		      const Teuchos::RCP<const OP>& M,
		      const std::string& label,
		      const Teuchos::RCP<const Teuchos::ParameterList>& params)
    {
#ifdef HAVE_BELOS_TSQR
      using Belos::TsqrOrthoManager;
#endif // HAVE_BELOS_TSQR
      using Teuchos::rcp;

#ifdef HAVE_BELOS_TSQR
      // TsqrMatOrthoManager has to store more things and do more work
      // than TsqrOrthoManager, in order for the former to be correct
      // for the case of a nondefault (non-Euclidean) inner product.
      // Thus, it's better to create a TsqrOrthoManager, when we know
      // the operator is the default operator (M is null).  Of course,
      // a MatOrthoManager is-an OrthoManager, so returning a
      // TsqrMatOrthoManager would still be correct; this is just an
      // optimization.
      if (ortho == "TSQR" && M.is_null())
	return rcp (new TsqrOrthoManager<Scalar, MV> (params, label));
#endif // HAVE_BELOS_TSQR

      if (ortho == "Simple")
	{
	  TEST_FOR_EXCEPTION(! M.is_null(), std::logic_error,
			     "SimpleOrthoManager is not yet supported "
			     "when the operator M is nontrivial (i.e., "
			     "M != null).");
	  return rcp (new SimpleOrthoManager<Scalar, MV> (label, params));
	}
      // A MatOrthoManager is-an OrthoManager.
      return makeMatOrthoManager (ortho, M, label, params);
    }
  };

} // namespace Belos

#endif // __Belos_OrthoManagerFactory_hpp

