//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
#include <BelosOutputManager.hpp>

#include <Teuchos_StandardCatchMacros.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept> // #include <string>
#include <vector>

//////////////////////////////////////////////////////////////////////

namespace Belos {

  /// \class OrthoManagerFactory
  /// \brief Enumeration of all valid Belos (Mat)OrthoManager classes.
  ///
  /// This factory class knows how to initialize any of Belos' \c
  /// MatOrthoManager subclasses, given a short name of the subclass
  /// (as would naturally belong in a SolverManager's parameter list).
  /// As such, it may be used by any of Belos' \c SolverManager
  /// subclasses that use a (Mat)OrthoManager subclass for
  /// orthogonalization.
  ///
  /// This class' template parameters are the same as those of \c
  /// MatOrthoManager: Scalar is the scalar type (of entries in the
  /// multivector), MV is the multivector type, and OP is the operator
  /// type.  For example: Scalar=double, MV=Epetra_MultiVector, and
  /// OP=Epetra_Operator.
  template<class Scalar, class MV, class OP>
  class OrthoManagerFactory {
  private:
    //! List of valid OrthoManager names.
    std::vector<std::string> theList_;

  public:
    /// \brief Number of MatOrthoManager subclasses this factory recognizes.
    static int numOrthoManagers () {
#ifdef HAVE_BELOS_TSQR
      return 5; 
#else
      return 4;
#endif // HAVE_BELOS_TSQR
    }

    /// \brief Is the given MatOrthoManager subclass rank-reealing?
    ///
    /// Return true if and only if the given MatOrthoManager name is
    /// that of a MatOrthoManager subclass with rank-revealing
    /// capability.
    static bool isRankRevealing (const std::string& name) {
#ifdef HAVE_BELOS_TSQR
      // Currently only TSQR has a full rank-revealing capability.
      return (name == "TSQR");
#else
      return false;
#endif // HAVE_BELOS_TSQR
    }

    //! Constructor.
    OrthoManagerFactory () : theList_ (numOrthoManagers())
    {
      int index = 0;
      theList_[index++] = "DGKS";
      theList_[index++] = "ICGS";
      theList_[index++] = "IMGS";
#ifdef HAVE_BELOS_TSQR
      theList_[index++] = "TSQR";
#endif // HAVE_BELOS_TSQR
      theList_[index++] = "Simple";
    }

    /// \brief List of MatOrthoManager subclasses this factory recognizes.
    /// 
    /// This is useful as a list of valid command-line parameter
    /// values for choosing a MatOrthoManager subclass to test.
    ///
    /// \note To implementers: Anasazi and Belos currently implement
    ///   different sets of (Mat)OrthoManagers.  This method returns a
    ///   valid list of Belos MatOrthoManager subclasses.
    const std::vector<std::string>& 
    validNames () const { return theList_; }

    //! Whether this factory recognizes the MatOrthoManager with the given name.
    bool
    isValidName (const std::string& name) const 
    {
      return (std::find (theList_.begin(), theList_.end(), name) != theList_.end());
    }

    //! Print all recognized MatOrthoManager names to the given ostream.
    std::ostream&
    printValidNames (std::ostream& out) const
    {
      const int numValid = numOrthoManagers();
      TEUCHOS_TEST_FOR_EXCEPTION(numValid <= 0, std::logic_error,
			 "Invalid number " << numValid << " of valid MatOrtho"
			 "Manager names.  Please report this bug to the Belos "
			 "developers." );
      if (numValid > 1) {
	for (int k = 0; k < numValid - 1; ++k)
	  out << "\"" << theList_[k] << "\", ";
	out << "or ";
      }
      out << "\"" << theList_[numValid-1] << "\"";
      return out;
    }

    /// \brief List (as a string) of recognized MatOrthoManager names.
    ///
    /// This is useful for generating help for command-line arguments,
    /// when writing a test that uses different orthogonalizations.
    std::string
    validNamesString () const
    {
      std::ostringstream os;
      (void) printValidNames (os);
      return os.str();
    }

    /// \brief Name of the "default" MatOrthoManager subclass.
    ///
    /// This is the name of the MatOrthoManager subclass that serves
    /// as a reasonable default for all Belos solvers that use general
    /// orthogonalizations.  It may not be the fastest or the most
    /// accurate, but it should be the most reasonable.
    const std::string& defaultName () const { return theList_[0]; }

    /// \brief Default parameters for the given MatOrthoManager subclass.
    ///
    /// \param name [in] MatOrthoManager subclass short name, for
    ///   which isValidName(name) returns true.
    ///
    /// \note This method does not cache its return value, so every
    ///   time you call this method with the same \c name argument, a
    ///   new parameter list will be created.  You can save the
    ///   parameter list yourself if you want to reuse it.
    Teuchos::RCP<const Teuchos::ParameterList> 
    getDefaultParameters (const std::string& name) const
    {
      if (name == "DGKS") {
	DGKSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getValidParameters ();
      }
#ifdef HAVE_BELOS_TSQR
      else if (name == "TSQR") {
	TsqrMatOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getValidParameters ();
      }
#endif // HAVE_BELOS_TSQR
      else if (name == "ICGS") {
	ICGSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getValidParameters ();
      }
      else if (name == "IMGS") {
	IMGSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getValidParameters ();
      }
      else if (name == "Simple") {
	IMGSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getValidParameters ();
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
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

    /// \brief "Fast" parameters for the given MatOrthoManager subclass.
    ///
    /// "Fast" usually means that accuracy and/or robustness (with
    /// respect to rank deficiency) might be compromised in order to
    /// improve performance.
    ///
    /// \param name [in] MatOrthoManager subclass short name, for
    ///   which isValidName(name) returns true.
    ///
    /// \note This method does not cache its return value, so every
    ///   time you call this method with the same \c name argument, a
    ///   new parameter list will be created.  You can save the
    ///   parameter list yourself if you want to reuse it.
    Teuchos::RCP<const Teuchos::ParameterList> 
    getFastParameters (const std::string& name) const
    {
      if (name == "DGKS") {
	DGKSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getFastParameters ();
      }
#ifdef HAVE_BELOS_TSQR
      else if (name == "TSQR") {
	TsqrMatOrthoManager<Scalar, MV, OP> orthoMan; 
	return orthoMan.getFastParameters ();
      }
#endif // HAVE_BELOS_TSQR
      else if (name == "ICGS") {
	ICGSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getFastParameters ();
      }
      else if (name == "IMGS") {
	IMGSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getFastParameters ();
      }
      else if (name == "Simple") {
	IMGSOrthoManager<Scalar, MV, OP> orthoMan;
	return orthoMan.getFastParameters ();
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
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

    /// \brief Return an instance of the specified MatOrthoManager subclass.
    ///
    /// \param ortho [in] Name of the MatOrthoManager subclass
    ///   instance to return.  The \c validNames() method returns a
    ///   list of the supported names.
    /// \param M [in] Inner product operator.  If Teuchos::null,
    ///   orthogonalize with respect to the standard Euclidean 
    ///   inner product.
    /// \param outMan [in/out] Output manager, which the OrthoManager
    ///   instance may use (but is not required to use) for various
    ///   kinds of status output.
    /// \param label [in] Label for Belos-specific timers, if Belos
    ///   timers were enabled at compile time.  Otherwise, this
    ///   parameter's value doesn't matter.
    /// \param params [in/out] Optional list of parameters for setting
    ///   up the specific MatOrthoManager subclass.  A default
    ///   parameter list with embedded documentation is available for
    ///   each MatOrthoManager subclass that this factory knows how to
    ///   make.
    Teuchos::RCP<Belos::MatOrthoManager<Scalar, MV, OP> >
    makeMatOrthoManager (const std::string& ortho, 
			 const Teuchos::RCP<const OP>& M,
			 const Teuchos::RCP<OutputManager<Scalar> >& outMan,
			 const std::string& label,
			 const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
#ifdef HAVE_BELOS_TSQR
      using Belos::TsqrMatOrthoManager;
#endif // HAVE_BELOS_TSQR
      using Belos::ICGSOrthoManager;
      using Belos::IMGSOrthoManager;
      using Belos::DGKSOrthoManager;
      using Belos::SimpleOrthoManager;
      using Teuchos::rcp;

      if (ortho == "DGKS") {
	typedef DGKSOrthoManager<Scalar, MV, OP> dgks_type;
	return rcp (new dgks_type (params, label, M));
      }
#ifdef HAVE_BELOS_TSQR
      else if (ortho == "TSQR") {
	typedef TsqrMatOrthoManager<Scalar, MV, OP> ortho_type;
	return rcp (new ortho_type (params, label, M));
      }
#endif // HAVE_BELOS_TSQR
      else if (ortho == "ICGS") {
	typedef ICGSOrthoManager<Scalar, MV, OP> ortho_type;
	return rcp (new ortho_type (params, label, M));
      }
      else if (ortho == "IMGS") {
	typedef IMGSOrthoManager<Scalar, MV, OP> ortho_type;
	return rcp (new ortho_type (params, label, M));
      } 
      else if (ortho == "Simple") {
	TEUCHOS_TEST_FOR_EXCEPTION(ortho == "Simple", std::logic_error,
				   "SimpleOrthoManager does not yet support "
				   "the MatOrthoManager interface");
      } 
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
			   "Invalid orthogonalization manager name: Valid names"
			   " are " << validNamesString() << ".  For many of "
			   "the test executables, the orthogonalization manager"
			   " name often corresponds to the \"ortho\" command-"
			   "line argument.");
      return Teuchos::null; // Guard to avoid compiler warnings.
    }

    /// \brief Return an instance of the specified OrthoManager subclass.
    ///
    /// \param ortho [in] Name of the OrthoManager subclass instance
    ///   to return.  The \c validNames() method returns a list of the
    ///   supported names.
    /// \param M [in] Inner product operator.  If Teuchos::null,
    ///   orthogonalize with respect to the standard Euclidean inner
    ///   product.
    /// \param outMan [in/out] Output manager, which the OrthoManager
    ///   instance may use (but is not required to use) for various
    ///   kinds of status output.
    /// \param label [in] Label for timers.
    /// \param params [in/out] Optional list of parameters for 
    ///   setting up the specific OrthoManager subclass.
    ///
    /// \return OrthoManager instance.
    Teuchos::RCP<Belos::OrthoManager<Scalar, MV> >
    makeOrthoManager (const std::string& ortho, 
		      const Teuchos::RCP<const OP>& M,
		      const Teuchos::RCP<OutputManager<Scalar> >& outMan,
		      const std::string& label,
		      const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
#ifdef HAVE_BELOS_TSQR
      using Belos::TsqrOrthoManager;
#endif // HAVE_BELOS_TSQR
      using Teuchos::rcp;

      if (ortho == "Simple") {
	TEUCHOS_TEST_FOR_EXCEPTION(! M.is_null(), std::logic_error,
				   "SimpleOrthoManager is not yet supported "
				   "when the operator M is nontrivial (i.e., "
				   "M != null).");
	return rcp (new SimpleOrthoManager<Scalar, MV> (outMan, label, params));
      }
#ifdef HAVE_BELOS_TSQR
      // TsqrMatOrthoManager has to store more things and do more work
      // than TsqrOrthoManager, in order for the former to be correct
      // for the case of a nondefault (non-Euclidean) inner product.
      // Thus, it's better to create a TsqrOrthoManager, when we know
      // the operator is the default operator (M is null).  Of course,
      // a MatOrthoManager is-an OrthoManager, so returning a
      // TsqrMatOrthoManager would still be correct; this is just an
      // optimization.
      else if (ortho == "TSQR" && M.is_null()) {
	return rcp (new TsqrOrthoManager<Scalar, MV> (params, label));
      }
#endif // HAVE_BELOS_TSQR
      else {
	// A MatOrthoManager is-an OrthoManager.
	return makeMatOrthoManager (ortho, M, outMan, label, params);
      }
    }
  };

} // namespace Belos

#endif // __Belos_OrthoManagerFactory_hpp

