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

#ifndef __Belos_AkxFactory_hpp
#define __Belos_AkxFactory_hpp

#include <BelosMonomialOpAkx.hpp>
#include <Teuchos_ParameterList.hpp>

#include <algorithm>
#include <cctype>

// TODO Include other OpAkx implementations, and other Akx
// implementations (that don't implement OpAkx) as they are completed.

namespace Belos {

  /// \class AkxFactory
  /// \brief Create and return Akx subclass instances.
  /// \author Mark Hoemmen
  ///
  /// This factory knows how to create various types of matrix powers
  /// kernel bases, all of which implement the Belos::Akx interface.
  /// As new matrix powers kernel implementations are added to
  /// Trilinos, this factory class will help you instantiate them and
  /// use them in solvers.
  ///
  template<class Scalar, class MV>
  class AkxFactory {
  public:
    /// \typedef akx_base_type
    ///
    /// Type of the matrix powers kernel abstract interface.
    typedef Akx<Scalar, MV> akx_base_type;    

    /// \brief Constructor.
    ///
    /// It doesn't do anything exciting.
    AkxFactory();

    /// \brief Make an OpAkx<Scalar, MV, OP> subclass instance.
    ///
    /// OpAkx implements the matrix powers kernel using
    /// straightforward applications of the operator A (and the left
    /// and/or right preconditioner, if applicable).  It's a good
    /// default choice for (Flexible) CA-GMRES, especially with
    /// reasonably long restart cycles (so that the main benefit of
    /// CA-GMRES would likely come from TSQR and Block Gram-Schmidt,
    /// not from the matrix powers kernel).
    ///
    /// \param A [in] The matrix of interest.
    /// \param M_left [in] If not null, the left preconditioner,
    ///   or split preconditioner if M_right is also not null.
    /// \param M_right [in] If not null, the right preconditioner,
    ///   or split preconditioner if M_left is also not null.
    /// \param params [in] List of parameters to help pick and
    ///   configure the matrix powers kernel implementation.  If null,
    ///   default parameters are used.
    ///   
    /// \return Instance of an OpAkx<Scalar, MV, OP> subclass.
    ///
    template<class OP>
    Teuchos::RCP<akx_base_type>
    makeOpAkx (const Teuchos::RCP<const OP>& A, 
	       const Teuchos::RCP<const OP>& M_left, 
	       const Teuchos::RCP<const OP>& M_right,
	       const Teuchos::RCP<const Teuchos::ParameterList>& params);

    //! Return the default matrix powers kernel parameter list.
    Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ();

    /// \brief Return name of the matrix powers kernel parameter list.
    ///
    /// This method is useful if you are looking for parameters for
    /// the matrix powers kernel as a sublist of a solver's parameter
    /// list.
    const std::string& parameterListName () const {
      return listName_;
    }

  private:
    /// \brief Validate and return canonical basis name.
    ///
    /// First, canonicalize basisName.  Return the canonical name if
    /// it names one of the one of the accepted types of matrix powers
    /// kernel basis.  Otherwise, throw an invalid_argument exception.
    ///
    /// \param basisName [in] Name of a matrix powers kernel basis.
    /// 
    /// \return Canonical matrix powers kernel basis name.
    std::string
    validateBasisName (const std::string& basisName) const;

    //! Name of the matrix powers kernel parameter list.
    const std::string listName_;

    //! List of valid matrix powers kernel basis names.
    std::vector<std::string> validBasisNames_;
    
    /// \brief The default parameter list.  
    ///
    /// This is initialized lazily in \c getDefaultParameters().
    Teuchos::RCP<const Teuchos::ParameterList> defaultParams_;
  };


  template<class Scalar, class MV>
  AkxFactory<Scalar, MV>::AkxFactory() : listName_ ("Akx") 
  {
    typedef std::vector<std::string>::size_type size_type;
    //
    // TODO (mfh 14 Feb 2011) Support more basis types later.
    //
    const size_type numValidBasisNames = 1;
    const char* validBasisNames[] = {"Monomial"};
    validBasisNames_.resize (numValidBasisNames);

    for (size_type k = 0; k < numValidBasisNames; ++k)
      validBasisNames_.push_back (validBasisNames[k]);
  }

  template<class Scalar, class MV>
  Teuchos::RCP<typename AkxFactory<Scalar, MV>::akx_base_type>
  AkxFactory<Scalar, MV>::
  makeOpAkx (const Teuchos::RCP<const OP>& A, 
	     const Teuchos::RCP<const OP>& M_left, 
	     const Teuchos::RCP<const OP>& M_right,
	     const Teuchos::RCP<const Teuchos::ParameterList>& params)
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::Exceptions::InvalidParameter;
    using Teuchos::Exceptions::InvalidParameterName;
    using Teuchos::Exceptions::InvalidParameterType;

    RCP<const ParameterList> plist = 
      params.is_null() ? getDefaultParameters() : params;
    
    std::string basisType;
    try {
      basisType = plist->get<std::string> ("Basis Type");
    } catch (InvalidParameter& e) {
      // FIXME (mfh 14 Feb 2011, 03 Mar 2011)
      //
      // Just rethrow for now.  Later we might wrap and rethrow, or
      // implement more tolerant behavior.
      throw e; 
    }
    // Canonicalize (case, etc.) the string and validate it.  Throw an
    // exception if invalid, else return the canonicalized name.
    const std::string validName = validateBasisName (basisType);

    // Recommended (at least initially) basis length.
    int basisLength = 5; // a reasonable default
    try {
      basisLength = plist->get<int> ("Basis Length");
    } catch (InvalidParameter& e) {
      // Do nothing; let the default stay
    }
    TEST_FOR_EXCEPTION(basisLength < 1, std::invalid_argument,
		       "The \"Basis Length\" parameter must be >= 1, "
		       "but its value here is " << basisLength << ".");
    if (validName == "Monomial")
      {
	typedef MonomialOpAkx<Scalar, MV, OP> akx_impl_type;
	RCP<akx_impl_type> akx (new akx_impl_type (A, M_left, M_right, basisLength));
	return akx;
      }
    else
      {
	TEST_FOR_EXCEPTION(validName != "Monomial", std::logic_error,
			   "We have not yet implemented the \"" 
			   << validName << "\" basis.");
	throw std::logic_error("Should never get here!");
      }
  }

  template<class Scalar, class MV>
  Teuchos::RCP<const Teuchos::ParameterList>
  AkxFactory<Scalar, MV>::getDefaultParameters ()
  {
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;

    if (defaultParams_.is_null())
      {
	RCP<ParameterList> plist = parameterList (parameterListName ());
	// Default basis type is monomial for now; this will likely
	// change to the Newton basis.
	const std::string defaultBasis ("Monomial");
	{
	  std::ostringstream os;
	  os << "Default basis type.  Valid choices include: {";
	  for (std::vector<std::string>::size_type k = 0; 
	       k < validBasisNames_.size(); ++k)
	    {
	      os << validBasisNames_[k];
	      if (k < validBasisNames_.size() - 1)
		os << ", ";
	    }
	  os << "}.";
	  plist->set ("Basis Type", defaultBasis, os.str());
	}
	// An initial basis length of 5 is a reasonable first guess.
	// Solvers should revise this dynamically, based on their
	// estimates of the condition number and rank of the generated
	// candidate basis.
	const int defaultBasisLength = 5;
	plist->set ("Basis Length", defaultBasisLength, 
		    "Default recommended basis length.");
	defaultParams_ = plist;
      }
    return defaultParams_;
  }

  template<class Scalar, class MV>
  std::string
  AkxFactory<Scalar, MV>::validateBasisName (const std::string& basisName) const
  {
    // Canonicalize the basis name.  First, strip whitespace.
    // Second, capitalize the first letter and lowercase the rest.
    TEST_FOR_EXCEPTION(basisName.empty(), std::invalid_argument,
		       "The matrix powers kernel basis name is an empty "
		       "string.");
    const size_t npos = std::string::npos;
    size_t firstNonSpacePos = basisName.find_first_not_of (" \t\n");
    size_t lastNonSpacePos = basisName.find_last_not_of (" \t\n");
    TEST_FOR_EXCEPTION(firstNonSpacePos == npos, std::invalid_argument,
		       "The matrix powers kernel basis name \"" << basisName 
		       << "\" contains only whitespace.");
    // canonName must have length at least one.
    std::string canonName = 
      basisName.substr (firstNonSpacePos, lastNonSpacePos-firstNonSpacePos+1);
    // Capitalize the first letter and lowercase the rest, in place.
    canonName[0] = toupper (canonName[0]);
    std::transform (canonName.begin()+1, canonName.end(), 
		    canonName.begin()+1, tolower);
    const bool foundIt = validBasisNames_.end() != 
      std::find (validBasisNames_.begin(), validBasisNames_.end(), canonName);
    TEST_FOR_EXCEPTION(! foundIt, std::invalid_argument,
		       "Invalid basis name \"" << basisName << "\".");
    return canonName;
  }


} // namespace Belos

#endif // __Belos_AkxFactory_hpp
