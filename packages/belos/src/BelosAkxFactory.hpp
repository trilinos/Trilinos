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

// TODO Include other OpAkx implementations, and other Akx
// implementations (that don't implement OpAkx) as they are completed.

namespace Belos {

  /// \class AkxFactory
  /// \brief Create and return Akx subclass instances.
  /// \author Mark Hoemmen
  ///
  template<class Scalar, class MV>
  class AkxFactory {
  public:
    typedef Akx<Scalar, MV> akx_base_type;    

    AkxFactory() {}

    /// \brief Make an OpAkx<Scalar, MV, OP> subclass instance.
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
    template<class OP>
    Teuchos::RCP<akx_base_type>
    makeOpAkx (const Teuchos::RCP<const OP>& A, 
	       const Teuchos::RCP<const OP>& M_left, 
	       const Teuchos::RCP<const OP>& M_right,
	       const Teuchos::RCP<const Teuchos::ParameterList>& params);

    //! Return the default parameter list.
    Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ();

  private:
    /// Throw an invalid_argument exception if basisType does not name
    /// one of the accepted matrix powers kernel basis types.
    void 
    validateBasisType (const std::string& basisType);

    //! The default parameter list.  Initialized lazily.
    Teuchos::RCP<const Teuchos::ParameterList> defaultParams_;
  };

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
    // FIXME (mfh 14 Feb 2011) Catch these exceptions.  For now we let
    // them pass through, since users aren't obligated to supply
    // parameters.
    try {
      basisType = plist->get<std::string> ("Basis Type");
    } catch (InvalidParameter& e) {
      // Just rethrow for now.  Later we might wrap and rethrow, or
      // implement more tolerant behavior.
      throw e; 
    }
    // TODO (mfh 14 Feb 2011) Canonicalize (case, etc.) the string
    // before validating it.
    validateBasisType (basisType);

    if (basisType == "Monomial")
      {
	typedef MonomialOpAkx<Scalar, MV, OP> akx_impl_type;
	return RCP<akx_base_type> (new akx_impl_type (A, M_left, M_right));
      }
    else
      throw std::logic_error("Should never get here!");
  }

  template<class Scalar, class MV>
  Teuchos::RCP<const Teuchos::ParameterList>
  AkxFactory<Scalar, MV>::getDefaultParameters ()
  {
    using Teuchos::RCP;
    using Teuchos::ParameterList;

    if (defaultParams_.is_null())
      {
	RCP<ParameterList> plist = Teuchos::parameterList();

	const std::string defaultBasis ("Monomial");
	plist->set ("Basis Type", defaultBasis, "Default basis type.");

	defaultParams_ = plist;
      }
    return defaultParams_;
  }

  template<class Scalar, class MV>
  void
  AkxFactory<Scalar, MV>::validateBasisType (const std::string& basisType)
  {
    // TODO (mfh 14 Feb 2011) Support more basis types later.
    const char* validValues[] = {"Monomial"};
    const int numValidValues = 1;
    const bool foundIt = 
      validValues+numValid == 
      std::find (validValues, validValues+numValid, basisType);
    TEST_FOR_EXCEPTION(! foundIt, std::invalid_argument,
		       "Invalid basis type \"" << basisType << "\".");
  }


} // namespace Belos

#endif // __Belos_AkxFactory_hpp
