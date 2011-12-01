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

#ifndef __Belos_GmresBaseFactory_hpp
#define __Belos_GmresBaseFactory_hpp

#include <BelosStandardGmres.hpp>

namespace Belos {

  /// \class GmresBaseFactory
  /// \author Mark Hoemmen
  /// \brief Factory for creating GmresBase subclass instances.
  ///
  /// \warning This is EXPERIMENTAL CODE.  DO NOT RELY ON THIS CODE.
  ///   The interface or implementation may change at any time.
  ///
  /// GmresBase describes a general GMRES or Flexible GMRES (FGMRES)
  /// iteration, and leaves the specifics up to the subclasses.  \c
  /// GmresBaseFactory is a factory class that generates subclass
  /// instances.  The factory is responsible for:
  /// - Picking the appropriate subclass of GmresBase (which may
  ///   depend upon the static type of OP (and MV), as well as what
  ///   information the caller provides in the parameter list
  /// - Reading and validating parameters from the given parameter
  ///   list, and filling in missing parameters with default values.
  ///   (The parameter list will indicate whether or not Flexible
  ///   GMRES is to be used, and also other information that the
  ///   factory may use to help pick the subclass.)
  /// - Constructing specialized kernels that depend on the
  ///   operator(s), such as the matrix powers kernel (if applicable).
  ///   This should actually be delegated to a matrix powers kernel
  ///   factory class that "does the right thing," depending on the
  ///   types and any additional information available.
  ///
  template<class Scalar, class MV, class OP>
  class GmresBaseFactory {
  public:
    typedef GmresBase<Scalar, MV, OP> base_type;
    typedef LinearProblem<Scalar, MV, OP> lp_type;
    typedef OrthoManager<Scalar, MV> ortho_type;

    /// \brief Instantiate the appropriate GmresBase subclass instance.
    ///
    /// \param lp [in/out] The linear problem to solve.
    /// \param ortho [in] The orthogonalization method that (F)GMRES
    ///   should use when solving the linear problem.
    /// \param params [in] Optional parameters; if null, defaults are
    ///   used.
    ///
    /// \return GmresBase subclass instance for solving the linear
    ///   system.
    static Teuchos::RCP<base_type> 
    create (const Teuchos::RCP<lp_type>& lp,
	    const Teuchos::RCP<const ortho_type>& ortho,
	    const Teuchos::RCP<OutputManager<Scalar> >& outMan,
	    const Teuchos::RCP<const Teuchos::ParameterList>& params)
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;

      const char prefix[] = "Belos::GmresBaseFactory::create: ";

      // Get the maximum number of iterations.
      //
      // TODO (mfh {16,21} Feb 2011) Pass in maxIterCount and flexible
      // through the parameter list.  They shouldn't be passed as
      // arguments to create(), since they are really specific to the
      // GmresBase implementation.  For example, one could imagine a
      // GmresBase implementation that doesn't do restarting, and that
      // allocates new basis vectors as necessary until convergence or
      // running out of memory.
      //
      // TODO (mfh 16 Feb 2011) Fill in other GMRES implementations.
      const int maxIterCount = 20;
      const bool flexible = false;
      TEUCHOS_TEST_FOR_EXCEPTION(lp.is_null(), std::logic_error, 
			 prefix << "LinearProblem instance is null.");
      TEUCHOS_TEST_FOR_EXCEPTION(ortho.is_null(), std::logic_error, 
			 prefix << "OrthoManager instance is null.");
      return rcp (new StandardGmres<Scalar, MV, OP> (lp, ortho, outMan, 
						     maxIterCount, flexible));
    }
  };

} // namespace Belos


#endif // __Belos_GmresBaseFactory_hpp
