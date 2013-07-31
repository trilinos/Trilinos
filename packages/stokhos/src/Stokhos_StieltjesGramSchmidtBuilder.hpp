// $Id: Stokhos_SGModelEvaluator.cpp,v 1.10 2009/10/06 16:51:22 agsalin Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_SGModelEvaluator.cpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_STIELTJES_GRAM_SCHMIDT_BUILDER_HPP
#define STOKHOS_STIELTJES_GRAM_SCHMIDT_BUILDER_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_GramSchmidtBasis.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

namespace Stokhos {

  /*! 
   * \brief Class for building a reduced-dimension basis and quadrature from
   * a given set of polynomial chaos expansions.  First generates 1-D
   * orthogonal bases using the discretized Stieltjes procedure, forms their
   * tensor product, and then orthogonalizes using Gram-Schmidt.
   */
  template <typename ordinal_type, typename value_type>
  class StieltjesGramSchmidtBuilder {
  public:
  
    //! Constructor
    StieltjesGramSchmidtBuilder(
      const Teuchos::RCP<const Quadrature<ordinal_type, value_type> >& quad,
      const Teuchos::Array< OrthogPolyApprox<ordinal_type, value_type> >& pces,
      ordinal_type new_order,  bool use_pce_qp, bool normalize);

    //! Destructor
    ~StieltjesGramSchmidtBuilder() {}

    //! Get reduced basis
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >
    getReducedBasis() const;

    //! Get reduced quadrature
    Teuchos::RCP<Quadrature<ordinal_type, value_type> >
    getReducedQuadrature() const;

    //! Get reduced PCEs
    void 
    computeReducedPCEs(
      const Teuchos::Array< OrthogPolyApprox<ordinal_type, value_type> >& pces,
      Teuchos::Array< OrthogPolyApprox<ordinal_type, value_type> >& new_pces);

  private:

    // Prohibit copying
    StieltjesGramSchmidtBuilder(const StieltjesGramSchmidtBuilder&);

    // Prohibit Assignment
    StieltjesGramSchmidtBuilder& operator=(const StieltjesGramSchmidtBuilder& b);
    
  protected:
    
    //! Quadrature object for original basis
    Teuchos::RCP<const Quadrature<ordinal_type, value_type> > quad;

    //! Reduced tensor basis
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> > tensor_basis;

    //! Reduced Gram-Schmidt basis
    Teuchos::RCP< GramSchmidtBasis<ordinal_type,value_type> > gs_basis;

    //! Reduced quadrature
    Teuchos::RCP< UserDefinedQuadrature<ordinal_type, value_type> > gs_quad;

  }; // class StieltjesGramSchmidtBuilder

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_StieltjesGramSchmidtBuilderImp.hpp"

#endif // STOKHOS_STIELTJES_GRAM_SCHMIDT_BUILDER_HPP
