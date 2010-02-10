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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
