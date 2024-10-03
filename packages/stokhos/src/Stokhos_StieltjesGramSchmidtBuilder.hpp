// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
