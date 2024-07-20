// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"

#include "Stokhos_BasisFactory.hpp"
#include "Stokhos_QuadratureFactory.hpp"
#include "Stokhos_AlgebraicOrthogPolyExpansion.hpp"
#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_ForUQTKOrthogPolyExpansion.hpp"
//#include "Stokhos_DerivOrthogPolyExpansion.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_PseudoSpectralOrthogPolyExpansion.hpp"
#include "Stokhos_PseudoSpectralOperatorFactory.hpp"

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OrthogPolyExpansion<ordinal_type, value_type> >
Stokhos::ExpansionFactory<ordinal_type, value_type>::
create(Teuchos::ParameterList& sgParams)
{
  // Check if expansion is already there
  Teuchos::ParameterList& expParams = sgParams.sublist("Expansion");
  Teuchos::RCP< Stokhos::OrthogPolyExpansion<ordinal_type,value_type> > expansion = expParams.template get< Teuchos::RCP< Stokhos::OrthogPolyExpansion<ordinal_type,value_type> > >("Stochastic Galerkin Expansion", Teuchos::null);
  if (expansion != Teuchos::null)
    return expansion;

  // Get basis
  Teuchos::ParameterList& basisParams = sgParams.sublist("Basis");
  Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis;
  if (basisParams.template isType< Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis"))
    basis = basisParams.template get< Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis");
  else {
    basis = Stokhos::BasisFactory<ordinal_type,value_type>::create(sgParams);
    basisParams.set("Stochastic Galerkin Basis", basis);
  }
    
  // Get 3-tensor
  Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type,value_type> > Cijk;
  if (sgParams.template isType<Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type,value_type> > >("Triple Product Tensor"))
    Cijk = sgParams.template get<Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type,value_type> > >("Triple Product Tensor");
  else {
    std::string tp_type = sgParams.get("Triple Product Size", "Full");
    TEUCHOS_TEST_FOR_EXCEPTION(
      tp_type != "Full" && tp_type != "Linear", 
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<  "Invalid triple product expansion type  " << tp_type <<
      std::endl);

    if (tp_type == "Full")
      Cijk = basis->computeTripleProductTensor();
    else
      Cijk = basis->computeLinearTripleProductTensor();
    
    sgParams.set("Triple Product Tensor", Cijk);
  }

  // Create expansion
  std::string exp_type = expParams.get("Type", "Algebraic");
  if (exp_type == "Algebraic")
    expansion = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type,value_type>(basis, Cijk, Teuchos::rcp(&expParams,false)));

  else if (exp_type == "Quadrature") {
    Teuchos::ParameterList& quadParams = sgParams.sublist("Quadrature");
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > quad;
    if (quadParams.template isType<Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > >("Stochastic Galerkin Quadrature"))
      quad = quadParams.template get<Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > >("Stochastic Galerkin Quadrature");
    else {
      quad = 
	Stokhos::QuadratureFactory<ordinal_type,value_type>::create(sgParams);
      quadParams.set("Stochastic Galerkin Quadrature", quad);
    }
    expansion = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<ordinal_type,value_type>(basis, Cijk, quad, Teuchos::rcp(&expParams,false)));
 }

  else if (exp_type == "For UQTK") {
#ifdef HAVE_STOKHOS_FORUQTK
    typename Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type,value_type>::EXPANSION_METHOD method = 
      expParams.get("ForUQTK Expansion Method", 
		    Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type,value_type>::TAYLOR);
    value_type rtol = expParams.get("ForUQTK Expansion Tolerance", 1e-12);
    expansion = 
      Teuchos::rcp(new Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type,value_type>(basis, Cijk, method, rtol));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Stokhos::ExpansionFactory::create():  " <<
		       "ForUQTK expansion requires ForUQTK!" << std::endl);
#endif
  }

  /*
  else if (exp_type == "Derivative") {
    Teuchos::RCP<const Stokhos::DerivBasis<ordinal_type,value_type> > deriv_basis = Teuchos::rcp_dynamic_cast<const Stokhos::DerivBasis<ordinal_type,value_type> >(basis, true);
    Teuchos::RCP<Teuchos::SerialDenseMatrix<ordinal_type,value_type> > Bij;
    if (sgParams.template isType<Teuchos::RCP<const Teuchos::SerialDenseMatrix<ordinal_type,value_type> > >("Derivative Double Product Tensor"))
      Bij = sgParams.template get<const Teuchos::RCP<Teuchos::SerialDenseMatrix<ordinal_type,value_type> > >("Derivative Double Product Tensor");
    else {
      Bij = deriv_basis->computeDerivDoubleProductTensor();
      sgParams.set("Derivative Double Product Tensor", Bij);
    }
    Teuchos::RCP<const Stokhos::Dense3Tensor<ordinal_type,value_type> > Dijk;
    if (sgParams.template isType<Teuchos::RCP<const Stokhos::Dense3Tensor<ordinal_type,value_type> > >("Derivative Triple Product Tensor"))
      Dijk = sgParams.template get<Teuchos::RCP<const Stokhos::Dense3Tensor<ordinal_type,value_type> > >("Derivative Triple Product Tensor");
    else {
      Dijk = deriv_basis->computeDerivTripleProductTensor(Bij, Cijk);
      sgParams.set("Derivative Triple Product Tensor", Dijk);
    }
    expansion = 
      Teuchos::rcp(new 
		   Stokhos::DerivOrthogPolyExpansion<ordinal_type,value_type>(
		     deriv_basis, Bij, Cijk, Dijk));
  }
  */

  else if (exp_type == "Pseudospectral") {
    typedef Stokhos::PseudoSpectralOperator<ordinal_type,value_type> psop_type;
    Teuchos::ParameterList& psopParams = 
      sgParams.sublist("Pseudospectral Operator");
    Teuchos::RCP<const psop_type> psop;
    if (psopParams.template isType<Teuchos::RCP<const psop_type> >(
	  "Stochastic Galerkin Pseudospectral Operator"))
      psop = psopParams.template get<Teuchos::RCP<const psop_type> >(
	"Stochastic Galerkin Pseudospectral Operator");
    else {
      psop = 
	Stokhos::PseudoSpectralOperatorFactory<ordinal_type,value_type>::create(sgParams);
      psopParams.set("Stochastic Galerkin Pseudospectral Operator", psop);
    }
    expansion = 
      Teuchos::rcp(new Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type,value_type>(basis, Cijk, psop, Teuchos::rcp(&expParams,false)));
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Stokhos::ExpansionFactory::create():  " <<
		       "Invalid expansion type  " << exp_type << std::endl);

  expParams.set("Stochastic Galerkin Expansion", expansion);
  return expansion;
}
