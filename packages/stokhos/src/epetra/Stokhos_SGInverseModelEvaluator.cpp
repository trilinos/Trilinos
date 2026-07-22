// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_SGInverseModelEvaluator.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "Epetra_Map.h"
#include "Teuchos_Assert.hpp"

Stokhos::SGInverseModelEvaluator::SGInverseModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::Array<int>& sg_p_index_map_,
  const Teuchos::Array<int>& sg_g_index_map_,
  const Teuchos::Array< Teuchos::RCP<const Epetra_Map> >& base_g_maps_) : 
  me(me_),
  sg_p_index_map(sg_p_index_map_),
  sg_g_index_map(sg_g_index_map_),
  base_g_maps(base_g_maps_),
  num_p(0),
  num_g(0),
  num_p_sg(sg_p_index_map.size()),
  num_g_sg(sg_g_index_map.size())
{
  InArgs me_inargs = me->createInArgs();
  OutArgs me_outargs = me->createOutArgs();
  num_p = me_inargs.Np() - num_p_sg;
  num_g = me_outargs.Ng();

  TEUCHOS_TEST_FOR_EXCEPTION(
    base_g_maps.size() != num_g_sg, std::logic_error,
    std::endl 
    << "Error!  Stokhos::SGInverseModelEvaluator::SGInverseModelEvaluator():"
    << "  Base response map array size does not match size of index array!");
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::SGInverseModelEvaluator::
get_x_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGInverseModelEvaluator::
get_f_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGInverseModelEvaluator::
get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::SGInverseModelEvaluator::get_p_map():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGInverseModelEvaluator::
get_g_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_g || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::SGInverseModelEvaluator::get_g_map():"
    << "  Invalid response index l = " << l << std::endl);
  return base_g_maps[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGInverseModelEvaluator::
get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::SGInverseModelEvaluator::get_p_names():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector> 
Stokhos::SGInverseModelEvaluator::
get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::SGInverseModelEvaluator::get_p_init():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::SGInverseModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  for (int i=0; i<num_p_sg; i++)
    inArgs.setSupports(IN_ARG_p_sg, sg_p_index_map[i], true);
  inArgs.setSupports(IN_ARG_sg_basis, me_inargs.supports(IN_ARG_sg_basis));
  inArgs.setSupports(IN_ARG_sg_quadrature, 
		     me_inargs.supports(IN_ARG_sg_quadrature));
  inArgs.setSupports(IN_ARG_sg_expansion, 
		     me_inargs.supports(IN_ARG_sg_expansion));

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::SGInverseModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());

  outArgs.set_Np_Ng(num_p, num_g);
  for (int i=0; i<num_g_sg; i++) {
    outArgs.setSupports(OUT_ARG_g_sg, sg_g_index_map[i], true);
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp_sg, sg_g_index_map[i], j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
  }
  
  return outArgs;
}

void 
Stokhos::SGInverseModelEvaluator::evalModel(const InArgs& inArgs,
					    const OutArgs& outArgs) const
{
  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));

  // Pass SG parameters
  for (int i=0; i<num_p_sg; i++) {
    InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(sg_p_index_map[i]);
    if (p_sg != Teuchos::null) {
      me_inargs.set_p(i+num_p, p_sg->getBlockVector());
    }
  }

  // Pass SG data
  if (me_inargs.supports(IN_ARG_sg_basis))
    me_inargs.set_sg_basis(inArgs.get_sg_basis());
  if (me_inargs.supports(IN_ARG_sg_quadrature))
    me_inargs.set_sg_quadrature(inArgs.get_sg_quadrature());
  if (me_inargs.supports(IN_ARG_sg_expansion))
    me_inargs.set_sg_expansion(inArgs.get_sg_expansion());

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // SG Responses
  for (int i=0; i<num_g_sg; i++) {
    int ii = sg_g_index_map[i];
    
    // g
    OutArgs::sg_vector_t g_sg = outArgs.get_g_sg(ii);
    if (g_sg != Teuchos::null) {
      me_outargs.set_g(i, Teuchos::rcp_dynamic_cast<Epetra_Vector>(g_sg->getBlockVector()));
    }

    // dg/dp
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_sg, ii, j).none()) {
	SGDerivative dgdp_sg = outArgs.get_DgDp_sg(ii,j);
	Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp_sg_mv = 
	  dgdp_sg.getMultiVector();
	Teuchos::RCP<Epetra_Operator> dgdp_sg_op =
	  dgdp_sg.getLinearOp();
	if (dgdp_sg_mv != Teuchos::null) {
	  me_outargs.set_DgDp(
	    i, j, Derivative(dgdp_sg_mv->getBlockMultiVector(),
			     dgdp_sg.getMultiVectorOrientation()));
	}
	else if (dgdp_sg_op != Teuchos::null) {
	  me_outargs.set_DgDp(i, j, Derivative(dgdp_sg_op));
	}
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

}
