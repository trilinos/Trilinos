// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_MPInverseModelEvaluator.hpp"

#include "Teuchos_Assert.hpp"
#include "Stokhos_ProductEpetraVector.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "Stokhos_ProductEpetraOperator.hpp"
#include "Epetra_Map.h"

Stokhos::MPInverseModelEvaluator::MPInverseModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::Array<int>& mp_p_index_map_,
  const Teuchos::Array<int>& mp_g_index_map_,
  const Teuchos::Array< Teuchos::RCP<const Epetra_Map> >& base_g_maps_) : 
  me(me_),
  mp_p_index_map(mp_p_index_map_),
  mp_g_index_map(mp_g_index_map_),
  base_g_maps(base_g_maps_),
  num_p(0),
  num_g(0),
  num_p_mp(mp_p_index_map.size()),
  num_g_mp(mp_g_index_map.size())
{
  InArgs me_inargs = me->createInArgs();
  OutArgs me_outargs = me->createOutArgs();
  num_p = me_inargs.Np() - num_p_mp;
  num_g = me_outargs.Ng();

  TEUCHOS_TEST_FOR_EXCEPTION(
    base_g_maps.size() != num_g_mp, std::logic_error,
    std::endl 
    << "Error!  Stokhos::MPInverseModelEvaluator::MPInverseModelEvaluator():"
    << "  Base response map array size does not match size of index array!");
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_x_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_f_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_map():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_g_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_g || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_g_map():"
    << "  Invalid response index l = " << l << std::endl);
  return base_g_maps[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::MPInverseModelEvaluator::
get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_names():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector> 
Stokhos::MPInverseModelEvaluator::
get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_init():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::MPInverseModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  for (int i=0; i<num_p_mp; i++)
    inArgs.setSupports(IN_ARG_p_mp, mp_p_index_map[i], true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::MPInverseModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());

  outArgs.set_Np_Ng(num_p, num_g);
  for (int i=0; i<num_g_mp; i++) {
    outArgs.setSupports(OUT_ARG_g_mp, mp_g_index_map[i], true);
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp_mp, mp_g_index_map[i], j, 
			  me_outargs.supports(OUT_ARG_DgDp,i,j));
  }
  
  return outArgs;
}

void 
Stokhos::MPInverseModelEvaluator::evalModel(const InArgs& inArgs,
					    const OutArgs& outArgs) const
{
  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));

  // Pass MP parameters
  for (int i=0; i<num_p_mp; i++) {
    mp_const_vector_t p_mp = inArgs.get_p_mp(mp_p_index_map[i]);
    if (p_mp != Teuchos::null) {
      me_inargs.set_p(i+num_p, p_mp->getBlockVector());
    }
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();


  // MP Responses
  for (int i=0; i<num_g_mp; i++) {
    int ii = mp_g_index_map[i];

    // g
    mp_vector_t g_mp = outArgs.get_g_mp(ii);
    if (g_mp != Teuchos::null) {
      me_outargs.set_g(i, Teuchos::rcp_dynamic_cast<Epetra_Vector>(g_mp->getBlockVector()));
    }

    // dg/dp
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_mp, ii, j).none()) {
	MPDerivative dgdp_mp = outArgs.get_DgDp_mp(ii,j);
	Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dgdp_mp_mv = 
	  dgdp_mp.getMultiVector();
	Teuchos::RCP<Epetra_Operator> dgdp_mp_op =
	  dgdp_mp.getLinearOp();
	if (dgdp_mp_mv != Teuchos::null) {
	  me_outargs.set_DgDp(
	    i, j, Derivative(dgdp_mp_mv->getBlockMultiVector(),
			     dgdp_mp.getMultiVectorOrientation()));
	}
	else if (dgdp_mp_op != Teuchos::null) {
	  me_outargs.set_DgDp(i, j, Derivative(dgdp_mp_op));
	}
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

}
