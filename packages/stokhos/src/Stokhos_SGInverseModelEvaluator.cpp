// $Id$ 
// $Source$ 
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

#include "Stokhos_SGInverseModelEvaluator.hpp"

#include "Teuchos_TestForException.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Epetra_Map.h"

Stokhos::SGInverseModelEvaluator::SGInverseModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::Array<int>& sg_p_index_,
  const Teuchos::Array<int>& sg_g_index_) 
  : me(me_),
    sg_p_index(sg_p_index_),
    sg_g_index(sg_g_index_),
    num_p_sg(sg_p_index.size()),
    num_g_sg(sg_g_index.size()),
    block_p(num_p_sg),
    block_g(num_g_sg),
    block_dgdp(num_g_sg)
{
  InArgs me_inargs = me->createInArgs();
  OutArgs me_outargs = me->createOutArgs();
  
  // Create parameter SG blocks
  for (int i=0; i<num_p_sg; i++) {
    Teuchos::RCP<const Epetra_Map> sg_p_map = me->get_p_map(sg_p_index[i]);
    block_p[i] = Teuchos::rcp(new Epetra_Vector(*sg_p_map));
  }

  // Create response SG blocks
  for (int i=0; i<num_g_sg; i++) {
    Teuchos::RCP<const Epetra_Map> sg_g_map = me->get_g_map(sg_g_index[i]);
    
    // Create g SG blocks
    block_g[i] = Teuchos::rcp(new Epetra_Vector(*sg_g_map));
    
    // Create dg/dp SG blocks
    block_dgdp[i].resize(me_inargs.Np());
    for (int j=0; j<me_inargs.Np(); j++) {
      Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(j);
      block_dgdp[i][j] = 
	Teuchos::rcp(new Epetra_MultiVector(*sg_g_map, p_map->NumMyElements()));
    }
  }
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
  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGInverseModelEvaluator::
get_g_map(int l) const
{
  return me->get_g_map(l);
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGInverseModelEvaluator::
get_p_names(int l) const
{
  return me->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector> 
Stokhos::SGInverseModelEvaluator::
get_p_init(int l) const
{
  return me->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::SGInverseModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(me_inargs.Np());
  inArgs.set_Np_sg(num_p_sg);  

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::SGInverseModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());

  outArgs.set_Np_Ng(me_outargs.Np(), me_outargs.Ng());
  for (int i=0; i<me_outargs.Ng(); i++)
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));

  outArgs.set_Np_Ng_sg(me_outargs.Np(), num_g_sg);
  for (int i=0; i<num_g_sg; i++)
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp_sg, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, sg_g_index[i], j));
  
  return outArgs;
}

void 
Stokhos::SGInverseModelEvaluator::evalModel(const InArgs& inArgs,
				     const OutArgs& outArgs) const
{
  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();

  // Pass parameters
  for (int i=0; i<me_inargs.Np(); i++)
    me_inargs.set_p(i, inArgs.get_p(i));

  // Pass SG parameters
  for (int i=0; i<num_p_sg; i++) {
    InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(i);
    if (p_sg != Teuchos::null) {
      p_sg->assignToBlockVector(*(block_p[i]));
      me_inargs.set_p(sg_p_index[i], block_p[i]);
    }
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // Responses
  for (int i=0; i<me_outargs.Ng(); i++) {
    // g
    me_outargs.set_g(i, outArgs.get_g(i));

    // dg/dp
    for (int j=0; j<outArgs.Np(); j++)
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	me_outargs.set_DgDp(i, j, outArgs.get_DgDp(i,j));
  }


  // SG Responses
  for (int i=0; i<num_g_sg; i++) {
    // g
    OutArgs::sg_vector_t g_sg = outArgs.get_g_sg(i);
    if (g_sg != Teuchos::null) {
      g_sg->assignToBlockVector(*(block_g[i]));
      me_outargs.set_g(sg_g_index[i], block_g[i]);
    }

    // dg/dp
    for (int j=0; j<me_outargs.Np(); j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_sg, i, j).none()) {
	SGDerivative dgdp_sg = outArgs.get_DgDp_sg(i,j);
	Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp_sg_mv
	  = dgdp_sg.getMultiVector();
	if (dgdp_sg_mv != Teuchos::null) {
	  dgdp_sg_mv->assignToBlockMultiVector(*(block_dgdp[i][j]));
	  me_outargs.set_DgDp(sg_g_index[i], j, 
			      Derivative(block_dgdp[i][j],
					 dgdp_sg.getMultiVectorOrientation()));
	}
	TEST_FOR_EXCEPTION(dgdp_sg.getLinearOp() != Teuchos::null, 
			   std::logic_error,
			   "Error!  Stokhos::SGInverseModelEvaluator::evalModel"
			   << " cannot handle operator form of dg/dp!");
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Copy block SG components
  for (int i=0; i<num_g_sg; i++) {
    // g
    OutArgs::sg_vector_t g_sg = outArgs.get_g_sg(i);
    if (g_sg != Teuchos::null) {
      g_sg->assignFromBlockVector(*(block_g[i]));
    }

    // dg/dp
    for (int j=0; j<me_outargs.Np(); j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_sg, i, j).none()) {
	SGDerivative dgdp_sg = outArgs.get_DgDp_sg(i,j);
	Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp_sg_mv
	  = dgdp_sg.getMultiVector();
	if (dgdp_sg_mv != Teuchos::null)
	  dgdp_sg_mv->assignFromBlockMultiVector(*(block_dgdp[i][j]));
      }
    }
  }

}
