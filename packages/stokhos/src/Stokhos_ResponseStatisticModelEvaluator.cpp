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

#include "Stokhos_ResponseStatisticModelEvaluator.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "Teuchos_TestForException.hpp"

Stokhos::ResponseStatisticModelEvaluator::ResponseStatisticModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::Array< Teuchos::RCP<const Epetra_Map> >& base_g_maps_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Epetra_BlockMap>& block_map_) : 
  me(me_),
  base_g_maps(base_g_maps_),
  sg_basis(sg_basis_),
  sg_comm(sg_comm_),
  block_map(block_map_),
  num_p(0),
  num_g(0)
{
  InArgs me_inargs = me->createInArgs();
  OutArgs me_outargs = me->createOutArgs();
  num_p = me_inargs.Np();
  num_g = me_outargs.Ng();
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::ResponseStatisticModelEvaluator::
get_x_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::ResponseStatisticModelEvaluator::
get_f_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::ResponseStatisticModelEvaluator::
get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << 
    "Error!  Stokhos::ResponseStatisticModelEvaluator::get_p_map():  " << 
    "Invalid parameter index l = " << l << std::endl);

  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::ResponseStatisticModelEvaluator::
get_g_map(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= 3*num_g || l < 0, std::logic_error,
    std::endl << 
    "Error!  Stokhos::ResponseStatisticModelEvaluator::get_g_map():  " << 
    "Invalid response index l = " << l << std::endl);

  if (l < num_g)
    return me->get_g_map(l);
  else if (l < 2*num_g)
    return base_g_maps[l-num_g];
  else
    return base_g_maps[l-2*num_g];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::ResponseStatisticModelEvaluator::
get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << 
    "Error!  Stokhos::ResponseStatisticModelEvaluator::get_p_names():  " << 
    "Invalid parameter index l = " << l << std::endl);

  return me->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector> 
Stokhos::ResponseStatisticModelEvaluator::
get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << 
    "Error!  Stokhos::ResponseStatisticModelEvaluator::get_p_init():  " << 
    "Invalid parameter index l = " << l << std::endl);

  return me->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::ResponseStatisticModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);

  // Support pass-through of sg data
  inArgs.setSupports(IN_ARG_sg_basis, me_inargs.supports(IN_ARG_sg_basis));
  inArgs.setSupports(IN_ARG_sg_quadrature, 
		     me_inargs.supports(IN_ARG_sg_quadrature));
  inArgs.setSupports(IN_ARG_sg_expansion, 
		     me_inargs.supports(IN_ARG_sg_expansion));

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::ResponseStatisticModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());

  outArgs.set_Np_Ng(num_p, 3*num_g);
  for (int i=0; i<num_g; i++) {
    for (int j=0; j<num_p; j++) {
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
      outArgs.setSupports(OUT_ARG_DgDp, i+num_g, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
      outArgs.setSupports(OUT_ARG_DgDp, i+2*num_g, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
    }
  }
  
  return outArgs;
}

void 
Stokhos::ResponseStatisticModelEvaluator::evalModel(const InArgs& inArgs,
					    const OutArgs& outArgs) const
{
  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));

  // Pass SG data
  if (me_inargs.supports(IN_ARG_sg_basis))
    me_inargs.set_sg_basis(inArgs.get_sg_basis());
  if (me_inargs.supports(IN_ARG_sg_quadrature))
    me_inargs.set_sg_quadrature(inArgs.get_sg_quadrature());
  if (me_inargs.supports(IN_ARG_sg_expansion))
    me_inargs.set_sg_expansion(inArgs.get_sg_expansion());

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  Teuchos::Array<bool> need_g_var(num_g);

  // Responses
  for (int i=0; i<num_g; i++) {

    // See if we need g_var for d/dp(g_var) below
    need_g_var[i] = false;
    for (int j=0; j<num_p; j++)
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	if (!outArgs.get_DgDp(i+2*num_g,j).isEmpty())
	  need_g_var[i] = true;
    
    // g
    Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(i);
    Teuchos::RCP<Epetra_Vector> g_mean = outArgs.get_g(i+num_g);
    Teuchos::RCP<Epetra_Vector> g_var = outArgs.get_g(i+2*num_g);
    if ((g_mean != Teuchos::null || g_var != Teuchos::null || need_g_var[i]) && 
	g == Teuchos::null) {
      g = Teuchos::rcp(new Epetra_Vector(*(me->get_g_map(i))));
    }
    me_outargs.set_g(i, g);

    // dg/dp
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	Derivative dgdp = outArgs.get_DgDp(i,j);
	Teuchos::RCP<Epetra_MultiVector> dgdp_mean = 
	  outArgs.get_DgDp(i+num_g,j).getMultiVector();
	Teuchos::RCP<Epetra_MultiVector> dgdp_var = 
	  outArgs.get_DgDp(i+2*num_g,j).getMultiVector();
	if ((dgdp_mean != Teuchos::null || dgdp_var != Teuchos::null) && 
	    dgdp.isEmpty()) {
	  Teuchos::RCP<const Epetra_Map> g_map = me->get_g_map(i);
	  Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(j);
	  DerivativeSupport ds = me_outargs.supports(OUT_ARG_DgDp,i,j);
	  if (ds.supports(DERIV_LINEAR_OP))
	    dgdp = Derivative(me->create_DgDp_op(i,j));
	  else if (ds.supports(DERIV_MV_BY_COL))
	    dgdp = 
	      Derivative(Teuchos::rcp(new Epetra_MultiVector(
					*g_map, 
					p_map->NumMyElements())),
			 DERIV_MV_BY_COL);
	  else
	    dgdp = 
	      Derivative(Teuchos::rcp(new Epetra_MultiVector(
					*p_map, 
					g_map->NumMyElements())),
			 DERIV_TRANS_MV_BY_ROW);
	}
	me_outargs.set_DgDp(i, j, dgdp);

	// We need g to compute d/dp(g_var)
	if (dgdp_var != Teuchos::null && g == Teuchos::null) {
	  g = Teuchos::rcp(new Epetra_Vector(*(me->get_g_map(i))));
	  me_outargs.set_g(i, g);
	}
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Compute statistics
  for (int i=0; i<num_g; i++) {

    // g
    Teuchos::RCP<Epetra_Vector> g = me_outargs.get_g(i);
    Teuchos::RCP<Epetra_Vector> g_mean = outArgs.get_g(i+num_g);
    Teuchos::RCP<Epetra_Vector> g_var = outArgs.get_g(i+2*num_g);
    if (need_g_var[i] && g_var == Teuchos::null)
      g_var = Teuchos::rcp(new Epetra_Vector(*(base_g_maps[i])));
    Teuchos::RCP<EpetraVectorOrthogPoly> g_sg;
    if (g_mean != Teuchos::null || g_var != Teuchos::null) {
      g_sg = Teuchos::rcp(new EpetraVectorOrthogPoly(
		   sg_basis, block_map, base_g_maps[i], me->get_g_map(i), 
		   sg_comm, View, *g));
    }

    if (g_mean != Teuchos::null)
      g_sg->computeMean(*g_mean);
    if (g_var != Teuchos::null)
      g_sg->computeVariance(*g_var);

    // dg/dp
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	Derivative dgdp = me_outargs.get_DgDp(i,j);
	Teuchos::RCP<Epetra_MultiVector> dgdp_mean = 
	  outArgs.get_DgDp(i+num_g,j).getMultiVector();
	Teuchos::RCP<Epetra_MultiVector> dgdp_var = 
	  outArgs.get_DgDp(i+2*num_g,j).getMultiVector();
	if (dgdp_mean != Teuchos::null || dgdp_var != Teuchos::null) {
	  if (dgdp.getLinearOp() != Teuchos::null) {
	    Teuchos::RCP<Epetra_Operator> dgdp_op = dgdp.getLinearOp();
	    int n = base_g_maps[i]->NumMyElements();
	    EpetraMultiVectorOrthogPoly X_sg(sg_basis, block_map, 
					     base_g_maps[i], sg_comm, n);
	    dgdp_op->SetUseTranspose(true);
	    if (dgdp_mean != Teuchos::null) {
	      X_sg.init(0.0);
	      for (int l=0; l<n; l++)
		X_sg[0][l][l] = (*g_sg)[0][l];
	      TEST_FOR_EXCEPTION(
		outArgs.get_DgDp(i+num_g,j).getMultiVectorOrientation() == DERIV_MV_BY_COL, 
		std::logic_error,
		"Error!  ResponseStatisticModelEvaluator does not support " <<
		" dg/dp DERIV_MV_BY_COL for underlying operator dg/dp");
	      dgdp_op->Apply(*(X_sg.getBlockMultiVector()), *dgdp_mean);
	    }
	    if (dgdp_var != Teuchos::null) {
	      X_sg.init(0.0);
	      for (int k=1; k<sg_basis->size(); k++)
		for (int l=0; l<n; l++)
		  X_sg[k][l][l] = 2.0*(*g_sg)[k][l]*sg_basis->norm_squared(k);
	      TEST_FOR_EXCEPTION(
		outArgs.get_DgDp(i+2*num_g,j).getMultiVectorOrientation() == DERIV_MV_BY_COL, 
		std::logic_error,
		"Error!  ResponseStatisticModelEvaluator does not support " <<
		" dg/dp DERIV_MV_BY_COL for underlying operator dg/dp");
	      dgdp_op->Apply(*(X_sg.getBlockMultiVector()), *dgdp_var);
	    }
	    
	  }
	  else if (dgdp.getMultiVector() != Teuchos::null) {
	    Teuchos::RCP<const Epetra_MultiVector> dgdp_mv = 
	      dgdp.getMultiVector();
	    Teuchos::RCP<const Epetra_Map> row_map;
	    if (dgdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	      row_map = me->get_g_map(i);
	    else
	      row_map = me->get_p_map(j);
	    EpetraMultiVectorOrthogPoly dgdp_sg(
	      sg_basis, block_map, row_map, Teuchos::rcpFromRef(dgdp_mv->Map()),
	      sg_comm, View, *dgdp_mv);
	    if (dgdp_mean != Teuchos::null)
	      dgdp_sg.computeMean(*dgdp_mean);
	    if (dgdp_var != Teuchos::null) {
	      dgdp_var->PutScalar(0.0);
	      for (int k=1; k<sg_basis->size(); k++)
		for (int m=0; m<dgdp_var->NumVectors(); m++)
		  for (int l=0; l<row_map->NumMyElements(); l++)
		    (*dgdp_var)[m][l] += 2.0*(*g_sg)[k][l]*dgdp_sg[k][m][l]*
		      sg_basis->norm_squared(k);
	    }
	  }
	}
      }
    }
  }
}
