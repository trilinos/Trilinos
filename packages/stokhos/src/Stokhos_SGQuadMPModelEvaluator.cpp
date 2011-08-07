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

#include "Stokhos_SGQuadMPModelEvaluator.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Stokhos_ProductEpetraVector.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TestForException.hpp"

Stokhos::SGQuadMPModelEvaluator::
SGQuadMPModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm_,
  const Teuchos::RCP<const Epetra_Map>& mp_block_map_) : 
  me(me_),
  mp_comm(mp_comm_),
  mp_block_map(mp_block_map_),
  num_p(0),
  num_g(0),
  x_dot_mp(),
  x_mp(),
  p_mp(),
  f_mp(),
  W_mp(),
  dfdp_mp(),
  g_mp(),
  dgdx_mp(),
  dgdx_dot_mp(),
  dgdp_mp()
{
  int num_qp = mp_block_map->NumMyElements();

  // Create storage for x_dot, x, and p at a quad point
  InArgs me_inargs = me->createInArgs();
  num_p = me_inargs.Np();
  num_p_mp = me_inargs.Np_mp();
  if (me_inargs.supports(IN_ARG_x_dot))
    x_dot_mp = Teuchos::rcp(new ProductEpetraVector(
			      mp_block_map, me->get_x_map(), mp_comm));
  if (me_inargs.supports(IN_ARG_x))
    x_mp = Teuchos::rcp(new ProductEpetraVector(
			mp_block_map, me->get_x_map(), mp_comm));
  p_mp.resize(num_p_mp);
  for (int i=0; i<num_p_mp; i++)
    p_mp[i] = Teuchos::rcp(new ProductEpetraVector(
			     mp_block_map, me->get_p_mp_map(i), mp_comm));

  // Create storage for f and W at a quad point
  OutArgs me_outargs = me->createOutArgs();
  num_g = me_outargs.Ng();
  num_g_mp = me_outargs.Ng_mp();

  // f
  if (me_outargs.supports(OUT_ARG_f))
    f_mp = Teuchos::rcp(new ProductEpetraVector(
			  mp_block_map, me->get_f_map(), mp_comm));

  // W
  if (me_outargs.supports(OUT_ARG_W)) {
    W_mp = Teuchos::rcp(new ProductContainer<Epetra_Operator>(mp_block_map));
    for (int i=0; i<num_qp; i++)
      W_mp->setCoeffPtr(i, me->create_W());
  }

  // df/dp
  dfdp_mp.resize(num_p);
  for (int i=0; i<num_p; i++) {
    if (me_outargs.supports(OUT_ARG_DfDp,i).supports(DERIV_MV_BY_COL))
      dfdp_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(
	Teuchos::rcp(new ProductEpetraMultiVector(
		       mp_block_map, me->get_f_map(), mp_comm,
		       me->get_p_map(i)->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DfDp,i).supports(DERIV_TRANS_MV_BY_ROW))
      dfdp_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(
	Teuchos::rcp(new ProductEpetraMultiVector(
		       mp_block_map, me->get_p_map(i), mp_comm,
		       me->get_f_map()->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DfDp,i).supports(DERIV_LINEAR_OP)) {
      Teuchos::RCP< ProductContainer<Epetra_Operator> > dfdp_mp_op =
	Teuchos::rcp(new ProductContainer<Epetra_Operator>(mp_block_map));
      for (int j=0; j<num_qp; j++)
	dfdp_mp_op->setCoeffPtr(j, me->create_DfDp_op(i));
      dfdp_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(dfdp_mp_op);
    }

  }

  g_mp.resize(num_g_mp);
  dgdx_mp.resize(num_g_mp);
  dgdx_dot_mp.resize(num_g_mp);
  dgdp_mp.resize(num_g_mp);
  for (int i=0; i<num_g_mp; i++) {

    // g
    g_mp[i] = 
      Teuchos::rcp(new ProductEpetraVector(
		     mp_block_map, me->get_g_mp_map(i), mp_comm));

    // dg/dx
    if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_TRANS_MV_BY_ROW))
      dgdx_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(
	Teuchos::rcp(new ProductEpetraMultiVector(
		       mp_block_map, me->get_x_map(), mp_comm,
		       me->get_g_mp_map(i)->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_MV_BY_COL))
      dgdx_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(
	Teuchos::rcp(new ProductEpetraMultiVector(
		       mp_block_map, me->get_g_mp_map(i), mp_comm,
		       me->get_x_map()->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP)) {
      Teuchos::RCP< ProductContainer<Epetra_Operator> > dgdx_mp_op =
	Teuchos::rcp(new ProductContainer<Epetra_Operator>(mp_block_map));
      for (int j=0; j<num_qp; j++)
	dgdx_mp_op->setCoeffPtr(j, me->create_DgDx_op(i));
      dgdx_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(dgdx_mp_op);
    }
    
    // dg/dx_dot
    if (me_outargs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_TRANS_MV_BY_ROW))
      dgdx_dot_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(
	Teuchos::rcp(new ProductEpetraMultiVector(
		       mp_block_map, me->get_x_map(), mp_comm,
		       me->get_g_mp_map(i)->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_MV_BY_COL))
      dgdx_dot_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(
	Teuchos::rcp(new ProductEpetraMultiVector(
		       mp_block_map, me->get_g_mp_map(i), mp_comm,
		       me->get_x_map()->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_LINEAR_OP)) {
      Teuchos::RCP< ProductContainer<Epetra_Operator> > dgdx_dot_mp_op =
	Teuchos::rcp(new ProductContainer<Epetra_Operator>(mp_block_map));
      for (int j=0; j<num_qp; j++)
	dgdx_dot_mp_op->setCoeffPtr(j, me->create_DgDx_dot_op(i));
      dgdx_dot_mp[i] = EpetraExt::ModelEvaluator::MPDerivative(dgdx_dot_mp_op);
    }

    // dg/dp
    dgdp_mp[i].resize(num_p);
    for (int j=0; j<num_p; j++) {
      if (me_outargs.supports(OUT_ARG_DgDp, i, j).supports(DERIV_TRANS_MV_BY_ROW))
	dgdp_mp[i][j] = EpetraExt::ModelEvaluator::MPDerivative(
	  Teuchos::rcp(new ProductEpetraMultiVector(
			 mp_block_map, me->get_p_map(j), mp_comm,
			 me->get_g_mp_map(i)->NumGlobalElements())));
      else if (me_outargs.supports(OUT_ARG_DgDp, i, j).supports(DERIV_MV_BY_COL))
	dgdp_mp[i][j] = EpetraExt::ModelEvaluator::MPDerivative(
	  Teuchos::rcp(new ProductEpetraMultiVector(
			 mp_block_map, me->get_g_mp_map(i), mp_comm,
			 me->get_p_map(j)->NumGlobalElements())));
      else if (me_outargs.supports(OUT_ARG_DgDp, i, j).supports(DERIV_LINEAR_OP)) {
	Teuchos::RCP< ProductContainer<Epetra_Operator> > dgdp_mp_op =
	  Teuchos::rcp(new ProductContainer<Epetra_Operator>(mp_block_map));
	for (int k=0; k<num_qp; k++)
	  dgdp_mp_op->setCoeffPtr(k, me->create_DgDp_op(i,j));
	dgdp_mp[i][j] = EpetraExt::ModelEvaluator::MPDerivative(dgdp_mp_op);
      }

    }

  }
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadMPModelEvaluator::
get_x_map() const
{
  return me->get_x_map();
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadMPModelEvaluator::
get_f_map() const
{
  return me->get_f_map();
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadMPModelEvaluator::
get_p_map(int l) const
{
  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadMPModelEvaluator::
get_p_sg_map(int l) const
{
  return me->get_p_mp_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadMPModelEvaluator::
get_g_map(int l) const
{
  return me->get_g_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadMPModelEvaluator::
get_g_sg_map(int l) const
{
  return me->get_g_mp_map(l);
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGQuadMPModelEvaluator::
get_p_names(int l) const
{
  return me->get_p_names(l);
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGQuadMPModelEvaluator::
get_p_sg_names(int l) const
{
  return me->get_p_mp_names(l);
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGQuadMPModelEvaluator::
get_x_init() const
{
  return me->get_x_init();
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGQuadMPModelEvaluator::
get_p_init(int l) const
{
  return me->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGQuadMPModelEvaluator::
create_W() const
{
  return me->create_W();
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::SGQuadMPModelEvaluator::
createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(me_inargs.Np()); 
  inArgs.setSupports(IN_ARG_x_dot, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_t, me_inargs.supports(IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha, me_inargs.supports(IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta, me_inargs.supports(IN_ARG_beta));

  inArgs.set_Np_sg(num_p_mp);
  inArgs.setSupports(IN_ARG_x_sg, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_x_dot_sg, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_sg_basis, true);
  inArgs.setSupports(IN_ARG_sg_quadrature, true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::SGQuadMPModelEvaluator::
createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_outargs.Np(), me_outargs.Ng());
  outArgs.setSupports(OUT_ARG_f, me_outargs.supports(OUT_ARG_f));
  outArgs.setSupports(OUT_ARG_W, me_outargs.supports(OUT_ARG_W));
  for (int j=0; j<me_outargs.Np(); j++)
    outArgs.setSupports(OUT_ARG_DfDp, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  for (int i=0; i<me_outargs.Ng(); i++) {
    outArgs.setSupports(OUT_ARG_DgDx, i, 
			me_outargs.supports(OUT_ARG_DgDx, i));
    outArgs.setSupports(OUT_ARG_DgDx_dot, i, 
			me_outargs.supports(OUT_ARG_DgDx_dot, i));
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
  }

  outArgs.setSupports(OUT_ARG_f_sg, me_outargs.supports(OUT_ARG_f));
  if (me_outargs.supports(OUT_ARG_W)) {
    outArgs.set_W_properties(me_outargs.get_W_properties());
    outArgs.setSupports(OUT_ARG_W_sg, true);
  }
  outArgs.set_Np_Ng_sg(num_p, num_g_mp);
  for (int j=0; j<me_outargs.Np(); j++)
    outArgs.setSupports(OUT_ARG_DfDp_sg, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  for (int i=0; i<num_g_mp; i++) {
    outArgs.setSupports(OUT_ARG_DgDx_sg, i, 
			me_outargs.supports(OUT_ARG_DgDx, i));
    outArgs.setSupports(OUT_ARG_DgDx_dot_sg, i, 
			me_outargs.supports(OUT_ARG_DgDx_dot, i));
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp_sg, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
  }
  
  return outArgs;
}

void 
Stokhos::SGQuadMPModelEvaluator::
evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();
  if (me_inargs.supports(IN_ARG_x))
    me_inargs.set_x(inArgs.get_x());
  if (me_inargs.supports(IN_ARG_x_dot))
    me_inargs.set_x_dot(inArgs.get_x_dot());
  if (me_inargs.supports(IN_ARG_alpha))
    me_inargs.set_alpha(inArgs.get_alpha());
  if (me_inargs.supports(IN_ARG_beta))
    me_inargs.set_beta(inArgs.get_beta());
  if (me_inargs.supports(IN_ARG_t))
    me_inargs.set_t(inArgs.get_t());
  for (int i=0; i<inArgs.Np(); i++)
    me_inargs.set_p(i, inArgs.get_p(i));

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();
  if (me_outargs.supports(OUT_ARG_f))
    me_outargs.set_f(outArgs.get_f());
  if (me_outargs.supports(OUT_ARG_W))
    me_outargs.set_W(outArgs.get_W());
  for (int j=0; j<outArgs.Np(); j++)
    if (!outArgs.supports(OUT_ARG_DfDp, j).none())
      me_outargs.set_DfDp(j, outArgs.get_DfDp(j));
  for (int i=0; i<outArgs.Ng(); i++) {
    me_outargs.set_g(i, outArgs.get_g(i));
    if (!outArgs.supports(OUT_ARG_DgDx, i).none())
	me_outargs.set_DgDx(i, outArgs.get_DgDx(i));
    if (!outArgs.supports(OUT_ARG_DgDx_dot, i).none())
	me_outargs.set_DgDx(i, outArgs.get_DgDx_dot(i));
    for (int j=0; j<outArgs.Np(); j++)
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	me_outargs.set_DgDp(i, j, outArgs.get_DgDp(i,j));
  }

  bool do_quad = false;
  InArgs::sg_const_vector_t x_sg;
  InArgs::sg_const_vector_t x_dot_sg;
  Teuchos::Array<InArgs::sg_const_vector_t> p_sg(inArgs.Np_sg());
  OutArgs::sg_vector_t f_sg;
  OutArgs::sg_operator_t W_sg;
  Teuchos::Array<SGDerivative> dfdp_sg(outArgs.Np());
  Teuchos::Array<OutArgs::sg_vector_t> g_sg(outArgs.Ng_sg());
  Teuchos::Array<SGDerivative> dgdx_sg(outArgs.Ng_sg());
  Teuchos::Array<SGDerivative> dgdx_dot_sg(outArgs.Ng_sg());
  Teuchos::Array< Teuchos::Array<SGDerivative> > dgdp_sg(outArgs.Ng_sg());
  TEST_FOR_EXCEPTION(inArgs.get_sg_basis() == Teuchos::null, 
		     std::logic_error,
		     "Error!  Stokhos::SGQuadModelEvaluator::evalModel():  " <<
		     "SG basis inArg cannot be null!");
  TEST_FOR_EXCEPTION(inArgs.get_sg_quadrature() == Teuchos::null, 
		     std::logic_error,
		     "Error!  Stokhos::SGQuadModelEvaluator::evalModel():  " <<
		     "SG quadrature inArg cannot be null!");
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    inArgs.get_sg_basis();
  Teuchos::RCP< const Stokhos::Quadrature<int,double> > quad = 
    inArgs.get_sg_quadrature();
  
  if (inArgs.supports(IN_ARG_x_sg)) {
    x_sg = inArgs.get_x_sg();
    if (x_sg != Teuchos::null) {
      do_quad = true;
    }
  }
  if (inArgs.supports(IN_ARG_x_dot_sg)) {
    x_dot_sg = inArgs.get_x_dot_sg();
    if (x_dot_sg != Teuchos::null) {
      do_quad = true;
    }
  }
  for (int i=0; i<inArgs.Np_sg(); i++) {
    p_sg[i] = inArgs.get_p_sg(i);
    if (p_sg[i] != Teuchos::null) {
      do_quad = true;
    }
  }
  if (outArgs.supports(OUT_ARG_f_sg)) {
    f_sg = outArgs.get_f_sg();
    if (f_sg != Teuchos::null)
      f_sg->init(0.0);
  }
  if (outArgs.supports(OUT_ARG_W_sg)) {
    W_sg = outArgs.get_W_sg();
    if (W_sg != Teuchos::null)
      W_sg->init(0.0);
  }
  for (int i=0; i<inArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp_sg, i).none()) {
      dfdp_sg[i] = outArgs.get_DfDp_sg(i);
      if (dfdp_sg[i].getMultiVector() != Teuchos::null)
	dfdp_sg[i].getMultiVector()->init(0.0);
      else if (dfdp_sg[i].getLinearOp() != Teuchos::null)
	dfdp_sg[i].getLinearOp()->init(0.0);
    }
  }
      
  for (int i=0; i<outArgs.Ng_sg(); i++) {
    g_sg[i] = outArgs.get_g_sg(i);
    if (g_sg[i] != Teuchos::null)
      g_sg[i]->init(0.0);
    
    if (!outArgs.supports(OUT_ARG_DgDx_sg, i).none()) {
      dgdx_sg[i] = outArgs.get_DgDx_sg(i);
      if (dgdx_sg[i].getMultiVector() != Teuchos::null)
	dgdx_sg[i].getMultiVector()->init(0.0);
      else if (dgdx_sg[i].getLinearOp() != Teuchos::null)
	dgdx_sg[i].getLinearOp()->init(0.0);
    }

    if (!outArgs.supports(OUT_ARG_DgDx_dot_sg, i).none()) {
      dgdx_dot_sg[i] = outArgs.get_DgDx_dot_sg(i);
      if (dgdx_dot_sg[i].getMultiVector() != Teuchos::null)
	dgdx_dot_sg[i].getMultiVector()->init(0.0);
      else if (dgdx_dot_sg[i].getLinearOp() != Teuchos::null)
	dgdx_dot_sg[i].getLinearOp()->init(0.0);
    }

    dgdp_sg[i].resize(outArgs.Np());
    for (int j=0; j<outArgs.Np(); j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_sg, i, j).none()) {
	dgdp_sg[i][j] = outArgs.get_DgDp_sg(i,j);
	if (dgdp_sg[i][j].getMultiVector() != Teuchos::null)
	  dgdp_sg[i][j].getMultiVector()->init(0.0);
	else if (dgdp_sg[i][j].getLinearOp() != Teuchos::null)
	  dgdp_sg[i][j].getLinearOp()->init(0.0);
      }
    }
  }

  if (do_quad) {
    // Get quadrature data
    const Teuchos::Array<double>& quad_weights = 
      quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> > & quad_values = 
      quad->getBasisAtQuadPoints();
    const Teuchos::Array<double>& basis_norms = basis->norm_squared();

    // Evaluate inputs at quadrature points
    int nqp = mp_block_map->NumMyElements();
    for (int qp=0; qp<nqp; qp++) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
      TEUCHOS_FUNC_TIME_MONITOR_DIFF(
	"Stokhos: SGQuadMPModelEvaluator -- Polynomial Evaluation",
	PolyEvaluation);
#endif

      int gqp = mp_block_map->GID(qp);
      
      if (x_sg != Teuchos::null) {
	x_sg->evaluate(quad_values[gqp], (*x_mp)[qp]);
	me_inargs.set_x_mp(x_mp);
      }
      if (x_dot_sg != Teuchos::null) {
	x_dot_sg->evaluate(quad_values[gqp], (*x_dot_mp)[qp]);
	me_inargs.set_x_dot_mp(x_mp);
      }
      for (int i=0; i<inArgs.Np_sg(); i++) {
	if (p_sg[i] != Teuchos::null) {
	  p_sg[i]->evaluate(quad_values[gqp], (*(p_mp[i]))[qp]);
	  me_inargs.set_p_mp(i, p_mp[i]);
	}
      }

    }

    // Set OutArgs
    if (f_sg != Teuchos::null)
      me_outargs.set_f_mp(f_mp);
    if (W_sg != Teuchos::null)
      me_outargs.set_W_mp(W_mp);
    for (int i=0; i<inArgs.Np_sg(); i++) {
      if (!dfdp_sg[i].isEmpty())
	me_outargs.set_DfDp_mp(i, dfdp_mp[i]);
    }
    for (int i=0; i<outArgs.Ng_sg(); i++) {
      if (g_sg[i] != Teuchos::null)
	me_outargs.set_g_mp(i, g_mp[i]);
      if (!dgdx_dot_sg[i].isEmpty())
	me_outargs.set_DgDx_dot_mp(i, dgdx_dot_mp[i]);
      if (!dgdx_sg[i].isEmpty())
	me_outargs.set_DgDx_mp(i, dgdx_mp[i]);
      for (int j=0; j<outArgs.Np(); j++)
	if (!dgdp_sg[i][j].isEmpty())
	  me_outargs.set_DgDp_mp(i, j, dgdp_mp[i][j]);
    }
    

    {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
      TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadMPModelEvaluator -- Model Evaluation");
#endif
      
      // Evaluate multi-point model at quadrature points
      me->evalModel(me_inargs, me_outargs);
      
    }

    // Perform integrations
    for (int qp=0; qp<nqp; qp++) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
      TEUCHOS_FUNC_TIME_MONITOR_DIFF(
	"Stokhos: SGQuadMPModelEvaluator -- Polynomial Integration", Integration);
#endif

      int gqp = mp_block_map->GID(qp);

      // Sum in results
      if (f_sg != Teuchos::null) {
	f_sg->sumIntoAllTerms(quad_weights[gqp], quad_values[gqp], basis_norms,
			      (*f_mp)[qp]);
      }
      if (W_sg != Teuchos::null) {
	W_sg->sumIntoAllTerms(quad_weights[gqp], quad_values[gqp], basis_norms,
			      (*W_mp)[qp]);
      }
      for (int j=0; j<outArgs.Np(); j++) {
	if (!dfdp_sg[j].isEmpty()) {
	  if (dfdp_sg[j].getMultiVector() != Teuchos::null) {
	    dfdp_sg[j].getMultiVector()->sumIntoAllTerms(
	      quad_weights[gqp], quad_values[gqp], basis_norms, 
	      (*(dfdp_mp[j].getMultiVector()))[qp]);
	  }
	  else if (dfdp_sg[j].getLinearOp() != Teuchos::null) {
	    dfdp_sg[j].getLinearOp()->sumIntoAllTerms(
	      quad_weights[gqp], quad_values[gqp], basis_norms, 
	      (*(dfdp_mp[j].getLinearOp()))[qp]);
	  }
	}
      }
      for (int i=0; i<outArgs.Ng_sg(); i++) {
	if (g_sg[i] != Teuchos::null) {
	  g_sg[i]->sumIntoAllTerms(quad_weights[gqp], quad_values[gqp], 
				   basis_norms, (*g_mp[i])[qp]);
	}
	if (!dgdx_dot_sg[i].isEmpty()) {
	  if (dgdx_dot_sg[i].getMultiVector() != Teuchos::null) {
	    dgdx_dot_sg[i].getMultiVector()->sumIntoAllTerms(
	      quad_weights[gqp], quad_values[gqp], basis_norms, 
	      (*(dgdx_dot_mp[i].getMultiVector()))[qp]);
	  }
	  else if (dgdx_dot_sg[i].getLinearOp() != Teuchos::null) {
	    dgdx_dot_sg[i].getLinearOp()->sumIntoAllTerms(
	      quad_weights[gqp], quad_values[gqp], basis_norms, 
	      (*(dgdx_dot_mp[i].getLinearOp()))[qp]);
	  }
	}
	if (!dgdx_sg[i].isEmpty()) {
	  if (dgdx_sg[i].getMultiVector() != Teuchos::null) {
	    dgdx_sg[i].getMultiVector()->sumIntoAllTerms(
	      quad_weights[gqp], quad_values[gqp], basis_norms, 
	      (*(dgdx_mp[i].getMultiVector()))[qp]);
	  }
	  else if (dgdx_sg[i].getLinearOp() != Teuchos::null) {
	    dgdx_sg[i].getLinearOp()->sumIntoAllTerms(
	      quad_weights[gqp], quad_values[gqp], basis_norms, 
	      (*(dgdx_mp[i].getLinearOp()))[qp]);
	  }
	}
	for (int j=0; j<outArgs.Np(); j++) {
	  if (!dgdp_sg[i][j].isEmpty()) {
	    if (dgdp_sg[i][j].getMultiVector() != Teuchos::null) {
	      dgdp_sg[i][j].getMultiVector()->sumIntoAllTerms(
		quad_weights[gqp], quad_values[gqp], basis_norms, 
		(*(dgdp_mp[i][j].getMultiVector()))[qp]);
	    }
	    else if (dgdp_sg[i][j].getLinearOp() != Teuchos::null) {
	      dgdp_sg[i][j].getLinearOp()->sumIntoAllTerms(
		quad_weights[gqp], quad_values[gqp], basis_norms, 
		(*(dgdp_mp[i][j].getLinearOp()))[qp]);
	    }
	  }
	}
      }
      
    }

    // Now communicate partial sums across processors
    if (mp_block_map->DistributedGlobal()) {
      for (int i=0; i<outArgs.Ng_sg(); i++) {
	if (g_sg[i] != Teuchos::null) {
	  g_sg[i]->sumAll();
	}
      }
    }
  }
  else {
    // Compute the non-SG functions
    me->evalModel(me_inargs, me_outargs);
  }
}
