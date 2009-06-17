// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

#include "Stokhos_SGQuadModelEvaluator.hpp"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Teuchos_TimeMonitor.hpp"

Stokhos::SGQuadModelEvaluator::
SGQuadModelEvaluator(
	    const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
	    const Teuchos::RCP< const Stokhos::Quadrature<int,double> >& quad_,
	    const std::vector<int>& sg_p_index_,
	    const std::vector<int>& sg_g_index_) : 
  me(me_),
  quad(quad_),
  sg_p_index(sg_p_index_),
  sg_g_index(sg_g_index_),
  num_p(sg_p_index.size()),
  num_g(sg_g_index.size()),
  x_dot_qp(),
  x_qp(),
  p_qp(),
  f_qp(),
  W_qp()
{
  // Create storage for x_dot, x, and p at a quad point
  InArgs me_inargs = me->createInArgs();
  if (me_inargs.supports(IN_ARG_x_dot))
    x_dot_qp = Teuchos::rcp(new Epetra_Vector(*(me->get_x_map())));
  if (me_inargs.supports(IN_ARG_x))
    x_qp = Teuchos::rcp(new Epetra_Vector((*me->get_x_map())));
  p_qp.resize(me_inargs.Np());
  for (int i=0; i<me_inargs.Np(); i++)
    p_qp[i] = Teuchos::rcp(new Epetra_Vector(*(me->get_p_map(i))));

  // Create storage for f and W at a quad point
  OutArgs me_outargs = me->createOutArgs();
  if (me_outargs.supports(OUT_ARG_f))
    f_qp = Teuchos::rcp(new Epetra_Vector(*(me->get_f_map())));
  if (me_outargs.supports(OUT_ARG_W))
    W_qp = me->create_W();
  g_qp.resize(me_outargs.Ng());
  for (int i=0; i<me_outargs.Ng(); i++)
    g_qp[i] = Teuchos::rcp(new Epetra_Vector(*(me->get_g_map(i))));
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadModelEvaluator::
get_x_map() const
{
  return me->get_x_map();
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadModelEvaluator::
get_f_map() const
{
  return me->get_f_map();
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadModelEvaluator::
get_p_map(int l) const
{
  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGQuadModelEvaluator::
get_g_map(int l) const
{
  return me->get_g_map(l);
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGQuadModelEvaluator::
get_p_names(int l) const
{
  return me->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGQuadModelEvaluator::
get_x_init() const
{
  return me->get_x_init();
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGQuadModelEvaluator::
get_p_init(int l) const
{
  return me->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGQuadModelEvaluator::
create_W() const
{
  return me->create_W();
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::SGQuadModelEvaluator::
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

  inArgs.set_Np_sg(num_p);
  if (me_inargs.supports(IN_ARG_x))
    inArgs.setSupports(IN_ARG_x_sg, true);
  if (me_inargs.supports(IN_ARG_x_dot))
    inArgs.setSupports(IN_ARG_x_dot_sg, true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::SGQuadModelEvaluator::
createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_outargs.Np(), me_outargs.Ng());
  outArgs.setSupports(OUT_ARG_f, me_outargs.supports(OUT_ARG_f));
  outArgs.setSupports(OUT_ARG_W, me_outargs.supports(OUT_ARG_W));

  outArgs.set_Ng_sg(num_g);
  if (me_outargs.supports(OUT_ARG_f))
    outArgs.setSupports(OUT_ARG_f_sg, true);
  if (me_outargs.supports(OUT_ARG_W)) {
    outArgs.set_W_properties(me_outargs.get_W_properties());
    outArgs.setSupports(OUT_ARG_W_sg, true);
  }
  
  return outArgs;
}

void 
Stokhos::SGQuadModelEvaluator::
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
  for (int i=0; i<outArgs.Ng(); i++)
    me_outargs.set_g(i, outArgs.get_g(i));

  bool do_quad = false;
  InArgs::sg_const_vector_t x_sg;
  InArgs::sg_const_vector_t x_dot_sg;
  Teuchos::Array<InArgs::sg_const_vector_t> p_sg(inArgs.Np_sg());
  OutArgs::sg_vector_t f_sg;
  OutArgs::sg_operator_t W_sg;
  Teuchos::Array<OutArgs::sg_vector_t> g_sg(outArgs.Ng_sg());
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
  if (inArgs.supports(IN_ARG_x_sg)) {
    x_sg = inArgs.get_x_sg();
    if (x_sg != Teuchos::null) {
      basis = x_sg->basis();
      do_quad = true;
    }
  }
  if (inArgs.supports(IN_ARG_x_dot_sg)) {
    x_dot_sg = inArgs.get_x_dot_sg();
    if (x_dot_sg != Teuchos::null) {
      basis = x_dot_sg->basis();
      do_quad = true;
    }
  }
  for (int i=0; i<inArgs.Np_sg(); i++) {
    p_sg[i] = inArgs.get_p_sg(i);
    if (p_sg[i] != Teuchos::null) {
      basis = p_sg[i]->basis();
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
  for (int i=0; i<outArgs.Ng_sg(); i++) {
    g_sg[i] = outArgs.get_g_sg(i);
    if (g_sg[i] != Teuchos::null)
      g_sg[i]->init(0.0);
  }

  if (do_quad) {
    // Get quadrature data
    const std::vector< std::vector<double> >& quad_points = 
      quad->getQuadPoints();
    const std::vector<double>& quad_weights = 
      quad->getQuadWeights();
    const std::vector< std::vector<double> > & quad_values = 
      quad->getBasisAtQuadPoints();
    const std::vector<double>& basis_norms = basis->norm_squared();

    // Perform integrations
    for (unsigned int qp=0; qp<quad_points.size(); qp++) {

      {
      TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- Polynomial Evaluation");

      // Evaluate inputs at quadrature points
      if (x_sg != Teuchos::null) {
	TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- X Evaluation");
	x_sg->evaluate(quad_values[qp], *x_qp);
	me_inargs.set_x(x_qp);
      }
      if (x_dot_sg != Teuchos::null) {
	TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- X_dot Evaluation");
	x_dot_sg->evaluate(quad_values[qp], *x_dot_qp);
	me_inargs.set_x_dot(x_qp);
      }
      for (int i=0; i<inArgs.Np_sg(); i++) {
	if (p_sg[i] != Teuchos::null) {
	  TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- P Evaluation");
	  p_sg[i]->evaluate(quad_values[qp], *(p_qp[i]));
	  me_inargs.set_p(sg_p_index[i], p_qp[i]);
	}
      }
      if (f_sg != Teuchos::null)
	me_outargs.set_f(f_qp);
      if (W_sg != Teuchos::null)
	me_outargs.set_W(W_qp);
      for (int i=0; i<outArgs.Ng_sg(); i++)
	me_outargs.set_g(sg_g_index[i], g_qp[i]);
      }

      {
      TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- Model Evaluation");

      // Evaluate model at quadrature points
      me->evalModel(me_inargs, me_outargs);

      }

      {
      TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- Polynomial Integration");

      // Sum in results
      if (f_sg != Teuchos::null) {
	TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- F Integration");
	f_sg->sumIntoAllTerms(quad_weights[qp], quad_values[qp], basis_norms,
			      *f_qp);
      }
      if (W_sg != Teuchos::null) {
	TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- W Integration");
	W_sg->sumIntoAllTerms(quad_weights[qp], quad_values[qp], basis_norms,
			      *W_qp);
      }
      for (int i=0; i<outArgs.Ng_sg(); i++)
	if (g_sg[i] != Teuchos::null) {
	  TEUCHOS_FUNC_TIME_MONITOR("SGQuadModelEvaluator -- G Integration");
	  g_sg[i]->sumIntoAllTerms(quad_weights[qp], quad_values[qp], 
				   basis_norms, *g_qp[i]);
	}
      }
    }
  }
  else {
    // Compute the non-SG functions
    me->evalModel(me_inargs, me_outargs);
  }
}
