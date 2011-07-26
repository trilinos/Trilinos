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

#include "Stokhos_MPModelEvaluatorAdapter.hpp"
#include "Stokhos_ProductEpetraVector.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "Stokhos_ProductEpetraOperator.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_Map.h"

Stokhos::MPModelEvaluatorAdapter::
MPModelEvaluatorAdapter(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_) : 
  me(me_)
{
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluatorAdapter::
get_x_map() const
{
  return me->get_x_map();
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluatorAdapter::
get_f_map() const
{
  return me->get_f_map();
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluatorAdapter::
get_p_map(int l) const
{
  return me->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluatorAdapter::
get_g_map(int l) const
{
  return me->get_g_map(l);
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::MPModelEvaluatorAdapter::
get_p_names(int l) const
{
  return me->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::MPModelEvaluatorAdapter::
get_x_init() const
{
  return me->get_x_init();
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::MPModelEvaluatorAdapter::
get_p_init(int l) const
{
  return me->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluatorAdapter::
create_W() const
{
  return me->create_W();
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::MPModelEvaluatorAdapter::
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

  for (int i=0; i<me_inargs.Np(); i++)
    inArgs.setSupports(IN_ARG_p_mp, i, true);
  inArgs.setSupports(IN_ARG_x_mp, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_x_dot_mp, me_inargs.supports(IN_ARG_x_dot));
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::MPModelEvaluatorAdapter::
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

  outArgs.setSupports(OUT_ARG_f_mp, me_outargs.supports(OUT_ARG_f));
  if (me_outargs.supports(OUT_ARG_W)) {
    outArgs.set_W_properties(me_outargs.get_W_properties());
    outArgs.setSupports(OUT_ARG_W_mp, true);
  }
  for (int j=0; j<me_outargs.Np(); j++)
    outArgs.setSupports(OUT_ARG_DfDp_mp, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  for (int i=0; i<me_outargs.Ng(); i++) {
    outArgs.setSupports(OUT_ARG_g_mp, i, true);
    outArgs.setSupports(OUT_ARG_DgDx_mp, i, 
			me_outargs.supports(OUT_ARG_DgDx, i));
    outArgs.setSupports(OUT_ARG_DgDx_dot_mp, i, 
			me_outargs.supports(OUT_ARG_DgDx_dot, i));
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp_mp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
  }
  
  return outArgs;
}

void 
Stokhos::MPModelEvaluatorAdapter::
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

  mp_const_vector_t x_mp;
  mp_const_vector_t x_dot_mp;
  Teuchos::Array<mp_const_vector_t> p_mp(inArgs.Np());
  mp_vector_t f_mp;
  mp_operator_t W_mp;
  Teuchos::Array<MPDerivative> dfdp_mp(outArgs.Np());
  Teuchos::Array<mp_vector_t> g_mp(outArgs.Ng());
  Teuchos::Array<MPDerivative> dgdx_mp(outArgs.Ng());
  Teuchos::Array<MPDerivative> dgdx_dot_mp(outArgs.Ng());
  Teuchos::Array< Teuchos::Array<MPDerivative> > dgdp_mp(outArgs.Ng());
  int num_mp;
  
  if (inArgs.supports(IN_ARG_x_mp)) {
    x_mp = inArgs.get_x_mp();
    if (x_mp != Teuchos::null) {
      num_mp = x_mp->size();
    }
  }
  if (inArgs.supports(IN_ARG_x_dot_mp)) {
    x_dot_mp = inArgs.get_x_dot_mp();
    if (x_dot_mp != Teuchos::null) {
      num_mp = x_dot_mp->size();
    }
  }
  for (int i=0; i<inArgs.Np(); i++) {
    p_mp[i] = inArgs.get_p_mp(i);
    if (p_mp[i] != Teuchos::null) {
      num_mp = p_mp[i]->size();
    }
  }
  if (outArgs.supports(OUT_ARG_f_mp)) {
    f_mp = outArgs.get_f_mp();
    if (f_mp != Teuchos::null)
      f_mp->init(0.0);
  }
  if (outArgs.supports(OUT_ARG_W_mp)) {
    W_mp = outArgs.get_W_mp();
    if (W_mp != Teuchos::null)
      W_mp->init(0.0);
  }
  for (int i=0; i<inArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp_mp, i).none()) {
      dfdp_mp[i] = outArgs.get_DfDp_mp(i);
      if (dfdp_mp[i].getMultiVector() != Teuchos::null)
	dfdp_mp[i].getMultiVector()->init(0.0);
      else if (dfdp_mp[i].getLinearOp() != Teuchos::null)
	dfdp_mp[i].getLinearOp()->init(0.0);
    }
  }
      
  for (int i=0; i<outArgs.Ng(); i++) {
    g_mp[i] = outArgs.get_g_mp(i);
    if (g_mp[i] != Teuchos::null)
      g_mp[i]->init(0.0);
    
    if (!outArgs.supports(OUT_ARG_DgDx_mp, i).none()) {
      dgdx_mp[i] = outArgs.get_DgDx_mp(i);
      if (dgdx_mp[i].getMultiVector() != Teuchos::null)
	dgdx_mp[i].getMultiVector()->init(0.0);
      else if (dgdx_mp[i].getLinearOp() != Teuchos::null)
	dgdx_mp[i].getLinearOp()->init(0.0);
    }

    if (!outArgs.supports(OUT_ARG_DgDx_dot_mp, i).none()) {
      dgdx_dot_mp[i] = outArgs.get_DgDx_dot_mp(i);
      if (dgdx_dot_mp[i].getMultiVector() != Teuchos::null)
	dgdx_dot_mp[i].getMultiVector()->init(0.0);
      else if (dgdx_dot_mp[i].getLinearOp() != Teuchos::null)
	dgdx_dot_mp[i].getLinearOp()->init(0.0);
    }

    dgdp_mp[i].resize(outArgs.Np());
    for (int j=0; j<outArgs.Np(); j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_mp, i, j).none()) {
	dgdp_mp[i][j] = outArgs.get_DgDp_mp(i,j);
	if (dgdp_mp[i][j].getMultiVector() != Teuchos::null)
	  dgdp_mp[i][j].getMultiVector()->init(0.0);
	else if (dgdp_mp[i][j].getLinearOp() != Teuchos::null)
	  dgdp_mp[i][j].getLinearOp()->init(0.0);
      }
    }
  }

  for (int qp=0; qp<num_mp; qp++) {
    
    // Set InArgs
    if (x_mp != Teuchos::null)
      me_inargs.set_x(x_mp->getCoeffPtr(qp));
    if (x_dot_mp != Teuchos::null)
      me_inargs.set_x_dot(x_dot_mp->getCoeffPtr(qp));
    for (int i=0; i<inArgs.Np(); i++)
      if (p_mp[i] != Teuchos::null)
	me_inargs.set_p(i, p_mp[i]->getCoeffPtr(qp));
    
    // Set OutArgs
    if (f_mp != Teuchos::null)
      me_outargs.set_f(f_mp->getCoeffPtr(qp));
    if (W_mp != Teuchos::null)
      me_outargs.set_W(W_mp->getCoeffPtr(qp));
    
    for (int i=0; i<inArgs.Np(); i++) {
      if (!dfdp_mp[i].isEmpty()) {
	if (dfdp_mp[i].getMultiVector() != Teuchos::null) {
	  Derivative deriv(dfdp_mp[i].getMultiVector()->getCoeffPtr(qp),
			   dfdp_mp[i].getMultiVectorOrientation());
	  me_outargs.set_DfDp(i, deriv);
	}
	else if (dfdp_mp[i].getLinearOp() != Teuchos::null) {
	  Derivative deriv(dfdp_mp[i].getLinearOp()->getCoeffPtr(qp));
	  me_outargs.set_DfDp(i, deriv);
	}
      }
    }
    
    for (int i=0; i<outArgs.Ng(); i++) {
      if (g_mp[i] != Teuchos::null)
	me_outargs.set_g(i, g_mp[i]->getCoeffPtr(qp));
      if (!dgdx_dot_mp[i].isEmpty()) {
	if (dgdx_dot_mp[i].getMultiVector() != Teuchos::null) {
	  Derivative deriv(dgdx_dot_mp[i].getMultiVector()->getCoeffPtr(qp),
			   dgdx_dot_mp[i].getMultiVectorOrientation());
	  me_outargs.set_DgDx_dot(i, deriv);
	}
	else if (dgdx_dot_mp[i].getLinearOp() != Teuchos::null) {
	  Derivative deriv(dgdx_dot_mp[i].getLinearOp()->getCoeffPtr(qp));
	  me_outargs.set_DgDx_dot(i, deriv);
	}
      }
      if (!dgdx_mp[i].isEmpty()) {
	if (dgdx_mp[i].getMultiVector() != Teuchos::null) {
	  Derivative deriv(dgdx_mp[i].getMultiVector()->getCoeffPtr(qp),
			   dgdx_mp[i].getMultiVectorOrientation());
	  me_outargs.set_DgDx(i, deriv);
	}
	else if (dgdx_mp[i].getLinearOp() != Teuchos::null) {
	  Derivative deriv(dgdx_mp[i].getLinearOp()->getCoeffPtr(qp));
	  me_outargs.set_DgDx(i, deriv);
	}
      }
      for (int j=0; j<outArgs.Np(); j++) {
	if (!dgdp_mp[i][j].isEmpty()) {
	  if (dgdp_mp[i][j].getMultiVector() != Teuchos::null) {
	    Derivative deriv(dgdp_mp[i][j].getMultiVector()->getCoeffPtr(qp),
			     dgdp_mp[i][j].getMultiVectorOrientation());
	    me_outargs.set_DgDp(i, j, deriv);
	  }
	  else if (dgdp_mp[i][j].getLinearOp() != Teuchos::null) {
	    Derivative deriv(dgdp_mp[i][j].getLinearOp()->getCoeffPtr(qp));
	    me_outargs.set_DgDp(i, j, deriv);
	  }
	}
      }
    }
    
    // Evaluate model
    me->evalModel(me_inargs, me_outargs);

  }

  // Evaluate single-point model
  if (num_mp == 0)
    me->evalModel(me_inargs, me_outargs);
}
