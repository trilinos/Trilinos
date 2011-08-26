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

#include "Stokhos_SGQuadModelEvaluator.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TestForException.hpp"

Stokhos::SGQuadModelEvaluator::
SGQuadModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_) : 
  me(me_),
  num_p(0),
  num_g(0),
  x_dot_qp(),
  x_qp(),
  p_qp(),
  f_qp(),
  W_qp(),
  dfdp_qp(),
  g_qp(),
  dgdx_qp(),
  dgdx_dot_qp(),
  dgdp_qp()
{
  // Create storage for x_dot, x, and p at a quad point
  InArgs me_inargs = me->createInArgs();
  num_p = me_inargs.Np();
  if (me_inargs.supports(IN_ARG_x_dot))
    x_dot_qp = Teuchos::rcp(new Epetra_Vector(*(me->get_x_map())));
  if (me_inargs.supports(IN_ARG_x))
    x_qp = Teuchos::rcp(new Epetra_Vector((*me->get_x_map())));
  p_qp.resize(num_p);
  for (int i=0; i<num_p; i++)
    p_qp[i] = Teuchos::rcp(new Epetra_Vector(*(me->get_p_map(i))));

  // Create storage for f and W at a quad point
  OutArgs me_outargs = me->createOutArgs();
  num_g = me_outargs.Ng();

  // f
  if (me_outargs.supports(OUT_ARG_f))
    f_qp = Teuchos::rcp(new Epetra_Vector(*(me->get_f_map())));

  // W
  if (me_outargs.supports(OUT_ARG_W))
    W_qp = me->create_W();

  // df/dp
  dfdp_qp.resize(num_p);
  for (int i=0; i<num_p; i++)
    if (me_outargs.supports(OUT_ARG_DfDp,i).supports(DERIV_MV_BY_COL))
      dfdp_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	Teuchos::rcp(new Epetra_MultiVector(
		       *(me->get_f_map()),
		       me->get_p_map(i)->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DfDp,i).supports(DERIV_TRANS_MV_BY_ROW))
      dfdp_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	Teuchos::rcp(new Epetra_MultiVector(
		       *(me->get_p_map(i)),
		       me->get_f_map()->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DfDp,i).supports(DERIV_LINEAR_OP))
      dfdp_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	me->create_DfDp_op(i));

  g_qp.resize(num_g);
  dgdx_qp.resize(num_g);
  dgdx_dot_qp.resize(num_g);
  dgdp_qp.resize(num_g);
  for (int i=0; i<num_g; i++) {

    // g
    g_qp[i] = 
      Teuchos::rcp(new Epetra_Vector(*(me->get_g_map(i))));

    // dg/dx
    if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_TRANS_MV_BY_ROW))
      dgdx_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	Teuchos::rcp(new Epetra_MultiVector(
		       *(me->get_x_map()),
		       me->get_g_map(i)->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_MV_BY_COL))
      dgdx_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	Teuchos::rcp(new Epetra_MultiVector(
		       *(me->get_g_map(i)),
		       me->get_x_map()->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP))
      dgdx_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	me->create_DgDx_op(i));
    
    // dg/dx_dot
    if (me_outargs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_TRANS_MV_BY_ROW))
      dgdx_dot_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	Teuchos::rcp(new Epetra_MultiVector(
		       *(me->get_x_map()),
		       me->get_g_map(i)->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_MV_BY_COL))
      dgdx_dot_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	Teuchos::rcp(new Epetra_MultiVector(
		       *(me->get_g_map(i)),
		       me->get_x_map()->NumGlobalElements())));
    else if (me_outargs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_LINEAR_OP))
      dgdx_dot_qp[i] = EpetraExt::ModelEvaluator::Derivative(
	me->create_DgDx_dot_op(i));

    // dg/dp
    dgdp_qp[i].resize(num_p);
    for (int j=0; j<num_p; j++)
      if (me_outargs.supports(OUT_ARG_DgDp, i, j).supports(DERIV_TRANS_MV_BY_ROW))
	dgdp_qp[i][j] = EpetraExt::ModelEvaluator::Derivative(
	  Teuchos::rcp(new Epetra_MultiVector(
			 *(me->get_p_map(j)),
			 me->get_g_map(i)->NumGlobalElements())));
      else if (me_outargs.supports(OUT_ARG_DgDp, i, j).supports(DERIV_MV_BY_COL))
	dgdp_qp[i][j] = EpetraExt::ModelEvaluator::Derivative(
	  Teuchos::rcp(new Epetra_MultiVector(
			 *(me->get_g_map(i)),
			 me->get_p_map(j)->NumGlobalElements())));
      else if (me_outargs.supports(OUT_ARG_DgDp, i, j).supports(DERIV_LINEAR_OP))
	dgdp_qp[i][j] = EpetraExt::ModelEvaluator::Derivative(
	  me->create_DgDp_op(i,j));
  }
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
  inArgs.set_Np(num_p); 
  inArgs.setSupports(IN_ARG_x_dot, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_t, me_inargs.supports(IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha, me_inargs.supports(IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta, me_inargs.supports(IN_ARG_beta));

  for (int i=0; i<num_p; i++)
    inArgs.setSupports(IN_ARG_p_sg, i, true);
  inArgs.setSupports(IN_ARG_x_sg, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_x_dot_sg, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_sg_basis, true);
  inArgs.setSupports(IN_ARG_sg_quadrature, true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::SGQuadModelEvaluator::
createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(num_p, num_g);
  outArgs.setSupports(OUT_ARG_f, me_outargs.supports(OUT_ARG_f));
  outArgs.setSupports(OUT_ARG_W, me_outargs.supports(OUT_ARG_W));
  for (int j=0; j<num_p; j++)
    outArgs.setSupports(OUT_ARG_DfDp, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  for (int i=0; i<num_g; i++) {
    outArgs.setSupports(OUT_ARG_DgDx, i, 
			me_outargs.supports(OUT_ARG_DgDx, i));
    outArgs.setSupports(OUT_ARG_DgDx_dot, i, 
			me_outargs.supports(OUT_ARG_DgDx_dot, i));
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
  }

  outArgs.setSupports(OUT_ARG_f_sg, me_outargs.supports(OUT_ARG_f));
  if (me_outargs.supports(OUT_ARG_W)) {
    outArgs.set_W_properties(me_outargs.get_W_properties());
    outArgs.setSupports(OUT_ARG_W_sg, true);
  }
  for (int j=0; j<num_p; j++) {
    outArgs.setSupports(OUT_ARG_DfDp_sg, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  }
  for (int i=0; i<num_g; i++) {
    outArgs.setSupports(OUT_ARG_g_sg, i, true);
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
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();
  if (me_outargs.supports(OUT_ARG_f))
    me_outargs.set_f(outArgs.get_f());
  if (me_outargs.supports(OUT_ARG_W))
    me_outargs.set_W(outArgs.get_W());
  for (int j=0; j<num_p; j++)
    if (!outArgs.supports(OUT_ARG_DfDp, j).none())
      me_outargs.set_DfDp(j, outArgs.get_DfDp(j));
  for (int i=0; i<num_g; i++) {
    me_outargs.set_g(i, outArgs.get_g(i));
    if (!outArgs.supports(OUT_ARG_DgDx, i).none())
	me_outargs.set_DgDx(i, outArgs.get_DgDx(i));
    if (!outArgs.supports(OUT_ARG_DgDx_dot, i).none())
	me_outargs.set_DgDx(i, outArgs.get_DgDx_dot(i));
    for (int j=0; j<num_p; j++)
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	me_outargs.set_DgDp(i, j, outArgs.get_DgDp(i,j));
  }

  bool do_quad = false;
  InArgs::sg_const_vector_t x_sg;
  InArgs::sg_const_vector_t x_dot_sg;
  Teuchos::Array<InArgs::sg_const_vector_t> p_sg(num_p);
  OutArgs::sg_vector_t f_sg;
  OutArgs::sg_operator_t W_sg;
  Teuchos::Array<SGDerivative> dfdp_sg(num_p);
  Teuchos::Array<OutArgs::sg_vector_t> g_sg(num_g);
  Teuchos::Array<SGDerivative> dgdx_sg(num_g);
  Teuchos::Array<SGDerivative> dgdx_dot_sg(num_g);
  Teuchos::Array< Teuchos::Array<SGDerivative> > dgdp_sg(num_g);
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
  for (int i=0; i<num_p; i++) {
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
  for (int i=0; i<num_p; i++) {
    if (!outArgs.supports(OUT_ARG_DfDp_sg, i).none()) {
      dfdp_sg[i] = outArgs.get_DfDp_sg(i);
      if (dfdp_sg[i].getMultiVector() != Teuchos::null)
	dfdp_sg[i].getMultiVector()->init(0.0);
      else if (dfdp_sg[i].getLinearOp() != Teuchos::null)
	dfdp_sg[i].getLinearOp()->init(0.0);
    }
  }
      
  for (int i=0; i<num_g; i++) {
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

    dgdp_sg[i].resize(num_p);
    for (int j=0; j<num_p; j++) {
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
    const Teuchos::Array< Teuchos::Array<double> >& quad_points = 
      quad->getQuadPoints();
    const Teuchos::Array<double>& quad_weights = 
      quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> > & quad_values = 
      quad->getBasisAtQuadPoints();
    const Teuchos::Array<double>& basis_norms = basis->norm_squared();

    // Perform integrations
    for (int qp=0; qp<quad_points.size(); qp++) {

      // StieltjesPCEBasis can introduce quadrature points with zero weight
      // Don't do those evaluations, since the model might not like the
      // quadrature points (i.e., zero)
      if (quad_weights[qp] == 0.0)
	continue;

      {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Stokhos: SGQuadModelEvaluator -- Polynomial Evaluation",
          PolyEvaluation);
#endif

        // Evaluate inputs at quadrature points
        if (x_sg != Teuchos::null) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
          TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- X Evaluation");
#endif
          x_sg->evaluate(quad_values[qp], *x_qp);
          me_inargs.set_x(x_qp);
        }
        if (x_dot_sg != Teuchos::null) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
          TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- X_dot Evaluation");
#endif
          x_dot_sg->evaluate(quad_values[qp], *x_dot_qp);
          me_inargs.set_x_dot(x_qp);
        }
        for (int i=0; i<num_p; i++) {
          if (p_sg[i] != Teuchos::null) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
            TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- P Evaluation");
#endif
            p_sg[i]->evaluate(quad_values[qp], *(p_qp[i]));
            me_inargs.set_p(i, p_qp[i]);
          }
        }
        if (f_sg != Teuchos::null)
          me_outargs.set_f(f_qp);
        if (W_sg != Teuchos::null)
          me_outargs.set_W(W_qp);
	for (int i=0; i<num_p; i++) {
	  if (!dfdp_sg[i].isEmpty())
	    me_outargs.set_DfDp(i, dfdp_qp[i]);
	}
        for (int i=0; i<num_g; i++) {
	  if (g_sg[i] != Teuchos::null)
	    me_outargs.set_g(i, g_qp[i]);
	  if (!dgdx_dot_sg[i].isEmpty())
	    me_outargs.set_DgDx_dot(i, dgdx_dot_qp[i]);
	  if (!dgdx_sg[i].isEmpty())
	    me_outargs.set_DgDx(i, dgdx_qp[i]);
          for (int j=0; j<num_p; j++)
            if (!dgdp_sg[i][j].isEmpty())
              me_outargs.set_DgDp(i, j, dgdp_qp[i][j]);
        }

      }

      {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
        TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- Model Evaluation");
#endif

        // Evaluate model at quadrature points
        me->evalModel(me_inargs, me_outargs);

      }

      {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
        TEUCHOS_FUNC_TIME_MONITOR_DIFF(
	  "Stokhos: SGQuadModelEvaluator -- Polynomial Integration", Integration);
#endif

        // Sum in results
        if (f_sg != Teuchos::null) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
          TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- F Integration");
#endif
          f_sg->sumIntoAllTerms(quad_weights[qp], quad_values[qp], basis_norms,
				*f_qp);
        }
        if (W_sg != Teuchos::null) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
          TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- W Integration");
#endif
          W_sg->sumIntoAllTerms(quad_weights[qp], quad_values[qp], basis_norms,
				*W_qp);
        }
	for (int j=0; j<num_p; j++) {
	  if (!dfdp_sg[j].isEmpty()) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	    TEUCHOS_FUNC_TIME_MONITOR(
	      "Stokhos: SGQuadModelEvaluator -- df/dp Integration");
#endif
	    if (dfdp_sg[j].getMultiVector() != Teuchos::null) {
	      dfdp_sg[j].getMultiVector()->sumIntoAllTerms(
		quad_weights[qp], quad_values[qp], basis_norms, 
		*(dfdp_qp[j].getMultiVector()));
	    }
	    else if (dfdp_sg[j].getLinearOp() != Teuchos::null) {
	      dfdp_sg[j].getLinearOp()->sumIntoAllTerms(
		quad_weights[qp], quad_values[qp], basis_norms, 
		*(dfdp_qp[j].getLinearOp()));
	    }
	  }
	}
        for (int i=0; i<num_g; i++) {
          if (g_sg[i] != Teuchos::null) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
            TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SGQuadModelEvaluator -- G Integration");
#endif
            g_sg[i]->sumIntoAllTerms(quad_weights[qp], quad_values[qp], 
				     basis_norms, *g_qp[i]);
          }
	  if (!dgdx_dot_sg[i].isEmpty()) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	    TEUCHOS_FUNC_TIME_MONITOR(
	      "Stokhos: SGQuadModelEvaluator -- dg/dx_dot Integration");
#endif
	    if (dgdx_dot_sg[i].getMultiVector() != Teuchos::null) {
	      dgdx_dot_sg[i].getMultiVector()->sumIntoAllTerms(
		quad_weights[qp], quad_values[qp], basis_norms, 
		*(dgdx_dot_qp[i].getMultiVector()));
	    }
	    else if (dgdx_dot_sg[i].getLinearOp() != Teuchos::null) {
	      dgdx_dot_sg[i].getLinearOp()->sumIntoAllTerms(
		quad_weights[qp], quad_values[qp], basis_norms, 
		*(dgdx_dot_qp[i].getLinearOp()));
	    }
	  }
	  if (!dgdx_sg[i].isEmpty()) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	    TEUCHOS_FUNC_TIME_MONITOR(
	      "Stokhos: SGQuadModelEvaluator -- dg/dx Integration");
#endif
	    if (dgdx_sg[i].getMultiVector() != Teuchos::null) {
	      dgdx_sg[i].getMultiVector()->sumIntoAllTerms(
		quad_weights[qp], quad_values[qp], basis_norms, 
		*(dgdx_qp[i].getMultiVector()));
	    }
	    else if (dgdx_sg[i].getLinearOp() != Teuchos::null) {
	      dgdx_sg[i].getLinearOp()->sumIntoAllTerms(
		quad_weights[qp], quad_values[qp], basis_norms, 
		*(dgdx_qp[i].getLinearOp()));
	    }
	  }
          for (int j=0; j<num_p; j++) {
	    if (!dgdp_sg[i][j].isEmpty()) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	      TEUCHOS_FUNC_TIME_MONITOR(
		"Stokhos: SGQuadModelEvaluator -- dg/dp Integration");
#endif
	      if (dgdp_sg[i][j].getMultiVector() != Teuchos::null) {
		dgdp_sg[i][j].getMultiVector()->sumIntoAllTerms(
		  quad_weights[qp], quad_values[qp], basis_norms, 
		  *(dgdp_qp[i][j].getMultiVector()));
	      }
	      else if (dgdp_sg[i][j].getLinearOp() != Teuchos::null) {
		dgdp_sg[i][j].getLinearOp()->sumIntoAllTerms(
		  quad_weights[qp], quad_values[qp], basis_norms, 
		  *(dgdp_qp[i][j].getLinearOp()));
	      }
	    }
          }
        }

      }
    }
  }
  else {
    // Compute the non-SG functions
    me->evalModel(me_inargs, me_outargs);
  }
}
