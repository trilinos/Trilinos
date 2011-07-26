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

#include "Stokhos_MPModelEvaluator.hpp"

#include <algorithm>
#include "Teuchos_TestForException.hpp"
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_MPPreconditionerFactory.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"
#include "Stokhos_ProductEpetraOperator.hpp"


Stokhos::MPModelEvaluator::MPModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm_,
  const Teuchos::RCP<const Epetra_Map>& mp_block_map_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) 
  : me(me_),
    num_mp_blocks(mp_block_map_->NumMyElements()),
    mp_comm(mp_comm_),
    mp_block_map(mp_block_map_),
    params(params_),
    supports_x(false),
    x_map(me->get_x_map()),
    f_map(me->get_f_map()),
    mp_x_map(),
    mp_f_map(),
    num_p(0),
    num_p_mp(0),
    mp_p_index_map(),
    mp_p_map(),
    mp_p_names(),
    num_g(0),
    num_g_mp(0),
    mp_g_index_map(),
    mp_g_map(),
    W_mp_blocks(),
    dgdx_dot_mp_blocks(),
    dgdx_mp_blocks(),
    mp_p_init(),
    my_W(),
    my_x()
{
  if (x_map != Teuchos::null)
    supports_x = true;
  
  if (supports_x) {

    // Create block MP x and f maps
    mp_x_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *x_map, *mp_block_map, *mp_comm));

    mp_f_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *f_map, *mp_block_map, *mp_comm));

    // Create default mp_x_init
    mp_x_init = Teuchos::rcp(new ProductEpetraVector(
			       mp_block_map, x_map, mp_x_map, mp_comm));
    for (unsigned int i=0; i<num_mp_blocks; i++)
      (*mp_x_init)[i] = *(me->get_x_init());

    // Preconditioner needs an x
    my_x = Teuchos::rcp(new Epetra_Vector(*mp_x_map));

    W_mp_blocks = 
      Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		     mp_block_map, x_map, f_map, mp_f_map, mp_comm));
    for (unsigned int i=0; i<num_mp_blocks; i++)
      W_mp_blocks->setCoeffPtr(i, me->create_W());
  }
    
  // Parameters -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the parameters

  InArgs me_inargs = me->createInArgs();
  num_p = me_inargs.Np();
  
  // Get the p_mp's supported and build index map
  for (int i=0; i<num_p; i++) {
    if (me_inargs.supports(IN_ARG_p_mp, i))
      mp_p_index_map.push_back(i);
  }
  num_p_mp = mp_p_index_map.size();

  mp_p_map.resize(num_p_mp);
  mp_p_names.resize(num_p_mp);
  mp_p_init.resize(num_p_mp);
  
  // Create parameter maps, names, and initial values
  for (int i=0; i<num_p_mp; i++) {
    Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(mp_p_index_map[i]);
    mp_p_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *p_map, *mp_block_map, *mp_comm));
    
    Teuchos::RCP<const Teuchos::Array<std::string> > p_names = 
      me->get_p_names(mp_p_index_map[i]);
    if (p_names != Teuchos::null) {
      mp_p_names[i] = 
	Teuchos::rcp(new Teuchos::Array<std::string>(num_mp_blocks*(p_names->size())));
      for (int j=0; j<p_names->size(); j++) {
	std::stringstream ss;
	ss << (*p_names)[j] << " -- MP Coefficient " << i;
	(*mp_p_names[i])[j] = ss.str();
      }
    }
    
    // Create default mp_p_init
    mp_p_init[i] = 
      Teuchos::rcp(new ProductEpetraVector(
		     mp_block_map, p_map, mp_p_map[i], mp_comm));
    mp_p_init[i]->init(0.0);
    // for (unsigned int j=0; j<num_mp_blocks; j++)
    //   (*mp_p_init[i])[j] = *(me->get_p_init(i));
    
  }

  // Responses -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the respones

  // Get number of MP parameters model supports derivatives w.r.t.
  OutArgs me_outargs = me->createOutArgs();
  num_g = me_outargs.Ng();
  
  // Get the g_mp's supported and build index map
  for (int i=0; i<num_g; i++) {
    if (me_outargs.supports(OUT_ARG_g_mp, i))
      mp_g_index_map.push_back(i);
  }
  num_g_mp = mp_g_index_map.size();

  mp_g_map.resize(num_g_mp);
  dgdx_dot_mp_blocks.resize(num_g_mp);
  dgdx_mp_blocks.resize(num_g_mp);

  // Create response maps
  for (int i=0; i<num_g_mp; i++) {
    Teuchos::RCP<const Epetra_Map> g_map = me->get_g_map(mp_g_index_map[i]);
    mp_g_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *g_map, *mp_block_map, *mp_comm));
    
    // Create dg/dxdot, dg/dx MP blocks
    dgdx_dot_mp_blocks[i] = 
      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(mp_block_map));
    dgdx_mp_blocks[i] = 
      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(mp_block_map));
  }
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_x_map() const
{
  return mp_x_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_f_map() const
{
  return mp_f_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l < 0 || l >= num_p + num_p_mp, std::logic_error,
		     "Error!  Invalid p map index " << l);
  if (l < num_p)
    return me->get_p_map(l);
  else 
    return mp_p_map[l-num_p];
  
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j < 0 || j >= num_g_mp, std::logic_error,
		     "Error!  Invalid g map index " << j);
  return mp_g_map[j];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::MPModelEvaluator::get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(l < 0 || l >= num_p + num_p_mp, std::logic_error,
		     "Error!  Invalid p map index " << l);
  if (l < num_p)
    return me->get_p_names(l);
  else
    return mp_p_names[l-num_p];
  
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::MPModelEvaluator::get_x_init() const
{
  return mp_x_init->getBlockVector();
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::MPModelEvaluator::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l < 0 || l >= num_p + num_p_mp, std::logic_error,
		     "Error!  Invalid p map index " << l);
  if (l < num_p)
    return me->get_p_init(l);
  else if (l < num_p + num_p_mp)
    return mp_p_init[l-num_p]->getBlockVector();
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_W() const
{
  if (supports_x) {
    my_W = 
      Teuchos::rcp(new Stokhos::BlockDiagonalOperator(mp_comm, num_mp_blocks, 
						      x_map, f_map, 
						      mp_x_map, mp_f_map));
    my_W->setupOperator(W_mp_blocks);
    
    return my_W;
  }
  
  return Teuchos::null;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Stokhos::MPModelEvaluator::create_WPrec() const
{
  if (supports_x) {
    Teuchos::RCP<Teuchos::ParameterList> mpPrecParams =
      Teuchos::rcp(&(params->sublist("MP Preconditioner")), false);
    Stokhos::MPPreconditionerFactory mp_prec_factory(mpPrecParams);
    Teuchos::RCP<Epetra_Operator> precOp = 
      mp_prec_factory.build(mp_comm, num_mp_blocks, x_map, mp_x_map);
    return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,
								      true));
  }
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_DgDx_op(int j) const
{
  TEST_FOR_EXCEPTION(j < 0 || j >= num_g_mp || !supports_x, 
		     std::logic_error,
		     "Error:  dg/dx index " << j << " is not supported!");
  
  int jj = mp_g_index_map[j];
  Teuchos::RCP<const Epetra_Map> g_map = me->get_g_map(jj);
  Teuchos::RCP< Stokhos::ProductEpetraOperator > mp_blocks =
    Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		   mp_block_map, x_map, g_map, mp_g_map[j], mp_comm));
  OutArgs me_outargs = me->createOutArgs();
  DerivativeSupport ds = me_outargs.supports(OUT_ARG_DgDx_mp, jj);
  if (ds.supports(DERIV_LINEAR_OP))
    for (unsigned int i=0; i<num_mp_blocks; i++)
      mp_blocks->setCoeffPtr(i, me->create_DgDx_op(jj));
  else if (ds.supports(DERIV_MV_BY_COL))
    for (unsigned int i=0; i<num_mp_blocks; i++) {
      Teuchos::RCP<Epetra_MultiVector> mv = 
	Teuchos::rcp(new Epetra_MultiVector(*g_map, x_map->NumMyElements()));
      Teuchos::RCP<Epetra_Operator> block = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
      mp_blocks->setCoeffPtr(i, block);
    }
  else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
    for (unsigned int i=0; i<num_mp_blocks; i++) {
      Teuchos::RCP<Epetra_MultiVector> mv = 
	Teuchos::rcp(new Epetra_MultiVector(*x_map, g_map->NumMyElements()));
      Teuchos::RCP<Epetra_Operator> block = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
      mp_blocks->setCoeffPtr(i, block);
    }
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  me_outargs.supports(OUT_ARG_DgDx_mp, " << j
		       << ").none() is true!");
  
  Teuchos::RCP<Stokhos::BlockDiagonalOperator> dgdx_mp = 
    Teuchos::rcp(new Stokhos::BlockDiagonalOperator(
		   mp_comm, num_mp_blocks, x_map, g_map,
		   mp_x_map, mp_g_map[j]));
  dgdx_mp->setupOperator(mp_blocks);
  return dgdx_mp;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_DgDx_dot_op(int j) const
{
  TEST_FOR_EXCEPTION(j < 0 || j >= num_g_mp || !supports_x, 
		     std::logic_error,
		     "Error:  dg/dx_dot index " << j << " is not supported!");
  
  int jj = mp_g_index_map[j];
  Teuchos::RCP<const Epetra_Map> g_map = me->get_g_map(jj);
  Teuchos::RCP< Stokhos::ProductEpetraOperator > mp_blocks =
    Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		   mp_block_map, x_map, g_map, mp_g_map[j], 
		   mp_comm));
  OutArgs me_outargs = me->createOutArgs();
  DerivativeSupport ds = me_outargs.supports(OUT_ARG_DgDx_dot_mp, jj);
  if (ds.supports(DERIV_LINEAR_OP))
    for (unsigned int i=0; i<num_mp_blocks; i++)
      mp_blocks->setCoeffPtr(i, me->create_DgDx_dot_op(jj));
  else if (ds.supports(DERIV_MV_BY_COL))
    for (unsigned int i=0; i<num_mp_blocks; i++) {
      Teuchos::RCP<Epetra_MultiVector> mv = 
	Teuchos::rcp(new Epetra_MultiVector(*g_map, x_map->NumMyElements()));
      Teuchos::RCP<Epetra_Operator> block = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
      mp_blocks->setCoeffPtr(i, block);
    }
  else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
    for (unsigned int i=0; i<num_mp_blocks; i++) {
      Teuchos::RCP<Epetra_MultiVector> mv = 
	Teuchos::rcp(new Epetra_MultiVector(*x_map, g_map->NumMyElements()));
      
      Teuchos::RCP<Epetra_Operator> block = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
      mp_blocks->setCoeffPtr(i, block);
    }
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  me_outargs.supports(OUT_ARG_DgDx_dot_mp, " 
		       << j << ").none() is true!");
  
  Teuchos::RCP<Stokhos::BlockDiagonalOperator> dgdx_dot_mp = 
    Teuchos::rcp(new Stokhos::BlockDiagonalOperator(
		   mp_comm, num_mp_blocks, x_map, g_map,
		   mp_x_map, mp_g_map[j]));
  dgdx_dot_mp->setupOperator(mp_blocks);
  return dgdx_dot_mp;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_DgDp_op(int j, int i) const
{
  TEST_FOR_EXCEPTION(
    j < 0 || j >= num_g_mp || i < 0 || i >= num_p+num_p_mp, 
    std::logic_error,
    "Error:  dg/dp index " << j << "," << i << " is not supported!");

  OutArgs me_outargs = me->createOutArgs();
  int jj = mp_g_index_map[j];
  Teuchos::RCP<const Epetra_Map> g_map = me->get_g_map(jj);
  if (i < num_p) {
    if (me_outargs.supports(OUT_ARG_DgDp_mp,jj,i).supports(DERIV_LINEAR_OP)) {
      Teuchos::RCP<Stokhos::ProductEpetraOperator> mp_blocks =
	Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		       mp_block_map, me->get_p_map(i), g_map, 
		       mp_g_map[j], mp_comm));
      for (unsigned int l=0; l<num_mp_blocks; l++)
	mp_blocks->setCoeffPtr(l, me->create_DgDp_op(i,j));
      return mp_blocks;
    }
    else
      TEST_FOR_EXCEPTION(
	true, std::logic_error,
	"Error:  Underlying model evaluator must support DERIV_LINER_OP " << 
	"to create operator for dg/dp index " << j << "," << i << "!");
  }
  else {
    int ii = mp_p_index_map[i-num_p];
    Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(ii);
    Teuchos::RCP< Stokhos::ProductEpetraOperator> mp_blocks =
      Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		     mp_block_map, p_map, g_map, 
		     mp_g_map[j], mp_comm));
    DerivativeSupport ds = me_outargs.supports(OUT_ARG_DgDp_mp,jj,ii);
    if (ds.supports(DERIV_LINEAR_OP))
      for (unsigned int l=0; l<num_mp_blocks; l++)
	mp_blocks->setCoeffPtr(l, me->create_DgDp_op(jj,ii));
    else if (ds.supports(DERIV_MV_BY_COL))
      for (unsigned int l=0; l<num_mp_blocks; l++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*g_map, p_map->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(l, block);
      }
    else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
      for (unsigned int l=0; l<num_mp_blocks; l++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*p_map, g_map->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(l, block);
      }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error!  me_outargs.supports(OUT_ARG_DgDp_mp, " << jj
			 << "," << ii << ").none() is true!");

    Teuchos::RCP<Stokhos::BlockDiagonalOperator> dgdp_mp = 
    Teuchos::rcp(new Stokhos::BlockDiagonalOperator(
		   mp_comm, num_mp_blocks, p_map, g_map,
		   mp_p_map[i-num_p], mp_g_map[j]));
    dgdp_mp->setupOperator(mp_blocks);
    return dgdp_mp;
  }
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_DfDp_op(int i) const
{
  TEST_FOR_EXCEPTION(i < 0 || i >= num_p+num_p_mp, 
                     std::logic_error,
                     "Error:  df/dp index " << i << " is not supported!");

  OutArgs me_outargs = me->createOutArgs();
  if (i < num_p) {
    if (me_outargs.supports(OUT_ARG_DfDp_mp,i).supports(DERIV_LINEAR_OP)) {
      Teuchos::RCP<Stokhos::ProductEpetraOperator> mp_blocks =
	Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		       mp_block_map, me->get_p_map(i), me->get_f_map(), 
		       mp_f_map, mp_comm));
      for (unsigned int l=0; l<num_mp_blocks; l++)
	mp_blocks->setCoeffPtr(l, me->create_DfDp_op(i));
      return mp_blocks;
    }
    else
      TEST_FOR_EXCEPTION(
	true, std::logic_error,
	"Error:  Underlying model evaluator must support DERIV_LINER_OP " << 
	"to create operator for df/dp index " << i << "!");
  }
  else {
    int ii = mp_p_index_map[i-num_p];
    Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(ii);
    Teuchos::RCP<Stokhos::ProductEpetraOperator> mp_blocks =
	Teuchos::rcp(new Stokhos::ProductEpetraOperator(
		       mp_block_map, me->get_p_map(ii), me->get_f_map(), 
		       mp_f_map, mp_comm));
    DerivativeSupport ds = me_outargs.supports(OUT_ARG_DfDp_mp,ii);
    if (ds.supports(DERIV_LINEAR_OP))
      for (unsigned int l=0; l<num_mp_blocks; l++)
	mp_blocks->setCoeffPtr(l, me->create_DfDp_op(ii));
    else if (ds.supports(DERIV_MV_BY_COL))
      for (unsigned int l=0; l<num_mp_blocks; l++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*f_map, p_map->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(l, block);
      }
    else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
      for (unsigned int l=0; l<num_mp_blocks; l++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*p_map, f_map->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(l, block);
      }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error!  me_outargs.supports(OUT_ARG_DfDp_mp, " << ii 
			 << ").none() is true!");

    
    Teuchos::RCP<Stokhos::BlockDiagonalOperator> dfdp_mp = 
      Teuchos::rcp(new Stokhos::BlockDiagonalOperator(
		     mp_comm, num_mp_blocks, 
		     p_map, f_map, mp_p_map[i-num_p], mp_f_map));
    dfdp_mp->setupOperator(mp_blocks);
    return dfdp_mp;
  }
  
  return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::MPModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p + num_p_mp); 
  inArgs.setSupports(IN_ARG_x_dot, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_t, me_inargs.supports(IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha, me_inargs.supports(IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta, me_inargs.supports(IN_ARG_beta));
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::MPModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(num_p + num_p_mp, num_g_mp);
  outArgs.setSupports(OUT_ARG_f, me_outargs.supports(OUT_ARG_f));
  outArgs.setSupports(OUT_ARG_W, me_outargs.supports(OUT_ARG_W));
  outArgs.setSupports(OUT_ARG_WPrec, me_outargs.supports(OUT_ARG_W));
  for (int j=0; j<num_p; j++)
    outArgs.setSupports(OUT_ARG_DfDp, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  for (int j=0; j<num_p_mp; j++)
    if (!me_outargs.supports(OUT_ARG_DfDp_mp, mp_p_index_map[j]).none())
      outArgs.setSupports(OUT_ARG_DfDp, j+num_p, DERIV_LINEAR_OP);
  for (int i=0; i<num_g_mp; i++) {
    int ii = mp_g_index_map[i];
    if (!me_outargs.supports(OUT_ARG_DgDx_dot_mp, ii).none())
      outArgs.setSupports(OUT_ARG_DgDx_dot, i,  DERIV_LINEAR_OP);
    if (!me_outargs.supports(OUT_ARG_DgDx_mp, ii).none())
      outArgs.setSupports(OUT_ARG_DgDx, i,  DERIV_LINEAR_OP);
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp_mp, ii, j));
    for (int j=0; j<num_p_mp; j++)
      if (!me_outargs.supports(OUT_ARG_DgDp_mp, ii, mp_p_index_map[j]).none())
	outArgs.setSupports(OUT_ARG_DgDp, i, j+num_p, DERIV_LINEAR_OP);
  }
  
  return outArgs;
}

void 
Stokhos::MPModelEvaluator::evalModel(const InArgs& inArgs,
				     const OutArgs& outArgs) const
{
  // Get the input arguments
  Teuchos::RCP<const Epetra_Vector> x;
  if (inArgs.supports(IN_ARG_x)) {
    x = inArgs.get_x();
    if (x != Teuchos::null)
      *my_x = *x;
  }
  Teuchos::RCP<const Epetra_Vector> x_dot;
  if (inArgs.supports(IN_ARG_x_dot))
    x_dot = inArgs.get_x_dot();

  // Get the output arguments
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f_out;
  if (outArgs.supports(OUT_ARG_f))
    f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out;
  if (outArgs.supports(OUT_ARG_W))
    W_out = outArgs.get_W();
  Teuchos::RCP<Epetra_Operator> WPrec_out;
  if (outArgs.supports(OUT_ARG_WPrec))
    WPrec_out = outArgs.get_WPrec();

  // Check if we are using the "matrix-free" method for W and we are 
  // computing a preconditioner.  
  bool eval_prec = (W_out == Teuchos::null && WPrec_out != Teuchos::null);

  // Here we are assuming a full W fill occurred previously which we can use
  // for the preconditioner.  Given the expense of computing the MP W blocks
  // this saves significant computational cost
  if (eval_prec) {
    Teuchos::RCP<Stokhos::MPPreconditioner> W_prec = 
      Teuchos::rcp_dynamic_cast<Stokhos::MPPreconditioner>(WPrec_out, true);
    W_prec->setupPreconditioner(my_W, *my_x);
    
    // We can now quit unless a fill of f, g, or dg/dp was also requested
    bool done = (f_out == Teuchos::null);
    for (int i=0; i<outArgs.Ng(); i++) {
      done = done && (outArgs.get_g(i) == Teuchos::null);
      for (int j=0; j<outArgs.Np(); j++)
	if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	done = done && (outArgs.get_DgDp(i,j).isEmpty());
    }
    if (done)
      return;
  }

  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();
  if (x != Teuchos::null) {
   Teuchos::RCP<Stokhos::ProductEpetraVector> x_mp = 
      create_x_mp(View, x.get());
    me_inargs.set_x_mp(x_mp);
  }
  if (x_dot != Teuchos::null) {
    Teuchos::RCP<Stokhos::ProductEpetraVector> x_dot_mp = 
      create_x_mp(View, x_dot.get());
    me_inargs.set_x_dot_mp(x_dot_mp);
  }
  if (me_inargs.supports(IN_ARG_alpha))
    me_inargs.set_alpha(inArgs.get_alpha());
  if (me_inargs.supports(IN_ARG_beta))
    me_inargs.set_beta(inArgs.get_beta());
  if (me_inargs.supports(IN_ARG_t))
    me_inargs.set_t(inArgs.get_t());

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));
  for (int i=0; i<num_p_mp; i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i+num_p);

    // We always need to pass in the MP parameters, so just use
    // the initial parameters if it is NULL
    if (p == Teuchos::null)
      p = mp_p_init[i]->getBlockVector();

    // Convert block p to MP polynomial
    Teuchos::RCP<Stokhos::ProductEpetraVector> p_mp = 
      create_p_mp(mp_p_index_map[i], View, p.get());
    me_inargs.set_p_mp(mp_p_index_map[i], p_mp);
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // f
  if (f_out != Teuchos::null) {
    Teuchos::RCP<Stokhos::ProductEpetraVector> f_mp = 
      create_f_mp(View, f_out.get());
    me_outargs.set_f_mp(f_mp);
  }

  // W
  if (W_out != Teuchos::null && !eval_prec)
     me_outargs.set_W_mp(W_mp_blocks);

  // df/dp -- scalar p
  for (int i=0; i<num_p; i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dfdp_mp; 
	if (dfdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	  dfdp_mp = 
	    Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			   mp_block_map, me->get_f_map(), mp_f_map, mp_comm, 
			   View, *(dfdp.getMultiVector())));
	else if (dfdp.getMultiVectorOrientation() == DERIV_TRANS_MV_BY_ROW)
	  dfdp_mp = 
	    Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			   mp_block_map, me->get_p_map(i), mp_p_map[i], mp_comm,
			   View, *(dfdp.getMultiVector())));
	me_outargs.set_DfDp_mp(i, 
			       MPDerivative(dfdp_mp,
					    dfdp.getMultiVectorOrientation()));
      }
      else if (dfdp.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::ProductEpetraOperator> dfdp_mp = 
	  Teuchos::rcp_dynamic_cast<Stokhos::ProductEpetraOperator>(dfdp.getLinearOp(), true);
	me_outargs.set_DfDp_mp(i, MPDerivative(dfdp_mp));
      }
    }
  }

  // dfdp -- block p.  Here we only support DERIV_LINEAR_OP
  for (int i=0; i<num_p_mp; i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i+num_p).none()) {
      Derivative dfdp = outArgs.get_DfDp(i+num_p);
      if (dfdp.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::BlockDiagonalOperator> dfdp_op =
	  Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(dfdp.getLinearOp(), true);
	Teuchos::RCP<Stokhos::ProductEpetraOperator> dfdp_op_mp = 
	  dfdp_op->getMPOps();
	int ii = mp_p_index_map[i];
	if (me_outargs.supports(OUT_ARG_DfDp_mp,ii).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DfDp_mp(ii, MPDerivative(dfdp_op_mp));
	else {
	  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dfdp_mp = 
	    Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(mp_block_map));
	  for (unsigned int l=0; l<num_mp_blocks; l++) {
	    Teuchos::RCP<Stokhos::EpetraMultiVectorOperator> mv_op =
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(dfdp_op_mp->getCoeffPtr(i), true);
	    dfdp_mp->setCoeffPtr(l, mv_op->getMultiVector());
	  }
	  if (me_outargs.supports(OUT_ARG_DfDp_mp,ii).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DfDp_mp(
	      ii, MPDerivative(dfdp_mp, DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DfDp_mp(
	      ii, MPDerivative(dfdp_mp, DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(
	dfdp.getLinearOp() == Teuchos::null && dfdp.isEmpty() == false, 
	std::logic_error,
	"Error!  Stokhos::MPModelEvaluator::evalModel: " << 
	"Operator form of df/dp(" << i+num_p << ") is required!");
    }
  }

  // Responses (g, dg/dx, dg/dp, ...)
  for (int i=0; i<num_g_mp; i++) {
    int ii = mp_g_index_map[i];

    // g
    Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(i);
    if (g != Teuchos::null) {
      Teuchos::RCP<Stokhos::ProductEpetraVector> g_mp =
	create_g_mp(ii, View, g.get());
      me_outargs.set_g_mp(ii, g_mp);
    }

    // dg/dxdot
    if (outArgs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx_dot = outArgs.get_DgDx_dot(i);
      if (dgdx_dot.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::BlockDiagonalOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(
	    dgdx_dot.getLinearOp(), true);
	Teuchos::RCP<Stokhos::ProductEpetraOperator> mp_blocks =
	  op->getMPOps();
	if (me_outargs.supports(OUT_ARG_DgDx, ii).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_dot_mp(ii, mp_blocks);
	else {
	  for (unsigned int k=0; k<num_mp_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		mp_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_dot_mp_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_dot_mp, ii).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_dot_mp(ii, MPDerivative(dgdx_dot_mp_blocks[i],
							DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_dot_mp(ii, MPDerivative(dgdx_dot_mp_blocks[i],
							DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(dgdx_dot.getLinearOp() == Teuchos::null &&
			 dgdx_dot.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::MPModelEvaluator::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dx
    if (outArgs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx = outArgs.get_DgDx(i);
      if (dgdx.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::BlockDiagonalOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(
	    dgdx.getLinearOp(), true);
	Teuchos::RCP<Stokhos::ProductEpetraOperator> mp_blocks =
	  op->getMPOps();
	if (me_outargs.supports(OUT_ARG_DgDx, ii).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_mp(ii, mp_blocks);
	else {
	  for (unsigned int k=0; k<num_mp_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		mp_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_mp_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_mp, ii).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_mp(ii, MPDerivative(dgdx_mp_blocks[i],
						    DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_mp(ii, MPDerivative(dgdx_mp_blocks[i],
						    DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(dgdx.getLinearOp() == Teuchos::null &&
			 dgdx.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::MPModelEvaluator::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dp -- scalar p
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	Derivative dgdp = outArgs.get_DgDp(i,j);
	if (dgdp.getMultiVector() != Teuchos::null) {
	  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dgdp_mp; 
	  if (dgdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	    dgdp_mp = 
	      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			     mp_block_map, me->get_g_map(ii), mp_g_map[i], 
			     mp_comm, View, *(dgdp.getMultiVector())));
	  else if (dgdp.getMultiVectorOrientation() == DERIV_TRANS_MV_BY_ROW)
	    dgdp_mp = 
	      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			     mp_block_map, me->get_p_map(j), mp_p_map[j], 
			     mp_comm, View, *(dgdp.getMultiVector())));
	  me_outargs.set_DgDp_mp(
	    ii, j, MPDerivative(dgdp_mp, dgdp.getMultiVectorOrientation()));
	}
	else if (dgdp.getLinearOp() != Teuchos::null) {
	  Teuchos::RCP<Stokhos::ProductEpetraOperator> dgdp_mp = 
	    Teuchos::rcp_dynamic_cast<Stokhos::ProductEpetraOperator>(dgdp.getLinearOp(), true);
	  me_outargs.set_DgDp_mp(ii, j, MPDerivative(dgdp_mp));
	}
      }
    }

    // dgdp -- block p.  Here we only support DERIV_LINEAR_OP
    for (int j=0; j<num_p_mp; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i, j+num_p).none()) {
	Derivative dgdp = outArgs.get_DgDp(i,j+num_p);
	if (dgdp.getLinearOp() != Teuchos::null) {
	  Teuchos::RCP<Stokhos::BlockDiagonalOperator> dgdp_op =
	    Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(dgdp.getLinearOp(), true);
	  Teuchos::RCP<Stokhos::ProductEpetraOperator> dgdp_op_mp =
	    dgdp_op->getMPOps();
	  int jj = mp_p_index_map[j];
	  if (me_outargs.supports(OUT_ARG_DgDp_mp,ii,jj).supports(DERIV_LINEAR_OP))
	    me_outargs.set_DgDp_mp(ii, jj, MPDerivative(dgdp_op_mp));
	  else {
	    Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dgdp_mp = 
	      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(mp_block_map));
	    for (unsigned int l=0; l<num_mp_blocks; l++) {
	      Teuchos::RCP<Stokhos::EpetraMultiVectorOperator> mv_op =
		Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(dgdp_op_mp->getCoeffPtr(i), true);
	      dgdp_mp->setCoeffPtr(l, mv_op->getMultiVector());
	    }
	    if (me_outargs.supports(OUT_ARG_DgDp_mp,ii,jj).supports(DERIV_MV_BY_COL))
	      me_outargs.set_DgDp_mp(
		ii, jj, MPDerivative(dgdp_mp, DERIV_MV_BY_COL));
	    else
	      me_outargs.set_DgDp_mp(
		ii, jj, MPDerivative(dgdp_mp, DERIV_TRANS_MV_BY_ROW));
	  }
	}
	TEST_FOR_EXCEPTION(
	  dgdp.getLinearOp() == Teuchos::null && dgdp.isEmpty() == false, 
	  std::logic_error,
	  "Error!  Stokhos::MPModelEvaluator::evalModel: " << 
	  "Operator form of dg/dp(" << i << "," << j+num_p << ") is required!");
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Copy block MP components for W
  if (W_out != Teuchos::null && !eval_prec) {
    Teuchos::RCP<Epetra_Operator> W;
    if (W_out != Teuchos::null)
      W = W_out;
    else
      W = my_W;
    Teuchos::RCP<Stokhos::BlockDiagonalOperator> W_mp = 
      Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(W, true);
    W_mp->setupOperator(W_mp_blocks);
      
    if (WPrec_out != Teuchos::null) {
      Teuchos::RCP<Stokhos::MPPreconditioner> W_prec = 
	Teuchos::rcp_dynamic_cast<Stokhos::MPPreconditioner>(WPrec_out, true);
      W_prec->setupPreconditioner(W_mp, *my_x);
    }
  }
}

void 
Stokhos::MPModelEvaluator::set_x_mp_init(
  const Stokhos::ProductEpetraVector& x_mp_in)
{
  *mp_x_init = x_mp_in;
}

void 
Stokhos::MPModelEvaluator::set_p_mp_init(
  int i, const Stokhos::ProductEpetraVector& p_mp_in)
{
  *mp_p_init[i] = p_mp_in;
}

Teuchos::Array<int> 
Stokhos::MPModelEvaluator::get_p_mp_map_indices() const
{
  return mp_p_index_map;
}

Teuchos::Array<int> 
Stokhos::MPModelEvaluator::get_g_mp_map_indices() const
{
  return mp_g_index_map;
}

Teuchos::Array< Teuchos::RCP<const Epetra_Map> > 
Stokhos::MPModelEvaluator::get_g_mp_base_maps() const
{
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_maps(num_g);
  for (int i=0; i<num_g; i++)
    base_maps[i] = me->get_g_map(i);
  return base_maps;
 }

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_x_mp(Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_x;
  if (v == NULL)
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, x_map, mp_x_map, mp_comm));
  else
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, x_map, mp_x_map, mp_comm,
			  CV, *v));
  return mp_x;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector>
Stokhos::MPModelEvaluator::create_x_mv_mp(int num_vecs, Epetra_DataAccess CV, 
					  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_x;
  if (v == NULL)
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, x_map, mp_x_map, mp_comm,
			  num_vecs));
  else
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, x_map, mp_x_map, mp_comm,
			  CV, *v));
  return mp_x;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_p_mp(int l, Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_p;
  Teuchos::Array<int>::const_iterator it = std::find(mp_p_index_map.begin(),
						     mp_p_index_map.end(), 
						     l);
  TEST_FOR_EXCEPTION(it == mp_p_index_map.end(), std::logic_error,
		     "Error!  Invalid p map index " << l);
  int ll = it - mp_p_index_map.begin();
  if (v == NULL)
    mp_p = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, me->get_p_map(l),
			  mp_p_map[ll], mp_comm));
  else
    mp_p = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, me->get_p_map(l),
			  mp_p_map[ll], mp_comm, CV, *v));
  return mp_p;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_f_mp(Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_f;
  if (v == NULL)
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, f_map, mp_f_map, mp_comm));
  else
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, f_map, mp_f_map, mp_comm,
			  CV, *v));
  return mp_f;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector>
Stokhos::MPModelEvaluator::create_f_mv_mp(
  int num_vecs, 
  Epetra_DataAccess CV, 
  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_f;
  if (v == NULL)
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, f_map, mp_f_map, mp_comm,
			  num_vecs));
  else
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, f_map, mp_f_map, mp_comm,
			  CV, *v));
  return mp_f;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_g_mp(int l, Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_g;
  Teuchos::Array<int>::const_iterator it = std::find(mp_g_index_map.begin(),
						     mp_g_index_map.end(), 
						     l);
  TEST_FOR_EXCEPTION(it == mp_g_index_map.end(), std::logic_error,
		     "Error!  Invalid g map index " << l);
  int ll = it - mp_g_index_map.begin();
  if (v == NULL)
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, 
			  me->get_g_map(l), 
			  mp_g_map[ll], mp_comm));
  else
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, 
			  me->get_g_map(l), 
			  mp_g_map[ll], mp_comm, CV, *v));
  return mp_g;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector>
Stokhos::MPModelEvaluator::create_g_mv_mp(int l, int num_vecs,
					  Epetra_DataAccess CV, 
					  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_g;
  Teuchos::Array<int>::const_iterator it = std::find(mp_g_index_map.begin(),
						     mp_g_index_map.end(), 
						     l);
  TEST_FOR_EXCEPTION(it == mp_g_index_map.end(), std::logic_error,
		     "Error!  Invalid g map index " << l);
  int ll = it - mp_g_index_map.begin();
  if (v == NULL)
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, 
			  me->get_g_map(l), 
			  mp_g_map[ll], mp_comm, num_vecs));
  else
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, 
			  me->get_g_map(l), 
			  mp_g_map[ll], mp_comm, CV, *v));
  return mp_g;
}
