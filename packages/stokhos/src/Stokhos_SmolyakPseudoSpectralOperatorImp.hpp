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
template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
template <typename coeff_compare_type, typename point_growth_rule_type>
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>:: 
SmolyakPseudoSpectralOperator(
  const SmolyakBasis<ordinal_type,value_type,coeff_compare_type>& smolyak_basis, 
  const point_growth_rule_type& point_growth_rule,
  bool use_smolyak_apply,
  bool use_pst,
  const point_compare_type& point_compare) :
  use_smolyak(use_smolyak_apply),
  points(point_compare) {

  typedef SmolyakBasis<ordinal_type,value_type,coeff_compare_type> smolyak_basis_type;
  typedef typename smolyak_basis_type::tensor_product_basis_type tensor_product_basis_type;

  // Generate sparse grid and tensor operators
  coeff_sz = smolyak_basis.size();
  ordinal_type dim = smolyak_basis.dimension();
  ordinal_type num_terms = smolyak_basis.getNumSmolyakTerms();
  for (ordinal_type i=0; i<num_terms; ++i) {
    
    // Get tensor product basis for given term
    Teuchos::RCP<const tensor_product_basis_type> tp_basis = 
      smolyak_basis.getTensorProductBasis(i);
    
    // Get coefficient multi-index defining basis orders
    multiindex_type coeff_index = tp_basis->getMaxOrders();
    
    // Apply growth rule to cofficient multi-index
    multiindex_type point_growth_index(dim);
    for (ordinal_type j=0; j<dim; ++j) {
      point_growth_index[j] = point_growth_rule[j](coeff_index[j]);
    }
    
    // Build tensor product operator for given index
    Teuchos::RCP<operator_type> op = 
      Teuchos::rcp(new operator_type(*tp_basis, use_pst,
				     point_growth_index));
    if (use_smolyak)
      operators.push_back(op);
    
    // Get smolyak cofficient for given index
    value_type c = smolyak_basis.getSmolyakCoefficient(i);
    if (use_smolyak)
      smolyak_coeffs.push_back(c);
    
    // Include points in union over all sets
    typename operator_type::set_iterator op_set_iterator = op->set_begin();
    typename operator_type::set_iterator op_set_end = op->set_end();
    for (; op_set_iterator != op_set_end; ++op_set_iterator) {
      const point_type& point = op_set_iterator->first;
      value_type w = op_set_iterator->second.first;
      set_iterator si = points.find(point);
      if (si == points.end())
	points[point] = std::make_pair(c*w,ordinal_type(0));
      else
	si->second.first += c*w;
    }
    
  }
  
  // Generate linear ordering of points
  ordinal_type idx = 0;
  point_map.resize(points.size());
  for (set_iterator si = points.begin(); si != points.end(); ++si) {
    si->second.second = idx;
    point_map[idx] = si->first;
    ++idx;
  }  
  
  if (use_smolyak) {
    
    // Build gather/scatter maps into global domain/range for each operator
    gather_maps.resize(operators.size());
    scatter_maps.resize(operators.size());
    for (ordinal_type i=0; i<operators.size(); i++) {
      Teuchos::RCP<operator_type> op = operators[i];
      
      gather_maps[i].reserve(op->point_size());
      typename operator_type::iterator op_iterator = op->begin();
      typename operator_type::iterator op_end = op->end();
      for (; op_iterator != op_end; ++op_iterator) {
	gather_maps[i].push_back(points[*op_iterator].second);
      }
      
      Teuchos::RCP<const tensor_product_basis_type> tp_basis = 
	smolyak_basis.getTensorProductBasis(i);
      ordinal_type op_coeff_sz = tp_basis->size();
      scatter_maps[i].reserve(op_coeff_sz);
      for (ordinal_type j=0; j<op_coeff_sz; ++j) {
	scatter_maps[i].push_back(smolyak_basis.index(tp_basis->term(j)));
      }
    }
  }
      
  //else {
    
    // Generate quadrature operator
    ordinal_type nqp = points.size();
    ordinal_type npc = coeff_sz;
    qp2pce.reshape(npc,nqp);
    pce2qp.reshape(nqp,npc);
    qp2pce.putScalar(1.0);
    pce2qp.putScalar(1.0);
    Teuchos::Array<value_type> vals(npc);
    for (set_iterator si = points.begin(); si != points.end(); ++si) {
      ordinal_type j = si->second.second;
      value_type w = si->second.first;
      point_type point = si->first;
      smolyak_basis.evaluateBases(point, vals);
      for (ordinal_type i=0; i<npc; ++i) {
	qp2pce(i,j) = w*vals[i] / smolyak_basis.norm_squared(i);
	pce2qp(j,i) = vals[i];
      }
    }
    //}
  
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
transformQP2PCE(
  const value_type& alpha, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
  const value_type& beta,
  bool trans) const {
  
  if (use_smolyak)
    transformQP2PCE_smolyak(alpha, input, result, beta, trans);
  else
    apply_direct(qp2pce, alpha, input, result, beta, trans);
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
transformPCE2QP(
  const value_type& alpha, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
  const value_type& beta,
  bool trans) const {

  // Currently we use the direct method for mapping PCE->QP because the
  // current implementation doesn't work.  Need to evaluate tensor bases
  // on all quad points, not just the quad points associated with that 
  // basis.
  
  // if (use_smolyak)
  //   transformPCE2QP_smolyak(alpha, input, result, beta, trans);
  // else
    apply_direct(pce2qp, alpha, input, result, beta, trans);
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
apply_direct(
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& A,
  const value_type& alpha, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
  const value_type& beta,
  bool trans) const {
  if (trans) {
    TEUCHOS_ASSERT(input.numCols() <= A.numCols());
    TEUCHOS_ASSERT(result.numCols() == A.numRows());
    TEUCHOS_ASSERT(result.numRows() == input.numRows());
    blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, input.numRows(), 
	      A.numRows(), input.numCols(), alpha, input.values(), 
	      input.stride(), A.values(), A.stride(), beta, 
	      result.values(), result.stride());
  }
  else {
    TEUCHOS_ASSERT(input.numRows() <= A.numCols());
    TEUCHOS_ASSERT(result.numRows() == A.numRows());
    TEUCHOS_ASSERT(result.numCols() == input.numCols());
    blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, A.numRows(), 
	      input.numCols(), input.numRows(), alpha, A.values(), 
	      A.stride(), input.values(), input.stride(), beta, 
	      result.values(), result.stride());
  }
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
transformQP2PCE_smolyak(
  const value_type& alpha, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
  const value_type& beta,
  bool trans) const {
  Teuchos::SerialDenseMatrix<ordinal_type,value_type> op_input, op_result;
  result.scale(beta);
  for (ordinal_type i=0; i<operators.size(); i++) {
    Teuchos::RCP<operator_type> op = operators[i];
    if (trans) {
      op_input.reshape(input.numRows(), op->point_size());
      op_result.reshape(result.numRows(), op->coeff_size());
    }
    else {
      op_input.reshape(op->point_size(), input.numCols());
      op_result.reshape(op->coeff_size(), result.numCols());
    }
    gather(gather_maps[i], input, trans, op_input);
    op->transformQP2PCE(smolyak_coeffs[i], op_input, op_result, 0.0, trans);
    scatter(scatter_maps[i], op_result, trans, result);
  }
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
transformPCE2QP_smolyak(
  const value_type& alpha, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
  const value_type& beta,
  bool trans) const {
  Teuchos::SerialDenseMatrix<ordinal_type,value_type> op_input, op_result;
  result.scale(beta);

  for (ordinal_type i=0; i<operators.size(); i++) {
    Teuchos::RCP<operator_type> op = operators[i];
    if (trans) {
      op_input.reshape(input.numRows(), op->coeff_size());
      op_result.reshape(result.numRows(), op->point_size());
    }
    else {
      op_input.reshape(op->coeff_size(), input.numCols());
      op_result.reshape(op->point_size(), result.numCols());
    }
    
    gather(scatter_maps[i], input, trans, op_input);
    op->transformPCE2QP(smolyak_coeffs[i], op_input, op_result, 0.0, trans);
    scatter(gather_maps[i], op_result, trans, result);
  }
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
gather(
  const Teuchos::Array<ordinal_type>& map, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input, 
  bool trans, 
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result) const {
  if (trans) {
    for (ordinal_type j=0; j<map.size(); j++)
      for (ordinal_type i=0; i<input.numRows(); i++)
	result(i,j) = input(i,map[j]);
  }
  else {
    for (ordinal_type j=0; j<input.numCols(); j++)
      for (ordinal_type i=0; i<map.size(); i++)
	result(i,j) = input(map[i],j);
  }
}

template <typename ordinal_type, typename value_type, 
	  typename point_compare_type>
void
Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type>::
scatter(
  const Teuchos::Array<ordinal_type>& map, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input, 
  bool trans, 
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result) const {
  if (trans) {
    for (ordinal_type j=0; j<map.size(); j++)
      for (ordinal_type i=0; i<input.numRows(); i++)
	result(i,map[j]) += input(i,j);
  }
  else {
    for (ordinal_type j=0; j<input.numCols(); j++)
      for (ordinal_type i=0; i<map.size(); i++)
	result(map[i],j) += input(i,j);
  }
}
