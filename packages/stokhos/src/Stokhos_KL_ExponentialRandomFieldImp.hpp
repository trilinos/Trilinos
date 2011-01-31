// $Id: Stokhos_Quadrature.hpp,v 1.4 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_Quadrature.hpp,v $ // @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <sstream>
#include <cmath>
#include <algorithm>
#include "Stokhos_KL_OneDExponentialCovarianceFunction.hpp"

template <typename value_type>
Stokhos::KL::ExponentialRandomField<value_type>::
ExponentialRandomField(Teuchos::ParameterList& solverParams)
{
  // Get required parameters
  num_KL = solverParams.get<int>("Number of KL Terms");
  mean = solverParams.get<double>("Mean");
  std_dev = solverParams.get<double>("Standard Deviation");

  Teuchos::Array<double> domain_upper_bound_double;
  Teuchos::Array<double> domain_lower_bound_double;
  Teuchos::Array<double> correlation_length_double;

  // Get Domain Upper Bounds
  if (solverParams.isType<std::string>("Domain Upper Bounds")) 
    domain_upper_bound_double =
      Teuchos::getArrayFromStringParameter<double>(
	solverParams, "Domain Upper Bounds");
  else
    domain_upper_bound_double = 
      solverParams.get< Teuchos::Array<double> >("Domain Upper Bounds");

  // (Convert each element from double to value_type)
  domain_upper_bound.resize(domain_upper_bound_double.size());
  for (int i=0; i<domain_upper_bound.size(); i++)
    domain_upper_bound[i]=domain_upper_bound_double[i];
    
  // Get Domain Lower Bounds
  if (solverParams.isType<std::string>("Domain Lower Bounds")) 
    domain_lower_bound_double =
      Teuchos::getArrayFromStringParameter<double>(
	solverParams, "Domain Lower Bounds");
  else
    domain_lower_bound_double = 
      solverParams.get< Teuchos::Array<double> >("Domain Lower Bounds");

  // (Convert each element from double to value_type)
  domain_lower_bound.resize(domain_lower_bound_double.size());
  for (int i=0; i<domain_lower_bound.size(); i++)
    domain_lower_bound[i]=domain_lower_bound_double[i];
    
  // Get Correlation Lengths
  if (solverParams.isType<std::string>("Correlation Lengths")) 
    correlation_length_double =
      Teuchos::getArrayFromStringParameter<double>(
	solverParams, "Correlation Lengths");
  else
    correlation_length_double = 
      solverParams.get< Teuchos::Array<double> >("Correlation Lengths");

  // (Convert each element from double to value_type)
  correlation_length.resize(correlation_length_double.size());
  for (int i=0; i<correlation_length.size(); i++)
    correlation_length[i]=correlation_length_double[i];
    
  // Compute 1-D eigenfunctions for each dimension
  dim = domain_upper_bound.size();
  Teuchos::Array< Teuchos::Array< OneDEigenPair<value_type> > > eig_pairs(dim);
  for (int i=0; i<dim; i++) {
    std::stringstream ss;
    ss << "x_" << i;
    eig_pairs[i].resize(num_KL);
    OneDExponentialCovarianceFunction<value_type> cov_func(
      num_KL, domain_lower_bound[i], domain_upper_bound[i],
      correlation_length[i], ss.str(), solverParams);
    eig_pairs[i] = cov_func.getEigenPairs();
  }

  // Compute all possible tensor product combinations of 1-D eigenfunctions
  int num_prod = static_cast<int>(std::pow(static_cast<double>(num_KL), 
					   static_cast<double>(dim)));
  product_eig_pairs.resize(num_prod);
  Teuchos::Array<int> index(dim, 0);
  int cnt = 0;
  Teuchos::Array<value_type> vals(dim);
  Teuchos::Array< OneDEigenPair<value_type> > eigs(dim);
  while (cnt < num_prod) {
    for (int i=0; i<dim; i++) {
      eigs[i] = eig_pairs[i][index[i]];
    }
    product_eig_pairs[cnt] = ProductEigenPair<value_type>(eigs);
    ++index[0];
    int j = 0;
    while (j < dim-1 && index[j] == num_KL) {
      index[j] = 0;
      ++j;
      ++index[j];
    }
    ++cnt;
  }

  // Sort product eigenfunctions based on product eigenvalue
  std::sort(product_eig_pairs.begin(), product_eig_pairs.end(),
	    ProductEigenPairGreater<value_type>());
}

template <typename value_type>
template <typename rvar_type>
typename Teuchos::PromotionTraits<rvar_type, value_type>::promote
Stokhos::KL::ExponentialRandomField<value_type>::
evaluate(const Teuchos::Array<value_type>& point,
	 const Teuchos::Array<rvar_type>& random_variables) const
{
  typedef typename Teuchos::PromotionTraits<rvar_type, value_type>::promote result_type;
  result_type result = 0.0;
  for (int i=0; i<num_KL; i++) {
    result += 
      random_variables[i]*(std::sqrt(product_eig_pairs[i].eig_val)*
			   product_eig_pairs[i].evalEigenfunction(point));
  }
  result = mean + std_dev*result;
  return result;
}

template <typename value_type>
value_type
Stokhos::KL::ExponentialRandomField<value_type>::
evaluate_eigenfunction(const Teuchos::Array<value_type>& point, int i) const
{
  return std_dev*std::sqrt(product_eig_pairs[i].eig_val)*
    product_eig_pairs[i].evalEigenfunction(point);
}

template <typename value_type>
const Teuchos::Array< Stokhos::KL::ProductEigenPair<value_type> >&
Stokhos::KL::ExponentialRandomField<value_type>::
getEigenPairs() const
{
  return product_eig_pairs;
}

template <typename value_type>
void
Stokhos::KL::ExponentialRandomField<value_type>::
print(std::ostream& os) const 
{
  os << "KL expansion using " << num_KL << " terms in " << dim 
     << " dimensions:" << std::endl;
  os << "\tDomain = ";
  for (int i=0; i<dim-1; i++)
    os << "[" << domain_lower_bound[i] << "," << domain_upper_bound[i] 
       << "] x ";
  os << "[" << domain_lower_bound[dim-1] << "," << domain_upper_bound[dim-1] 
     << "]" << std::endl;
  os << "\tCorrelation lengths = ";
  for (int i=0; i<dim-1; i++)
    os << correlation_length[i] << ", ";
  os << correlation_length[dim-1] << std::endl;
  os << "\tMean = " << mean << ", standard deviation = " << std_dev 
     << std::endl;
  os << "\tEigenvalues, Eigenfunctions:" << std::endl;
  for (int i=0; i<num_KL; i++)
    os << "\t\t" << product_eig_pairs[i] << std::endl;
}


