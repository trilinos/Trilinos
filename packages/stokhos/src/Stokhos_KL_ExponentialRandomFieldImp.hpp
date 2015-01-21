// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include <cmath>
#include <algorithm>
#include "Teuchos_Array.hpp"

#include "Stokhos_KL_OneDExponentialCovarianceFunction.hpp"

template <typename value_type, typename device_type>
Stokhos::KL::ExponentialRandomField<value_type,device_type>::
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
  Teuchos::Array<value_type> domain_upper_bound(domain_upper_bound_double.size());
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
  Teuchos::Array<value_type> domain_lower_bound(domain_lower_bound_double.size());
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
  Teuchos::Array<value_type> correlation_length(correlation_length_double.size());
  for (int i=0; i<correlation_length.size(); i++)
    correlation_length[i]=correlation_length_double[i];

  // Compute 1-D eigenfunctions for each dimension
  dim = domain_upper_bound.size();
  Teuchos::Array< Teuchos::Array< one_d_eigen_pair_type > > eig_pairs(dim);
  for (int i=0; i<dim; i++) {
    eig_pairs[i].resize(num_KL);
    OneDExponentialCovarianceFunction<value_type> cov_func(
      num_KL, domain_lower_bound[i], domain_upper_bound[i],
      correlation_length[i], i, solverParams);
    eig_pairs[i] = cov_func.getEigenPairs();
  }

  // Compute all possible tensor product combinations of 1-D eigenfunctions
  int num_prod = static_cast<int>(std::pow(static_cast<double>(num_KL),
                                           static_cast<double>(dim)));
  Teuchos::Array<product_eigen_pair_type> product_eig_pairs(num_prod);
  Teuchos::Array<int> index(dim, 0);
  int cnt = 0;
  Teuchos::Array<value_type> vals(dim);
  Teuchos::Array<one_d_eigen_pair_type> eigs(dim);
  while (cnt < num_prod) {
    for (int i=0; i<dim; i++) {
      eigs[i] = eig_pairs[i][index[i]];
    }
    product_eig_pairs[cnt].set(eigs);
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
            ProductEigenPairGreater<one_d_eigen_func_type, device_type>());

  // Copy eigenpairs into view
  product_eigen_funcs =
    eigen_func_array_type("product eigen functions", num_prod, dim);
  product_eigen_values =
    eigen_value_array_type("product eigen vvalues", num_prod);
  typename eigen_func_array_type::HostMirror host_product_eigen_funcs =
    ::Kokkos::create_mirror_view(product_eigen_funcs);
  typename eigen_value_array_type::HostMirror host_product_eigen_values =
    ::Kokkos::create_mirror_view(product_eigen_values);
  for (int i=0; i<num_prod; ++i) {
    host_product_eigen_values(i) = 1.0;
    for (int j=0; j<dim; ++j) {
      host_product_eigen_values(i) *= product_eig_pairs[i].eig_pairs[j].eig_val;
      host_product_eigen_funcs(i,j) =
        product_eig_pairs[i].eig_pairs[j].eig_func;
    }
  }
  ::Kokkos::deep_copy(product_eigen_funcs, host_product_eigen_funcs);
  ::Kokkos::deep_copy(product_eigen_values, host_product_eigen_values);
}

template <typename value_type, typename device_type>
template <typename point_type, typename rv_type>
KOKKOS_INLINE_FUNCTION
typename Teuchos::PromotionTraits<typename rv_type::value_type, value_type>::promote
Stokhos::KL::ExponentialRandomField<value_type,device_type>::
evaluate(const point_type& point, const rv_type& random_variables) const
{
  typedef typename Teuchos::PromotionTraits<typename rv_type::value_type, value_type>::promote result_type;
  result_type result = mean;
  for (int i=0; i<num_KL; ++i) {
    value_type t = std_dev*std::sqrt(product_eigen_values(i));
    for (int j=0; j<dim; ++j)
      t *= product_eigen_funcs(i,j).evaluate(point[j]);
    result += random_variables[i]*t;
  }
  return result;
}

template <typename value_type, typename device_type>
template <typename point_type>
KOKKOS_INLINE_FUNCTION
value_type
Stokhos::KL::ExponentialRandomField<value_type,device_type>::
evaluate_standard_deviation(const point_type& point) const
{
  value_type result = 0.0;
  for (int i=0; i<num_KL; i++) {
    value_type t = 1.0;
    for (int j=0; j<dim; ++j)
      t *= product_eigen_funcs(i,j).evaluate(point[j]);
    result += product_eigen_values(i).eig_val*t*t;
  }
  return std::sqrt(result);
}

template <typename value_type, typename device_type>
template <typename point_type>
KOKKOS_INLINE_FUNCTION
value_type
Stokhos::KL::ExponentialRandomField<value_type,device_type>::
evaluate_eigenfunction(const point_type& point, int i) const
{
  value_type t = std_dev*std::sqrt(product_eigen_values(i));
  for (int j=0; j<dim; ++j)
    t *= product_eigen_funcs(i,j).evaluate(point[j]);
  return t;
}

template <typename value_type, typename device_type>
void
Stokhos::KL::ExponentialRandomField<value_type,device_type>::
print(std::ostream& os) const
{
  os << "KL expansion using " << num_KL << " terms in " << dim
     << " dimensions:" << std::endl;
  os << "\tMean = " << mean << ", standard deviation = " << std_dev
     << std::endl;
  os << "\tEigenvalues, Eigenfunctions:" << std::endl;
  for (int i=0; i<num_KL; i++) {
    os << "\t\t" << product_eigen_values(i) << ", ";
    for (int j=0; j<dim-1; ++j) {
      os << "(";
      product_eigen_funcs(i,j).print(os);
      os << ") * ";
    }
    os << "(";
    product_eigen_funcs(i,dim-1).print(os);
    os << ")" << std::endl;
  }
}
