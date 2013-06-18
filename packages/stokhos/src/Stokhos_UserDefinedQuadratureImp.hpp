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

template <typename ordinal_type, typename value_type>
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
UserDefinedQuadrature(
const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
     const Teuchos::RCP<const Teuchos::Array<value_type> >& weights)
  : quad_points(points),
    quad_weights(weights)
{
  ordinal_type nqp = points->size();
  Teuchos::RCP<Teuchos::Array< Teuchos::Array<value_type> > > qv = 
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(nqp));
  for (ordinal_type i=0; i<nqp; i++) {
    (*qv)[i].resize(basis->size());
    basis->evaluateBases((*points)[i], (*qv)[i]);
  }
  quad_values = qv;
}

template <typename ordinal_type, typename value_type>
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
UserDefinedQuadrature(
const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
const Teuchos::RCP<const Teuchos::Array<value_type> >& weights,
const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& values)
  : quad_points(points),
    quad_weights(weights),
    quad_values(values)
{
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return *quad_points;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return *quad_weights;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return *quad_values;
}

template <typename ordinal_type, typename value_type>
std::ostream& 
Stokhos::UserDefinedQuadrature<ordinal_type,value_type>::
print(std::ostream& os) const
{
  ordinal_type nqp = quad_weights->size();
  os << "Sparse Grid Quadrature with " << nqp << " points:"
     << std::endl << "Weight : Points" << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    os << i << ": " << (*quad_weights)[i] << " : ";
    for (ordinal_type j=0; 
	 j<static_cast<ordinal_type>((*quad_points)[i].size()); 
	 j++)
      os << (*quad_points)[i][j] << " ";
    os << std::endl;
  }
  os << "Basis values at quadrature points:" << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    os << i << " " << ": ";
    for (ordinal_type j=0; 
	 j<static_cast<ordinal_type>((*quad_values)[i].size()); 
	 j++)
      os << (*quad_values)[i][j] << " ";
    os << std::endl;
  }

  return os;
}
