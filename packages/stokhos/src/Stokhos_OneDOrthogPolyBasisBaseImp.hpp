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

template <typename ordinal_type, typename value_type>
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
OneDOrthogPolyBasisBase(const std::string& name, ordinal_type p_) :
  name(name),
  p(p_),
  basis(p+1, Polynomial<value_type>(p)),
  double_basis(2*p+1, Polynomial<value_type>(2*p)),
  norms(p+1, value_type(0.0)),
  Cijk(),
  Bij()
{
}

template <typename ordinal_type, typename value_type>
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
~OneDOrthogPolyBasisBase()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
dimension() const
{
  return 1;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
size() const
{
  return p+1;
}

template <typename ordinal_type, typename value_type>
const std::vector<value_type>&
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
getTripleProductTensor() const
{
  ordinal_type sz = size();

  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  if (Cijk == Teuchos::null) {
    Cijk = Teuchos::rcp(new Dense3Tensor<ordinal_type, value_type>(sz));
    std::vector<value_type> a(2*sz);
    for (ordinal_type i=0; i<sz; i++) {
      for (ordinal_type j=0; j<sz; j++) {
	projectProduct(i, j, a);
	for (ordinal_type k=0; k<sz; k++) {
	  (*Cijk)(i,j,k) = a[k];
	}
      }
    }
  }

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
getDerivDoubleProductTensor() const
{
  ordinal_type sz = size();

  // Compute Bij = < \Psi_i \Psi_j' >
  if (Bij == Teuchos::null) {
    std::vector<value_type> b(sz);
    Bij = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type>(sz,sz));
    for (ordinal_type i=0; i<sz; i++) {
      projectDerivative(i, b);
      for (ordinal_type j=0; j<sz; j++)
	(*Bij)(i,j) = b[j]*norms[j];
    }
  }

  return Bij;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
projectProduct(ordinal_type i, ordinal_type j, std::vector<value_type>& coeffs) const
{
  // Important note:  To get the correct value for <\Psi_i \Psi_j \Psi_k>,
  // we have to expand \Psi_i*\Psi_j in the full order_i+order_j basis, not
  // just the order p basis.  There for the basis needs to be of size 2*p
  Polynomial<value_type> pij(2*p);
  std::vector<value_type> a(2*(p+1));

  // Multiply ith and jth basis polynomial
  pij.multiply(1.0, basis[i], basis[j], 0.0);

  // Project onto basis
  projectPoly(pij, a);
  for (ordinal_type k=0; k<=p; k++)
    coeffs[k] = a[k]*norms[k];
}

template <typename ordinal_type, typename value_type>
Stokhos::Polynomial<value_type>
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
toStandardBasis(const value_type coeffs[], ordinal_type n) const
{
  ordinal_type px = p;
  if (n < p+1)
    px = n-1;
  Polynomial<value_type> x(px);

  for (ordinal_type i=0; i<=px; i++)
    x.add(coeffs[i], basis[i], 1.0);

  return x;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  return basis[i].coeff(0);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << name << " basis of order " << p << ".  Basis vectors:\n";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(basis.size()); i++)
    os << "\t" << basis[i];
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(norms.size()); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
std::vector<ordinal_type>
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  std::vector<ordinal_type> t(1);
  t[0] = i;
  return t;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
getIndex(const std::vector<ordinal_type>& term) const
{
  return term[0];
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::OneDOrthogPolyBasisBase<ordinal_type, value_type>::
getName() const
{
  return name;
}
