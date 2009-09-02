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
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox() :
  basis_(),
  coeff_(1, value_type(0.0))
{
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis) :
  basis_(basis),
  coeff_()
{
  if (basis != Teuchos::null)
    coeff_.resize(basis->size(), value_type(0.0));
  else
    coeff_.resize(1, value_type(0.0));
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
		 ordinal_type sz) :
  basis_(basis),
  coeff_(sz, value_type(0.0))
{
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& x) :
  basis_(x.basis_),
  coeff_(x.coeff_) 
{
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
~OrthogPolyApprox()
{
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
operator=(const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& x) 
{
  basis_ = x.basis_;
  if (this != &x) {
    coeff_ = x.coeff_;
  }

  return *this;
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
init(const value_type& v) 
{
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
    coeff_[i] = v;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
basis() const
{
  return basis_;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
reset(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& new_basis)
{
  basis_ = new_basis;
  coeff_.resize(basis_->size());
}

template <typename ordinal_type, typename value_type> 
ordinal_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
size() const 
{ 
  return coeff_.size(); 
}

template <typename ordinal_type, typename value_type> 
value_type*
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
coeff() 
{ 
  return &coeff_[0]; 
}

template <typename ordinal_type, typename value_type> 
const value_type*
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
coeff() const 
{ 
  return &coeff_[0]; 
}

template <typename ordinal_type, typename value_type> 
value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
operator[](ordinal_type i) 
{ 
  return coeff_[i]; 
}

template <typename ordinal_type, typename value_type> 
const value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
operator[](ordinal_type i) const
{ 
  return coeff_[i]; 
}

template <typename ordinal_type, typename value_type> 
value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
term(ordinal_type dimension, ordinal_type order)
{
  ordinal_type d = basis_->dimension();
  Teuchos::Array<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = basis_->getIndex(term);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type> 
const value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
term(ordinal_type dimension, ordinal_type order) const
{
  ordinal_type d = basis_->dimension();
  Teuchos::Array<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = basis_->getIndex(term);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type> 
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
evaluate(const Teuchos::Array<value_type>& point) const
{
  Teuchos::Array<value_type> basis_vals(basis_->size()); 
  basis_->evaluateBases(point, basis_vals);
  return evaluate(point, basis_vals);
}

template <typename ordinal_type, typename value_type> 
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
evaluate(const Teuchos::Array<value_type>& point,
	 const Teuchos::Array<value_type>& basis_vals) const
{
  value_type val = value_type(0.0);
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
    val += coeff_[i]*basis_vals[i];

  return val;
}

template <typename ordinal_type, typename value_type> 
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
mean() const
{
  return coeff_[0];
}

template <typename ordinal_type, typename value_type> 
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
standard_deviation() const
{
  value_type std_dev = 0.0;
  for (ordinal_type i=1; i<static_cast<ordinal_type>(coeff_.size()); i++) {
    std_dev += coeff_[i]*coeff_[i]*basis_->norm_squared(i);
  }
  std_dev = std::sqrt(std_dev);
  return std_dev;
}

template <typename ordinal_type, typename value_type> 
std::ostream&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
print(std::ostream& os) const
{
  Teuchos::Array<ordinal_type> trm;
  os << "Stokhos::OrthogPolyApprox of size " << coeff_.size() << " in basis "
     << "\n\t" << basis_->getName() << ":" << std::endl;
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++) {
    trm = basis_->getTerm(i);
    os << "\t\t(";
    for (ordinal_type j=0; j<static_cast<ordinal_type>(trm.size())-1; j++)
      os << trm[j] << ", ";
    os << trm[trm.size()-1] << ") = " << coeff_[i] << std::endl;
  }

  return os;
}

template <typename ordinal_type, typename value_type>
std::ostream& 
Stokhos::operator << (std::ostream& os, 
		      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  os << "[ ";
      
  for (ordinal_type i=0; i<static_cast<ordinal_type>(a.size()); i++) {
    os << a[i] << " ";
  }

  os << "]\n";
  return os;
}
