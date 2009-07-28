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
  coeff_(1)
{
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(const value_type& x) :
  coeff_(1)
{
  coeff_[0] = x;
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(ordinal_type sz, const value_type& x) :
  coeff_(sz)
{
  coeff_[0] = x;
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(ordinal_type sz) :
  coeff_(sz)
{
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(ordinal_type sz, ordinal_type l) :
  coeff_(l)
{
  coeff_.resize(sz);
}

template <typename ordinal_type, typename value_type> 
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
OrthogPolyApprox(const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& x) :
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
  if (this != &x) {
    coeff_ = x.coeff_;
  }

  return *this;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
resize(ordinal_type sz)
{
  coeff_.resize(sz);
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
reserve(ordinal_type sz)
{
  coeff_.reserve(sz);
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
term(const BasisT& basis, 
     int i0, int i1, int i2, int i3, int i4,
     int i5, int i6, int i7, int i8, int i9)
{
  std::vector<ordinal_type> trm;
  ordinal_type d = basis.dimension();
  if (i0 >= 0 && d >= 1)
    trm.push_back(i0);
  if (i1 >= 0 && d >= 2)
    trm.push_back(i1);
  if (i2 >= 0 && d >= 3)
    trm.push_back(i2);
  if (i3 >= 0 && d >= 4)
    trm.push_back(i3);
  if (i4 >= 0 && d >= 5)
    trm.push_back(i4);
  if (i5 >= 0 && d >= 6)
    trm.push_back(i5);
  if (i6 >= 0 && d >= 7)
    trm.push_back(i6);
  if (i7 >= 0 && d >= 8)
    trm.push_back(i7);
  if (i8 >= 0 && d >= 9)
    trm.push_back(i8);
  if (i9 >= 0 && d >= 10)
    trm.push_back(i9);

  ordinal_type index = basis.getIndex(trm);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
const value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
term(const BasisT& basis, 
     int i0, int i1, int i2, int i3, int i4,
     int i5, int i6, int i7, int i8, int i9) const
{
  std::vector<ordinal_type> trm;
  ordinal_type d = basis.dimension();
  if (i0 >= 0 && d >= 1)
    trm.push_back(i0);
  if (i1 >= 0 && d >= 2)
    trm.push_back(i1);
  if (i2 >= 0 && d >= 3)
    trm.push_back(i2);
  if (i3 >= 0 && d >= 4)
    trm.push_back(i3);
  if (i4 >= 0 && d >= 5)
    trm.push_back(i4);
  if (i5 >= 0 && d >= 6)
    trm.push_back(i5);
  if (i6 >= 0 && d >= 7)
    trm.push_back(i6);
  if (i7 >= 0 && d >= 8)
    trm.push_back(i7);
  if (i8 >= 0 && d >= 9)
    trm.push_back(i8);
  if (i9 >= 0 && d >= 10)
    trm.push_back(i9);

  ordinal_type index = basis.getIndex(trm);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
term2(const BasisT& basis, ordinal_type dimension, ordinal_type order)
{
  ordinal_type d = basis.dimension();
  std::vector<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = basis.getIndex(term);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
const value_type&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
term2(const BasisT& basis, ordinal_type dimension, ordinal_type order) const
{
  ordinal_type d = basis.dimension();
  std::vector<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = basis.getIndex(term);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
Stokhos::Polynomial<value_type>
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
toStandardBasis(const BasisT& basis) const
{
  return basis.toStandardBasis(&coeff_[0], coeff_.size());
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
evaluate(const BasisT& basis, const std::vector<value_type>& point) const
{
  const std::vector<value_type>& basis_pts = basis.evaluateBases(point);
  value_type val = value_type(0.0);
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
    val += coeff_[i]*basis_pts[i];

  return val;
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
evaluate(const BasisT& basis, const std::vector<value_type>& point,
         const std::vector<value_type>& basis_vals) const
{
  value_type val = value_type(0.0);
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
    val += coeff_[i]*basis_vals[i];

  return val;
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
std::ostream&
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
print(const BasisT& basis, std::ostream& os) const
{
  std::vector<ordinal_type> trm;
  os << "Stokhos::OrthogPolyApprox of size " << coeff_.size() << " in basis "
     << "\n\t" << basis.getName() << ":" << std::endl;
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++) {
    trm = basis.getTerm(i);
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

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
mean(const BasisT& basis) const
{
  return coeff_[0];
}

template <typename ordinal_type, typename value_type> 
template <typename BasisT>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type>::
standard_deviation(const BasisT& basis) const
{
  value_type std_dev = 0.0;
  for (ordinal_type i=1; i<static_cast<ordinal_type>(coeff_.size()); i++) {
    std_dev += coeff_[i]*coeff_[i]*basis.norm_squared(i);
  }
  std_dev = std::sqrt(std_dev);
  return std_dev;
}
