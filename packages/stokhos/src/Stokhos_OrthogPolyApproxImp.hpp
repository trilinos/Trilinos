// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ProductBasis.hpp"

template <typename ordinal_type, typename value_type, typename storage_type>
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
OrthogPolyApprox(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& theBasis,
  ordinal_type sz,
  const value_type* vals) :
  basis_(theBasis),
  coeff_(0)
{
  if (sz != 0)
    coeff_.resize(sz);
  else if (theBasis != Teuchos::null)
    coeff_.resize(theBasis->size());
  else
    coeff_.resize(ordinal_type(1));
  if (vals != NULL)
    coeff_.init(vals);
}

template <typename ordinal_type, typename value_type, typename storage_type>
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
OrthogPolyApprox(const Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>& x) :
  basis_(x.basis_),
  coeff_(x.coeff_)
{
}

template <typename ordinal_type, typename value_type, typename storage_type>
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
~OrthogPolyApprox()
{
}

template <typename ordinal_type, typename value_type, typename storage_type>
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>&
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
operator=(const Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>& x)
{
  if (this != &x) {
    basis_ = x.basis_;
    coeff_ = x.coeff_;
  }

  return *this;
}

template <typename ordinal_type, typename value_type, typename storage_type>
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>&
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
operator=(const value_type& v)
{
  coeff_.init(value_type(0));
  coeff_.init(&v, 1);

  return *this;
}

template <typename ordinal_type, typename value_type, typename storage_type>
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
init(const value_type& v)
{
  coeff_.init(v);
}

template <typename ordinal_type, typename value_type, typename storage_type>
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
init(const value_type* v)
{
  coeff_.init(v);
}

template <typename ordinal_type, typename value_type, typename storage_type>
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
load(value_type* v)
{
  coeff_.load(v);
}

template <typename ordinal_type, typename value_type, typename storage_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
basis() const
{
  return basis_;
}

template <typename ordinal_type, typename value_type, typename storage_type>
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
reset(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& new_basis, ordinal_type sz)
{
  basis_ = new_basis;
  if (sz != 0)
    resize(sz);
  else if (basis_ != Teuchos::null)
    resize(basis_->size());
  else
    resize(1);
}

template <typename ordinal_type, typename value_type, typename storage_type>
void
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
resize(ordinal_type sz)
{
  coeff_.resize(sz);
}

template <typename ordinal_type, typename value_type, typename storage_type>
ordinal_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
size() const
{
  return coeff_.size();
}

template <typename ordinal_type, typename value_type, typename storage_type>
typename Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::pointer
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
coeff()
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(coeff_.size() == 0, std::logic_error,
                     "Stokhos::OrthogPolyApprox::coeff():  " <<
                     "Coefficient array is empty!");
#endif
  return coeff_.coeff();
}

template <typename ordinal_type, typename value_type, typename storage_type>
typename Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::const_pointer
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
coeff() const
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(coeff_.size() == 0, std::logic_error,
                     "Stokhos::OrthogPolyApprox::coeff():  " <<
                     "Coefficient array is empty!");
#endif
  return coeff_.coeff();
}

template <typename ordinal_type, typename value_type, typename storage_type>
typename Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::reference
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
operator[](ordinal_type i)
{
  return coeff_[i];
}

template <typename ordinal_type, typename value_type, typename storage_type>
typename Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::const_reference
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
operator[](ordinal_type i) const
{
  return coeff_[i];
}

template <typename ordinal_type, typename value_type, typename storage_type>
typename Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::reference
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
term(ordinal_type dimension, ordinal_type theOrder)
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> >
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  ordinal_type d = basis_->dimension();
  MultiIndex<ordinal_type> theTerm(d);
  theTerm[dimension] = theOrder;
  ordinal_type index = product_basis->index(theTerm);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type, typename storage_type>
typename Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::const_reference
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
term(ordinal_type dimension, ordinal_type theOrder) const
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> >
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  ordinal_type d = basis_->dimension();
  MultiIndex<ordinal_type> theTerm(d);
  theTerm[dimension] = theOrder;
  ordinal_type index = product_basis->index(theTerm);
  return coeff_[index];
}

template <typename ordinal_type, typename value_type, typename storage_type>
const Stokhos::MultiIndex<ordinal_type>&
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
order(ordinal_type theTerm) const
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> >
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  return product_basis->term(theTerm);
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
evaluate(const Teuchos::Array<value_type>& point) const
{
  Teuchos::Array<value_type> basis_vals(basis_->size());
  basis_->evaluateBases(point, basis_vals);
  return evaluate(point, basis_vals);
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
evaluate(const Teuchos::Array<value_type>& point,
         const Teuchos::Array<value_type>& basis_vals) const
{
  value_type val = value_type(0.0);
  for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
    val += coeff_[i]*basis_vals[i];

  return val;
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
mean() const
{
  return coeff_[0];
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
standard_deviation() const
{
  value_type std_dev = 0.0;
  for (ordinal_type i=1; i<static_cast<ordinal_type>(coeff_.size()); i++) {
    std_dev += coeff_[i]*coeff_[i]*basis_->norm_squared(i);
  }
  std_dev = std::sqrt(std_dev);
  return std_dev;
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
two_norm() const
{
  return std::sqrt(this->two_norm_squared());
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
two_norm_squared() const
{
  value_type nrm = 0.0;
  if (basis_ == Teuchos::null) { // Check for special case of constants
    TEUCHOS_TEST_FOR_EXCEPTION(
      coeff_.size() != 1, std::logic_error,
      "basis_ == null && coeff_.size() > 1");
    nrm = coeff_[0]*coeff_[0];
  }
  else {
    for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
      nrm += coeff_[i]*coeff_[i]*basis_->norm_squared(i);
  }
  return nrm;
}

template <typename ordinal_type, typename value_type, typename storage_type>
value_type
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
inner_product(const Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>& b) const
{
  // Check a and b are compatible
  TEUCHOS_TEST_FOR_EXCEPTION(
    basis_ == Teuchos::null && coeff_.size() != 1, std::logic_error,
    "basis_ == null && coeff_.size() > 1");
  TEUCHOS_TEST_FOR_EXCEPTION(
    b.basis_ == Teuchos::null && b.coeff_.size() != 1, std::logic_error,
      "b.basis_ == null && b.coeff_.size() > 1");
  TEUCHOS_TEST_FOR_EXCEPTION(
    coeff_.size() != b.coeff_.size() &&
    coeff_.size() != 1 && b.coeff_.size() != 1, std::logic_error,
    "Coefficient array sizes do not match");

  value_type v = 0.0;
  if (coeff_.size() == 1 || b.coeff_.size() == 1)
    v = coeff_[0]*b.coeff_[0];
  else
    for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++)
      v += coeff_[i]*b.coeff_[i]*basis_->norm_squared(i);

  return v;
}

template <typename ordinal_type, typename value_type, typename storage_type>
std::ostream&
Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>::
print(std::ostream& os) const
{
  os << "Stokhos::OrthogPolyApprox of size " << coeff_.size() << " in basis "
     << "\n\t" << basis_->getName() << ":" << std::endl;

  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> >
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_);

  if (product_basis != Teuchos::null) {
    for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++) {
      const Stokhos::MultiIndex<ordinal_type>& trm = product_basis->term(i);
      os << "\t\t(";
      for (ordinal_type j=0; j<static_cast<ordinal_type>(trm.size())-1; j++)
        os << trm[j] << ", ";
      os << trm[trm.size()-1] << ") = " << coeff_[i] << std::endl;
    }
  }
  else {
    os << "[ ";

    for (ordinal_type i=0; i<static_cast<ordinal_type>(coeff_.size()); i++) {
      os << coeff_[i] << " ";
    }

    os << "]\n";
  }

  return os;
}

template <typename ordinal_type, typename value_type, typename storage_type>
std::ostream&
Stokhos::operator << (
  std::ostream& os,
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, storage_type>& a)
{
  os << "[ ";

  for (ordinal_type i=0; i<static_cast<ordinal_type>(a.size()); i++) {
    os << a[i] << " ";
  }

  os << "]\n";
  return os;
}
