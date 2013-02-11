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

#include "Sacado_DynamicArrayTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_ConstantOrthogPolyExpansion.hpp"

namespace Sacado {
namespace ETPCE {

template <int NN, typename ExprT>
struct ExprQuadFuncWrapper {
  static const int N = NN;
  typedef typename ExprT::value_type value_type;
  const ExprT& expr_;
  ExprQuadFuncWrapper(const ExprT& expr) : expr_(expr) {}
  template <typename tuple_type>
  KERNEL_PREFIX value_type operator() (tuple_type x) const {
    return expr_.template eval_sample<0>(x);
  }
};

template <typename T, typename Storage> 
template <typename S>
void
OrthogPolyImpl<T,Storage>::
expressionCopy(const Expr<S>& x)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("ETPCE ExpressionCopy(" << x.name() << ")");
#endif
  int p = x.order();
  if (p == 0) {
    (*th_)[0] = x.val();
  }
  else if (p <= 2) {
    int sz = th_->size();
    approx_type* tc = th_.get();
    bool on_rhs = false;

    // Check if *this is on the right-hand-side of the expression
    if (p == 2) {
      const int N = Expr<S>::num_args;
      for (int i=0; i<N; i++) {
	const approx_type* opa = &(x.getArg(i));
	if (opa == &(*th_)) {
	  on_rhs = true;
	  break;
	}
      }
    }

    // If we are on the RHS, we have to put the results in a temporary
    if (on_rhs)
      tc = new approx_type(expansion_->getBasis(), sz);
    
    (*tc)[0] = x.val();
    if (x.has_fast_access(sz)) {
      for (int i=1; i<sz; i++)
  	(*tc)[i] = x.fast_higher_order_coeff(i);
    }
    else {
      for (int i=1; i<sz; i++)
  	(*tc)[i] = x.higher_order_coeff(i);
    }

    // Set underlying OPA if we created a temporary
    if (on_rhs) {
      th_ = Sacado::Handle<approx_type>(tc);
    }
  }
  else {
    const int N = Expr<S>::num_args;
    const approx_type* opas[N];
    for (int i=0; i<N; i++)
      opas[i] = &(x.getArg(i));
    ExprQuadFuncWrapper< N, Expr<S> > func(x);
    quad_expansion_->nary_op(func, *th_, opas);
  }
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>::
OrthogPolyImpl() :
  expansion_(const_expansion_),
  quad_expansion_(),
  th_(new Stokhos::OrthogPolyApprox<int,value_type,Storage>)
{
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>::
OrthogPolyImpl(const typename OrthogPolyImpl<T,Storage>::value_type& x) :
  expansion_(const_expansion_),
  quad_expansion_(),
  th_(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(Teuchos::null, 1, &x))
{
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>::
OrthogPolyImpl(const Teuchos::RCP<expansion_type>& expansion) :
  expansion_(expansion),
  quad_expansion_(Teuchos::rcp_dynamic_cast<quad_expansion_type>(expansion_)),
  th_(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(expansion_->getBasis()))
{
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>::
OrthogPolyImpl(const Teuchos::RCP<expansion_type>& expansion,
	   ordinal_type sz) :
  expansion_(expansion),
  quad_expansion_(Teuchos::rcp_dynamic_cast<quad_expansion_type>(expansion_)),
  th_(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(expansion_->getBasis(), sz))
{
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>::
OrthogPolyImpl(const OrthogPolyImpl<T,Storage>& x) :
  expansion_(x.expansion_),
  quad_expansion_(x.quad_expansion_),
  th_(x.th_)
{
}

template <typename T, typename Storage> 
template <typename S>
OrthogPolyImpl<T,Storage>::
OrthogPolyImpl(const Expr<S>& x) :
  expansion_(x.expansion()),
  quad_expansion_(Teuchos::rcp_dynamic_cast<quad_expansion_type>(expansion_)),
  th_(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(expansion_->getBasis(), x.size()))
{
  expressionCopy(x);
}

template <typename T, typename Storage> 
void
OrthogPolyImpl<T,Storage>::
reset(const Teuchos::RCP<expansion_type>& expansion)
{
  expansion_ = expansion;
  quad_expansion_ = Teuchos::rcp_dynamic_cast<quad_expansion_type>(expansion_);
  th_->reset(expansion_->getBasis());
}

template <typename T, typename Storage> 
void
OrthogPolyImpl<T,Storage>::
reset(const Teuchos::RCP<expansion_type>& expansion, ordinal_type sz)
{
  expansion_ = expansion;
  quad_expansion_ = Teuchos::rcp_dynamic_cast<quad_expansion_type>(expansion_);
  th_->reset(expansion_->getBasis(), sz);
}

template <typename T, typename Storage> 
typename OrthogPolyImpl<T,Storage>::value_type
OrthogPolyImpl<T,Storage>::
evaluate(const Teuchos::Array<typename OrthogPolyImpl<T,Storage>::value_type>& point) const
{
  return th_->evaluate(point);
}

template <typename T, typename Storage> 
typename OrthogPolyImpl<T,Storage>::value_type
OrthogPolyImpl<T,Storage>::
evaluate(
  const Teuchos::Array<typename OrthogPolyImpl<T,Storage>::value_type>& point,
  const Teuchos::Array<typename OrthogPolyImpl<T,Storage>::value_type>& bvals) const
{
  return th_->evaluate(point, bvals);
}

template <typename T, typename Storage> 
template <typename S>
bool 
OrthogPolyImpl<T,Storage>::
isEqualTo(const Expr<S>& x) const {
  typedef IsEqual<value_type> IE;
  if (x.size() != this->size()) return false;
  // Allow expansions to be different if their size is 1 and one
  // of them is a constant expansion
  if (expansion_ != x.expansion_) {
    if (x.size() != 1)
      return false;
    if ((expansion_ != const_expansion_) && 
	(x.expansion_ != const_expansion_))
      return false;
  }
  bool eq = true;
  for (int i=0; i<this->size(); i++)
    eq = eq && IE::eval(x.coeff(i), this->coeff(i));
  return eq;
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator=(const typename OrthogPolyImpl<T,Storage>::value_type& v) 
{
  th_.makeOwnCopy();
  *th_ = v;
  return *this;
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator=(const OrthogPolyImpl<T,Storage>& x) 
{
  expansion_ = x.expansion_;
  quad_expansion_ = x.quad_expansion_;
  th_ = x.th_;
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator=(const Expr<S>& x) 
{
  th_.makeOwnCopy();
  expansion_ = x.expansion();
  quad_expansion_ = Teuchos::rcp_dynamic_cast<quad_expansion_type>(expansion_);
  th_->reset(expansion_->getBasis(), x.size());
  expressionCopy(x);
  return *this;
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator+=(const typename OrthogPolyImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  expansion_->plusEqual(*th_, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator-=(const typename OrthogPolyImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  expansion_->minusEqual(*th_, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator*=(const typename OrthogPolyImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  expansion_->timesEqual(*th_, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPolyImpl<T,Storage>& 
OrthogPolyImpl<T,Storage>::
operator/=(const typename OrthogPolyImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  expansion_->divideEqual(*th_, v);
  return *this;
}

template <typename T, typename Storage>
std::ostream& 
operator << (std::ostream& os, const OrthogPoly<T,Storage>& a)
{
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;

  os << "[ ";
      
  for (ordinal_type i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

template <typename T, typename Storage>
std::istream& 
operator >> (std::istream& is, OrthogPoly<T,Storage>& a)
{
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;

  // Read in the opening "["
  char bracket;
  is >> bracket;
      
  for (ordinal_type i=0; i<a.size(); i++) {
    is >> a.fastAccessCoeff(i);
  }

  // Read in the closing "]"

  is >> bracket;
  return is;
}

} // namespace PCE
} // namespace Sacado
