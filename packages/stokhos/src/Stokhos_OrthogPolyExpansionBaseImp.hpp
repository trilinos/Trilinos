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

#include "Stokhos_DenseDirectDivisionExpansionStrategy.hpp"
#include "Stokhos_SPDDenseDirectDivisionExpansionStrategy.hpp"
#include "Stokhos_MeanBasedDivisionExpansionStrategy.hpp"
#include "Stokhos_GMRESDivisionExpansionStrategy.hpp"
#include "Stokhos_CGDivisionExpansionStrategy.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_DynamicArrayTraits.hpp"

template <typename ordinal_type, typename value_type, typename node_type> 
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
OrthogPolyExpansionBase(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  basis(basis_),
  Cijk(Cijk_),
  params(params_),
  sz(basis->size())
{
  if (params == Teuchos::null)
    params = Teuchos::rcp(new Teuchos::ParameterList);

  // Create division strategy
  std::string name = params->get("Division Strategy","Dense Direct");
  double tol = params->get("Division Tolerance", 1e-6);
  int prec_iter = params->get("prec_iter", 1);
  int max_it = params->get("max_it_div", 50);
  std::string prec = params->get("Prec Strategy","None");
  int PrecNum;
  if (prec == "None")
     PrecNum=0;
  else if (prec == "Diag")
     PrecNum=1;
  else if (prec == "Jacobi")
     PrecNum=2;
  else if (prec == "GS")
     PrecNum=3;
  else if (prec == "Schur")
     PrecNum=4;
  else 
     PrecNum=-1;
  std::string linopt = params->get("Prec option", "whole");
  int linear = 0;
  if (linopt == "linear")
	linear = 1; 

  std::string schuropt = params->get("Schur option", "diag");
  int diag = 0;
  if (schuropt == "full")
        diag = 1;

  int equil = params->get("Equilibrate", 0);

  
  if (name == "Dense Direct")
    division_strategy = 
      Teuchos::rcp(new DenseDirectDivisionExpansionStrategy<ordinal_type,value_type,node_type>(this->basis, this->Cijk));
  else if (name == "SPD Dense Direct")
    division_strategy = 
      Teuchos::rcp(new SPDDenseDirectDivisionExpansionStrategy<ordinal_type,value_type,node_type>(this->basis, this->Cijk));
  else if (name == "Mean-Based")
    division_strategy =
      Teuchos::rcp(new MeanBasedDivisionExpansionStrategy<ordinal_type,value_type,node_type>()); 
  else if (name == "GMRES")
    division_strategy =
      Teuchos::rcp(new GMRESDivisionExpansionStrategy<ordinal_type,value_type,node_type>(this->basis, this->Cijk, prec_iter, tol, PrecNum, max_it, linear, diag, equil));
  else if (name == "CG")
    division_strategy =
       Teuchos::rcp(new CGDivisionExpansionStrategy<ordinal_type,value_type,node_type>(this->basis, this->Cijk, prec_iter, tol, PrecNum, max_it, linear, diag, equil));

//    TEUCHOS_TEST_FOR_EXCEPTION(
//      true, std::logic_error,
//      "Invalid division strategy name" << name);
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
unaryMinus(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::unaryMinus(OPA)");
#endif
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* ca = a.coeff();

  for (ordinal_type i=0; i<pc; i++)
  {
    cc[i] = -ca[i];
  }
}
template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
plusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	  const value_type& val)
{
  c[0] += val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
minusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  c[0] -= val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::timesEqual(const)");
#endif
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] *= val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	    const value_type& val)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::divideEqual(const)");
#endif
  ordinal_type pc = c.size();
  value_type* cc=c.coeff();
  for (ordinal_type i=0; i<pc; i++){
    cc[i] /= val;
}
}
template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
plusEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::plusEqual(OPA)");
#endif
  ordinal_type xp = x.size();
  if (c.size() < xp)
    c.resize(xp);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (ordinal_type i=0; i<xp; i++)
    cc[i] += xc[i];
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
minusEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::minusEqual(OPA)");
#endif
  ordinal_type xp = x.size();
  if (c.size() < xp)
    c.resize(xp);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (ordinal_type i=0; i<xp; i++)
    cc[i] -= xc[i];
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
timesEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::timesEqual(OPA)");
#endif
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
  ordinal_type pc;
  if (p > 1 && xp > 1)
    pc = sz;
  else
    pc = p*xp;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::OrthogPolyExpansionBase::timesEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (p > 1 && xp > 1) {
    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);

    typename Cijk_type::i_iterator i_begin = Cijk->i_begin();
    typename Cijk_type::i_iterator i_end = Cijk->i_end();
    if (pc < Cijk->num_i())
      i_end = Cijk->find_i(pc);
    ordinal_type k_lim = p;
    ordinal_type j_lim = xp;
    const value_type* kc = tc;
    const value_type* jc = xc;
    if (xp < p) {
      k_lim = xp;
      j_lim = p;
      kc = xc;
      jc = tc;
    }

    value_type tmp, cijk;
    ordinal_type i,j,k;
    for (typename Cijk_type::i_iterator i_it=i_begin; i_it!=i_end; ++i_it) {
      i = index(i_it);
      tmp = value_type(0.0);
      for (typename Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it); 
	   k_it != Cijk->k_end(i_it); ++k_it) {
	k = index(k_it);
	if (k < k_lim) {
	  for (typename Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    j = index(j_it);
	    cijk = value(j_it);
	    if (j < j_lim)
	      tmp += cijk*kc[k]*jc[j];
	  }
	}
      }
      cc[i] = tmp / basis->norm_squared(i);
    }

    Stokhos::ds_array<value_type>::destroy_and_release(tc, p);
  }
  else if (p > 1) {
    for (ordinal_type i=0; i<p; i++)
      cc[i] *= xc[0];
  }
  else if (xp > 1) {
    for (ordinal_type i=1; i<xp; i++)
      cc[i] = cc[0]*xc[i];
    cc[0] *= xc[0];
  }
  else {
    cc[0] *= xc[0];
  }
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
divideEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type >& x)
{
  division_strategy->divide(c, 1.0, c, x, 0.0);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::plus(OPA,OPA)");
#endif
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > pb) {
    for (ordinal_type i=0; i<pb; i++)
      cc[i] = ca[i] + cb[i];
    for (ordinal_type i=pb; i<pc; i++)
      cc[i] = ca[i];
  }
  else {
    for (ordinal_type i=0; i<pa; i++)
      cc[i] = ca[i] + cb[i];
    for (ordinal_type i=pa; i<pc; i++)
      cc[i] = cb[i];
  }
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const value_type& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::plus(const,OPA)");
#endif
  ordinal_type pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a + cb[0];
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = cb[i];
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
     const value_type& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::plus(OPA,const)");
#endif
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::minus(OPA,OPA)");
#endif
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > pb) {
    for (ordinal_type i=0; i<pb; i++)
      cc[i] = ca[i] - cb[i];
    for (ordinal_type i=pb; i<pc; i++)
      cc[i] = ca[i];
  }
  else {
    for (ordinal_type i=0; i<pa; i++)
      cc[i] = ca[i] - cb[i];
    for (ordinal_type i=pa; i<pc; i++)
      cc[i] = -cb[i];
  }
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::minus(const,OPA)");
#endif
  ordinal_type pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a - cb[0];
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = -cb[i];
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::minus(OPA,const)");
#endif
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::times(OPA,OPA)");
#endif
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pa > 1 && pb > 1)
    pc = sz;
  else
    pc = pa*pb;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::OrthogPolyExpansionBase::times()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    typename Cijk_type::i_iterator i_begin = Cijk->i_begin();
    typename Cijk_type::i_iterator i_end = Cijk->i_end();
    if (pc < Cijk->num_i())
      i_end = Cijk->find_i(pc);
    ordinal_type k_lim = pa;
    ordinal_type j_lim = pb;
    const value_type* kc = ca;
    const value_type* jc = cb;
    if (pb < pa) {
      k_lim = pb;
      j_lim = pa;
      kc = cb;
      jc = ca;
    }
    value_type tmp, cijk;
    ordinal_type i,j,k;
    for (typename Cijk_type::i_iterator i_it=i_begin; i_it!=i_end; ++i_it) {
      i = index(i_it);
      tmp = value_type(0.0);
      for (typename Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it); 
	   k_it != Cijk->k_end(i_it); ++k_it) {
	k = index(k_it);
	if (k < k_lim) {
	  for (typename Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    j = index(j_it);
	    cijk = value(j_it);
	    if (j < j_lim)
	      tmp += cijk*kc[k]*jc[j];
	  }
	}
      }
      cc[i] = tmp / basis->norm_squared(i);
    }
  }
  else if (pa > 1) {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (pb > 1) {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[0]*cb[i];
  }
  else {
    cc[0] = ca[0]*cb[0];
  }
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::times(const,OPA)");
#endif
  ordinal_type pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = a*cb[i];
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::times(OPA,const)");
#endif
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = ca[i]*b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  division_strategy->divide(c, 1.0, a, b, 0.0);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  Stokhos::OrthogPolyApprox<ordinal_type,value_type> aa(b.basis(), 1, &a);
  division_strategy->divide(c, 1.0, aa, b, 0.0);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const value_type& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::divide(OPA,const)");
#endif
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = ca[i]/b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
fabs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::fabs(OPA)");
#endif
  c.init(0.0);
  c[0] = a.two_norm();
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
abs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::abs(OPA)");
#endif
  c.init(0.0);
  c[0] = a.two_norm();
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::max(OPA,OPA)");
#endif
  if (a.two_norm() >= b.two_norm())
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::max(const,OPA)");
#endif
  if (a >= b.two_norm()) {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = a;
  }
  else
    c = b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::max(OPA,const)");
#endif
  if (a.two_norm() >= b)
    c = a;
  else {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = b;
  }
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::min(OPA,OPA)");
#endif
  if (a.two_norm() <= b.two_norm())
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::min(const,OPA)");
#endif
  if (a <= b.two_norm()) {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = a;
  }
  else
    c = b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::min(OPA,const)");
#endif
  if (a.two_norm() <= b)
    c = a;
  else {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = b;
  }
}
