// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


namespace ROL {

/* ======================================================================
                             Kronecker Products
   ====================================================================== */ 
template <class Real>
Quadrature<Real> kron_prod(Quadrature<Real> & rule1, Quadrature<Real> & rule2) {
  /* 
    Compute the Kronecker Product of a Tensor Product rule and a 1D rule.
  */
  if (rule1.getNumPoints()==0) {
    Quadrature<Real> rule(rule2);
    return rule;
  }
  else {
    // Initialize Arrays Containing Updated Nodes and Weights
    int dim1 = rule1.getDimension();
    int dim2 = rule2.getDimension();
    Quadrature<Real> rule(dim1+dim2); 

    Real weight = 0.0;
    std::vector<Real> node2(dim2,0.0);

    // Perform Kronecker Products
    // Compute Kronecker Product of Nodes
    typename std::map<std::vector<Real>,int>::iterator it = rule.begin();
    typename std::map<std::vector<Real>,int>::iterator it_i;
    typename std::map<std::vector<Real>,int>::iterator it_j;
    for (it_i=rule1.begin(); it_i!=rule1.end(); it_i++) {
      for (it_j=rule2.begin(); it_j!=rule2.end(); it_j++) {
	std::vector<Real> node = rule1.getNode(it_i);
	node2   = rule2.getNode(it_j);
	weight  = rule1.getWeight(node)*rule2.getWeight(node2);
	//node.push_back(node2[0]);
        node.insert(node.end(),node2.begin(),node2.end());
	rule.insert(it,node,weight);
	it = rule.end();
      }
    }
    return rule;
  }
}

int growthRule1D(int index, EROLGrowth growth, EROLBurkardt rule) {
  //
  //  Compute the growth sequence for 1D quadrature rules according to growth.
  //  For more information on growth rules, see 
  //  
  //  J. Burkardt. 1D Quadrature Rules For Sparse Grids.
  //  http://people.sc.fsu.edu/~jburkardt/presentations/sgmga_1d_rules.pdf.
  //  
  //  J. Burkardt. SGMGA: Sparse Grid Mixed Growth Anisotropic Rules.
  //  http://people.sc.fsu.edu/~jburkardt/cpp_src/sgmga/sgmga.html.
  //  
  //  Drew P. Kouri
  //  Sandia National Laboratories - CSRI
  //  May 27, 2011
  //

  int level = index-1;
  //int level = index;
  if (rule==BURK_CLENSHAWCURTIS) { // Clenshaw-Curtis
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      if (level==0) {
	return 1;
      }
      else { 
	int o = 2;
	while(o<2*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      if (level==0) {
	return 1;
      }
      else {
	int o = 2;
	while (o<4*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_FULLEXP) {
      if (level==0) {
	return 1;
      }
      else {
	return (int)pow(2.0,(double)level)+1;
      }
    }
  }
  else if (rule==BURK_FEJER2) { // Fejer Type 2
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (o<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      int o = 1;
      while (o<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_PATTERSON) { // Gauss-Patterson
    if (growth==GROWTH_SLOWLIN||
	growth==GROWTH_SLOWLINODD||
	growth==GROWTH_MODLIN) {
      std::cout << "Specified Growth Rule Not Allowed!\n";
      return 0;
    }
    else if (growth==GROWTH_SLOWEXP) {
      if (level==0) {
	return 1;
      }
      else {
	int p = 5;
	int o = 3;
	while (p<2*level+1) {
	  p = 2*p+1;
	  o = 2*o+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      if (level==0) {
	return 1;
      }
      else {
	int p = 5;
	int o = 3;
	while (p<4*level+1) {
	  p = 2*p+1;
	  o = 2*o+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_LEGENDRE) { // Gauss-Legendre
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_HERMITE) { // Gauss-Hermite
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }
  
  else if (rule==BURK_LAGUERRE) { // Gauss-Laguerre
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }


  else if (rule==BURK_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }
  
  else if (rule==BURK_GENZKEISTER) { // Hermite-Genz-Keister  
    static int o_hgk[5] = { 1, 3, 9, 19, 35 };
    static int p_hgk[5] = { 1, 5, 15, 29, 51 };
    if (growth==GROWTH_SLOWLIN||
	growth==GROWTH_SLOWLINODD||
	growth==GROWTH_MODLIN) {
      std::cout << "Specified Growth Rule Not Allowed!\n";
      return 0;
    }
    else if (growth==GROWTH_SLOWEXP) { 
      int l = 0, p = p_hgk[l], o = o_hgk[l];
      while (p<2*level+1 && l<4) {
	l++;
	p = p_hgk[l];
	o = o_hgk[l];
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      int l = 0, p = p_hgk[l], o = o_hgk[l];
      while (p<4*level+1 && l<4) {
	l++;
	p = p_hgk[l];
	o = o_hgk[l];
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      int l = level; l = std::max(l,0); l = std::min(l,4);
      return o_hgk[l];
    }
  }  

  else if (rule==BURK_TRAPEZOIDAL) { // Trapezoidal
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      if (level==0) {
	return 1;
      }
      else { 
	int o = 2;
	while(o<2*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      if (level==0) {
	return 1;
      }
      else {
	int o = 2;
	while (o<4*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_FULLEXP) {
      if (level==0) {
	return 1;
      }
      else {
	return (int)pow(2.0,(double)level)+1;
      }
    }
  }
  return 0;
} // end growthRule1D
} // end ROL namespace
