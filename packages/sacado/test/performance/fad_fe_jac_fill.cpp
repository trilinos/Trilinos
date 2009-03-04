// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Sacado_Random.hpp"
#include "Sacado.hpp"
#include "Sacado_Fad_SimpleFad.hpp"

#include "Fad/fad.h"
#include "TinyFadET/tfad.h"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// ADOL-C includes
#ifdef HAVE_ADOLC
#ifdef PACKAGE
#undef PACKAGE
#endif
#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif
#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif
#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif
#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif
#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif
#ifdef VERSION
#undef VERSION
#endif
#define ADOLC_TAPELESS
#define NUMBER_DIRECTIONS 100
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/interfaces.h"
#endif

// A performance test that computes a finite-element-like Jacobian using
// several Fad variants

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;

void FAD::error(const char *msg) {
  std::cout << msg << std::endl;
}

struct ElemData {
  static const unsigned int nqp = 2;
  static const unsigned int nnode = 2;
  double w[nqp], jac[nqp], phi[nqp][nnode], dphi[nqp][nnode];
  unsigned int gid[nnode];

  ElemData(double mesh_spacing) {
    // Quadrature points
    double xi[nqp];
    xi[0] = -1.0 / std::sqrt(3.0);
    xi[1] =  1.0 / std::sqrt(3.0);
    
    for (unsigned int i=0; i<nqp; i++) {
      // Weights
      w[i] = 1.0;
      
      // Element transformation Jacobian
      jac[i] = 0.5*mesh_spacing;
      
      // Shape functions & derivatives
      phi[i][0] = 0.5*(1.0 - xi[i]);
      phi[i][1] = 0.5*(1.0 + xi[i]);
      dphi[i][0] = -0.5;
      dphi[i][1] = 0.5;
    }
  }
};

template <class FadType>
void fad_init_fill(const ElemData& e,
		   unsigned int neqn,
		   const std::vector<double>& x, 
		   std::vector<FadType>& x_fad) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_fad[node*neqn+eqn] = FadType(e.nnode*neqn, node*neqn+eqn, 
				     x[e.gid[node]*neqn+eqn]);
}

#ifdef HAVE_ADOLC
#ifndef ADOLC_TAPELESS
void adolc_init_fill(bool retape,
		     int tag,
		     const ElemData& e,
		     unsigned int neqn,
		     const std::vector<double>& x, 
		     std::vector<double>& x_local,
		     std::vector<adouble>& x_ad) {
  if (retape)
    trace_on(tag);
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      x_local[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];
      if (retape)
	x_ad[node*neqn+eqn] <<= x_local[node*neqn+eqn];
    }
}

#else

void adolc_tapeless_init_fill(const ElemData& e,
			      unsigned int neqn,
			      const std::vector<double>& x, 
			      std::vector<adtl::adouble>& x_ad) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      x_ad[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];
      for (unsigned int i=0; i<neqn*e.nnode; i++)
	x_ad[node*neqn+eqn].setADValue(i, 0.0);
      x_ad[node*neqn+eqn].setADValue(node*neqn+eqn, 1.0);
    }
  
}
#endif
#endif

void analytic_init_fill(const ElemData& e,
			unsigned int neqn,
			const std::vector<double>& x, 
			std::vector<double>& x_local) {
  for (unsigned int node=0; node<e.nnode; node++) 
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_local[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];
}

template <class T>
void template_element_fill(const ElemData& e, 
			   unsigned int neqn,
			   const std::vector<T>& x, 
			   std::vector<T>& u, 
			   std::vector<T>& du, 
			   std::vector<T>& f) {
  // Construct element solution, derivative
  for (unsigned int qp=0; qp<e.nqp; qp++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      u[qp*neqn+eqn] = 0.0;
      du[qp*neqn+eqn] = 0.0;
      for (unsigned int node=0; node<e.nnode; node++) {
	u[qp*neqn+eqn] += x[node*neqn+eqn] * e.phi[qp][node];
	du[qp*neqn+eqn] += x[node*neqn+eqn] * e.dphi[qp][node];
      }
    }
  }

  // Compute sum of equations for coupling
  std::vector<T> s(e.nqp*neqn);
  for (unsigned int qp=0; qp<e.nqp; qp++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      s[qp*neqn+eqn] = 0.0;
      for (unsigned int j=0; j<neqn; j++) {
      	if (j != eqn)
      	  s[qp*neqn+eqn] += u[qp*neqn+j]; 
      }
    }
  }

  // Evaluate element residual
  for (unsigned int node=0; node<e.nnode; node++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      unsigned int row = node*neqn+eqn;
      f[row] = 0.0;
      for (unsigned int qp=0; qp<e.nqp; qp++) {
	f[row] += 
	  e.w[qp]*e.jac[qp]*(-(1.0/(e.jac[qp]*e.jac[qp]))*
			     du[qp*neqn+eqn]*e.dphi[qp][node] + 
			     e.phi[qp][node]*(exp(u[qp*neqn+eqn]) + 
					      u[qp*neqn+eqn]*s[qp*neqn+eqn]));
      }
    }
  }
}

void analytic_element_fill(const ElemData& e, 
			   unsigned int neqn,
			   const std::vector<double>& x, 
			   std::vector<double>& u, 
			   std::vector<double>& du, 
			   std::vector<double>& f,
			   std::vector< std::vector<double> >& jac) {
  // Construct element solution, derivative
  for (unsigned int qp=0; qp<e.nqp; qp++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      u[qp*neqn+eqn] = 0.0;
      du[qp*neqn+eqn] = 0.0;
      for (unsigned int node=0; node<e.nnode; node++) {
	u[qp*neqn+eqn] += x[node*neqn+eqn] * e.phi[qp][node];
	du[qp*neqn+eqn] += x[node*neqn+eqn] * e.dphi[qp][node];
      }
    }
  }

  // Compute sum of equations for coupling
  std::vector<double> s(e.nqp*neqn);
  for (unsigned int qp=0; qp<e.nqp; qp++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      s[qp*neqn+eqn] = 0.0;
      for (unsigned int j=0; j<neqn; j++) {
	if (j != eqn)
	  s[qp*neqn+eqn] += u[qp*neqn+j]; 
      }
    }
  }

  // Evaluate element residual
  for (unsigned int node_row=0; node_row<e.nnode; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
      unsigned int row = node_row*neqn+eqn_row;
      f[row] = 0.0;
      for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
	for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	  unsigned int col = node_col*neqn+eqn_col;
	  jac[row][col] = 0.0;
	}
      }
      for (unsigned int qp=0; qp<e.nqp; qp++) {
	double c1 = e.w[qp]*e.jac[qp];
	double c2 = -(1.0/(e.jac[qp]*e.jac[qp]))*e.dphi[qp][node_row];
	double c3 = e.phi[qp][node_row]*(exp(u[qp*neqn+eqn_row]));
	double c4 = e.phi[qp][node_row]*u[qp*neqn+eqn_row];
	double c5 = e.phi[qp][node_row]*s[qp*neqn+eqn_row];
	double c35 = c3+c5;
	double c14 = c1*c4;
	f[row] += c1*(c2*du[qp*neqn+eqn_row] + c3 + c4*s[qp*neqn+eqn_row]);
	for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
	  jac[row][node_col*neqn+eqn_row] += 
	    c1*(c2*e.dphi[qp][node_col] + c35*e.phi[qp][node_col]);
	  for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	    if (eqn_col != eqn_row)
	      jac[row][node_col*neqn+eqn_col] += c14*e.phi[qp][node_col];
	  }
	}
      }
    }
  }
}

template <class FadType>
void fad_process_fill(const ElemData& e,
		      unsigned int neqn,
		      const std::vector<FadType>& f_fad, 
		      std::vector<double>& f,
		      std::vector< std::vector<double> >& jac) {
  for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
    f[e.gid[0]*neqn+eqn_row] += f_fad[0*neqn+eqn_row].val();
    f[e.gid[1]*neqn+eqn_row] += f_fad[1*neqn+eqn_row].val();
    for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
      for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	unsigned int col = node_col*neqn+eqn_col;
	unsigned int next_col = (node_col+1)*neqn+eqn_col;
	jac[e.gid[0]*neqn+eqn_row][next_col] += 
	  f_fad[0*neqn+eqn_row].fastAccessDx(col);
	jac[e.gid[1]*neqn+eqn_row][col] += 
	  f_fad[1*neqn+eqn_row].fastAccessDx(col);
      }
    }
  }
}

#ifdef HAVE_ADOLC
#ifndef ADOLC_TAPELESS
void adolc_process_fill(bool retape,
			int tag, 
			const ElemData& e,
			unsigned int neqn,
			std::vector<double>& x_local,
			std::vector<adouble>& f_ad, 
			std::vector<double>& f,
			std::vector<double>& f_local,
			std::vector< std::vector<double> >& jac,
			double **seed,
			double **jac_local) {
  if (retape) {
    for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++)
      f_ad[0*neqn+eqn_row] >>= f_local[0*neqn+eqn_row];
    for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++)
      f_ad[1*neqn+eqn_row] >>= f_local[1*neqn+eqn_row];
    trace_off();
  }

  //jacobian(tag, neqn*e.nnode, neqn*e.nnode, &x_local[0], jac_local);
  forward(tag, neqn*e.nnode, neqn*e.nnode, neqn*e.nnode, &x_local[0],
	  seed, &f_local[0], jac_local);
  
  for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
    f[e.gid[0]*neqn+eqn_row] += f_local[0*neqn+eqn_row];
    f[e.gid[1]*neqn+eqn_row] += f_local[1*neqn+eqn_row];
    for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
      for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	unsigned int col = node_col*neqn+eqn_col;
	unsigned int next_col = (node_col+1)*neqn+eqn_col;
	jac[e.gid[0]*neqn+eqn_row][next_col] += jac_local[0*neqn+eqn_row][col];
	jac[e.gid[1]*neqn+eqn_row][col] += jac_local[1*neqn+eqn_row][col];
      }
    }
  }
}

#else

void adolc_tapeless_process_fill(const ElemData& e,
				 unsigned int neqn,
				 const std::vector<adtl::adouble>& f_ad, 
				 std::vector<double>& f,
				 std::vector< std::vector<double> >& jac) {
  for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
    f[e.gid[0]*neqn+eqn_row] += f_ad[0*neqn+eqn_row].getValue();
    f[e.gid[1]*neqn+eqn_row] += f_ad[1*neqn+eqn_row].getValue();
    for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
      for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	unsigned int col = node_col*neqn+eqn_col;
	unsigned int next_col = (node_col+1)*neqn+eqn_col;
	jac[e.gid[0]*neqn+eqn_row][next_col] += 
	  f_ad[0*neqn+eqn_row].getADValue(col);
	jac[e.gid[1]*neqn+eqn_row][col] += 
	  f_ad[1*neqn+eqn_row].getADValue(col);
      }
    }
  }
}
#endif
#endif

void analytic_process_fill(const ElemData& e,
			   unsigned int neqn,
			   const std::vector<double>& f_local, 
			   const std::vector< std::vector<double> >& jac_local, 
			   std::vector<double>& f,
			   std::vector< std::vector<double> >& jac) {
  for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
    f[e.gid[0]*neqn+eqn_row] += f_local[0*neqn+eqn_row];
    f[e.gid[1]*neqn+eqn_row] += f_local[1*neqn+eqn_row];
    for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
      for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	unsigned int col = node_col*neqn+eqn_col;
	unsigned int next_col = (node_col+1)*neqn+eqn_col;
	jac[e.gid[0]*neqn+eqn_row][next_col] += jac_local[0*neqn+eqn_row][col];
	jac[e.gid[1]*neqn+eqn_row][col] += jac_local[1*neqn+eqn_row][col];
      }
    }
  }
}

void residual_process_fill(const ElemData& e,
			   unsigned int neqn,
			   const std::vector<double>& f_local, 
			   std::vector<double>& f) {
  for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
    f[e.gid[0]*neqn+eqn_row] += f_local[0*neqn+eqn_row];
    f[e.gid[1]*neqn+eqn_row] += f_local[1*neqn+eqn_row];
  }
}

template <class FadType>
double fad_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
		    double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > jac(num_nodes*num_eqns);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
      jac[row] = std::vector<double>((e.nnode+1)*num_eqns);
      for (unsigned int node_col=0; node_col<e.nnode+1; node_col++) {
	for (unsigned int eqn_col=0; eqn_col<num_eqns; eqn_col++) { 
	  unsigned int col = node_col*num_eqns + eqn_col;
	  jac[row][col] = 0.0;
	}
      }
    }
  }

  Teuchos::Time timer("FE Fad Jacobian Fill", false);
  timer.start(true);
  std::vector<FadType> x_fad(e.nnode*num_eqns), f_fad(e.nnode*num_eqns);
  std::vector<FadType> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    fad_init_fill(e, num_eqns, x, x_fad);
    template_element_fill(e, num_eqns, x_fad, u, du, f_fad);
    fad_process_fill(e, num_eqns, f_fad, f, jac);
  }
  timer.stop();

  // std::cout << "Fad Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Fad Jacobian = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

#ifdef HAVE_ADOLC
#ifndef ADOLC_TAPELESS
double adolc_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
		      double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > jac(num_nodes*num_eqns);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
      jac[row] = std::vector<double>((e.nnode+1)*num_eqns);
      for (unsigned int node_col=0; node_col<e.nnode+1; node_col++) {
	for (unsigned int eqn_col=0; eqn_col<num_eqns; eqn_col++) { 
	  unsigned int col = node_col*num_eqns + eqn_col;
	  jac[row][col] = 0.0;
	}
      }
    }
  }

  Teuchos::Time timer("FE ADOL-C Jacobian Fill", false);
  timer.start(true);
  std::vector<adouble> x_ad(e.nnode*num_eqns), f_ad(e.nnode*num_eqns);
  std::vector<adouble> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  std::vector<double> x_local(e.nnode*num_eqns);
  std::vector<double> f_local(e.nnode*num_eqns);
  double **jac_local = new double*[e.nnode*num_eqns];
  double **seed = new double*[e.nnode*num_eqns];
  for (unsigned int i=0; i<e.nnode*num_eqns; i++) {
    jac_local[i] = new double[e.nnode*num_eqns];
    seed[i] = new double[e.nnode*num_eqns];
    for (unsigned int j=0; j<e.nnode*num_eqns; j++)
      seed[i][j] = 0.0;
    seed[i][i] = 1.0;
  }
  
  // Tape first element
  e.gid[0] = 0;
  e.gid[1] = 1;
  adolc_init_fill(true, 0, e, num_eqns, x, x_local, x_ad);
  template_element_fill(e, num_eqns, x_ad, u, du, f_ad);
  adolc_process_fill(true, 0, e, num_eqns, x_local, f_ad, f, f_local,
		     jac, seed, jac_local);

  // Now do remaining fills reusing tape
  for (unsigned int i=1; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    adolc_init_fill(false, 0, e, num_eqns, x, x_local, x_ad);
    adolc_process_fill(false, 0, e, num_eqns, x_local, f_ad, f, f_local,
		       jac, seed, jac_local);
  }
  for (unsigned int i=0; i<e.nnode*num_eqns; i++) {
    delete [] jac_local[i];
    delete [] seed[i];
  }
  delete [] jac_local;
  delete [] seed;
  timer.stop();

  // std::cout << "ADOL-C Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "ADOL-C Jacobian = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

double adolc_retape_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			     double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > jac(num_nodes*num_eqns);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
      jac[row] = std::vector<double>((e.nnode+1)*num_eqns);
      for (unsigned int node_col=0; node_col<e.nnode+1; node_col++) {
	for (unsigned int eqn_col=0; eqn_col<num_eqns; eqn_col++) { 
	  unsigned int col = node_col*num_eqns + eqn_col;
	  jac[row][col] = 0.0;
	}
      }
    }
  }

  Teuchos::Time timer("FE ADOL-C Retape Jacobian Fill", false);
  timer.start(true);
  std::vector<adouble> x_ad(e.nnode*num_eqns), f_ad(e.nnode*num_eqns);
  std::vector<adouble> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  std::vector<double> x_local(e.nnode*num_eqns);
  std::vector<double> f_local(e.nnode*num_eqns);
  double **jac_local = new double*[e.nnode*num_eqns];
  double **seed = new double*[e.nnode*num_eqns];
  for (unsigned int i=0; i<e.nnode*num_eqns; i++) {
    jac_local[i] = new double[e.nnode*num_eqns];
    seed[i] = new double[e.nnode*num_eqns];
    for (unsigned int j=0; j<e.nnode*num_eqns; j++)
      seed[i][j] = 0.0;
    seed[i][i] = 1.0;
  }
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    adolc_init_fill(true, 1, e, num_eqns, x, x_local, x_ad);
    template_element_fill(e, num_eqns, x_ad, u, du, f_ad);
    adolc_process_fill(true, 1, e, num_eqns, x_local, f_ad, f, f_local,
		       jac, seed, jac_local);
  }
  for (unsigned int i=0; i<e.nnode*num_eqns; i++) {
    delete [] jac_local[i];
    delete [] seed[i];
  }
  delete [] jac_local;
  delete [] seed;
  timer.stop();

  // std::cout << "ADOL-C Residual (retaped) = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "ADOL-C Jacobian (retaped) = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

#else

double adolc_tapeless_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			       double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > jac(num_nodes*num_eqns);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
      jac[row] = std::vector<double>((e.nnode+1)*num_eqns);
      for (unsigned int node_col=0; node_col<e.nnode+1; node_col++) {
	for (unsigned int eqn_col=0; eqn_col<num_eqns; eqn_col++) { 
	  unsigned int col = node_col*num_eqns + eqn_col;
	  jac[row][col] = 0.0;
	}
      }
    }
  }

  Teuchos::Time timer("FE Tapeless ADOL-C Jacobian Fill", false);
  timer.start(true);
  std::vector<adtl::adouble> x_ad(e.nnode*num_eqns), f_ad(e.nnode*num_eqns);
  std::vector<adtl::adouble> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    adolc_tapeless_init_fill(e, num_eqns, x, x_ad);
    template_element_fill(e, num_eqns, x_ad, u, du, f_ad);
    adolc_tapeless_process_fill(e, num_eqns, f_ad, f, jac);
  }
  timer.stop();

  // std::cout << "Tapeless ADOL-C Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Tapeless ADOL-C Jacobian = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}
#endif
#endif

double analytic_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			 double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > jac(num_nodes*num_eqns);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
      jac[row] = std::vector<double>((e.nnode+1)*num_eqns);
      for (unsigned int node_col=0; node_col<e.nnode+1; node_col++) {
	for (unsigned int eqn_col=0; eqn_col<num_eqns; eqn_col++) { 
	  unsigned int col = node_col*num_eqns + eqn_col;
	  jac[row][col] = 0.0;
	}
      }
    }
  }

  Teuchos::Time timer("FE Analytic Jacobian Fill", false);
  timer.start(true);
  std::vector<double> x_local(e.nnode*num_eqns), f_local(e.nnode*num_eqns);
  std::vector< std::vector<double> > jac_local(e.nnode*num_eqns);
  std::vector<double> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<e.nnode*num_eqns; i++)
    jac_local[i] = std::vector<double>(e.nnode*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    analytic_init_fill(e, num_eqns, x, x_local);
    analytic_element_fill(e, num_eqns, x_local, u, du, f_local, jac_local);
    analytic_process_fill(e, num_eqns, f_local, jac_local, f, jac);
  }
  timer.stop();

  // std::cout.precision(8);
  // std::cout << "Analytic Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Analytic Jacobian = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

double residual_fill(unsigned int num_nodes, unsigned int num_eqns,
			 double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
    }
  }

  Teuchos::Time timer("FE Residual Fill", false);
  timer.start(true);
  std::vector<double> x_local(e.nnode*num_eqns), f_local(e.nnode*num_eqns);
  std::vector<double> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    analytic_init_fill(e, num_eqns, x, x_local);
    template_element_fill(e, num_eqns, x_local, u, du, f_local);
    residual_process_fill(e, num_eqns, f_local, f);
  }
  timer.stop();

  // std::cout.precision(8);
  // std::cout << "Analytic Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  return timer.totalElapsedTime();
}

int main(int argc, char* argv[]) {
  int ierr = 0;

  try {
    double t, ta, tr;
    int p = 2;
    int w = p+7;

    // Maximum number of derivative components for SLFad
    const int slfad_max = 100;

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a finite-element-like Jacobian fill");
    int num_nodes = 100000;
    int num_eqns = 2;
    int rt = 0;
    clp.setOption("n", &num_nodes, "Number of nodes");
    clp.setOption("p", &num_eqns, "Number of equations");
    clp.setOption("rt", &rt, "Include ADOL-C retaping test");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    double mesh_spacing = 1.0 / static_cast<double>(num_nodes - 1);

    // Memory pool & manager
    Sacado::Fad::MemPoolManager<double> poolManager(num_nodes*num_eqns);
    Sacado::Fad::MemPool* pool = poolManager.getMemoryPool(num_nodes*num_eqns);
    Sacado::Fad::DMFad<double>::setDefaultPool(pool);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "num_nodes =  " << num_nodes 
	      << ", num_eqns = " << num_eqns << ":  " << std::endl
	      << "           " << "   Time" << "\t"<< "Time/Analytic" << "\t"
	      << "Time/Residual" << std::endl;

    ta = 1.0;
    tr = 1.0;

    tr = residual_fill(num_nodes, num_eqns, mesh_spacing);

    ta = analytic_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "Analytic:  " << std::setw(w) << ta << "\t" << std::setw(w) << ta/ta << "\t" << std::setw(w) << ta/tr << std::endl;

#ifdef HAVE_ADOLC
#ifndef ADOLC_TAPELESS
    t = adolc_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ADOL-C:    " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    if (rt != 0) {
      t = adolc_retape_jac_fill(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ADOL-C(rt):" << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }

#else
    t = adolc_tapeless_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ADOL-C(tl):" << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
#endif
#endif

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< FAD::TFad<4,double> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "TFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }
    else if (num_eqns*2 == 40) {
      t = fad_jac_fill< FAD::TFad<40,double> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "TFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }

    t = fad_jac_fill< FAD::Fad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "Fad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< Sacado::Fad::SFad<double,4> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }
    else if (num_eqns*2 == 40) {
      t = fad_jac_fill< Sacado::Fad::SFad<double,40> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }

    if (num_eqns*2 < slfad_max) {
      t = fad_jac_fill< Sacado::Fad::SLFad<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SLFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }
    
    t = fad_jac_fill< Sacado::Fad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "DFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    t = fad_jac_fill< Sacado::Fad::SimpleFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "SimpleFad: " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    t = fad_jac_fill< Sacado::Fad::DMFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "DMFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl; 

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< Sacado::ELRFad::SFad<double,4> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }
    else if (num_eqns*2 == 40) {
      t = fad_jac_fill< Sacado::ELRFad::SFad<double,40> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }

    if (num_eqns*2 < slfad_max) {
      t = fad_jac_fill< Sacado::ELRFad::SLFad<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSLFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }

    t = fad_jac_fill< Sacado::ELRFad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ELRDFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    
    t = fad_jac_fill< Sacado::CacheFad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "CacheFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    t = fad_jac_fill< Sacado::Fad::DVFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "DVFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    
  }
  catch (std::exception& e) {
    cout << e.what() << endl;
    ierr = 1;
  }
  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

  return ierr;
}
