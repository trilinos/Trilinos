// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef FE_JAC_FILL_FUNCS_HPP
#define FE_JAC_FILL_FUNCS_HPP

#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Fad_SimpleFad.hpp"

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
//#define ADOLC_TAPELESS
#define NUMBER_DIRECTIONS 100
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/interfaces.h"
#include "adolc/taping.h"
#endif

#ifdef HAVE_ADIC
// We have included an ADIC differentiated version of the element fill 
// routine to compare the speed of operator overloading to source
// transformation.  To run this code, all that is necessary is to turn
// on the ADIC TPL.  However to modify the code, it is necessary to
// re-run the ADIC source transformation tool.  To do so, first update
// the changes to adic_element_fill.c.  Then set the following environment
// variables:
//     export ADIC_ARCH=linux
//     export ADIC=/home/etphipp/AD_libs/adic
// Next run ADIC via in the tests/performance source directory:
//     ${ADIC}/bin/linux/adiC -vd gradient -i adic_element_fill.init
// Finally, copy the resulting differentiated function in adic_element_fill.ad.c
// back into this file.  The function will need to be edited by changing
// the allocation of s to a std::vector<DERIV_TYPE> (the compiler 
// doesn't seem to like malloc), and commenting out the g_filenum lines.
#define ad_GRAD_MAX 130
#include "ad_deriv.h"
#endif


// A performance test that computes a finite-element-like Jacobian using
// several Fad variants

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
  std::vector<T> s(e.nqp);
  for (unsigned int qp=0; qp<e.nqp; qp++) {
    s[qp] = 0.0;
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      s[qp] += u[qp*neqn+eqn]*u[qp*neqn+eqn];
  }

  // Evaluate element residual
  for (unsigned int node=0; node<e.nnode; node++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      unsigned int row = node*neqn+eqn;
      f[row] = 0.0;
      for (unsigned int qp=0; qp<e.nqp; qp++) {
	double c1 = e.w[qp]*e.jac[qp];
	double c2 = -e.dphi[qp][node]/(e.jac[qp]*e.jac[qp]);
	f[row] += 
	  c1*(c2*du[qp*neqn+eqn] + e.phi[qp][node]*s[qp]*exp(u[qp*neqn+eqn]));
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

double analytic_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			 double mesh_spacing);
double residual_fill(unsigned int num_nodes, unsigned int num_eqns,
		     double mesh_spacing);

#ifdef HAVE_ADOLC
#ifndef ADOLC_TAPELESS
void adolc_init_fill(bool retape,
		     int tag,
		     const ElemData& e,
		     unsigned int neqn,
		     const std::vector<double>& x, 
		     std::vector<double>& x_local,
		     std::vector<adouble>& x_ad);

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
			double **jac_local);

double adolc_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
		      double mesh_spacing);

double adolc_retape_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			     double mesh_spacing);

#else

void adolc_tapeless_init_fill(const ElemData& e,
			      unsigned int neqn,
			      const std::vector<double>& x, 
			      std::vector<adtl::adouble>& x_ad);

void adolc_tapeless_process_fill(const ElemData& e,
				 unsigned int neqn,
				 const std::vector<adtl::adouble>& f_ad, 
				 std::vector<double>& f,
				 std::vector< std::vector<double> >& jac);

double adolc_tapeless_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			       double mesh_spacing);

#endif
#endif

#ifdef HAVE_ADIC
void adic_init_fill(const ElemData& e,
		    unsigned int neqn,
		    const std::vector<double>& x, 
		    std::vector<DERIV_TYPE>& x_fad);
void   adic_element_fill(ElemData  *e,unsigned int  neqn,DERIV_TYPE  *x,DERIV_TYPE  *u,DERIV_TYPE  *du,DERIV_TYPE  *f);
void adic_process_fill(const ElemData& e,
		       unsigned int neqn,
		       const std::vector<DERIV_TYPE>& f_fad, 
		       std::vector<double>& f,
		       std::vector< std::vector<double> >& jac);
double adic_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
		     double mesh_spacing);
inline void   AD_Init(int  arg0) {
  ad_AD_GradInit(arg0);  
}
inline void   AD_Final() {
  ad_AD_GradFinal();
}
#endif

#endif
