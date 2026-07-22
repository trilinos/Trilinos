// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_No_Kokkos.hpp"
#include "Sacado_tradvec.hpp"

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
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/interfaces.h"
#include "adolc/taping.h"
#endif

// A performance test that computes a finite-element-like adjoint using
// Rad and ADOL-C

typedef Sacado::Rad::ADvar<double> RadType;
typedef Sacado::RadVec::ADvar<double> VRadType;

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
		   const std::vector< std::vector<double> >& w,
		   std::vector<FadType>& x_fad,
		   std::vector< std::vector<double> >& w_local) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_fad[node*neqn+eqn] = FadType(e.nnode*neqn, node*neqn+eqn, 
				     x[e.gid[node]*neqn+eqn]);

  for (unsigned int col=0; col<w.size(); col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	w_local[col][node*neqn+eqn] = w[col][e.gid[node]*neqn+eqn];
}

void rad_init_fill(const ElemData& e,
		   unsigned int neqn,
		   const std::vector<double>& x,
		   const std::vector< std::vector<double> >& w,
		   std::vector<RadType>& x_rad,
		   std::vector< std::vector<double> >& w_local) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_rad[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];

  for (unsigned int col=0; col<w.size(); col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	w_local[col][node*neqn+eqn] = w[col][e.gid[node]*neqn+eqn];
}

void vrad_init_fill(const ElemData& e,
		   unsigned int neqn,
		   const std::vector<double>& x,
		   const std::vector< std::vector<double> >& w,
		   std::vector<VRadType>& x_rad,
		   double** w_local) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_rad[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];

  for (unsigned int col=0; col<w.size(); col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	w_local[col][node*neqn+eqn] = w[col][e.gid[node]*neqn+eqn];
}

#ifdef HAVE_ADOLC
void adolc_init_fill(bool retape,
		     int tag,
		     const ElemData& e,
		     unsigned int neqn,
		     const std::vector<double>& x,
		     const std::vector< std::vector<double> >& w,
		     std::vector<double>& x_local,
		     std::vector<adouble>& x_ad,
		     double** w_local) {
  if (retape)
    trace_on(tag, 1);
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      x_local[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];
      if (retape)
	x_ad[node*neqn+eqn] <<= x_local[node*neqn+eqn];
    }

  for (unsigned int col=0; col<w.size(); col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	w_local[col][node*neqn+eqn] = w[col][e.gid[node]*neqn+eqn];
}
#endif

void analytic_init_fill(const ElemData& e,
			unsigned int neqn,
			const std::vector<double>& x,
			const std::vector< std::vector<double> >& w,
			std::vector<double>& x_local,
			std::vector< std::vector<double> >& w_local) {
  for (unsigned int node=0; node<e.nnode; node++) 
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_local[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];

  for (unsigned int node=0; node<e.nnode; node++) 
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      for (unsigned int col=0; col<w.size(); col++)
	w_local[col][node*neqn+eqn] = w[col][e.gid[node]*neqn+eqn];
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
			   std::vector< std::vector<double> >& w,
			   std::vector<double>& u, 
			   std::vector<double>& du, 
			   std::vector<double>& f,
			   std::vector< std::vector<double> >& adj) {
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

  for (unsigned int col=0; col<w.size(); col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	adj[col][node*neqn+eqn] = 0.0;

  // Evaluate element residual
  for (unsigned int node=0; node<e.nnode; node++) {
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      unsigned int row = node*neqn+eqn;
      f[row] = 0.0;
      for (unsigned int qp=0; qp<e.nqp; qp++) {
	double c1 = e.w[qp]*e.jac[qp];
	double c2 = -(1.0/(e.jac[qp]*e.jac[qp]))*e.dphi[qp][node];
	double c3 = e.phi[qp][node]*(exp(u[qp*neqn+eqn]));
	double c4 = e.phi[qp][node]*u[qp*neqn+eqn];
	double c5 = e.phi[qp][node]*s[qp*neqn+eqn];
	double c35 = c3+c5;
	double c14 = c1*c4;
	f[row] += c1*(c2*du[qp*neqn+eqn] + c3 + c4*s[qp*neqn+eqn]);
	for (unsigned int col=0; col<w.size(); col++) {
	  for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
	    adj[col][node_col*neqn+eqn] += 
	      c1*(c2*e.dphi[qp][node_col] + c35*e.phi[qp][node_col])*w[col][row];
	    for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	      if (eqn_col != eqn)
		adj[col][node_col*neqn+eqn_col] += c14*e.phi[qp][node_col]*w[col][row];
	    }
	  }
	}
      }
    }
  }
}

template <class FadType>
void fad_process_fill(const ElemData& e,
		      unsigned int neqn,
		      const std::vector< std::vector<double> >& w_local,
		      const std::vector<FadType>& f_fad, 
		      std::vector<double>& f,
		      std::vector< std::vector<double> >& adj) {
  // Get residual
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      f[e.gid[node]*neqn+eqn] += f_fad[node*neqn+eqn].val();

  // Get adjoint for each adjoint direction
  for (unsigned int col=0; col<w_local.size(); col++) {
    FadType z = 0.0;
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++) 
	z += f_fad[node*neqn+eqn]*w_local[col][node*neqn+eqn];

    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++) 
	adj[col][e.gid[node]*neqn+eqn] += z.fastAccessDx(node*neqn+eqn);
  }
}

void rad_process_fill(const ElemData& e,
		      unsigned int neqn,
		      std::vector< std::vector<double> >& w_local,
		      std::vector<RadType>& x_rad,
		      std::vector<RadType>& f_rad, 
		      std::vector<RadType*>& ff_rad,
		      std::vector<double>& f,
		      std::vector< std::vector<double> >& adj) {
  // Get residual
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      f[e.gid[node]*neqn+eqn] += f_rad[node*neqn+eqn].val();

  // Get adjoint for each adjoint direction
  for (unsigned int col=0; col<w_local.size(); col++) {
    Sacado::Rad::ADcontext<double>::Weighted_Gradcomp(neqn*e.nnode, 
						      &ff_rad[0],
						      &w_local[col][0]);
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++) 
	adj[col][e.gid[node]*neqn+eqn] += x_rad[node*neqn+eqn].adj();
  }
}

void vrad_process_fill(const ElemData& e,
		       unsigned int neqn,
		       unsigned int ncol,
		       std::vector<std::size_t>& vn,
		       double** w_local,
		       std::vector<VRadType>& x_rad,
		       std::vector<VRadType>& f_rad, 
		       VRadType*** vf_rad,
		       std::vector<double>& f,
		       std::vector< std::vector<double> >& adj) {
  // Get residual
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      f[e.gid[node]*neqn+eqn] += f_rad[node*neqn+eqn].val();

  // Get adjoint for each adjoint direction
  Sacado::RadVec::ADcontext<double>::Weighted_GradcompVec(ncol,
							  &vn[0], 
							  vf_rad,
							  w_local);
  for (unsigned int col=0; col<ncol; col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++) 
	adj[col][e.gid[node]*neqn+eqn] += x_rad[node*neqn+eqn].adj(col);
}

#ifdef HAVE_ADOLC
void adolc_process_fill(bool retape,
			int tag, 
			const ElemData& e,
			unsigned int neqn,
			unsigned int ncol,
			std::vector<double>& x_local,
			double **w_local,
			std::vector<adouble>& f_ad,
			std::vector<double>& f_local,
			std::vector<double>& f,
			double **adj_local,
			std::vector< std::vector<double> >& adj) {
  if (retape) {
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	f_ad[node*neqn+eqn] >>= f_local[node*neqn+eqn];
    trace_off();
  }
  else
    zos_forward(tag, neqn*e.nnode, neqn*e.nnode, 1, &x_local[0], &f_local[0]);

  fov_reverse(tag, neqn*e.nnode, neqn*e.nnode, ncol, w_local, adj_local);
  
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      f[e.gid[node]*neqn+eqn] += f_local[node*neqn+eqn];

  for (unsigned int col=0; col<ncol; col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	adj[col][e.gid[node]*neqn+eqn] += adj_local[col][node*neqn+eqn];
}
#endif

void analytic_process_fill(const ElemData& e,
			   unsigned int neqn,
			   const std::vector<double>& f_local, 
			   const std::vector< std::vector<double> >& adj_local, 
			   std::vector<double>& f,
			   std::vector< std::vector<double> >& adj) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      f[e.gid[node]*neqn+eqn] += f_local[node*neqn+eqn];

  for (unsigned int col=0; col<adj.size(); col++)
    for (unsigned int node=0; node<e.nnode; node++)
      for (unsigned int eqn=0; eqn<neqn; eqn++)
	adj[col][e.gid[node]*neqn+eqn] += adj_local[col][node*neqn+eqn];
}

void residual_process_fill(const ElemData& e,
			   unsigned int neqn,
			   const std::vector<double>& f_local, 
			   std::vector<double>& f) {
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      f[e.gid[node]*neqn+eqn] += f_local[node*neqn+eqn];
}

template <typename FadType>
double fad_adj_fill(unsigned int num_nodes, unsigned int num_eqns,
		    unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, adjoint
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols), w_local(num_cols);
  std::vector< std::vector<double> > adj(num_cols);
  for (unsigned int node=0; node<num_nodes; node++)
    for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
      x[node*num_eqns+eqn] = 
	(mesh_spacing*node - 0.5)*(mesh_spacing*node - 0.5);
      f[node*num_eqns+eqn] = 0.0;
    }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    w_local[col] = std::vector<double>(e.nnode*num_eqns);
    adj[col] = std::vector<double>(num_nodes*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
	adj[col][node*num_eqns+eqn] = 0.0;
      }
  }

  Teuchos::Time timer("FE Fad Adjoint Fill", false);
  timer.start(true);
  std::vector<FadType> x_fad(e.nnode*num_eqns), f_fad(e.nnode*num_eqns);
  std::vector<FadType> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    fad_init_fill(e, num_eqns, x, w, x_fad, w_local);
    template_element_fill(e, num_eqns, x_fad, u, du, f_fad);
    fad_process_fill(e, num_eqns, w_local, f_fad, f, adj);
  }
  timer.stop();

  // std::cout << "Fad Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Fad Adjoint = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<num_cols; j++)
  //     std::cout << adj[j][i] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

double rad_adj_fill(unsigned int num_nodes, unsigned int num_eqns,
		    unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, adjoint
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols), w_local(num_cols);
  std::vector< std::vector<double> > adj(num_cols);
  for (unsigned int node=0; node<num_nodes; node++)
    for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
      x[node*num_eqns+eqn] = 
	(mesh_spacing*node - 0.5)*(mesh_spacing*node - 0.5);
      f[node*num_eqns+eqn] = 0.0;
    }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    w_local[col] = std::vector<double>(e.nnode*num_eqns);
    adj[col] = std::vector<double>(num_nodes*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
	adj[col][node*num_eqns+eqn] = 0.0;
      }
  }

  Teuchos::Time timer("FE Rad Adjoint Fill", false);
  timer.start(true);
  std::vector<RadType> x_rad(e.nnode*num_eqns), f_rad(e.nnode*num_eqns);
  std::vector<RadType> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  std::vector<RadType*> ff_rad(f_rad.size());
  for (unsigned int i=0; i<f_rad.size(); i++)
    ff_rad[i] = &f_rad[i];
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    rad_init_fill(e, num_eqns, x, w, x_rad, w_local);
    template_element_fill(e, num_eqns, x_rad, u, du, f_rad);
    rad_process_fill(e, num_eqns, w_local, x_rad, f_rad, ff_rad, f, adj);
  }
  timer.stop();

  // std::cout << "Rad Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Rad Adjoint = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<num_cols; j++)
  //     std::cout << adj[j][i] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

double vrad_adj_fill(unsigned int num_nodes, unsigned int num_eqns,
		     unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, adjoint
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols);
  double **w_local = new double*[num_cols];
  std::vector< std::vector<double> > adj(num_cols);
  for (unsigned int node=0; node<num_nodes; node++)
    for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
      x[node*num_eqns+eqn] = 
	(mesh_spacing*node - 0.5)*(mesh_spacing*node - 0.5);
      f[node*num_eqns+eqn] = 0.0;
    }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    w_local[col] = new double[e.nnode*num_eqns];
    adj[col] = std::vector<double>(num_nodes*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
	adj[col][node*num_eqns+eqn] = 0.0;
      }
  }

  Teuchos::Time timer("FE Vector Rad Adjoint Fill", false);
  timer.start(true);
  std::vector<VRadType> x_rad(e.nnode*num_eqns), f_rad(e.nnode*num_eqns);
  std::vector<VRadType> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  VRadType ***vf_rad = new VRadType**[num_cols];
  std::vector<std::size_t> vn(num_cols);
  for (unsigned int i=0; i<num_cols; i++) {
    vf_rad[i] = new VRadType*[num_eqns*e.nnode];
    vn[i] = num_eqns*e.nnode;
    for (unsigned int j=0; j<num_eqns*e.nnode; j++)
      vf_rad[i][j] = &f_rad[j];
  }
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    vrad_init_fill(e, num_eqns, x, w, x_rad, w_local);
    template_element_fill(e, num_eqns, x_rad, u, du, f_rad);
    vrad_process_fill(e, num_eqns, num_cols, vn, w_local, x_rad, f_rad, vf_rad,
		      f, adj);
  }
  for (unsigned int col=0; col<num_cols; col++) {
    delete [] w_local[col];
    delete [] vf_rad[col];
  }
  delete [] w_local;
  delete [] vf_rad;
  timer.stop();

  // std::cout << "Vector Rad Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Vector Rad Adjoint = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<num_cols; j++)
  //     std::cout << adj[j][i] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

#ifdef HAVE_ADOLC
double adolc_adj_fill(unsigned int num_nodes, unsigned int num_eqns,
		      unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, adjoint
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols);
  std::vector< std::vector<double> > adj(num_cols);
  for (unsigned int node=0; node<num_nodes; node++)
    for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
      x[node*num_eqns+eqn] = 
	(mesh_spacing*node - 0.5)*(mesh_spacing*node - 0.5);
      f[node*num_eqns+eqn] = 0.0;
    }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    adj[col] = std::vector<double>(num_nodes*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
	adj[col][node*num_eqns+eqn] = 0.0;
      }
  }

  Teuchos::Time timer("FE ADOL-C Adjoint Fill", false);
  timer.start(true);
  std::vector<adouble> x_ad(e.nnode*num_eqns), f_ad(e.nnode*num_eqns);
  std::vector<adouble> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  std::vector<double> x_local(e.nnode*num_eqns);
  std::vector<double> f_local(e.nnode*num_eqns);
  double **adj_local = new double*[num_cols];
  double **w_local = new double*[num_cols];
  for (unsigned int i=0; i<num_cols; i++) {
    adj_local[i] = new double[e.nnode*num_eqns];
    w_local[i] = new double[e.nnode*num_eqns];
  }
  
  // Tape first element
  e.gid[0] = 0;
  e.gid[1] = 1;
  adolc_init_fill(true, 0, e, num_eqns, x, w, x_local, x_ad, w_local);
  template_element_fill(e, num_eqns, x_ad, u, du, f_ad);
  adolc_process_fill(true, 0, e, num_eqns, num_cols, x_local, w_local, f_ad, 
		     f_local, f, adj_local, adj);

  // Now do remaining fills reusing tape
  for (unsigned int i=1; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    adolc_init_fill(false, 0, e, num_eqns, x, w, x_local, x_ad, w_local);
    adolc_process_fill(false, 0, e, num_eqns, num_cols, x_local, w_local, f_ad, 
		       f_local, f, adj_local, adj);
  }
  for (unsigned int i=0; i<num_cols; i++) {
    delete [] adj_local[i];
    delete [] w_local[i];
  }
  delete [] adj_local;
  delete [] w_local;
  timer.stop();

  // std::cout << "ADOL-C Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "ADOL-C Adjoint = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<num_cols; j++)
  //     std::cout << adj[j][i] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

double adolc_retape_adj_fill(unsigned int num_nodes, unsigned int num_eqns,
			     unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, adjoint
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols);
  std::vector< std::vector<double> > adj(num_cols);
  for (unsigned int node=0; node<num_nodes; node++)
    for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
      x[node*num_eqns+eqn] = 
	(mesh_spacing*node - 0.5)*(mesh_spacing*node - 0.5);
      f[node*num_eqns+eqn] = 0.0;
    }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    adj[col] = std::vector<double>(num_nodes*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
	adj[col][node*num_eqns+eqn] = 0.0;
      }
  }

  Teuchos::Time timer("FE ADOL-C Retape Adjoint Fill", false);
  timer.start(true);
  std::vector<adouble> x_ad(e.nnode*num_eqns), f_ad(e.nnode*num_eqns);
  std::vector<adouble> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  std::vector<double> x_local(e.nnode*num_eqns);
  std::vector<double> f_local(e.nnode*num_eqns);
  double **adj_local = new double*[num_cols];
  double **w_local = new double*[num_cols];
  for (unsigned int i=0; i<num_cols; i++) {
    adj_local[i] = new double[e.nnode*num_eqns];
    w_local[i] = new double[e.nnode*num_eqns];
  }

  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    adolc_init_fill(true, 1, e, num_eqns, x, w, x_local, x_ad, w_local);
    template_element_fill(e, num_eqns, x_ad, u, du, f_ad);
    adolc_process_fill(true, 1, e, num_eqns, num_cols, x_local, w_local, f_ad, 
		       f_local, f, adj_local, adj);
  }
  for (unsigned int i=0; i<num_cols; i++) {
    delete [] adj_local[i];
    delete [] w_local[i];
  }
  delete [] adj_local;
  delete [] w_local;
  timer.stop();

  // std::cout << "ADOL-C Residual (retaped) = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "ADOL-C Adjoint (retaped) = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<num_cols; j++)
  //     std::cout << adj[j][i] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}
#endif

double analytic_adj_fill(unsigned int num_nodes, unsigned int num_eqns,
			 unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, adjoint
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols), w_local(num_cols);
  std::vector< std::vector<double> > adj(num_cols);
  for (unsigned int node=0; node<num_nodes; node++)
    for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
      x[node*num_eqns+eqn] = 
	(mesh_spacing*node - 0.5)*(mesh_spacing*node - 0.5);
      f[node*num_eqns+eqn] = 0.0;
    }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    w_local[col] = std::vector<double>(e.nnode*num_eqns);
    adj[col] = std::vector<double>(num_nodes*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
	adj[col][node*num_eqns+eqn] = 0.0;
      }
  }

  Teuchos::Time timer("FE Analytic Adjoint Fill", false);
  timer.start(true);
  std::vector<double> x_local(e.nnode*num_eqns), f_local(e.nnode*num_eqns);
  std::vector< std::vector<double> > adj_local(num_cols);
  std::vector<double> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_cols; i++)
    adj_local[i] = std::vector<double>(e.nnode*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    analytic_init_fill(e, num_eqns, x, w, x_local, w_local);
    analytic_element_fill(e, num_eqns, x_local, w_local, u, du, f_local, 
     			  adj_local);
    analytic_process_fill(e, num_eqns, f_local, adj_local, f, adj);
  }
  timer.stop();

  // std::cout.precision(8);
  // std::cout << "Analytic Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "Analytic Adjoint = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<num_cols; j++)
  //     std::cout << adj[j][i] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}

double residual_fill(unsigned int num_nodes, unsigned int num_eqns,
		     unsigned int num_cols, double mesh_spacing) {
  ElemData e(mesh_spacing);

  // Solution vector, residual, jacobian
  std::vector<double> x(num_nodes*num_eqns), f(num_nodes*num_eqns);
  std::vector< std::vector<double> > w(num_cols), w_local(num_cols);
  for (unsigned int node_row=0; node_row<num_nodes; node_row++) {
    for (unsigned int eqn_row=0; eqn_row<num_eqns; eqn_row++) { 
      unsigned int row = node_row*num_eqns + eqn_row;
      x[row] = (mesh_spacing*node_row - 0.5)*(mesh_spacing*node_row - 0.5);
      f[row] = 0.0;
    }
  }

  for (unsigned int col=0; col<num_cols; col++) {
    w[col] = std::vector<double>(num_nodes*num_eqns);
    w_local[col] = std::vector<double>(e.nnode*num_eqns);
    for (unsigned int node=0; node<num_nodes; node++)
      for (unsigned int eqn=0; eqn<num_eqns; eqn++) {
	w[col][node*num_eqns+eqn] = x[node*num_eqns+eqn];
      }
  }

  Teuchos::Time timer("FE Residual Fill", false);
  timer.start(true);
  std::vector<double> x_local(e.nnode*num_eqns), f_local(e.nnode*num_eqns);
  std::vector<double> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    analytic_init_fill(e, num_eqns, x, w, x_local, w_local);
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

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a finite-element-like Jacobian fill");
    int num_nodes = 100000;
    int num_eqns = 2;
    int num_cols = 4;
    int rt = 0;
    clp.setOption("n", &num_nodes, "Number of nodes");
    clp.setOption("p", &num_eqns, "Number of equations");
    clp.setOption("q", &num_cols, "Number of adjoint directions");
    clp.setOption("rt", &rt, "Include ADOL-C retaping test");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    double mesh_spacing = 1.0 / static_cast<double>(num_nodes - 1);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "num_nodes = " << num_nodes 
	      << ", num_eqns = " << num_eqns 
	      << ", num_cols = " << num_cols << ":  " << std::endl
	      << "           " << "   Time" << "\t"<< "Time/Analytic" << "\t"
	      << "Time/Residual" << std::endl;

    ta = 1.0;
    tr = 1.0;

    tr = residual_fill(num_nodes, num_eqns, num_cols, mesh_spacing);

    ta = analytic_adj_fill(num_nodes, num_eqns, num_cols, mesh_spacing);
    std::cout << "Analytic:  " << std::setw(w) << ta << "\t" << std::setw(w) << ta/ta << "\t" << std::setw(w) << ta/tr << std::endl;

#ifdef HAVE_ADOLC
    t = adolc_adj_fill(num_nodes, num_eqns, num_cols, mesh_spacing);
    std::cout << "ADOL-C:    " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    if (rt != 0) {
      t = adolc_retape_adj_fill(num_nodes, num_eqns, num_cols, mesh_spacing);
      std::cout << "ADOL-C(rt):" << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;
    }
#endif
    
    t = fad_adj_fill< Sacado::Fad::DFad<double> >(num_nodes, num_eqns, num_cols, mesh_spacing);
    std::cout << "DFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    t = fad_adj_fill< Sacado::ELRFad::DFad<double> >(num_nodes, num_eqns, num_cols, mesh_spacing);
    std::cout << "ELRDFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    t = rad_adj_fill(num_nodes, num_eqns, num_cols, mesh_spacing);
    std::cout << "Rad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

    t = vrad_adj_fill(num_nodes, num_eqns, num_cols, mesh_spacing);
    std::cout << "Vector Rad:" << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/tr << std::endl;

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  return ierr;
}
