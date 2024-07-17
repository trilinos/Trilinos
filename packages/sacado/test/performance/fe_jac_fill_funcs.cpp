// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "fe_jac_fill_funcs.hpp"

void analytic_init_fill(const ElemData& e,
			unsigned int neqn,
			const std::vector<double>& x, 
			std::vector<double>& x_local) {
  for (unsigned int node=0; node<e.nnode; node++) 
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      x_local[node*neqn+eqn] = x[e.gid[node]*neqn+eqn];
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
  std::vector<double> s(e.nqp);
  for (unsigned int qp=0; qp<e.nqp; qp++) {
    s[qp] = 0.0;
    for (unsigned int eqn=0; eqn<neqn; eqn++)
      s[qp] += u[qp*neqn+eqn]*u[qp*neqn+eqn];
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
	double c2 = -e.dphi[qp][node_row]/(e.jac[qp]*e.jac[qp]);
	double c3 = exp(u[qp*neqn+eqn_row]);
	double c4 = e.phi[qp][node_row]*s[qp]*c3;
	f[row] += c1*(c2*du[qp*neqn+eqn_row] + c4);
	for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
	  jac[row][node_col*neqn+eqn_row] += 
	    c1*(c2*e.dphi[qp][node_col] + c4*e.phi[qp][node_col]);
	  for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++)
	      jac[row][node_col*neqn+eqn_col] += 
	  	2.0*c1*e.phi[qp][node_row]*u[qp*neqn+eqn_col]*e.phi[qp][node_col]*c3;
	}
      }
    }
  }
}

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
  // std::cout.setf(std::ios::scientific);
  // std::cout << "Analytic Jacobian = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
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

#ifdef HAVE_ADIC
void adic_init_fill(const ElemData& e,
		    unsigned int neqn,
		    const std::vector<double>& x, 
		    std::vector<DERIV_TYPE>& x_fad) {
  static bool first = true;
  for (unsigned int node=0; node<e.nnode; node++)
    for (unsigned int eqn=0; eqn<neqn; eqn++) {
      x_fad[node*neqn+eqn].value = x[e.gid[node]*neqn+eqn];
      if (first)
	ad_AD_SetIndep(x_fad[node*neqn+eqn]);
    }
  if (first) {
    ad_AD_SetIndepDone();
    first = false;
  }
}

/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 04/13/12 11:10:49 by the version of   */
/*   ADIC 1.2.3 compiled on  04/14/09 12:39:01                        */
/*                                                                    */
/*   ADIC was prepared as an account of work sponsored by an          */
/*   agency of the United States Government and the University of     */
/*   Chicago.  NEITHER THE AUTHOR(S), THE UNITED STATES GOVERNMENT    */
/*   NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, INCLUDING */
/*   ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS  */
/*   OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR */
/*   THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR  */
/*   PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE */
/*   PRIVATELY OWNED RIGHTS.                                          */
/*                                                                    */
/**********************************************************************/
void   adic_element_fill(ElemData  *e,unsigned int  neqn, DERIV_TYPE  *x,DERIV_TYPE  *u,DERIV_TYPE  *du,DERIV_TYPE  *f) {
unsigned int  var_0, var_1, var_2, var_3, var_4, var_5, var_6, var_7;
DERIV_TYPE  var_8;
double  adji_0;
    double  loc_0;
    double  loc_1;
    double  loc_2;
    double  loc_3;
    double  loc_4;
    double  loc_5;
    double  loc_6;
    double  loc_7;
    double  loc_8;
    double  loc_9;
    double  adj_0;
    double  adj_1;
    double  adj_2;
    double  adj_3;

        // static int g_filenum = 0;
        // if (g_filenum == 0) {
        //     adintr_ehsfid(&g_filenum, __FILE__, "adic_element_fill");
        // }
            for (unsigned int  qp = 0;     qp < e->nqp;     )    {
        for (unsigned int  eqn = 0;         eqn < neqn;         )        {
            {
                ad_grad_axpy_0(&(u[qp * neqn + eqn]));
                DERIV_val(u[qp * neqn + eqn]) = 0.0;
            }
            {
                ad_grad_axpy_0(&(du[qp * neqn + eqn]));
                DERIV_val(du[qp * neqn + eqn]) = 0.0;
            }
            for (unsigned int  node = 0;             node < e->nnode;             )            {
                {
                    loc_0 = DERIV_val(x[node * neqn + eqn]) * e->phi[qp][node];
                    loc_1 = DERIV_val(u[qp * neqn + eqn]) + loc_0;
                    ad_grad_axpy_2(&(u[qp * neqn + eqn]), 1.000000000000000e+00, &(u[qp * neqn + eqn]), e->phi[qp][node], &(x[node * neqn + eqn]));
                    DERIV_val(u[qp * neqn + eqn]) = loc_1;
                }
                {
                    loc_0 = DERIV_val(x[node * neqn + eqn]) * e->dphi[qp][node];
                    loc_1 = DERIV_val(du[qp * neqn + eqn]) + loc_0;
                    ad_grad_axpy_2(&(du[qp * neqn + eqn]), 1.000000000000000e+00, &(du[qp * neqn + eqn]), e->dphi[qp][node], &(x[node * neqn + eqn]));
                    DERIV_val(du[qp * neqn + eqn]) = loc_1;
                }
                var_2 = node++;
            }
            var_1 = eqn++;
        }
        var_0 = qp++;
    }
//DERIV_TYPE  *s = malloc(e->nqp *  sizeof (DERIV_TYPE ));
	    std::vector<DERIV_TYPE> s(e->nqp);
    for (unsigned int  qp = 0;     qp < e->nqp;     )    {
        {
            ad_grad_axpy_0(&(s[qp]));
            DERIV_val(s[qp]) = 0.0;
        }
        for (unsigned int  eqn = 0;         eqn < neqn;         )        {
            {
                loc_0 = DERIV_val(u[qp * neqn + eqn]) * DERIV_val(u[qp * neqn + eqn]);
                loc_1 = DERIV_val(s[qp]) + loc_0;
                ad_grad_axpy_3(&(s[qp]), 1.000000000000000e+00, &(s[qp]), DERIV_val(u[qp * neqn + eqn]), &(u[qp * neqn + eqn]), DERIV_val(u[qp * neqn + eqn]), &(u[qp * neqn + eqn]));
                DERIV_val(s[qp]) = loc_1;
            }
            var_4 = eqn++;
        }
        var_3 = qp++;
    }
    for (unsigned int  node = 0;     node < e->nnode;     )    {
        for (unsigned int  eqn = 0;         eqn < neqn;         )        {
unsigned int  row = node * neqn + eqn;
            {
                ad_grad_axpy_0(&(f[row]));
                DERIV_val(f[row]) = 0.0;
            }
            for (unsigned int  qp = 0;             qp < e->nqp;             )            {
     DERIV_val(var_8) = exp(( DERIV_val(u[qp * neqn + eqn])));
      adji_0 = DERIV_val(var_8);
                {
                    ad_grad_axpy_1(&(var_8), adji_0, &(u[qp * neqn + eqn]));
                }
                {
                    loc_0 = e->w[qp] * e->jac[qp];
                    loc_1 =  -e->dphi[qp][node];
                    loc_2 = e->jac[qp] * e->jac[qp];
                    loc_3 = loc_1 / loc_2;
                    loc_4 = loc_3 * DERIV_val(du[qp * neqn + eqn]);
                    loc_5 = e->phi[qp][node] * DERIV_val(s[qp]);
                    loc_6 = loc_5 * DERIV_val(var_8);
                    loc_7 = loc_4 + loc_6;
                    loc_8 = loc_0 * loc_7;
                    loc_9 = DERIV_val(f[row]) + loc_8;
                    adj_0 = loc_5 * loc_0;
                    adj_1 = DERIV_val(var_8) * loc_0;
                    adj_2 = e->phi[qp][node] * adj_1;
                    adj_3 = loc_3 * loc_0;
                    ad_grad_axpy_4(&(f[row]), 1.000000000000000e+00, &(f[row]), adj_3, &(du[qp * neqn + eqn]), adj_2, &(s[qp]), adj_0, &(var_8));
                    DERIV_val(f[row]) = loc_9;
                }
                var_7 = qp++;
            }
            var_6 = eqn++;
        }
        var_5 = node++;
    }
    //free(s);
}

void adic_process_fill(const ElemData& e,
		       unsigned int neqn,
		       const std::vector<DERIV_TYPE>& f_fad, 
		       std::vector<double>& f,
		       std::vector< std::vector<double> >& jac) {
  for (unsigned int eqn_row=0; eqn_row<neqn; eqn_row++) {
    f[e.gid[0]*neqn+eqn_row] += f_fad[0*neqn+eqn_row].value;
    f[e.gid[1]*neqn+eqn_row] += f_fad[1*neqn+eqn_row].value;
    for (unsigned int node_col=0; node_col<e.nnode; node_col++) {
      for (unsigned int eqn_col=0; eqn_col<neqn; eqn_col++) {
	unsigned int col = node_col*neqn+eqn_col;
	unsigned int next_col = (node_col+1)*neqn+eqn_col;
	jac[e.gid[0]*neqn+eqn_row][next_col] += 
	  f_fad[0*neqn+eqn_row].grad[col];
	jac[e.gid[1]*neqn+eqn_row][col] += 
	  f_fad[1*neqn+eqn_row].grad[col];
      }
    }
  }
}

double adic_jac_fill(unsigned int num_nodes, unsigned int num_eqns,
			 double mesh_spacing) {
  AD_Init(0);
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

  Teuchos::Time timer("FE ADIC Jacobian Fill", false);
  timer.start(true);
  std::vector<DERIV_TYPE> x_fad(e.nnode*num_eqns), f_fad(e.nnode*num_eqns);
  std::vector<DERIV_TYPE> u(e.nqp*num_eqns), du(e.nqp*num_eqns);
  for (unsigned int i=0; i<num_nodes-1; i++) {
    e.gid[0] = i;
    e.gid[1] = i+1;
    
    adic_init_fill(e, num_eqns, x, x_fad);
    adic_element_fill(&e, num_eqns, &x_fad[0], &u[0], &du[0], &f_fad[0]);
    adic_process_fill(e, num_eqns, f_fad, f, jac);
  }
  timer.stop();
  AD_Final();

  // std::cout.precision(8);
  // std::cout << "ADIC Residual = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++)
  //   std::cout << "\t" << f[i] << std::endl;

  // std::cout.precision(8);
  // std::cout << "ADIC Jacobian = " << std::endl;
  // for (unsigned int i=0; i<num_nodes*num_eqns; i++) {
  //   std::cout << "\t";
  //   for (unsigned int j=0; j<(e.nnode+1)*num_eqns; j++)
  //     std::cout << jac[i][j] << "\t";
  //   std::cout << std::endl;
  // }

  return timer.totalElapsedTime();
}
#endif

