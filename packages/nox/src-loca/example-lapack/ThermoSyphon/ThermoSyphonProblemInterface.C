// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "ThermoSyphonProblemInterface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Matrix.H"

// Declarations of extern thermosyphon fortran functions

extern "C" {

#define THERMINIT_F77 F77_FUNC(therminit,THERMINIT)
#define THERMFILL_F77 F77_FUNC(thermfill,THERMFILL)

void THERMINIT_F77(int *res, int *jac, int *mass, int *m, int *mm,
		   double *phi, double *psi, double *w, double *matR, 
		   double *matB, double *matI, double *matB2, double *matB2R, 
		   double *part1,
		   const int *lu_in, double *prandtl, double *rayleigh,
		   double *bval1, double *bval2, const char *fname1, 
		   const char *fname2);
void THERMFILL_F77(int *flag, const int *m_max, int *mm, double *u,
		   double *w, double *phi, double *psi, double *convw,
		   double *convphi, double *convpsi, double *matB2R, 
		   double *rayleigh, double *prandtl, double *block12,
		   double *block13, double *block23, double *block31,
		   double *part1, double *lapphi, double *lappsi, 
		   double *lapw, double *phiw, double *psiw, double *matB,
		   double *matR, double *matI, double *rhs_phi, 
		   double *rhs_psi, double *rhs_w, double *bval1, 
		   double *bval2, int *m, const int *buff,
		   double *J, double *g, double *M);
		   
}

ThermoSyphonProblemInterface::ThermoSyphonProblemInterface(const string& file1,
							   const string& file2,
							   int M_max)  : 
  // Initial guess (read from file)
  initialGuess(),

  // File names
  fname1(file1),
  fname2(file2),

  // Dimensions (m, mm read from file)
  m(0),
  mm(0),
  m_max(M_max),

  // Continuation/Bifurcation parameters (read from file)
  prandtl(0.0),
  rayleigh(0.0),

  // Work arrays
  w(2*m_max+1),
  phi(2*m_max+1),
  psi(2*m_max+1),
  rhs_w(m_max+1),
  rhs_phi(m_max+1),
  rhs_psi(m_max+1),
  lapphi(m_max+1),
  lappsi(m_max+1),
  lapw(m_max+1),
  phiw(m_max+1),
  psiw(m_max+1),
  convw(m_max+1,m_max+1),
  convpsi(m_max+1,m_max+1),
  convphi(m_max+1,m_max+1),
  matB2R(m_max+1,m_max+1),
  matI(m_max+1,m_max+1),
  matB(m_max+1,m_max+1),
  matB2(m_max+1,m_max+1),
  matR(m_max+1,m_max+1),
  block12(m_max+1,m_max+1),
  block13(m_max+1,m_max+1),
  block23(m_max+1,m_max+1),
  block31(m_max+1,m_max+1),
  part1(m_max+1,m_max+1),

  // Boundary values (read from file)
  bval1(0.0),
  bval2(0.0),

  // Constant
  buff(36),
  lu_in(17),

  // Flags (set in init)
  res(0),
  jac(0),
  mass(0)
{
  THERMINIT_F77(&res, &jac, &mass, &m, &mm, &phi(0), &psi(0), &w(0), 
		&matR(0,0), &matB(0,0), &matI(0,0), &matB2(0,0), &matB2R(0,0), 
		&part1(0,0), &lu_in, &prandtl, &rayleigh, &bval1, &bval2, 
		fname1.c_str(), fname2.c_str());

  initialGuess = NOX::LAPACK::Vector(3*mm+3, 3*m_max+3);

  for (int i=0; i<=mm; i++) {
    initialGuess(i) = phi(2*i);
    initialGuess(i+mm+1) = psi(2*i);
    initialGuess(i+2*mm+2) = w(2*i);
  }
 
}

const NOX::LAPACK::Vector&
ThermoSyphonProblemInterface::getInitialGuess()
{
  return initialGuess;
}

bool
ThermoSyphonProblemInterface::computeF(NOX::LAPACK::Vector& f, 
			       const NOX::LAPACK::Vector &x)
  {

  double *J = NULL;
  double *M = NULL;

  NOX::LAPACK::Vector& xx = const_cast<NOX::LAPACK::Vector&>(x);

  THERMFILL_F77(&res, &m_max, &mm, &xx(0), &w(0), &phi(0), &psi(0), 
		&convw(0,0), &convphi(0,0), &convpsi(0,0), &matB2R(0,0), 
		&rayleigh, &prandtl, &block12(0,0), &block13(0,0), 
		&block23(0,0), &block31(0,0), &part1(0,0), &lapphi(0), 
		&lappsi(0), &lapw(0), &phiw(0), &psiw(0), &matB(0,0), 
		&matR(0,0), &matI(0,0), &rhs_phi(0), &rhs_psi(0), &rhs_w(0), 
		&bval1, &bval2, &m, &buff, J, &f(0), M);
  
  return true;
}

bool
ThermoSyphonProblemInterface::computeJacobian(NOX::LAPACK::Matrix& J, 
					      const NOX::LAPACK::Vector & x)
{

  double *f = NULL;
  double *M = NULL;

  NOX::LAPACK::Vector& xx = const_cast<NOX::LAPACK::Vector&>(x);

  THERMFILL_F77(&jac, &m_max, &mm, &xx(0), &w(0), &phi(0), &psi(0), 
		&convw(0,0), &convphi(0,0), &convpsi(0,0), &matB2R(0,0), 
		&rayleigh, &prandtl, &block12(0,0), &block13(0,0), 
		&block23(0,0), &block31(0,0), &part1(0,0), &lapphi(0), 
		&lappsi(0), &lapw(0), &phiw(0), &psiw(0), &matB(0,0), 
		&matR(0,0), &matI(0,0), &rhs_phi(0), &rhs_psi(0), &rhs_w(0), 
		&bval1, &bval2, &m, &buff, &J(0,0), f, M);
 
  return true;
}

bool
ThermoSyphonProblemInterface::computeMass(NOX::LAPACK::Matrix& M, 
					  const NOX::LAPACK::Vector & x)
{

  double *f = NULL;
  double *J = NULL;

  NOX::LAPACK::Vector& xx = const_cast<NOX::LAPACK::Vector&>(x);
 
  THERMFILL_F77(&mass, &m_max, &mm, &xx(0), &w(0), &phi(0), &psi(0), 
		&convw(0,0), &convphi(0,0), &convpsi(0,0), &matB2R(0,0), 
		&rayleigh, &prandtl, &block12(0,0), &block13(0,0), 
		&block23(0,0), &block31(0,0), &part1(0,0), &lapphi(0), 
		&lappsi(0), &lapw(0), &phiw(0), &psiw(0), &matB(0,0), 
		&matR(0,0), &matI(0,0), &rhs_phi(0), &rhs_psi(0), &rhs_w(0), 
		&bval1, &bval2, &m, &buff, J, f, &M(0,0));
 
  return true;
}

void
ThermoSyphonProblemInterface::setParams(const LOCA::ParameterVector& p) {
  prandtl = p.getValue("prandtl");
  rayleigh = p.getValue("rayleigh");
}

LOCA::ParameterVector
ThermoSyphonProblemInterface::getParams() const {
  LOCA::ParameterVector p;
  p.addParameter("prandtl",prandtl);
  p.addParameter("rayleigh",rayleigh);

  return p;
}

int 
ThermoSyphonProblemInterface::getSize() const {
  return initialGuess.length();
}

void
ThermoSyphonProblemInterface::printSolution(const NOX::LAPACK::Vector &x,
					    const double conParam)
{
  int n = x.length();

  cout << "At parameter value: " << conParam << "   the solution vector is\n";
   
  if (n < 8) {
    for (int i=0; i<n; i++)  cout << " " << x(i);
  }
  else {
    for (int i=0; i<6; i++)  cout << " " << x(i);
    cout << " ...";
    for (int i=n-2; i<n; i++)  cout << " " << x(i);
  }
  cout << endl;

}
