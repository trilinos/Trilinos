// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/*! \file AnasaziBlockDavidson.hpp
 *  \brief Implementation of the block Jacobi-Davidson method
 *
 *  \fixme: STILL TO BE COMPLETED!
 */

#ifndef ANASAZI_BLOCK_JACOBI_DAVIDSON_HPP
#define ANASAZI_BLOCK_JACOBI_DAVIDSON_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziProjectedEigenproblem.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

/*!
\class Anasazi::BlockJacobiDavidson

\brief This class implements the block Jacobi-Davidson method, an iterative
method for solving eigenvalue problems.

The following discussion of the Jacobi-Davidson algorithm is derived from <I>
O. Chinellato, The Complex-Symmetric Jacobi-Davidson Algorithm and its
Application to the Computation of some Resonance Frequencies of Anisotropic
Lossy Axisymmetric Cavities, PhD Thesis, ETHZ, 2005</I>, Chapter 6.4.

The <I>Jacobi-Davidson Algorithm</I> was introduced by Sleijpen and van der
Vorst in 1996, and it is a typical subspace solver, designed to compute a
modest number of eigenpair approximations close to a given target \f$\tau\f$.
The eigenproblem to be solved is: Find <TT>nev</TT> eigenvectors of
\f[
A U = \lambda B U
\f]
that are close to \f$\tau\f$.

\fixme: something on A, B???

In an abstract fashion, the algorithm can be written as follows:
\f[
\begin{tabular}{p{1cm}{l}}
{\tt  1} &  Procedure JacobiDavidson$(A, B, \tau, U, \epsilon_{conv})$: \\
{\tt  2} &  Do: \\
{\tt  3} &  $\quad\quad$ $[\Theta, Z]$ := Extract$(A, B, U, \tau, \sigma)$ \\
{\tt  4} &  $\quad\quad$ $[\Lambda, Q, U, \Theta, Z, r]$ := Convergence/Deflation$(\Theta, Z, \epsilon_{conv})$ \\
{\tt  5} &  $\quad\quad$ $\tau$ := TargetSwitch$(r, \tau, \Theta)$ \\
{\tt  6} &  $\quad\quad$ Solve $(A - \sigma B) c \approx - r$ \\
{\tt  7} &  $\quad\quad$ $[U, \Theta, Z]$ := Restart($U, \Theta, Z, s_{min}, s_{max})$ \\
{\tt  8} &  $\quad\quad$ $U$ := Orthogonalize$(B, Q, U, C)$ \\
{\tt  9} &  EndProcedure 
\end{tabular}
\f]
This template is composed of several building blocks. We will now give some details on the implementation of each of them.

<P><B>Step 3: Subspace Extraction</B>. This step
finds eigenpairs residing in a given subspace \f$\mathcal{U}\f$,
which satify certain optimility conditions. Given the subspaces
\f$\mathcal{U}\f$, \f$A \mathcal{U}\f$ and \f$B \mathcal{U}\f$, several
strategies exist to extractthe (approximated) eigenpairs. Perhaps the best
known is the <I>Ritz-Galerkin approach</I>, which first defines the reduced 
operators:
\f[
\begin{tabular}{lcr}
$A_{\mathcal{U}}$ & = & $\mathcal{U}^H A \mathcal{U}$ \\
$B_{\mathcal{U}}$ & = & $\mathcal{U}^H B \mathcal{U}$, \\
\end{tabular}
\f]
then it computes the eigenpairs \f$(theta_i, z_i)\f$ of the reduced eigenproblem
\f[
A_{\mathcal{U}} z_i = \theta_i B_{\mathcal{U}} z_i.
\f]
Since this problem is dense, we can take advantage of LAPACK routines. This is done by the Anasazi::ProjectedEigenproblem class.

Given \f$z_i\f$, one can recover the eigenvalue corresponding to \f$\theta_i\f$
as
\f[
U_i = \mathcal{U} z_i.
\f]
The associated residuals can be expressed as \f$r(\theta_i, \mathcal{U} z_i) = A \mathcal{U} z_i - \theta_i B \mathcal{U} z_i\f$.


<P><B>Step 4: Convergence, Deflation</B>. Each eigenpair approximation obtained from the extract component is associated witha residual \f$r(\theta_i, \mathcal{U} z_i)\f$. If the corresponding <I>residual norm</I> is smaller than a required tolerence \f$\epsilon_{conv}\f$, i.e. if
\f[
\frac{
\| (\theta_i, \mathcal{U} z_i) \| 
}{
\| \mathcal{U} z_i  \|
} = 
\frac{
\| A \mathcal{U} z_i - \theta_i B \mathcal{U} z_i  \|
}{
\| \mathcal{U} z_i \|
} \leq \epsilon_{conv}
\f]
then the eigenpair is considered to be accurate enough and hence added to the
space of already converged eigenpairs \f$(\Lambda, Q)\f$. The Convergence/Deflation component iterates over the extracted eigenpair approximations one at a time and checks for thir convergence. Each accepted eigenpari is appended to the space of the already copnvernged ones and the subspace \f$\mathcal{U}\f$ is adapted appropriately. If the overlall number of converged eigenpairs reached the desired number of eigenpairs <TT>nev</TT>, the process terminates. 
 When an eigenpair is accepted, the search spaces are purged of all the components pointing in the directions of this pair. This guarantees that the same eigenvectors will not be found again, and it is called Deflation.
In class
Anasazi::BlockJacobiDavidson, \f$\Lambda\f$ and \f$Q\f$ are defined within the
input Anasazi::Eigenproblem class. 


<P><B>Step 5: Target Switch</B>. The Jacobi-Davidson algorithm is a typical shift-and-invert algorithm; however, in contrast to Krylov space methods, the shift is allowed to vary from iteration to iteration. Choosing the shift close to an actual eigenvalue of the matrix pencil \f$(A, B)\f$ typically yeilds good subspace extensions in the direction of the eigenvector associated with this eigenvalue, which in turn will speed up the convergence towards this eigenpair. A possible strategy to define the target switch is as follows. As soon as the relative norm falls below the requested threshold, the shift is chosen to equal the approximation eigenvalue \f$\theta_i\f$.


<P><B>Step 6: Correction Equation</B>. The Jacobi-Davidson algorithm requires,
  at each step, the expansion of the search space \f$\mathcal{U}\f$. Clearly,
  the quality of the approximations depends on the procedure adopted to expand
  this search space. A very elegant way to obtain good additional candidates
  consists in solving the so called <I>correction equation</I>.  This equation
  can be derived by considering the application of an eigenpair correction
  \f$(\delta, c)\f$ to the eigenpair approximation 
  \f$(\theta_1, \mathcal{U} z_1)\f$ in such a way
  that the new residual \f$r(\theta_1 + \delta, \mathcal{U} z_1 + c)\f$ will vanish. By expanding the "corrected" residual we have
\f[
A(\mathcal{U} z_1 + c) - (\theta_1 + \delta) B (\mathcal{U} z_1 + c) = 0.
\f]
We can replace \f$(\theta_1 + \delta)\f$ with either \f$theta_1\f$ 
(if the algorithm is close to converge to the desired eigenvalue, and
 therefore \f$\delta\f$ will be small) or \f$\tau\f$ 
(if the current approximation is too far from an actual eigenvalue). Therefore,
  the above equation can be rewritten as
\f[
(A - \chi B)c = - r(\theta_1, \mathcal{U} z_1).
\f]
with an appropriate \f$\chi\f$. Obviously, \f$c = \mathcal{U} z_1\f$ satisfies
this equation, yet without yeilding an additional direction. Therefore,
  the space where the additional correction vector is searched must be restricted. According to Sleijpen and van der Vorst (1996), the space being B-orthogonal to the matrix \f$\hat{Q} = (Q, q)\f$, where \f$q\f$ is the actual eigenvector approximation, represents a reasonable choice. Hence, the correction equation to be solved reads
\f[
(A - \chi B) c =  r(\theta_1, \mathcal{U} z_1) \mbox{ s.t. }
\hat{Q}^T B c = 0.
\f]
Note that this constraint offers the further advantage of inhibiting a potentially rapid condition number growth induced by the \f$\sigma's\f$ convergence towards a true eigenvalue. 


Note that the correction equation is the major advantage of the
Jacovi-Davidson algorthm over other iterative eigenvalue solvers. In contract
to other algorithms where the subspace is augmented following strict rules,
  the correction equation, which reflects the actual approximation status, is
  solved in each iteration. In this way the algorithm is guarantee to adapt to
  the actual situation which in turn promises a good performance.
            

<P><B>Step 7: Restart</B> 
In order to keep the storage requirements low, the dimension of the
subspace is reduced as soon as its dimension reaches an upper limit,
\f$s_{max}\f$. When doing so, the keep the eigenvectors approximations
associated with the \f$s_{min} -1\f$ best eigenpair approximations contained
in the subspace, where \f$s_{min}\f$ denotes the lower dimension limit of
the search space \f$\mathcal{U}\f$. This way, the algorithm maintains some
of the approximations gained during the previous iterations. In contrast to Krylov-based solvers, for which a restart constritutes an involved task, this step can be carried out easily and efficiently, and it is implemented as follows:
\f[
\begin{tabular}{p{1cm}{l}}
{\tt  1} &  Procedure Restart$(\mathcal{U}, \Theta, Z, s_{min}, s_{max})$ \\
{\tt  2} &  Do: \\
{\tt  3} &  $\quad\quad$If ($s = s_{max}$) Then \\
{\tt  4} &  $\quad\quad\quad\quad$\mathcal{U} := \mathcal{U} Z(1:s_{max}, 1:s_{min -1)$ \\
{\tt  5} &  $\quad\quad$End \\
{\tt  6} &  $\quad\quad$Return (\mathcal{U}, \Theta, Z)$ \\
{\tt  7} &  EndProcedure  \\
\end{tabular}
\f]

\author Oscar Chinellato (ETHZ/ICOS), Sabine Emch (ETHZ), Marzio Sala (ETHZ/COLAB)

\date Last updated on 01-Nov-05
*/
namespace Anasazi {

template <class ScalarType, class MV, class OP>
class GMRES {
protected:
	typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType magnitudeType;	
	
	int _LSitmax; //maximum of iterations
	int _LSitRestart; //number of iterations at which restarts occur
	int _LSit; //number of performed iterations
	int _restarts; //number of performed restarts
	magnitudeType _epslin; //tolerance
	
	Teuchos::BLAS<int,ScalarType> blas; //used to access Teuchos blas functions
	std::string _errorString; //error string (stores error information, if funcition solve returns Anasazi::Failed)
	
	static void UCosSin(const ScalarType x1, const ScalarType x2, magnitudeType &c, ScalarType &s); //calculates cos und sin for Givens Rotations
	static void GivensU(const magnitudeType c, const ScalarType &s, ScalarType &x1, ScalarType &x2); //applies a Givens rotation to 2 specific numbers
	
	static double conjugate(const double t) {return t;} //calculates the complex conjugate
	static std::complex<double> conjugate(const std::complex<double> &t) {return std::complex<double>(t.real(),-t.imag());} //calculates the complex conjugate
	
public:
	GMRES() : _LSitmax(500), _epslin(1e-6), _LSitRestart(20), _LSit(0), _restarts(0)  {} //contructor
	GMRES(const int itmax, const double tol, const int restart) : _LSitmax(itmax), _epslin(tol), _LSitRestart(restart) {
		if (_LSitmax<0) _LSitmax=500;
		if (_LSitRestart<1) _LSitRestart=20;	
	}	
	
	void setMaxIter (const int i) {if (i>=0) _LSitmax=i;} //sets the maximum of iterations
	void setTolerance(const typename Teuchos::ScalarTraits<ScalarType>::magnitudeType e) {_epslin=e;} //sets the tolerance
	void setRestart(const int i) {if (i>=1) _LSitRestart=i;} //sets the number of iterations at which restart occurs
	std::string getErrString() {return _errorString;} //returns a failure string, which is only set, if the function solve returns Anasazi::Failed (undefined content in other cases)
	
	int getIterations() {return _LSit;}
	int getRestarts() {return _restarts;}
	
	virtual Anasazi::ReturnType solve(OP &A, OP &Kinv, MV &b, MV &sol); //solves a specific equation A*sol=b with preconditioner Kinv
};

//----------------------------------------------------------------------------------------------
//                                 *** IMPLEMENTATION ***
//----------------------------------------------------------------------------------------------

#if 0
//solves a specific equation A*sol=b with preconditioner Kinv
template <class ScalarType, class MV, class OP>
Anasazi::ReturnType GMRES<ScalarType, MV, OP>::solve(OP &A, OP &Kinv, MV &b, MV &sol) {
	int k=0, i, LSit;
	magnitudeType rho;
	
	std::vector<magnitudeType> vm(1);
	std::vector<ScalarType> vs(1);
	std::vector<int> vi(1), vi2;
	
	MV *r, *ua, *uK;
	std::vector<MV*> u(_LSitRestart+1);
	Teuchos::SerialDenseMatrix<int,ScalarType> e, *z, R, s;
	Teuchos::SerialDenseMatrix<int,magnitudeType> c;
	Anasazi::ReturnType ok;
	
	_LSit=0;
	_restarts=0;
	//FIXME: eventuell auf mehr als einen Eingabe-Vektor erweitern
	if (b.GetNumberVecs()!=1) {
		_errorString="GMRES: only implemented for paramter b with 1 column";
		return Anasazi::Failed;
	}
	if (sol.GetNumberVecs()!=1) {
		_errorString="GMRES: only implemented for parameter sol with 1 column";
		return Anasazi::Failed;
	}
	if (b.GetVecLength()!=sol.GetVecLength()) {
		_errorString="GMRES: only implemented for b.GetLength()==sol.GetLength()";
		return Anasazi::Failed;
	}
	
	vi[0]=0;

	r=b.Clone(1); //(*) 1 col
	ok=A.Apply(sol,*r);
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator A to r";
		delete(r);
		r=0;
		return Anasazi::Failed;
	}
	r->MvAddMv((ScalarType)1,b,(ScalarType)-1,*r); //r:=b-A*x

	ok=Kinv.Apply(*r,*r); //r:=K^(-1)*(b-A*x)
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator Kinv to r";
		delete(r);
		r=0;
		return Anasazi::Failed;
	}
	r->MvNorm(&vm);
	rho=vm[0];
	e.shape(_LSitRestart+1,1); //alle Werte mit 0 initialisiert
	e(0,0)=rho;
	
	for (i=0; i<=_LSitRestart; i++) u[i]=r->Clone(1);
	if (rho!=(magnitudeType)0) u[0]->MvAddMv((ScalarType)1/rho,*r,(ScalarType)0,*r); //u[0]=r/rho
	else u[0]->MvAddMv((ScalarType)1,*r,(ScalarType)0,*(u[0]));
	
	ua=r->Clone(1); //(*) 1 col
	uK=r->Clone(1); //(*) 1 col
	R.shape(_LSitRestart+1,_LSitRestart+1); //(m) (_LSitRestart+1 x _LSitRestart+1)
	c.shapeUninitialized(_LSitRestart,1); //(m) (_LSitRestart x 1)
	s.shapeUninitialized(_LSitRestart,1); //(m) (_LSitRestart x 1)

	k=0;
	for (LSit=1; LSit<=_LSitmax; LSit++) {
		ok=A.Apply(*(u[k]),*ua); //ua:=A*u[k]
		assert(ok==Anasazi::Ok);
		ok=Kinv.Apply(*ua,*uK); //uK:=Kinv*ua
		assert(ok==Anasazi::Ok);
		
		for (i=0; i<=k; i++) {
			u[i]->MvDot(*uK,&vs);
			R(i,k)=vs[0]; //R[i,k]:=u[i]^(H)*u[k]
			uK->MvAddMv((ScalarType)1,*uK,-R(i,k),*(u[i])); //uK:=uK-R[i,k]*u[i]
		}
	
		uK->MvNorm(&vm);
		R(k+1,k)=vm[0]; //R[k+1,k]=norm(u[k+1]);
		if (vm[0]!=0) u[k+1]->MvAddMv((ScalarType)1/R(k+1,k),*uK,(ScalarType)0,*uK); //u[k+1]:=uK/R[k+1,k]
		else u[k+1]->MvAddMv((ScalarType)1,*uK,(ScalarType)0,*uK);
		
		for (i=0; i<k; i++) {
			GivensU(c(i,0),s(i,0),R(i,k),R(i+1,k));
		}
		UCosSin(R(k,k),R(k+1,k),c(k,0),s(k,0));
		GivensU(c(k,0),s(k,0), R(k,k), R(k+1,k));
		
		//R(k+1,k)=(ScalarType)0; //Da ich weiss, dass dies 0 werden soll //unnoetig, wenn ich bei linSolve Dreiecksmatrix verwende
		
		GivensU(c(k,0),s(k,0), e(k,0),e(k+1,0));
			
		if (Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1,0))<=rho*_epslin) break;
		
		if (k+1==_LSitRestart) {
			_restarts++;
			
			//printf("GMRES: RESTART\n");
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, k+1, 1, (ScalarType)1, R.values(), R.stride(), e.values(), e.stride());
		
			for (i=0; i<=k; i++) sol.MvAddMv((ScalarType)1,sol,e(i,0),*(u[i]));  //x:=x+U[:,0:k]*z
			
			ok=A.Apply(sol,*r);
			assert(ok==Anasazi::Ok);
			r->MvAddMv((ScalarType)1,b,(ScalarType)-1,*r); //r:=b-A*x
			ok=Kinv.Apply(*r,*r); //r:=Kinv*r
			assert(ok==Anasazi::Ok);
			r->MvNorm(&vm); //eta:=vm[0]:=norm(r);
			
			assert(vm[0]!=(magnitudeType)0);
			e.shape(_LSitRestart+1,1); //alle Werte mit 0 initialisiert
			e(0,0)=vm[0];
			
			
			u[0]->MvAddMv((ScalarType)1/vm[0],*r,(ScalarType)0,*(u[0])); //u[0]:=r/eta
			
			k=-1;
		}
		
		k++;
	}
	_LSit=LSit;
	if ((LSit<=_LSitmax) && (rho!=0)) {
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, k+1, 1, (ScalarType)1, R.values(), R.stride(), e.values(), e.stride());
		
		for (i=0; i<=k; i++) sol.MvAddMv((ScalarType)1,sol,e(i,0),*(u[i])); //x:=x+U[:,0:k]*z
	}
	
	delete(r); //(#)
	r=0;
	for (i=0; i<=_LSitRestart; i++) delete(u[i]);
	delete(ua); //(#)
	ua=0;
	delete(uK); //(#)
	uK=0;
	
	if (LSit<=_LSitmax) {
		return Anasazi::Ok;
	} else {
		return Anasazi::Unconverged;
	}
}
#else
template <class ScalarType, class MV, class OP>
Anasazi::ReturnType GMRES<ScalarType, MV, OP>::solve(OP &A, OP &Kinv, MV &b, MV &x) {

  int k=0, i, j, info;
  ScalarType work[2*_LSitRestart];
  magnitudeType a, rho, eta;
	
  std::vector<magnitudeType> vm(1);
  std::vector<ScalarType> vs(1);
  std::vector<int> vi(1), vi2;
	
  MV *r, *u, *uview, *uA, *uK, *y;
  Teuchos::SerialDenseMatrix<int,ScalarType> e, *z, R, *eview, *Rview, *Rinv, s;
  Teuchos::SerialDenseMatrix<int,magnitudeType> c;
  Anasazi::ReturnType ok;
	
  //FIXME: eventuell auf mehr als einen Eingabe-Vektor erweitern
  assert(b.GetNumberVecs() == 1);
  assert(x.GetNumberVecs() == 1);
  assert(b.GetVecLength() == x.GetVecLength());

  u = b.Clone(1);                                       // u = b
  ok = A.Apply(x, *u);                                  // u = A*x
  assert(ok == Anasazi::Ok);
  u->MvAddMv((ScalarType)1, b, (ScalarType)-1, *u);     // u = b - u = b - A*x

  r = u->Clone(1);                                      
  ok = Kinv.Apply(*u, *r);                            // r = Kinv*u
  assert(ok == Anasazi::Ok);
  r->MvNorm(&vm);                            
  assert(vm[0] != 0.0);
  rho=vm[0];                                         // rho = ||r||
  delete(u);

  e.shape(_LSitRestart+1,1);                         // e = zeros(_LSitRestart+1, 1)
  e(0,0)=rho;                                        // e(1) = rho

  u=r->Clone(_LSitRestart+1);                        // u = zeros(x.GetVecLength(),_LSitRestart+1)
  vi[0]=0;
  uview=u->CloneView(vi);                                         // uview ~ u(:,1)
  uview->MvAddMv((ScalarType)1.0/rho, *r, (ScalarType)0.0, *uview);  // uview = uview/rho = u(:,1)/rho
  delete(uview);
	
  uA=r->Clone(1);                                    // uA ~ zeros(DIM, 1)
  uK=r->Clone(1);                                    // uK ~ zeros(DIM, 1)
  R.shape(_LSitRestart+1,_LSitRestart+1);            //  R ~ zeros(_LSitRestart+1, _LSitRestart+1)
  c.shape(_LSitRestart,1);                           //  c ~ zeros(_LSitRestart+1, 1)
  s.shape(_LSitRestart,1);                           //  s ~ zeros(_LSitRestart+1, 1)
	

  _restarts = 0;
  k=0;
  for (_LSit=1; _LSit<=_LSitmax; _LSit++) {
    vi[0] = k;
    uview = u->CloneView(vi); 
    ok = A.Apply(*uview, *uA);                          // uA = A*u[k]
    assert(ok==Anasazi::Ok);
    ok = Kinv.Apply(*uA, *uK);                          // uK = Kinv*uA
    assert(ok==Anasazi::Ok);
    delete(uview); 
		
    for (i=0; i<=k; i++) {
      vi[0] = i;
      uview = u->CloneView(vi); //(*)
      uK->MvDot((*uview), &vm);
      R(i,k) = vm[0];                                   //R[i+1,k] = u[i]^(H)*uK
      uK->MvAddMv((ScalarType)1,*uK,-R(i,k),*uview);    //u[k] = u[k] - R[i+1,k]*u[i]
      delete(uview); 
    }
		
    vi[0] = k+1;
    uview = u->CloneView(vi); 
    uK->MvNorm(&vm);
    assert(vm[0] != 0.0);
    R(k+1,k) = vm[0];                                                     //R[k+1,k] = ||uK||
    uview->MvAddMv((ScalarType)1/R(k+1,k), *uK, (ScalarType)0, *uview);   //u[k+1] = uK/R[k+1,k]
    delete(uview);
		
    for (i=0; i<k; i++) {
      GivensU(c(i,0), s(i,0), R(i,k), R(i+1,k));
    }
    UCosSin(R(k,k), R(k+1,k), c(k,0), s(k,0));
    GivensU(c(k,0), s(k,0), R(k,k), R(k+1,k));
    // R(k+1,k)=(ScalarType)0;
    GivensU(c(k,0), s(k,0), e(k,0), e(k+1,0));
		
    /// printf("%20.18g\n",Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1,0))/rho);

    if (Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1,0)) <= rho*_epslin) break;
		
    if ((k+1) == _LSitRestart) { //im Komplexen Fall auch fuer k=0 noetig

      //printf("Restart\n");			

      // FIXME: add underscore
      blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
		 k+1, k+1, (ScalarType)1.0, R.values(), _LSitRestart + 1, e.values(), _LSitRestart + 1);         // e = R\e
      
      vi.resize(k+1);
      for (i=0; i<=k; i++) vi[i]=i;
      uview = u->CloneView(vi);
      x.MvTimesMatAddMv((ScalarType)1, *uview, e, (ScalarType)1);         // x = x + U[:,1:k]*e
      delete(uview);
      vi.resize(1);
      
      ok = A.Apply(x,*r);                                  // r = A*x
      assert(ok == Anasazi::Ok);
      r->MvAddMv((ScalarType)1,b,(ScalarType)-1,*r);       // r = b - r = b - A*x
      
      vi[0]=0;
      uview = u->CloneView(vi);
      ok = Kinv.Apply(*r,*uview);                          // u(:,1) = Kinv*r
      assert(ok == Anasazi::Ok);
      uview->MvNorm(&vm);                            
      assert(vm[0] != 0.0);
      eta=vm[0];                                           // eta = ||r||
      uview->MvAddMv((ScalarType)1.0/eta, *uview, (ScalarType)0.0, *uview);  // uview = uview/rho = u(:,1)/rho
      delete(uview);
      
      e.shape(_LSitRestart+1,1);                         // e = zeros(_LSitRestart+1, 1)
      e(0,0)=eta;                                        // e(1) = rho
 
      _restarts++;
      k=-1;
    }
		
    k++;
  }

  if ((k != 0) || (_LSit != _LSitmax)) {
    
    blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
	       k+1, 1, (ScalarType)1.0, R.values(), _LSitRestart + 1, e.values(), _LSitRestart + 1);         // e = R\e

    vi.resize(1);
    for(i=0; i<=k; ++i){
      vi[0] = i;
      uview = u->CloneView(vi);
      x.MvAddMv((ScalarType)1.0, x, e(i,0), *uview);                 // x = x + U[i,1:k]*e(i)
      delete(uview);
    }
  }

	
  delete(r); //(#)
  delete(u); //(#)
  delete(uA); //(#)
  delete(uK); //(#)

  return(Anasazi::Ok);
}
#endif

//calculates cos und sin for Givens Rotations
template <class ScalarType, class MV, class OP>
void GMRES<ScalarType,MV,OP>::UCosSin(const ScalarType x1, const ScalarType x2, magnitudeType &c, ScalarType &s) {
	ScalarType r;
	magnitudeType l, l1, l2;
	
	l1=Teuchos::ScalarTraits<ScalarType>::magnitude(x1);
	
	if (l1 == (magnitudeType)0){
		c = (magnitudeType)0;
		s = (ScalarType)1;
	} else {
		l2=Teuchos::ScalarTraits<ScalarType>::magnitude(x2);
		l=sqrt(l1*l1+l2*l2); //FIXME: koennte Problem wegen Overflow geben
	
		c = l1/l;
		s = (x2/x1)*c;
	}
}

//applies a Givens rotation to 2 specific numbers
template <class ScalarType, class MV, class OP>
void GMRES<ScalarType, MV, OP>::GivensU(const magnitudeType c, const ScalarType &s, ScalarType &x1, ScalarType &x2) {
	ScalarType h=x1;
	
	x1=c*x1+conjugate(s)*x2;
	x2=c*x2-s*h;
}


  
template <class ScalarType, class MV>
class JacDavOp1 : public virtual Anasazi::Operator<ScalarType> {
protected:
	MV& _Q,& _Qb;
public:
	JacDavOp1(MV &Q, MV &Qb) : _Q(Q), _Qb(Qb) {}

  //
  // mv2 = (I - Q*Qb^T)*mv1
  //
  virtual Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType> &mv1, Anasazi::MultiVec<ScalarType> &mv2) const {
    int a=mv1.GetVecLength();
    Teuchos::SerialDenseMatrix<int,ScalarType> H;
    
    if ((mv1.GetNumberVecs()!=mv1.GetNumberVecs()) || (a!=mv2.GetVecLength()) || (a!=_Qb.GetVecLength()) || (a!=_Q.GetVecLength()) || (_Q.GetNumberVecs()!=_Qb.GetNumberVecs())) return Anasazi::Failed;
    
    H.shapeUninitialized(_Qb.GetNumberVecs(), mv1.GetNumberVecs());
    mv1.MvTransMv((ScalarType)1,_Qb,H);
    mv2.MvAddMv((ScalarType)1,mv1,(ScalarType)0,mv2);
    mv2.MvTimesMatAddMv((ScalarType)-1,_Q,H,(ScalarType)1);
    
    return Anasazi::Ok;
  }
};

template <class ScalarType, class MV, class OP>
class JacDavOp2 : public virtual Anasazi::Operator<ScalarType> {
protected:
	MV& _Qb,& _Q;
	Teuchos::RefCountPtr<OP> _A, _B;
	ScalarType _sigma;
public:
	JacDavOp2(MV &Qb, MV &Q, Teuchos::RefCountPtr<OP> &A, ScalarType sigma, Teuchos::RefCountPtr<OP> &B) : _Q(Q), _Qb(Qb) {
		_A=A;
		_B=B;
		_sigma=sigma;
	}
	
	virtual Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType> &mv1, Anasazi::MultiVec<ScalarType> &mv2) const {
		Anasazi::ReturnType info;
		Teuchos::SerialDenseMatrix<int,ScalarType> H;
		Anasazi::MultiVec<ScalarType> *h1;
		
		//mv2:=(A-sigma*B)*mv1
		info=_A->Apply(mv1,mv2);
		if (info!=Anasazi::Ok) return info;
		
		h1=mv2.Clone(mv2.GetNumberVecs());
		info=_B->Apply(mv1,(*h1));
		if (info!=Anasazi::Ok) return info;
		
		mv2.MvAddMv((ScalarType)1,mv2,-_sigma,(*h1));
		
		//H:=Q^(T)*(A-sigma*B)*mv1
		H.shapeUninitialized(_Q.GetNumberVecs(), mv1.GetNumberVecs());
		mv2.MvTransMv((ScalarType)1,_Q,H);
				
		//mv2:=(I-Qb*Q^(T))*(A-sigma*B)*mv1
		mv2.MvTimesMatAddMv((ScalarType)-1,_Qb,H,(ScalarType)1);
		
		return Anasazi::Ok;
	}
};

template <class ScalarType, class MV, class OP>
class JacDavOp4 : public virtual Anasazi::Operator<ScalarType> {
protected:
	MV& _Qk,& _Qb;
	Teuchos::SerialDenseMatrix<int,ScalarType> _LkTrans, _Rk;
	Teuchos::RefCountPtr<OP> _Kinv;
	Teuchos::BLAS<int,ScalarType> blas;
public:
	JacDavOp4(MV &Qk, Teuchos::SerialDenseMatrix<int,ScalarType> LkTrans, Teuchos::SerialDenseMatrix<int,ScalarType> Rk, MV &Qb, Teuchos::RefCountPtr<OP> &Kinv)  : _Qk(Qk), _Qb(Qb), _LkTrans(LkTrans), _Rk(Rk) {
		_Kinv=Kinv;
	}
	
	virtual Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType> &mv1, Anasazi::MultiVec<ScalarType> &mv2) const {
		Anasazi::ReturnType info;
		Teuchos::SerialDenseMatrix<int,ScalarType> H;
		
		//mv2:=Kinv*mv1
		info=_Kinv->Apply(mv1,mv2);
		if (info!=Anasazi::Ok) return info;
		
		//H:=Qb^(T)*Kinv*mv1
		H.shapeUninitialized(_Qb.GetNumberVecs(), mv1.GetNumberVecs());
		mv2.MvTransMv((ScalarType)1,_Qb,H);
		
		//H=Lk^(-1)*Qb^(T)*Kinv*mv1 -> solve Lk*x=H
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, H.numRows(), H.numCols(), (ScalarType)1, _LkTrans.values(), _LkTrans.stride(), H.values(), H.stride());
		
		//H=Rk^(-1)*Lk^(-1)*Qb^(T)*Kinv*mv1 -> solve Rk*x=H
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, H.numRows(), H.numCols(), (ScalarType)1, _Rk.values(), _Rk.stride(), H.values(), H.stride());
		
		//mv2:=mv1+Qk*Rk^(-1)*Lk^(-1)*Qb^(T)*Kinv*mv1
		mv2.MvTimesMatAddMv((ScalarType)-1, _Qk, H, (ScalarType)1);
		
		return Anasazi::Ok;
	}
};


#define MAX(x,y) (x)>=(y)?(x):(y)

  template<class ScalarType, class MagnitudeType,  class MV, class OP>
  class BlockJacobiDavidson : public Eigensolver<ScalarType,MV,OP> 
  {
  public:
    //@{ \name Constructor/Destructor.

    //! %Anasazi::Block Jacobi-Davidson constructor.
    BlockJacobiDavidson (const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                   const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                   Teuchos::ParameterList &pl
                  );

    //! %Anasazi::BlockDavidson destructor.
    virtual ~BlockJacobiDavidson() {};
    //@}

    //@{ \name Solver application methods.
    /*! 
     * \brief This method uses iterate to compute approximate solutions to the
     * original problem.  It may return without converging if it has
     * taken the maximum number of iterations or numerical breakdown is
     * observed.
     */
    ReturnType solve();
    //@}

    //@{ \name Solver status methods.

    //! Get the current restart count of the iteration method.
    int GetNumRestarts() const 
    {
      return(_numRestarts); 
    }

    //! Get the number of iterations.
    int GetNumIters() const 
    {
      return(_iters); 
    }

    /*! \brief Get a constant reference to the current linear problem, 
      which may include a current solution.
    */
    Eigenproblem<ScalarType,MV,OP>& GetEigenproblem() const { return(*_problem); };
    
    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int GetBlockSize() const 
    { 
      return(_blockSize); 
    }

    //@}
    //@{ \name Output methods.

    //! This method requests that the solver print out its current status to screen.
    void currentStatus();

    //@}
  private:

    /*! \brief These methods will not be defined.     */
    BlockJacobiDavidson(const BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP> &method);
    BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>& operator=
      (const BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP> &method);

    // 
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    //! Problem to be solved.
    Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
    //! Output manager to be used for any output operation.
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _om;
    //! Sorting manager to be used to sort the eigenvalues.
    Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sm;
    //! Parameter list containing options for this solver.
    Teuchos::ParameterList _pl;
    //! Output stream, reference from the output manager
    std::ostream& _os;
    
    // @}
    // @{ \name Information obtained from the eigenproblem
    
    //! Operator A
    Teuchos::RefCountPtr<OP> _A;
    //! Operator B
    Teuchos::RefCountPtr<OP> _B;
    //! Preconditioner
    Teuchos::RefCountPtr<OP> _Prec;
    //! Multi-vector that will contain the computed eigenvectors.
    Teuchos::RefCountPtr<MV> _evecs;
    //! Will contain the computed eigenvalues.
    Teuchos::RefCountPtr<std::vector<ScalarType> > _evals;
    
    // @}
    // @{ \name Internal data
    
    //! Tolerance for convergence
    MagnitudeType _residual_tolerance;
    //! Target, the solver will compute the eigenvalues close to this value.
    ScalarType _TARGET;
    //! Maximum number of allowed iterations.
    int _maxIters;
    //! Block size.
    int _blockSize;
    int _SMIN;
    int _SMAX;
    //! Number of requested eigenvalues.
    int _nev;
    //! Number of performed restarts.
    int _numRestarts;
    //! Number of performed iterations.
    int _iters;
    //! Number of computed and accepted eigenvalues.
    int _knownEV;
    //! Residual norm of computed and accepted eigenvalues.
    std::vector<MagnitudeType> _normR;
    bool _haveKinv;
    int _LSIterMax, _LSRestart, _gamma;
    double _elin, _eswitch;
    
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
  };
  
  //
  // Implementation
  //
  template <class ScalarType, class MagnitudeType, class MV, class OP>
  BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>::
  BlockJacobiDavidson(const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                        const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                        const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                        Teuchos::ParameterList &pl
                       ):
    _problem(problem), 
    _om(om),
    _sm(sm),
    _pl(pl),
    _os(_om->GetOStream()),
    _A(_problem->GetOperator()),
    _B(_problem->GetM()),
    _Prec(_problem->GetPrec()),
    _evecs(_problem->GetEvecs()),
    _evals(problem->GetEvals()), 
    _nev(problem->GetNEV()),
    _maxIters(_pl.get("Max Iters", 100)),
    _blockSize(_pl.get("Block Size", 1)),
    _residual_tolerance(_pl.get("Tol", 10e-8)),
    _numRestarts(0),
    _iters(0),
    _TARGET(_pl.get("Target", 0.0)),
    _SMIN(_pl.get("SMIN", 20)), 
    _SMAX(_pl.get("SMAX", 30)),
    _knownEV(0),
    _LSIterMax(_pl.get("KrylovSolver: MaxIters", 100)),
    _LSRestart(_pl.get("KrylovSolver: Restart",100)),
    _elin(_pl.get("KrylovSolver: Tolerance", 10e-8)),
    _gamma(_pl.get("gamma", 2)),
    _eswitch(_pl.get("Switch", 0.1))
  {
  }

  template <class ScalarType, class MagnitudeType, class MV, class OP>
  ReturnType BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>::solve()
  {
    // ============= //
    // Sanity checks //
    // ============= //
    
    // Check the Anasazi::Eigenproblem was set by user, if not, return failed.
    if (!_problem->IsProblemSet()) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error))
        _os << "ERROR : Anasazi::Eigenproblem was not set, call Anasazi::Eigenproblem::SetProblem() before calling solve"<< endl;
      return Failed;
    }
    
    // Symmetric or Hermitian, but not general
    // Check the Anasazi::Eigenproblem is symmetric, if not, return failed.
    if (!_problem->IsSymmetric()) 
    {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Anasazi::Eigenproblem is not symmetric" << endl;
      return Failed;
    }
    
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    Teuchos::RefCountPtr<MV> iVec = _problem->GetInitVec();
    
    if (iVec.get() == 0) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : Initial vector is not specified, set initial vector in eigenproblem "<<endl;
      return Failed;
    }
  
    if (iVec->GetNumberVecs() != _blockSize) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : Search space has incorrect dimension" << endl;
      return Failed;
    }
    
    // Check that the maximum number of blocks for the eigensolver is a positive number
    if (_blockSize <= 0) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : blockSize = "<< _blockSize <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    
    // Check that the maximum number of iterations is a positive number
    if (_maxIters <= 0) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : maxIter = "<< _maxIters <<" [ should be positive number ] " << endl;
      return Failed;
    } 

    // Check that the eigenvalue vector is allocated with correct size
    if (_evals->size() < _nev) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : evals->size() is less than nev" << endl;
      return Failed;
    } 

    // =============================== //
    // templated values of 1's and 0's //
    // =============================== //

    ScalarType ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();

    // ===================== //
    // Start with allocation //
    // ===================== //
    
    // Allocate the subspace U (and AU, BU), the Ritz-Galerkin matrix and the Ritz values ...
    int SearchSpaceSize = 0, outer=0;
    std::vector<int> sublist;

    Teuchos::RefCountPtr<MV> U, AU, BU;
    U  = MVT::Clone(*iVec, _SMAX);
    AU = MVT::Clone(*iVec, _SMAX);
    BU = MVT::Clone(*iVec, _SMAX);

    // Convenience multivectors, used to extract views of U, AU, BU
    MV *Utmp;
    MV *AUtmp;
    MV *BUtmp;
    MV *Rtmp;

    _normR.resize(_nev); // container for norms
    
    // working array
    std::vector<ScalarType> theta(_SMAX), vc(1);
    std::vector<int> perm(_SMAX), iperm(_SMAX), vQ2;
    Teuchos::SerialDenseMatrix<int, ScalarType> Z(_SMAX, _SMAX), LkTrans(_nev,_nev), Rk(_nev,_nev);
    Teuchos::SerialDenseMatrix<int, ScalarType> Ztmp, *lk, *rk, *Ztmp2;
    Anasazi::Operator<ScalarType> *Ksys, *Asys;

    // Allocate space for the residuals ...
    Teuchos::RefCountPtr<MV> R;
    R = MVT::CloneCopy(*iVec);

    // Allocate space for the result ...
    int conv = 0;   // number of converged eigenpairs per iteration
    int purged = 0; // number of deflated eigenpairs per iteration
    std::vector<ScalarType> sigma(_blockSize);
    MV *Qtmp, *Qtmp2;
    MV *BQtmp, *BQtmp2;
    MV *KBQtmp, *KBQtmp2, *h;

    Teuchos::RefCountPtr<MV> BQ, KBQ;
    Teuchos::BLAS<int,ScalarType> blas;
    Anasazi::ReturnType ok;
    
    BQ = MVT::Clone(*iVec, _nev);
    KBQ = MVT::Clone(*iVec, _nev);
    
    // Instantiate the ProjectedEigenproblem component and set some of the values ...
    // FIXME: for Hermitian problems...
    ProjectedEigenproblem<int, ScalarType, MV> PE("Symmetric", _SMAX);
    PE.SetTheta(&theta);
    PE.SetZ(&Z);

    _knownEV = 0;
    SearchSpaceSize = 0;
    _iters = 0;

    if (_om->doPrint()) 
    {
      _os << "[Starting Solver]" << endl;
    }

    // ================================================ //
    // Start the (Block-) Jacobi/Davidson iteration ... //
    // ================================================ //
    
    while(_iters < _maxIters)
    {
      // ================================================== //
      // Add the (preconditioned) residuals to the search   //
      // space. To this end, orthogonalise w.r.t. foregoing //
      // subspace vectors and w.r.t. found eigenvectors Q   //
      // ================================================== //
      
      sublist.resize(1);
      std::vector<ScalarType> eta(1);
      std::vector<MagnitudeType> nrm(1);
      std::vector<MagnitudeType> nrm_q(1);

      for(int i=0; i < _blockSize; i++)
      {
        sublist[0] = i;
        Rtmp = R->CloneView(sublist); //(*)

        // Gram-Schmidt reduce the vector ...
        for(int j=0; j < SearchSpaceSize ; j++){
          sublist[0] = j;
          Utmp  = U->CloneView(sublist); //(*)
          BUtmp = BU->CloneView(sublist); //(*)
          Rtmp->MvPseudoDot((*BUtmp), &eta);
          Rtmp->MvAddMv(ScalarOne, (*Rtmp), -eta[0], (*Utmp));
          delete Utmp, BUtmp; //(#)(#)
        }

        for(int j=0; j < _knownEV; j++){
          sublist[0] = j;
          Qtmp = _evecs->CloneView(sublist);
          BQtmp = BQ->CloneView(sublist); //(*)
          Rtmp->MvPseudoDot((*BQtmp), &eta); 
          Rtmp->MvAddMv(ScalarOne, (*Rtmp), -eta[0], (*Qtmp));	      
          delete Qtmp, BQtmp; //(#)(#)
        }
	
        // Now B-normalise the vector ...
        sublist[0] = _knownEV;
	Qtmp = _evecs->CloneView(sublist); //(*)
        if (_B.get() != 0)
          _B->Apply((*Rtmp), (*Qtmp));
        else
          Qtmp->SetBlock(*Rtmp, sublist);
	
	Rtmp->MvPseudoDot((*Qtmp), &eta);
        eta[0] = Teuchos::ScalarTraits<ScalarType>::squareroot(eta[0]);
        // scaling of Rtmp
        Rtmp->MvAddMv(ScalarOne / eta[0], (*Rtmp), ScalarZero, (*Rtmp));
	
        sublist[0] = SearchSpaceSize;
        U->SetBlock((*Rtmp),sublist);
	
        BUtmp = BU->CloneView(sublist);
        if (_B.get() != 0)
          _B->Apply((*Rtmp), (*BUtmp));
        else
          BUtmp->SetBlock((*Rtmp), sublist);

        delete Rtmp, Qtmp, BUtmp; //(#)(#)
        SearchSpaceSize++;
      }
      
      // Update AU
     
      sublist.resize(_blockSize);
      for(int i=0; i<_blockSize; ++i) sublist[i]=(SearchSpaceSize - _blockSize) + i;

      Utmp = U->CloneView(sublist); //(*)
      AUtmp = AU->CloneView(sublist); //(*)
      BUtmp = BU->CloneView(sublist); //(*)
      _A->Apply((*Utmp), (*AUtmp));
      
      // Update the ProjectedEigenproblem component by telling it
      // about the space increase
      PE.Add((*Utmp), (*AUtmp), (*BUtmp));
      delete Utmp; //(#)
      delete AUtmp; //(#)
      delete BUtmp; //(#)

      // ====================================================== //
      // Extract eigenpair approximations from the space, then  //
      // Sort eigenpairs with respect to TARGET and compute the //
      // inverse permutation.                                   //
      // CAUTION: At return of the sorting manager, the thetas  //
      // are sorted!!                                           //
      // ====================================================== //

      PE.Extract();

      for(int i=0; i<SearchSpaceSize ; ++i) {theta[i] -= _TARGET;}
      _sm->sort((Eigensolver<ScalarType,MV,OP>*)NULL, SearchSpaceSize , &theta[0], &perm);
      
      for(int i=0; i<SearchSpaceSize ; ++i) {theta[i] += _TARGET;}
      for(int i=0; i<SearchSpaceSize ; ++i) { iperm[perm[i]]=i;}

      // Permute the Z entries according to perm
      Ztmp.shape(SearchSpaceSize ,1);

      for(int j=0; j<SearchSpaceSize ; ++j){
        if (perm[j] != j){
          for(int i=0; i<SearchSpaceSize ; ++i) {Ztmp[0][i] = Z[j][i];}
          for(int i=0; i<SearchSpaceSize ; ++i) {Z[j][i] = Z[perm[j]][i];}
          for(int i=0; i<SearchSpaceSize ; ++i) {Z[perm[j]][i] = Ztmp[0][i];}

          perm[iperm[j]] = perm[j];	  
        }
      }
      
      // Check for convergence of the approximative eigenpairs. If
      // some of them are accurate enough, add them to the space _evecs
      Ztmp.shape(SearchSpaceSize, 1);

      sublist.resize(SearchSpaceSize);
      for(int i=0; i<SearchSpaceSize ; ++i) sublist[i]=i;	  

      Utmp = U->CloneView(sublist); //(*)
      AUtmp = AU->CloneView(sublist); //(*)
      BUtmp = BU->CloneView(sublist); //(*)
      
      sublist.resize(1);
      sublist[0] = 0;
      Rtmp = R->CloneView(sublist); //(*)
     
      for(int i=0; i<SearchSpaceSize ; ++i) {Ztmp[0][i] = Z[0][i];}
      
      // compute BQtmp
      sublist[0] = _knownEV;
      BQtmp = BQ->CloneView(sublist);
      BQtmp->MvTimesMatAddMv(ScalarOne, *BUtmp, Ztmp, ScalarZero);

      // compute K^{-1} BQtmp
      KBQtmp = KBQ->CloneView(sublist);
      _Prec->Apply(*BQtmp, *KBQtmp); //Qk=[_Qk _qk]

      // Now compute residual
      // Rtmp = AU*Z(:,1)
      Rtmp->MvTimesMatAddMv (ScalarOne, (*AUtmp), Ztmp, ScalarZero);
      // Rtmp = Rtmp - theta[0]*AU*Z(:,1) = Rtmp - theta[0] * BQtmp
      Rtmp->MvAddMv (ScalarOne, *Rtmp, -theta[0], *BQtmp);
      // nrm  = ||Rtmp||
      Rtmp->MvNorm(&nrm);
      Qtmp = _evecs->CloneView(sublist);
      // Qtmp = U*Z(:,1)
      Qtmp->MvTimesMatAddMv (ScalarOne, (*Utmp), Ztmp, ScalarZero);

      // nrm  = ||Rtmp||
      Qtmp->MvNorm(&nrm_q);
      nrm[0] /= nrm_q[0];

      if (_om->doPrint()) 
      {
        _os << "It: " << _iters << ", knownEV: " << _knownEV;
       // _os << ", s: " << s << ", rel_normr = " << nrm[0];
        _os << ", SearchSpaceSize : " << SearchSpaceSize  << ", rel_normr = " << nrm[0]/nrm_q[0];
        _os << ", theta_i = ";
        _os << "(" << theta[0] << ")";
        int j = ANASAZI_MIN(SearchSpaceSize ,5);
        for(int i=1; i<j; ++i){
          _os << ", (" << theta[i] << ")";
        }
        _os << endl;
      }
      
      conv = 0;
    
      while (nrm[0] <= _residual_tolerance) 
      {
        // Qtmp is already in the proper state...
	
        (*_evals)[_knownEV] = theta[conv];
	
        _normR[_knownEV] = nrm[0]; // relative norm
        _knownEV++;

        if (_knownEV == _nev) 
        {
          break;
        }
        conv++;
	
	outer=0;
	
	if (_Prec.get()) 
        {
		if (_knownEV==1) {
			LkTrans(0,0)=ScalarOne;
			BQtmp->MvPseudoDot(*KBQtmp, &vc);
			Rk(0,0)=vc[0];
		} else {
			lk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, LkTrans[_knownEV-1], LkTrans.stride(), _knownEV-1, 1); //(*)
			vQ2.resize(_knownEV-1);
			for (int i=0; i<_knownEV-1; i++) vQ2[i]=i;
			h=KBQ->CloneView(vQ2); //(*)
			BQtmp->MvTransMv(ScalarOne, *h, *lk);
			delete(h); //(#)
						
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, lk->numRows(), lk->numCols(), ScalarOne, Rk.values(), Rk.stride(), lk->values(), lk->stride());
			
			rk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Rk[_knownEV-1], Rk.stride(), _knownEV-1, 1); //(*)
			h=KBQ->CloneView(vQ2); //(*)
			KBQtmp->MvTransMv(ScalarOne, *h, *rk);
			delete(h); //(#)
			
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, rk->numRows(), rk->numCols(), ScalarOne, LkTrans.values(), LkTrans.stride(), rk->values(), rk->stride());
			
			BQtmp->MvPseudoDot(*KBQtmp, &vc);
			ScalarType rho=vc[0];
			Teuchos::SerialDenseMatrix<int,ScalarType> *H=new Teuchos::SerialDenseMatrix<int,ScalarType> (1,1); //(*)
			H->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, ScalarOne, *lk, *rk, ScalarZero); //rho+=lk*rk
			rho-=(*H)(0,0);
			delete(H); //(#)
			delete(lk); //(#)
			delete(rk); //(#)
			
			LkTrans(_knownEV-1,_knownEV-1)=Teuchos::ScalarTraits<ScalarType>::one(); //_Lk:=[Lk 0; lk^T 1]			
			Rk(_knownEV-1,_knownEV-1)=rho; //Rk:=[Rk rk; 0 rho]
		}
	}
	
        for(int i=0; i<SearchSpaceSize ; ++i) {Ztmp[0][i] = Z[conv][i];}

        // compute BQtmp
        sublist[0] = _knownEV;
        BQtmp = BQ->CloneView(sublist);
        BQtmp->MvTimesMatAddMv(ScalarOne, *BUtmp, Ztmp, ScalarZero);

        // compute K^{-1} BQtmp
        KBQtmp = KBQ->CloneView(sublist);
        _Prec->Apply(*BQtmp, *KBQtmp); //Qk=[_Qk _qk]

        // Now compute residual
        // Rtmp = AU*Z(:,1)
        Rtmp->MvTimesMatAddMv (ScalarOne, (*AUtmp), Ztmp, ScalarZero);
        // Rtmp = Rtmp - theta[conv]*AU*Z(:,1) = Rtmp - theta[conv] * BQtmp
        Rtmp->MvAddMv (ScalarOne, *Rtmp, -theta[conv], *BQtmp);
        // nrm  = ||Rtmp||
        Rtmp->MvNorm(&nrm);
        Qtmp = _evecs->CloneView(sublist);
        // Qtmp = U*Z(:,1)
        Qtmp->MvTimesMatAddMv (ScalarOne, (*Utmp), Ztmp, ScalarZero);

        // nrm  = ||Rtmp||
        Qtmp->MvNorm(&nrm_q);
        nrm[0] /= nrm_q[0];

      }

      delete Utmp, AUtmp, BUtmp;
      delete Rtmp, Qtmp;

      if (_knownEV == _nev) {
      	break;
      }

      // ========================================================= //
      // conv of the approximations are saved since they are       //
      // accurate enough Perform a deflation (so far only on the Z //
      // matrix. The U, AU and BU matrices are adapted ONLY in the //
      // end)                                                      //
      // ========================================================= //
      
      if (conv > 0)
      {
        if (_om->doPrint()) 
        {
          _os << "[Converged (" << conv << ")]" << endl;
        }
	
	Ztmp2=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Z[conv], Z.stride(), SearchSpaceSize, SearchSpaceSize - conv);
	
	sublist.resize(SearchSpaceSize);
        for(int i=0; i<SearchSpaceSize; ++i) sublist[i] = i;

        Utmp  = U->CloneView(sublist); //(*)
        AUtmp = AU->CloneView(sublist); //(*)
        BUtmp = BU->CloneView(sublist); //(*)
	

        Utmp->MvTimesMatAddMv(ScalarOne, (*Utmp), *Ztmp2, ScalarZero);
        AUtmp->MvTimesMatAddMv(ScalarOne, (*AUtmp), *Ztmp2, ScalarZero);
        BUtmp->MvTimesMatAddMv(ScalarOne, (*BUtmp), *Ztmp2, ScalarZero);

        PE.Rotate(*Ztmp2);

	Ztmp2->putScalar(ScalarZero);	
	for(int i=0; i<(SearchSpaceSize-conv); ++i){
	  theta[i] = theta[i+conv];
	  Z(i,i) = ScalarOne;
        }

        delete Utmp, AUtmp, BUtmp, Ztmp2; //(#)(#)(#)(#)(#)(#)(#)

        // Finally, compute the new search space size
        SearchSpaceSize = SearchSpaceSize - conv;
      }	  


      
      //TARGET SWITCH 
      for(int i=0; i<_blockSize; ++i){
 	if (nrm[0]<nrm_q[0]*_eswitch) sigma[i] = theta[0];
 	else sigma[i] = _TARGET;
      }


      // ====================================================== //
      // Compute the residuals of the best blockSize eigenpair  //
      // approximations and precondition them, i.e. compute the //
      // new search space directions                            //
      // ====================================================== //

      sublist.resize(SearchSpaceSize);
      for(int i=0; i<SearchSpaceSize; ++i) sublist[i]=i;

      Utmp  = U->CloneView(sublist); //(*)
      AUtmp = AU->CloneView(sublist); //(*)
      BUtmp = BU->CloneView(sublist); //(*)

      sublist.resize(1);
      for (int i=0; i<_blockSize ; ++i){
 	sublist[0] = i;
 	Rtmp = R->CloneView(sublist);
 	Ztmp2 = new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Z[i], Z.stride(), SearchSpaceSize, 1); //(*)	
 	Rtmp->MvTimesMatAddMv (ScalarOne, (*AUtmp), (*Ztmp2), ScalarZero);	
	
 	sublist[0] = _knownEV;
	Qtmp = _evecs->CloneView(sublist);
 	Qtmp->MvTimesMatAddMv (ScalarOne, (*BUtmp), (*Ztmp2), ScalarZero);	
	
	Rtmp->MvAddMv (ScalarOne, *Rtmp, -theta[i], *Qtmp);
#if 0
	// nrm = ||Rtmp||
        Rtmp->MvNorm(&nrm);

 	sublist[0] = _knownEV;
 	Qtmp->MvTimesMatAddMv (ScalarOne, (*Utmp), (*Ztmp2), ScalarZero);	
        // nrm  = ||Qtmp||	
        Qtmp->MvNorm(&nrm_q);

	//TARGET SWITCH
	if (nrm[0]<nrm_q[0]*_eswitch) sigma[i]=theta[i];
	else sigma[i]=_TARGET;
#endif
	
 	delete Rtmp, Qtmp, Ztmp2; //(#)
      }

      sublist.resize(1);
      sublist[0] = _knownEV;

      Qtmp = _evecs->CloneView(sublist); //(*)
      BQtmp = BQ->CloneView(sublist); //(*)
      KBQtmp = KBQ->CloneView(sublist); //(*)
      
      vQ2.resize(_knownEV+1);
      for (int i=0; i<=_knownEV; i++) vQ2[i]=i;
      Qtmp2 = _evecs->CloneView(vQ2); //(*)
      BQtmp2=BQ->CloneView(vQ2); //(*)
      KBQtmp2=KBQ->CloneView(vQ2); //(*)
      
      for(int j=0; j<_blockSize; ++j){

        sublist[0] = j;
        Rtmp = R->CloneView(sublist); //(*)

	Ztmp2 = new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Z[j], Z.stride(), SearchSpaceSize, 1); //(*)	
	Qtmp->MvTimesMatAddMv (ScalarOne, *Utmp, *Ztmp2, ScalarZero); //Qtmp:=U*Z(:,conv+j)
	delete(Ztmp2); //(#)
	
	if (_B.get()) _B->Apply((*Qtmp), (*BQtmp)); //Qb(:,_knownEV)=BQtmp:=B*Q(;,_knownEV)
        else BQtmp->SetBlock(*Qtmp, sublist); //Qb(:,_knownEV)=BQtmp:=Q(;,_knownEV)

	//CORRECTION
      	sublist[0]=_knownEV;
		
	if (!_Prec.get()) {
		Ksys = new JacDavOp1<ScalarType,MV>(*Qtmp2, *BQtmp2); //Ksys:=I-Q*Qb^(T) //(*)
	} else {
		sublist[0]=_knownEV;
		_Prec->Apply(*BQtmp, *KBQtmp); //Qk(:,_knownEV)=KBQtmp:=Prec*Qb(:,_knownEV)
		
		if (_knownEV==0) {
			LkTrans(0,0)=ScalarOne; //Lk(0,0):=1
			BQtmp->MvPseudoDot(*KBQtmp, &vc); //vc[0]:=BQtmp*KBQtmp //FIXME: Transpose in complex case
			Rk(0,0)=vc[0];//Rk(0,0)=vc[0]
		} else {
			lk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, LkTrans[_knownEV], LkTrans.stride(), _knownEV, 1);
			vQ2.resize(_knownEV);
			h=KBQ->CloneView(vQ2); //(*)
			BQtmp->MvTransMv(ScalarOne, *h, *lk);
			delete(h); //(#)
						
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, lk->numRows(), lk->numCols(), ScalarOne, Rk.values(), Rk.stride(), lk->values(), lk->stride());
			
			rk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Rk[_knownEV], Rk.stride(), _knownEV, 1);
			h=BQ->CloneView(vQ2); //(*)
			KBQtmp->MvTransMv(ScalarOne, *h, *rk);
			delete(h); //(#)
		
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, rk->numRows(), rk->numCols(), ScalarOne, LkTrans.values(), LkTrans.stride(), rk->values(), rk->stride());
			
			BQtmp->MvPseudoDot(*KBQtmp, &vc);
			ScalarType rho=vc[0];
			Teuchos::SerialDenseMatrix<int,ScalarType> *H=new Teuchos::SerialDenseMatrix<int,ScalarType> (1,1);
			H->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, ScalarOne, *lk, *rk, ScalarZero); //rho+=lk*rk
			rho-=(*H)(0,0);
			delete(H); //(#)
			delete(rk); //(#)
			delete(lk); //(#)

			LkTrans(_knownEV,_knownEV)=Teuchos::ScalarTraits<ScalarType>::one(); //_Lk:=[Lk 0; lk^T 1]
			Rk(_knownEV,_knownEV)=rho; //Rk:=[Rk rk; 0 rho]
		}
		
		Teuchos::SerialDenseMatrix<int,ScalarType> *LkTransView=new Teuchos::SerialDenseMatrix<int,ScalarType> (Teuchos::View, LkTrans.values(), LkTrans.stride(), _knownEV+1, _knownEV+1);
		
		Teuchos::SerialDenseMatrix<int,ScalarType> *RkView=new Teuchos::SerialDenseMatrix<int,ScalarType> (Teuchos::View, Rk.values(), Rk.stride(), _knownEV+1, _knownEV+1);
	
		Ksys=new JacDavOp4<ScalarType,MV,OP>(*KBQtmp2, *LkTransView, *RkView, *BQtmp2, _Prec); //Ksys:=(I-Qk*(Lk*Rk)^(-1)*Qb^(T))*K^(-1)
	}

	Asys=new JacDavOp2<ScalarType, MV, OP>(*BQtmp2, *Qtmp2, _A, sigma[j], _B); //Asys:=(I-Qb*Q^(T))*(_A-_sigma*_B) //(*)
	
	Teuchos::SerialDenseMatrix<int,ScalarType> *H=new Teuchos::SerialDenseMatrix<int,ScalarType>(); //(*)
	H->shapeUninitialized(_knownEV+1,1);

	Rtmp->MvTransMv(ScalarOne, *Qtmp2, *H);
	MV *r=Rtmp->CloneCopy(); //(*)
	r->MvTimesMatAddMv(ScalarOne, *BQtmp2, *H, -ScalarOne); //r:=-(I-Qb*Q^(T))*_r
	delete H; //(#)
	
	GMRES<ScalarType,MV,Anasazi::Operator<ScalarType> > gmres2;
	
	gmres2.setMaxIter(_LSIterMax);
	//printf("outer=%d, gamma=%d, tol=%f\n",outer,_gamma,MAX(pow(_gamma,-(double)outer),_elin));
	gmres2.setTolerance(MAX(pow(_gamma,-(double)outer),_elin));
	gmres2.setRestart(_LSRestart);
	
	Rtmp->MvInit(ScalarZero);
	
	ok=gmres2.solve(*Asys, *Ksys, *r, *Rtmp); 
	
	delete Rtmp; //(#)
	delete r; //(#)
	delete Asys; //(#)
	delete Ksys; //(#)
	
	if (ok!=Anasazi::Ok) {
		printf("GMRES failed\n");
		if (ok==Anasazi::Unconverged) printf("GMRES: not converged\n");
		else printf("numerical failure\n");
	}
	//assert(ok==Anasazi::Ok); //ignore failure
	
	printf("GMRES (it=%d, restarts=%d)\n",gmres2.getIterations(), gmres2.getRestarts());
      }

      delete Utmp, AUtmp, BUtmp, Qtmp, BQtmp, KBQtmp, Qtmp2, BQtmp2, KBQtmp2; //(#)(#)(#)(#)(#)(#)(#)(#)(#)

      // ========================================================= //
      // Check the actual search space size and perform a restart, //
      // if adding the new directions would lead to an oversized   //
      // subspace                                                  //
      // ========================================================= //
      
      purged = 0;
      if ((SearchSpaceSize + _blockSize)>_SMAX){
        purged = SearchSpaceSize - (_SMIN - _blockSize);

        if (_om->doPrint()) 
        {
          _os << "[Restart (" << SearchSpaceSize  << " -> "
              << SearchSpaceSize - purged << ")]" << endl;
        }

        // increment restart counter
        ++_numRestarts;
	
	// Reshape the Z matrix according to the restart behaviour	
	Ztmp2=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Z[0], Z.stride(), SearchSpaceSize, SearchSpaceSize-purged);
	
	sublist.resize(SearchSpaceSize);
        for(int i=0; i<SearchSpaceSize; ++i) sublist[i] = i;

        Utmp = U->CloneView(sublist); //(*)
        AUtmp = AU->CloneView(sublist); //(*)
        BUtmp = BU->CloneView(sublist); //(*)
	
	sublist.resize(SearchSpaceSize-purged);
	
        Utmp->MvTimesMatAddMv(ScalarOne, (*Utmp), *Ztmp2, ScalarZero);
        AUtmp->MvTimesMatAddMv(ScalarOne, (*AUtmp), *Ztmp2, ScalarZero);
        BUtmp->MvTimesMatAddMv(ScalarOne, (*BUtmp), *Ztmp2, ScalarZero);

        PE.Rotate(*Ztmp2);

        delete Utmp, AUtmp, BUtmp, Ztmp2; //(#)(#)(#)(#)(#)(#)(#)

        // Finally, compute the new search space size
        SearchSpaceSize = SearchSpaceSize - purged;
      }	  
    

      // Update the iteration step counter
      _iters++;
      
      outer++;
    }

    //return(Ok);
    if (_knownEV==_nev) return Ok;
    else return Unconverged;
    
  } // solve()

  // ==========================================================================
  template <class ScalarType, class MagnitudeType, class MV, class OP>
  void BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>::
  currentStatus()
  {
    if (_om->doPrint()) 
    {
      _os.setf(ios::scientific, ios::floatfield);
      _os.precision(6);
      _os <<endl;
      _os << "******************* CURRENT STATUS *******************" << endl;
      _os << "Anasazi Block Jacobi-Davidson solver" << endl;
      _os << "The number of iterations performed is " <<_iters << endl;
      _os << "The number of restarts performed is " << _numRestarts << endl;
      _os << "The block size is " << _blockSize << endl;
      _os << "The number of eigenvalues requested is " << _nev << endl;
      _os << "The number of eigenvalues computed is "<< _knownEV << endl;
      _os << "The target is " << _TARGET << endl;
      _os << "The requested residual tolerance is "<< _residual_tolerance << endl;
      _os << endl;
      _os << "COMPUTED EIGENVALUES                 "<<endl;
      _os << std::setw(16) << std::right << "Eigenvalue" 
          << std::setw(24) << std::right << "Relative Residual"
          << endl;
      _os << "------------------------------------------------------"<<endl;

      if (_knownEV > 0) 
      {
        for (int i = 0; i < _knownEV; i++) 
        {
          _os << std::setw(16) << std::right << (*_evals)[i]
              << std::setw(24) << std::right << _normR[i]
              << endl;
        }
      } 
      else 
      {
        _os << "[none computed]"<<endl;
      }
      _os << "******************************************************" << endl;  
      _os << endl; 
    }
  }
} // namespace Anasazi

#endif
