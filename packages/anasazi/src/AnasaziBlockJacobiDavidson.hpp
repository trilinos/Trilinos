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

#include "AnasaziOperator.hpp"
template <class ScalarType, class MV>
class JacDavOp1 : public virtual Anasazi::Operator<ScalarType> {
protected:
	MV& _Q;
        MV& _Qb;
public:
	JacDavOp1(MV &Q, MV &Qb) : _Q(Q), _Qb(Qb) {}
	
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
	MV& _Qb;
        MV& _Q;
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
	MV& _Qk;
        MV& _Qb;
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

#include "AnasaziConfigDefs.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziProjectedEigenproblem.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"

namespace Anasazi {

// single-vector GMRES class, several minor fixes still needed

template <class ScalarType, class MV, class OP>
class GMRES {
protected:
	typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType magnitudeType;	
	
	int _LSitmax; //maximum of iterations
	int _LSitRestart; //number of iterations at which restarts occur
	magnitudeType _epslin; //tolerance
	
	Teuchos::BLAS<int,ScalarType> blas; //used to access Teuchos blas functions
	std::string _errorString; //error string (stores error information, if funcition solve returns Anasazi::Failed)
	
	static void UCosSin(const ScalarType x1, const ScalarType x2, magnitudeType &c, ScalarType &s); //calculates cos und sin for Givens Rotations
	#if 0
	void GivensU(const int i, const int j, const magnitudeType c, const ScalarType &s, Teuchos::SerialDenseMatrix<int, ScalarType> &m); //applies a Givens rotation to a Teuchos-Matrix
	#endif
	static void GivensU(const magnitudeType c, const ScalarType &s, ScalarType &x1, ScalarType &x2); //applies a Givens rotation to 2 specific numbers
	
	static double conjugate(const double t) {return t;} //calculates the complex conjugate
	static std::complex<double> conjugate(const std::complex<double> &t) {return std::complex<double>(t.real(),-t.imag());} //calculates the complex conjugate
	
public:
	GMRES() : _LSitmax(500), _epslin(1e-6), _LSitRestart(20)  {} //contructor
	GMRES(const int itmax, const double tol, const int restart) : _LSitmax(itmax), _epslin(tol), _LSitRestart(restart) {
		if (_LSitmax<0) _LSitmax=500;
		if (_LSitRestart<1) _LSitRestart=20;	
	}
		
	
	void setMaxIter (const int i) {if (i>=0) _LSitmax=i;} //sets the maximum of iterations
	void setTolerance(const typename Teuchos::ScalarTraits<ScalarType>::magnitudeType e) {_epslin=e;} //sets the tolerance
	void setRestart(const int i) {if (i>=1) _LSitRestart=i;} //sets the number of iterations at which restart occurs
	std::string getErrString() {return _errorString;} //returns a failure string, which is only set, if the function solve returns Anasazi::Failed (undefined content in other cases)
	
	virtual Anasazi::ReturnType solve(OP &A, OP &Kinv, MV &b, MV &sol); //solves a specific equation A*sol=b with preconditioner Kinv
};

//----------------------------------------------------------------------------------------------
//                                 *** IMPLEMENTATION ***
//----------------------------------------------------------------------------------------------

//solves a specific equation A*sol=b with preconditioner Kinv
template <class ScalarType, class MV, class OP>
Anasazi::ReturnType GMRES<ScalarType, MV, OP>::solve(OP &A, OP &Kinv, MV &b, MV &sol) {
	int k=0, i, LSit;
	magnitudeType rho;
	
	std::vector<magnitudeType> vn(1);
	std::vector<ScalarType> vs(1);
	std::vector<int> vi(1), vi2;
	
	MV *r, *ua, *uK;
	std::vector<MV*> u(_LSitRestart+1);
	Teuchos::SerialDenseMatrix<int,ScalarType> e, *z, R, s;
	Teuchos::SerialDenseMatrix<int,magnitudeType> c;
	Anasazi::ReturnType ok;
	
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
	r->MvAddMv((ScalarType)1,b,(ScalarType)-1,*r); //r=b-A*x
	r->MvNorm(&vn);

	ok=Kinv.Apply(*r,*r); //r=K^(-1)*(b-A*x)
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator Kinv to r";
		delete(r);
		r=0;
		return Anasazi::Failed;
	}
	r->MvNorm(&vn);
	rho=vn[0];

	e.shape(_LSitRestart+1,1); //alle Werte mit 0 initialisiert
	e(0,0)=rho;
	
	for (i=0; i<=_LSitRestart; i++) u[i]=r->Clone(1);
	if (rho!=(ScalarType)0) u[0]->MvAddMv((ScalarType)1/rho,*r,(ScalarType)0,*r); //u[0]=r/rho
	else u[0]->MvAddMv((ScalarType)1,*r,(ScalarType)0,*(u[0]));
	
	ua=r->Clone(1); //(*) 1 col
	uK=r->Clone(1); //(*) 1 col
	R.shape(_LSitRestart+1,_LSitRestart+1); //(m) (_LSitRestart+1 x _LSitRestart+1)
	c.shapeUninitialized(_LSitRestart,1); //(m) (_LSitRestart x 1)
	s.shapeUninitialized(_LSitRestart,1); //(m) (_LSitRestart x 1)

	k=0;
	for (LSit=1; LSit<=_LSitmax; LSit++) {
		ok=A.Apply(*(u[k]),*ua); //ua=A*u[k]
		assert(ok==Anasazi::Ok);
		ok=Kinv.Apply(*ua,*uK); //uK=A*u[k]
		assert(ok==Anasazi::Ok);
		
		for (i=0; i<=k; i++) {
			u[i]->MvDot(*uK,&vs);
			R(i,k)=vs[0]; //R[i,k]=u[i]^(H)*u[k]
			uK->MvAddMv((ScalarType)1,*uK,-R(i,k),*(u[i])); //uK=uK-R[i,k]*u[i]
		}
	
		uK->MvNorm(&vn);
		R(k+1,k)=vn[0]; //R[k+1,k]=norm(u[k+1]);
		if (vn[0]!=0) u[k+1]->MvAddMv((ScalarType)1/R(k+1,k),*uK,(ScalarType)0,*uK); //u[k+1]=uK/R[k+1,k]
		else u[k+1]->MvAddMv((ScalarType)1,*uK,(ScalarType)0,*uK);
		
		for (i=0; i<k; i++) {
			GivensU(c(i,0),s(i,0),R(i,k),R(i+1,k));
		}
		UCosSin(R(k,k),R(k+1,k),c(k,0),s(k,0));
		GivensU(c(k,0),s(k,0), R(k,k), R(k+1,k));
		
		//R(k+1,k)=(ScalarType)0; //Da ich weiss, dass dies 0 werden soll //unnoetig, wenn ich bei linSolve Dreiecksmatrix verwende
		
		GivensU(c(k,0),s(k,0), e(k,0),e(k+1,0));
			
                if (LSit % 10 == 1)
                  cout << "iteration = " << LSit << ", res = " << Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1, 0)) << endl;

		if (Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1,0))<=rho*_epslin) break;
		
		if (k+1==_LSitRestart) {
			//printf("GMRES: RESTART\n");
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, k+1, 1, (ScalarType)1, R.values(), R.stride(), e.values(), e.stride());
		
			for (i=0; i<=k; i++) sol.MvAddMv((ScalarType)1,sol,e(i,0),*(u[i]));  //x:=x+U[:,0:k]*z
			
			ok=A.Apply(sol,*r);
			assert(ok==Anasazi::Ok);
			r->MvAddMv((ScalarType)1,b,(ScalarType)-1,*r); //r=b-A*x
			ok=Kinv.Apply(*r,*r); //r=Kinv*r
			assert(ok==Anasazi::Ok);
			r->MvNorm(&vn); //eta=vn[0]=norm(r);
			
			assert(vn[0]!=(magnitudeType)0);
			e.shape(_LSitRestart+1,1); //alle Werte mit 0 initialisiert
			e(0,0)=vn[0];
			
			
			u[0]->MvAddMv((ScalarType)1/vn[0],*r,(ScalarType)0,*(u[0])); //u[0]=r/eta
			
			k=-1;
		}
		
		k++;
	}
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

#if 0
//applies a Givens rotation to a Teuchos-Matrix
template <class ScalarType, class MV, class OP>
void GMRES<ScalarType, MV, OP>::GivensU(const int i, const int j, const magnitudeType c, const ScalarType &s, Teuchos::SerialDenseMatrix<int, ScalarType> &m) {
	Teuchos::SerialDenseMatrix<int, ScalarType> *tmp;
	
	tmp=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::Copy,m.values()+i,m.stride(),1,m.numCols());
	
	this->blas.SCAL(m.numCols(), (ScalarType) c, m.values()+i, m.stride());
	this->blas.AXPY(m.numCols(), conjugate(s), m.values()+j, m.stride(), m.values()+i, m.stride()); //m[:,i]=c*m[i]+conj(s)*m[:,j]
	
	this->blas.SCAL(m.numCols(), c, m.values()+j, m.stride());
	this->blas.AXPY(m.numCols(), -s, tmp->values(), tmp->stride(), m.values()+j, m.stride()); //m[:,j]=-s*zi+c*z[:,j]
	delete(tmp);
	tmp=0;
}
#endif

//applies a Givens rotation to 2 specific numbers
template <class ScalarType, class MV, class OP>
void GMRES<ScalarType, MV, OP>::GivensU(const magnitudeType c, const ScalarType &s, ScalarType &x1, ScalarType &x2) {
	ScalarType h=x1;
	
	x1=c*x1+conjugate(s)*x2;
	x2=c*x2-s*h;
}


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
{\tt  4} &  $\quad\quad\quad\quad$\mathcal{U} := \mathcal{U} Z(1:s_{max}, 1:s_{min -1}$ \\
{\tt  5} &  $\quad\quad$End \\
{\tt  6} &  $\quad\quad$Return (\mathcal{U}, \Theta, Z)$ \\
{\tt  7} &  EndProcedure  \\
\end{tabular}
\f]






\author Oscar Chinellato (ETHZ/ICOS), Sabine Emch (ETHZ), Marzio Sala (ETHZ/COLAB)

\date Last updated on 01-Nov-05
*/

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
    int _LSIterMax;
    double _elin, _gamma;
    
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
    _maxIters(_pl.get("Max Iters", 300)),
    _blockSize(_pl.get("Block Size", 1)),
    _residual_tolerance(_pl.get("Tol", Teuchos::ScalarTraits<MagnitudeType>::eps())),
    _numRestarts(0),
    _iters(0),
    _TARGET(_pl.get("Target", Teuchos::ScalarTraits<ScalarType>::zero())),
    _SMIN(_pl.get("SMIN", 20)), 
    _SMAX(_pl.get("SMAX", 30)),
    _knownEV(0),
    _haveKinv(false), //FIXME, don't needed
    _LSIterMax(_pl.get("KrylovSolver: MaxIters", 1550)),
    _elin(_pl.get("KrylovSolver: Tolerance", 10e5 * Teuchos::ScalarTraits<MagnitudeType>::eps())), 
    _gamma(_pl.get("gamma", 0.5)) // FIXME
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
  
    if (iVec->GetNumberVecs() < 1) 
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

    // Check that the maximum number of iterations is a positive number
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
    int SearchSpaceDim;
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
    std::vector<int> perm(_SMAX), iperm(_SMAX), vQ(2), vQ2;
    Teuchos::SerialDenseMatrix<int, ScalarType> Z(_SMAX, _SMAX), LkTrans(_nev,_nev), Rk(_nev,_nev);
    Teuchos::SerialDenseMatrix<int, ScalarType> Ztmp, *lk, *rk;
    Anasazi::Operator<ScalarType> *Ksys, *Asys;

    // Allocate space for the residuals ...
    Teuchos::RefCountPtr<MV> R;
    R = MVT::CloneCopy(*iVec);
    ///R = MVT::Clone(*iVec, _blockSize);
    ///R->MvRandom();

    // Allocate space for the result ...
    int conv;   // number of converged eigenpairs per iteration
    int purged; // number of deflated eigenpairs per iteration
    ScalarType sigma;
    MV *Qtmp, *Qtmp2;
    MV *BQtmp, *BQtmp2;
    MV *KBQtmp, *KBQtmp2, *h;

    Teuchos::RefCountPtr<MV> Q, BQ, KBQ;
    Teuchos::BLAS<int,ScalarType> blas;
    Anasazi::ReturnType ok;

    // FIXME: why nev + 1???
    Q = MVT::Clone(*iVec, _nev+1);
    BQ = MVT::Clone(*iVec, _nev); 
    KBQ = MVT::Clone(*iVec, _nev+1);

    // Instantiate the ProjectedEigenproblem component and set some of the values ...
    ProjectedEigenproblem<int, ScalarType, MV> PE("Symmetric", _SMAX);
    PE.SetTheta(&theta);
    PE.SetZ(&Z);

    _knownEV = 0;
    SearchSpaceDim = 0;
    _iters = 0;

    if (_om->doPrint()) 
    {
      _os << "[Starting Solver]" << endl;
    }

    // ========================================= //
    // Start the (Block-) Davidson iteration ... //
    // ========================================= //
    
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
        Rtmp = R->CloneView(sublist);

        // Gram-Schmidt reduce the vector ...
        for(int j=0; j < SearchSpaceDim; j++){
          sublist[0] = j;
          Utmp  = U->CloneView(sublist);
          BUtmp = BU->CloneView(sublist);
          Rtmp->MvDot((*BUtmp), &eta);
          Rtmp->MvAddMv(ScalarOne, (*Rtmp), -eta[0], (*Utmp));
          delete Utmp, BUtmp;
        }

        for(int j=0; j < _knownEV; j++){
          sublist[0] = j;
          // FIXME: _evecs or Q???
          Qtmp = _evecs->CloneView(sublist);
          BQtmp = BQ->CloneView(sublist);
          Rtmp->MvDot((*BQtmp), &eta);
          Rtmp->MvAddMv(ScalarOne, (*Rtmp), -eta[0], (*Qtmp));	      
          delete Qtmp, BQtmp;
        }

        // Now B-normalise the vector ...
        sublist[0] = _knownEV;
        // FIXME: evecs or Q???
        Qtmp = _evecs->CloneView(sublist);
        if (_B.get() != 0)
          _B->Apply((*Rtmp), (*Qtmp));
        else
          Qtmp->SetBlock(*Rtmp, sublist);

        // FIXME: Oscar, va bene? ho cambiato nrm con nrm_q
        std::vector<ScalarType> nrm_2(1);

        Rtmp->MvDot((*Qtmp), &nrm_2);
        nrm_2[0] = Teuchos::ScalarTraits<ScalarType>::squareroot(nrm_2[0]);
        nrm[0] = Teuchos::ScalarTraits<ScalarType>::magnitude(nrm_2[0]);
        Rtmp->MvAddMv(ScalarOne / nrm_2[0], (*Rtmp), ScalarZero, (*Rtmp)); 

        sublist[0] = SearchSpaceDim;
        U->SetBlock((*Rtmp),sublist);

        BUtmp = BU->CloneView(sublist);
        if (_B.get() != 0)
          _B->Apply((*Rtmp), (*BUtmp));
        else
          BUtmp->SetBlock((*Rtmp), sublist);

        delete Rtmp, Qtmp;
        SearchSpaceDim++;
      }

      // Update AU and BU ...
      sublist.resize(_blockSize);
      for(int i=0; i<_blockSize; ++i) sublist[i]=(SearchSpaceDim-_blockSize) + i;

      Utmp = U->CloneView(sublist);

      AUtmp = AU->CloneView(sublist);
      _A->Apply((*Utmp), (*AUtmp));

      BUtmp = BU->CloneView(sublist);

      // Update the ProjectedEigenproblem component by telling it
      // about the space increase
      PE.Add((*Utmp), (*AUtmp), (*BUtmp));
      delete Utmp;
      delete AUtmp;
      delete BUtmp;

      // ====================================================== //
      // Extract eigenpair approximations from the space, then  //
      // Sort eigenpairs with respect to TARGET and compute the //
      // inverse permutation.                                   //
      // CAUTION: At return of the sorting manager, the thetas  //
      // are sorted!!                                           //
      // ====================================================== //
      
      PE.Extract();

      for(int i=0; i<SearchSpaceDim; ++i) {theta[i] -= _TARGET;}
      _sm->sort((Eigensolver<ScalarType,MV,OP>*)NULL, SearchSpaceDim, &theta[0], &perm);
      for(int i=0; i<SearchSpaceDim; ++i) {theta[i] += _TARGET;}
      for(int i=0; i<SearchSpaceDim; ++i) { iperm[perm[i]]=i;}

      // Permute the Z entries according to perm
      Ztmp.shape(SearchSpaceDim,1);

      for(int j=0; j<SearchSpaceDim; ++j){
        if (perm[j] != j){
          for(int i=0; i<SearchSpaceDim; ++i) {Ztmp[0][i] = Z[j][i];}
          for(int i=0; i<SearchSpaceDim; ++i) {Z[j][i] = Z[perm[j]][i];}
          for(int i=0; i<SearchSpaceDim; ++i) {Z[perm[j]][i] = Ztmp[0][i];}

          perm[iperm[j]] = perm[j];	  
        }  
      }

      // Check for convergence of the approximative eigenpairs. If
      // some of them are accurate enough, add them to the space _evecs
      Ztmp.shape(SearchSpaceDim,1);

      sublist.resize(SearchSpaceDim);
      for(int i=0; i<SearchSpaceDim; ++i) sublist[i]=i;	  

      Utmp = U->CloneView(sublist);
      AUtmp = AU->CloneView(sublist);
      BUtmp = BU->CloneView(sublist);

      sublist.resize(1);
      sublist[0] = 0;
      Rtmp = R->CloneView(sublist);

      for(int i=0; i<SearchSpaceDim; ++i) {Ztmp[0][i] = Z[0][i];}
      // Rtmp = AU*Z(:,1)
      Rtmp->MvTimesMatAddMv (ScalarOne, (*AUtmp), Ztmp, ScalarZero); 
      // Rtmp = Rtmp - theta[0]*AU*Z(:,1)
      Rtmp->MvTimesMatAddMv (-theta[0], (*BUtmp), Ztmp, ScalarOne);
      // nrm  = ||Rtmp|| 
      Rtmp->MvNorm(&nrm);
      // 
      sublist[0] = _knownEV;
      Qtmp = _evecs->CloneView(sublist);
      // Qtmp = U*Z(:,1)
      Qtmp->MvTimesMatAddMv (ScalarOne, (*Utmp), Ztmp, ScalarOne);
      // nrm  = ||Rtmp|| 
      Qtmp->MvNorm(&nrm_q);

      nrm[0] /= nrm_q[0];

#if 1
      if (_om->doPrint()) 
      {
        _os << "It: " << _iters << ", knownEV: " << _knownEV;
        _os << ", SearchSpaceDim: " << SearchSpaceDim << ", rel_normr = " << nrm[0];
        _os << ", theta_i = ";
        _os << "(" << theta[0] << ")";
        int j = ANASAZI_MIN(SearchSpaceDim,5);
        for(int i=1; i<j; ++i){
          _os << ", (" << theta[i] << ")";
        }
        _os << endl;
      }
#endif

      conv = 0;
      while(nrm[0] < _residual_tolerance)
      {
        sublist.resize(1);
        sublist[0] = _knownEV;
        // FIXME: evecs or Q??
        Qtmp = _evecs->CloneView(sublist);
        BQtmp = BQ->CloneView(sublist);
        KBQtmp = KBQ->CloneView(sublist);

        // Qtmp = U*Z(:,1)
        Qtmp->MvTimesMatAddMv (ScalarOne, (*Utmp), Ztmp, ScalarZero);
        // BQtmp = B*Qtmp
        if (_B.get() != 0)
          _B->Apply((*Qtmp), (*BQtmp));
        else
          BQtmp->SetBlock(*Qtmp, sublist);

        (*_evals)[_knownEV] = theta[conv];
        _normR[_knownEV] = nrm[0];
        _knownEV++;

        if (_knownEV == _nev) break;
        conv++;

        if (_haveKinv) { //Kinv!=I //Annahme _Kinv!=I, wenn gesetzt
          sublist[0]=_knownEV-1;
          KBQ->SetBlock(*KBQtmp, sublist); //_Qk=[_Qk _qk]

          if (_knownEV==1) {
            LkTrans(0,0)=(ScalarType)1;
            BQtmp->MvDot(*KBQtmp, &vc);
            Rk(0,0)=vc[0];
          } else {
            lk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, LkTrans[_knownEV-1], LkTrans.stride(), _knownEV-1, 1); //(*)
            vQ2.resize(_knownEV-1);
            for (int i=0; i<_knownEV-1; i++) vQ2[i]=i;
            h=KBQ->CloneView(vQ2); //(*)
            BQtmp->MvTransMv((ScalarType)1, *h, *lk);
            delete(h); //(#)

            blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, lk->numRows(), lk->numCols(), (ScalarType)1, Rk.values(), Rk.stride(), lk->values(), lk->stride());

            rk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Rk[_knownEV-1], Rk.stride(), _knownEV-1, 1); //(*)
            h=KBQ->CloneView(vQ2); //(*)
            KBQtmp->MvTransMv((ScalarType)1, *h, *rk);
            delete(h); //(#)

            blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, rk->numRows(), rk->numCols(), (ScalarType)1, LkTrans.values(), LkTrans.stride(), rk->values(), rk->stride());

            BQtmp->MvDot(*KBQtmp, &vc);
            ScalarType rho=vc[0];
            Teuchos::SerialDenseMatrix<int,ScalarType> *H=new Teuchos::SerialDenseMatrix<int,ScalarType> (1,1);
            H->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, (ScalarType)1, *lk, *rk, (ScalarType)0); //rho+=lk*rk
            rho-=(*H)(0,0);
            delete(H);
            delete(lk);
            delete(rk);

            LkTrans(_knownEV-1,_knownEV-1)=Teuchos::ScalarTraits<ScalarType>::one(); //_Lk:=[Lk 0; lk^T 1]

            Rk(_knownEV-1,_knownEV-1)=rho; //Rk:=[Rk rk; 0 rho]
          }
        }
        // Ztmp = Z(:,conv)
        for(int i=0; i<SearchSpaceDim; ++i) {Ztmp[0][i] = Z[conv][i];}
        // Rtmp = AU*Z(:,conv)
        Rtmp->MvTimesMatAddMv (ScalarOne, (*AUtmp), Ztmp, ScalarZero);
        // Rtmp = Rtmp - theta[conv]*AU*Z(:,conv)
        Rtmp->MvTimesMatAddMv (-theta[conv], (*BUtmp), Ztmp, 1.0);
        // nrm  = ||Rtmp||
        Rtmp->MvNorm(&nrm);

        sublist[0] = _knownEV;
        Qtmp = _evecs->CloneView(sublist);
        // Qtmp = U*Z(:,1)
        Qtmp->MvTimesMatAddMv (ScalarOne, (*Utmp), Ztmp, ScalarOne);
        // nrm  = ||Rtmp|| 
        Qtmp->MvNorm(&nrm_q);

        nrm[0] /= nrm_q[0];

        delete Qtmp;
      }

      delete Utmp, AUtmp, BUtmp;
      if (_knownEV == _nev) break;

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
      }

      // ====================================================== //
      // Compute the residuals of the best blockSize eigenpair  //
      // approximations and precondition them, i.e. compute the //
      // new search space directions                            //
      // ====================================================== //
      
      // - projector
      // - implicit matrix A - \sigma B
      // - Krylov solver
      // - Preconditioner (???)
      
      Ztmp.shape(SearchSpaceDim,1);

      sublist.resize(SearchSpaceDim);
      for(int i=0; i<SearchSpaceDim; ++i) sublist[i]=i;	  

      AUtmp = AU->CloneView(sublist);
      BUtmp = BU->CloneView(sublist);

      sublist.resize(1);
      sublist[0]=_knownEV;	  
      Qtmp = _evecs->CloneView(sublist);
      Qtmp = Q->CloneView(sublist);
      BQtmp = BQ->CloneView(sublist);
      KBQtmp = KBQ->CloneView(sublist);

      for(int j=0; j<_blockSize; ++j){
        sublist[0] = j;
        Rtmp = R->CloneView(sublist);

        for(int i=0; i<SearchSpaceDim; ++i) {Ztmp[0][i] = Z[conv+j][i];}
        // Qtmp = AU*Z(:,1)
        Qtmp->MvTimesMatAddMv (ScalarOne, (*AUtmp), Ztmp, ScalarZero);
        // Qtmp = Qtmp - theta[0]*AU*Z(:,1)
        Qtmp->MvTimesMatAddMv (-theta[conv+j], (*BUtmp), Ztmp, ScalarOne);

        if (_Prec.get() != 0)
          _Prec->Apply((*Qtmp), (*Rtmp));
        else
          Rtmp->SetBlock(*Qtmp, sublist);

	  //CORRECTION
      	sublist[0]=_knownEV;
	vQ2.resize(_knownEV+1);
	for (int i=0; i<=_knownEV; i++) vQ2[i]=i;
	//Qtmp2=_evecs->CloneView(vQ2);
	Qtmp2 = Q->CloneView(sublist);
	BQtmp2=BQ->CloneView(vQ2);
	
        // FIXME: Ksys is never deleted!
	if (!_haveKinv) 
        {
		Ksys=new JacDavOp1<ScalarType,MV>(*Qtmp2, *BQtmp2); //Ksys:=I-Q*Qb^(T)
        } else 
        {
          KBQtmp2=KBQ->CloneView(vQ2); //(*)

          if (_knownEV==0) {
            LkTrans(0,0)=(ScalarType)1;
            BQtmp->MvDot(*KBQtmp, &vc);
            Rk(0,0)=vc[0];
          } else {
            lk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, LkTrans[_knownEV], LkTrans.stride(), _knownEV, 1);
            vQ2.resize(_knownEV);
            h=KBQ->CloneView(vQ2); //(*)
            BQtmp->MvTransMv((ScalarType)1, *h, *lk);
            delete(h); //(#)

            blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, lk->numRows(), lk->numCols(), (ScalarType)1, Rk.values(), Rk.stride(), lk->values(), lk->stride());

            rk=new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, Rk[_knownEV], Rk.stride(), _knownEV, 1); //.shapeUninitialized(n-1,1);
            h=BQ->CloneView(vQ2); //(*)
            KBQtmp->MvTransMv((ScalarType)1, *h, *rk);
            delete(h); //(#)

            blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, rk->numRows(), rk->numCols(), (ScalarType)1, LkTrans.values(), LkTrans.stride(), rk->values(), rk->stride());

            BQtmp->MvDot(*KBQtmp, &vc);
            ScalarType rho=vc[0];
            Teuchos::SerialDenseMatrix<int,ScalarType> *H=new Teuchos::SerialDenseMatrix<int,ScalarType> (1,1);
            H->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, (ScalarType)1, *lk, *rk, (ScalarType)0); //rho+=lk*rk
            rho-=(*H)(0,0);
            delete(H);

            //FIXME: kann ich eventuell wegrationalisieren, wenn ich Matrix am Anfang auf 0 setze		
            LkTrans(_knownEV,_knownEV)=Teuchos::ScalarTraits<ScalarType>::one(); //_Lk:=[Lk 0; lk^T 1]
            Rk(_knownEV,_knownEV)=rho; //Rk:=[Rk rk; 0 rho]
          }

          Teuchos::SerialDenseMatrix<int,ScalarType> *LkTransView=new Teuchos::SerialDenseMatrix<int,ScalarType> (Teuchos::View, LkTrans.values(), LkTrans.stride(), _knownEV+1, _knownEV+1);

          Teuchos::SerialDenseMatrix<int,ScalarType> *RkView=new Teuchos::SerialDenseMatrix<int,ScalarType> (Teuchos::View, Rk.values(), Rk.stride(), _knownEV+1, _knownEV+1);

          //FIXME: entweder LkRkinv berechnen oder Operator der Gleichungssysteme loest		
          Ksys=new JacDavOp4<ScalarType,MV,OP>(*KBQtmp2, *LkTransView, *RkView, *BQtmp2, _Prec); //Ksys:=(I-Qk*(Lk*Rk)^(-1)*Qb^(T))*K^(-1)
        }

	Asys=new JacDavOp2<ScalarType, MV, OP>(*BQtmp2, *Qtmp2, _A, sigma, _B); //Asys:=(I-Qb*Q^(T))*(_A-_sigma*_B)
	
	Teuchos::SerialDenseMatrix<int,ScalarType> *H=new Teuchos::SerialDenseMatrix<int,ScalarType>();
	H->shapeUninitialized(_knownEV + 1,1);

	Rtmp->MvTransMv(ScalarOne, *Qtmp2, *H);
	MV * r=Rtmp->CloneCopy();
	r->MvTimesMatAddMv(ScalarOne, *BQtmp2, *H, (ScalarType)-1); //r:=-(I-Qb*Q^(T))*_r
	
	//c:=KrylovSpaceSolver(Asys, r, c, max(_gamma^(-it), _elin, Ksys, _LSIterMax)
	GMRES<ScalarType,MV,Anasazi::Operator<ScalarType> > gmres;
	
	gmres.setMaxIter(_LSIterMax);
	gmres.setTolerance(_elin); //FIXME: max(gamma^it,_elin)
	gmres.setRestart(50); //FIXME
	
	Rtmp->MvInit((ScalarType)0);

        // FIXME
        Rtmp->MvRandom();

	ok=gmres.solve(*Asys, *Ksys, *r, *Rtmp); 
	assert(ok==Anasazi::Ok); //FIXME: Fehlermeldung oder Anasazi::Failed
      
	
	
        delete Rtmp;
      }

      delete AUtmp, BUtmp, Qtmp;

      // ========================================================= //
      // Check the actual search space size and perform a restart, //
      // if adding the new directions would lead to an oversized   //
      // subspace                                                  //
      // ========================================================= //
      
      purged = 0;
      if ((SearchSpaceDim-conv+_blockSize)>_SMAX){
        // increment restart counter
        ++_numRestarts;

        // 
        // FIXME: Why -1 in Sabine's code??
        purged = (SearchSpaceDim-conv) - (_SMIN-_blockSize);

        if (_om->doPrint()) 
        {
          _os << "[Restart (" << SearchSpaceDim - conv  << " -> "
              << SearchSpaceDim - conv - purged << ")]" << endl;
        }
      } 

      // Reshape the Z matrix according to the convergence
      // and/restart behaviour
      if ((conv + purged)>0)
      {
        sublist.resize(SearchSpaceDim);
        for(int i=0; i<SearchSpaceDim; ++i) sublist[i] = i;

        Utmp = U->CloneView(sublist);
        AUtmp = AU->CloneView(sublist);
        BUtmp = BU->CloneView(sublist);

        sublist.resize(SearchSpaceDim-conv-purged);
        for(int i=0; i<(SearchSpaceDim-conv-purged); ++i) sublist[i] = conv+i;
        Ztmp.shapeUninitialized(SearchSpaceDim,SearchSpaceDim-conv-purged);
        for(int i=0; i<SearchSpaceDim; ++i)
          for(int j=0; j<(SearchSpaceDim-conv-purged); ++j)
            Ztmp(i,j) = Z(i,j+conv);

        Utmp->MvTimesMatAddMv(ScalarOne, (*Utmp), Ztmp, ScalarZero);
        AUtmp->MvTimesMatAddMv(ScalarOne, (*AUtmp), Ztmp, ScalarZero);
        BUtmp->MvTimesMatAddMv(ScalarOne, (*BUtmp), Ztmp, ScalarZero);

        PE.Rotate(Ztmp);

        for(int i=0; i<SearchSpaceDim; ++i){
          for(int j=0; j<(SearchSpaceDim-conv-purged); ++j){
            Z(i,j) = Z(i,j+conv);
          }
        }

        delete Utmp, AUtmp, BUtmp, Ztmp;

        // Finally, compute the new search space size
        SearchSpaceDim = SearchSpaceDim - conv - purged;
      }

      // Update the iteration step counter
      _iters++;
    }

    if (_knownEV==_nev) return Ok;
    else return(Unconverged);
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
