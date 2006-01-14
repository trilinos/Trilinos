#ifndef ANASAZI_UTILS_H
#define ANASAZI_UTILS_H

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

//#define _SAB_CONJ_REVERSED

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
	static int version();
	
	virtual Anasazi::ReturnType solve(OP &A, OP &Kinv, MV &b, MV &sol); //solves a specific equation A*sol=b with preconditioner Kinv
};

//----------------------------------------------------------------------------------------------
//                                 *** IMPLEMENTATION ***
//----------------------------------------------------------------------------------------------

#if 1
template <class ScalarType, class MV, class OP>
int GMRES<ScalarType, MV, OP>::version() {return 1;}

//solves a specific equation A*sol=b with preconditioner Kinv
template <class ScalarType, class MV, class OP>
Anasazi::ReturnType GMRES<ScalarType, MV, OP>::solve(OP &A, OP &Kinv, MV &b, MV &sol) {
	int k=0, i, LSit;
	magnitudeType rho;
	
	std::vector<magnitudeType> vm(1);
	std::vector<ScalarType> vs(1);
	
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

	ua=b.Clone(1); //(*) 1 col
	ok=A.Apply(sol,*ua);
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator A to r";
		delete(ua);
		return Anasazi::Failed;
	}
	
	ua->MvAddMv((ScalarType)1,b,(ScalarType)-1,*ua); //r:=b-A*x
	r=b.Clone(1); //(*) 1 col
	ok=Kinv.Apply(*ua,*r); //r:=K^(-1)*(b-A*x)
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator Kinv to r";
		delete(ua);
		delete(r);
		return Anasazi::Failed;
	}
	r->MvNorm(&vm);
	rho=vm[0];
	if (rho==(magnitudeType)0) {
		//sol is solution
		delete(ua);
		delete(r);
		return Anasazi::Ok;
	}
	
	e.shape(_LSitRestart+1,1); //alle Werte mit 0 initialisiert
	e(0,0)=rho;
	
	//for (i=0; i<=_LSitRestart; i++) u[i]=r->Clone(1);
	for (i=0; i<=_LSitRestart; i++) u[i]=0;
	u[0]=r->Clone(1);
	u[0]->MvAddMv((ScalarType)1/rho,*r,(ScalarType)0,*r); //u[0]=r/rho
	
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
			uK->MvDot(*(u[i]),&vs);
			R(i,k)=vs[0]; //R[i,k]:=u[i]^(H)*u[k]
			uK->MvAddMv((ScalarType)1,*uK,-R(i,k),*(u[i])); //uK:=uK-R[i,k]*u[i]
		}
	
		uK->MvNorm(&vm);
		R(k+1,k)=vm[0]; //R[k+1,k]=norm(u[k+1]);
		if (u[k+1]==0) u[k+1]=r->Clone(1);
		if (vm[0]!=0) u[k+1]->MvAddMv((ScalarType)1/R(k+1,k),*uK,(ScalarType)0,*uK); //u[k+1]:=uK/R[k+1,k]
		else u[k+1]->MvAddMv((ScalarType)1,*uK,(ScalarType)0,*uK);
		
		for (i=0; i<k; i++) {
			GivensU(c(i,0),s(i,0),R(i,k),R(i+1,k));
		}
		UCosSin(R(k,k),R(k+1,k),c(k,0),s(k,0));
		GivensU(c(k,0),s(k,0), R(k,k), R(k+1,k));
		
		//R(k+1,k)=(ScalarType)0; //Da ich weiss, dass dies 0 werden soll //unnoetig, wenn ich bei linSolve Dreiecksmatrix verwende
		
		GivensU(c(k,0),s(k,0),e(k,0),e(k+1,0));
			
		if (Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1,0))<=rho*_epslin) break;
		
		if (k+1==_LSitRestart) {
			_restarts++;

			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, k+1, 1, (ScalarType)1, R.values(), R.stride(), e.values(), e.stride());
		
			for (i=0; i<=k; i++) sol.MvAddMv((ScalarType)1,sol,e(i,0),*(u[i]));  //x:=x+U[:,0:k]*z			
			
			ok=A.Apply(sol,*ua);
			assert(ok==Anasazi::Ok);
			ua->MvAddMv((ScalarType)1,b,(ScalarType)-1,*ua); //r:=b-A*x
			ok=Kinv.Apply(*ua,*r); //r:=Kinv*r
			assert(ok==Anasazi::Ok);
			r->MvNorm(&vm); //eta:=vm[0]:=norm(r);
			
			assert(vm[0]!=(magnitudeType)0);
			e.shape(_LSitRestart+1,1); //alle Werte mit 0 initialisiert
			e(0,0)=vm[0];
			
			assert(vm[0]!=0);
			u[0]->MvAddMv((ScalarType)1/vm[0],*r,(ScalarType)0,*(u[0])); //u[0]:=r/eta
			
			k=-1;
		}
		
		k++;
	}
	_LSit=LSit;
	if (LSit<=_LSitmax) {
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, k+1, 1, (ScalarType)1, R.values(), R.stride(), e.values(), e.stride());
		
		for (i=0; i<=k; i++) sol.MvAddMv((ScalarType)1,sol,e(i,0),*(u[i])); //x:=x+U[:,0:k]*z
	}
	
	delete(r); //(#)
	//for (i=0; i<=_LSitRestart; i++) delete(u[i]);
	for (i=0; i<=_LSitRestart; i++) if (u[i]!=0) delete(u[i]);
	delete(ua); //(#)
	delete(uK); //(#)
	
	if (LSit<=_LSitmax) {
		return Anasazi::Ok;
	} else {
		return Anasazi::Unconverged;
	}
}
#else
template <class ScalarType, class MV, class OP>
int GMRES<ScalarType, MV, OP>::version() {return 2;}

template <class ScalarType, class MV, class OP>
Anasazi::ReturnType GMRES<ScalarType, MV, OP>::solve(OP &A, OP &Kinv, MV &b, MV &x) {

	int k=0, i, j;
	ScalarType work[2*_LSitRestart];
	magnitudeType a, rho, eta;
		
	std::vector<magnitudeType> vm(1);
	std::vector<ScalarType> vs(1);
	std::vector<int> vi(1), vi2;
		
	MV *r, *u, *uview, *uA, *uK, *y;
	Teuchos::SerialDenseMatrix<int,ScalarType> e, *z, R, *eview, *Rview, *Rinv, s;
	Teuchos::SerialDenseMatrix<int,magnitudeType> c;
	Anasazi::ReturnType ok;
	
	_LSit=0;
	_restarts = 0;
	
	//FIXME: eventuell auf mehr als einen Eingabe-Vektor erweitern
	if (b.GetNumberVecs()!=1) {
		_errorString="GMRES: only implemented for paramter b with 1 column";
		return Anasazi::Failed;
	}
	if (x.GetNumberVecs()!=1) {
		_errorString="GMRES: only implemented for parameter sol with 1 column";
		return Anasazi::Failed;
	}
	if (b.GetVecLength()!=x.GetVecLength()) {
		_errorString="GMRES: only implemented for b.GetLength()==sol.GetLength()";
		return Anasazi::Failed;
	}
	
	u = b.Clone(1);                                       // u = b
	ok = A.Apply(x, *u);                                  // u = A*x
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator A to r";
		delete(u);
		return Anasazi::Failed;
	}
	
	u->MvAddMv((ScalarType)1, b, (ScalarType)-1, *u);     // u = b - u = b - A*x
	
	r = u->Clone(1);                                      
	ok = Kinv.Apply(*u, *r);                            // r = Kinv*u
	if (ok!=Anasazi::Ok) {
		_errorString="GMRES: error applying operator Kinv to r";
		delete(u);
		delete(r);
		return Anasazi::Failed;
	}

	r->MvNorm(&vm);                            
	//assert(vm[0] != 0.0); //Quatsch: ->wenn die Norm 0 ist, ist es bereits konvergiert
	rho=vm[0];                                         // rho = ||r||
	delete(u);
	if (rho==(magnitudeType)0) {
		//sol is solution
		delete(r);
		return Anasazi::Ok;
	}
	
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
	c.shapeUninitialized(_LSitRestart,1);                           //  c ~ zeros(_LSitRestart+1, 1)
	s.shapeUninitialized(_LSitRestart,1);                           //  s ~ zeros(_LSitRestart+1, 1)

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
			uK->MvDot((*uview), &vs);
			R(i,k) = vs[0];                                   //R[i+1,k] = u[i]^(H)*uK
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
		
		if (Teuchos::ScalarTraits<ScalarType>::magnitude(e(k+1,0)) <= rho*_epslin) break;
				
		if (k+1 == _LSitRestart) {
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, k+1, 1, (ScalarType)1.0, R.values(), R.stride(), e.values(), e.stride());         // e = R\e
			
			vi.resize(1);
			for(i=0; i<=k; ++i){
				vi[0] = i;
				uview = u->CloneView(vi);
				x.MvAddMv((ScalarType)1.0, x, e(i,0), *uview);                 // x = x + U[i,1:k]*e(i)
				delete(uview);
			}
				
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

	if (_LSit != _LSitmax) {
	
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
	
	if (_LSit<=_LSitmax) {
		return Anasazi::Ok;
	} else {
		return Anasazi::Unconverged;
	}
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
		
		#ifdef _SAB_CONJ_REVERSED
		s = conjugate((x2/x1)*c);
		#else
		s = (x2/x1)*c;
		#endif
	}
}

//applies a Givens rotation to 2 specific numbers
template <class ScalarType, class MV, class OP>
void GMRES<ScalarType, MV, OP>::GivensU(const magnitudeType c, const ScalarType &s, ScalarType &x1, ScalarType &x2) {
	ScalarType h=x1;
	
	#ifdef _SAB_CONJ_REVERSED
	x1=c*x1+s*x2;
	x2=c*x2-conjugate(s)*h;
	#else
	x1=c*x1+conjugate(s)*x2;
	x2=c*x2-s*h;
	#endif
}

template <class ScalarType, class MV>
class OpBase {
protected:
	OpBase() {}
public:
	virtual Anasazi::ReturnType Apply(const MV &mv1, MV &mv2) const =0;
};

template <class ScalarType, class MV> 
class OperatorTraits <ScalarType, MV, OpBase<ScalarType,MV> >
{
public:

static ReturnType Apply (const OpBase<ScalarType,MV> &Op, const MV& x, MV& y )
{ return Op.Apply( x, y ); }
};

template <class ScalarType, class MV>
class OpKsysSimple : public virtual OpBase<ScalarType,MV> {
protected:
	Teuchos::RefCountPtr<MV> _Q, _Qb;
public:
	OpKsysSimple(Teuchos::RefCountPtr<MV> Q, Teuchos::RefCountPtr<MV> Qb) : _Q(Q), _Qb(Qb) {}

  //
  // mv2 = (I - Q*Qb^T)*mv1
  //
  virtual Anasazi::ReturnType Apply(const MV &mv1, MV &mv2) const {
    int a=mv1.GetVecLength();
    Teuchos::SerialDenseMatrix<int,ScalarType> H;
    
    if ((mv1.GetNumberVecs()!=mv1.GetNumberVecs()) || (a!=mv2.GetVecLength()) || (a!=_Qb->GetVecLength()) || (a!=_Q->GetVecLength()) || (_Q->GetNumberVecs()!=_Qb->GetNumberVecs())) return Anasazi::Failed;
    
    H.shapeUninitialized(_Qb->GetNumberVecs(), mv1.GetNumberVecs());
    mv1.MvTransMv((ScalarType)1,*_Qb,H,Anasazi::NO_CONJ);
    mv2.MvAddMv((ScalarType)1,mv1,(ScalarType)0,mv2);
    mv2.MvTimesMatAddMv((ScalarType)-1,*_Q,H,(ScalarType)1);
    
    return Anasazi::Ok;
  }
};

template <class ScalarType, class MV, class OP>
class OpAsys : public virtual OpBase<ScalarType,MV> {
protected:
	Teuchos::RefCountPtr<MV> _Qb, _Q;
	Teuchos::RefCountPtr<OP> _A, _B;
	ScalarType _sigma;
public:
	OpAsys(Teuchos::RefCountPtr<MV> &Qb, Teuchos::RefCountPtr<MV> &Q, Teuchos::RefCountPtr<OP> &A, ScalarType sigma, Teuchos::RefCountPtr<OP> &B) : _Q(Q), _Qb(Qb) {
		_A=A;
		_B=B;
		_sigma=sigma;
	}
	
	virtual Anasazi::ReturnType Apply(const MV &mv1, MV &mv2) const {
		Anasazi::ReturnType info;
		Teuchos::SerialDenseMatrix<int,ScalarType> H;
		MV *h1;
		
		//mv2:=(A-sigma*B)*mv1
		info=_A->Apply(mv1,mv2);
		if (info!=Anasazi::Ok) return info;
		
		h1=mv2.Clone(mv2.GetNumberVecs()); //(*) !!!
		info=_B->Apply(mv1,(*h1));
		if (info!=Anasazi::Ok) {
			delete(h1); //(#) !!!
			return info;
		}
		
		mv2.MvAddMv((ScalarType)1,mv2,-_sigma,(*h1));
		delete(h1); //(#) !!!
		
		//H:=Q^(T)*(A-sigma*B)*mv1
		H.shapeUninitialized(_Q->GetNumberVecs(), mv1.GetNumberVecs());
		mv2.MvTransMv((ScalarType)1,*_Q,H, Anasazi::NO_CONJ);
				
		//mv2:=(I-Qb*Q^(T))*(A-sigma*B)*mv1
		mv2.MvTimesMatAddMv((ScalarType)-1,*_Qb,H,(ScalarType)1);
		
		return Anasazi::Ok;
	}
};

template <class ScalarType, class MV, class OP>
class OpKsysSoph : public virtual OpBase<ScalarType,MV> {
protected:
	Teuchos::RefCountPtr<MV> _Qk, _Qb;
	Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > _LkTrans, _Rk;
	Teuchos::RefCountPtr<OP> _Kinv;
	Teuchos::BLAS<int,ScalarType> blas;
public:
	OpKsysSoph(Teuchos::RefCountPtr<MV> &Qk, Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > LkTrans, Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > Rk, Teuchos::RefCountPtr<MV> &Qb, Teuchos::RefCountPtr<OP> &Kinv)  : _Qk(Qk), _Qb(Qb), _LkTrans(LkTrans), _Rk(Rk) {
		_Kinv=Kinv;
	}
	
	virtual Anasazi::ReturnType Apply(const MV &mv1, MV &mv2) const {
		Anasazi::ReturnType info;
		Teuchos::SerialDenseMatrix<int,ScalarType> H;
		
		//mv2:=Kinv*mv1
		info=_Kinv->Apply(mv1,mv2);
		if (info!=Anasazi::Ok) return info;
		
		//H:=Qb^(T)*Kinv*mv1
		H.shapeUninitialized(_Qb->GetNumberVecs(), mv1.GetNumberVecs());
		mv2.MvTransMv((ScalarType)1,*_Qb,H,Anasazi::NO_CONJ);
		
		//H=Lk^(-1)*Qb^(T)*Kinv*mv1 -> solve Lk*x=H
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, H.numRows(), H.numCols(), (ScalarType)1, _LkTrans->values(), _LkTrans->stride(), H.values(), H.stride());
		
		//H=Rk^(-1)*Lk^(-1)*Qb^(T)*Kinv*mv1 -> solve Rk*x=H
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, H.numRows(), H.numCols(), (ScalarType)1, _Rk->values(), _Rk->stride(), H.values(), H.stride());
		
		//mv2:=mv1+Qk*Rk^(-1)*Lk^(-1)*Qb^(T)*Kinv*mv1
		mv2.MvTimesMatAddMv((ScalarType)-1, *_Qk, H, (ScalarType)1);
		
		return Anasazi::Ok;
	}
};

} // namespace Anasazi
#endif
