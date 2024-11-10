// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// Extension of the RAD package (Reverse Automatic Differentiation) --
// a package specialized for function and gradient evaluations -- to
// Hessian-vector products.
// Written in 2007 by David M. Gay at Sandia National Labs, Albuquerque, NM.

#ifndef SACADO_TRAD2_H
#define SACADO_TRAD2_H

#include "Sacado_ConfigDefs.h"
#include "Sacado_trad2_Traits.hpp"

#include <stddef.h>
#include <Sacado_cmath.hpp>

#if defined(RAD_DEBUG_BLOCKKEEP) && !defined(HAVE_SACADO_UNINIT)
#undef RAD_DEBUG_BLOCKKEEP
#endif

#ifdef RAD_Const_WARN	// ==> RAD_AUTO_AD_Const and RAD_DEBUG
#ifndef RAD_AUTO_AD_Const
#define RAD_AUTO_AD_Const
#endif
#ifndef RAD_DEBUG
#define RAD_DEBUG
#endif
extern "C" int RAD_Const_Warn(const void*);// outside any namespace for
					// ease in setting breakpoints
#endif // RAD_Const_WARN

#ifdef RAD_DEBUG
#include <cstdio>
#include <stdlib.h>
#endif

#ifndef RAD_AUTO_AD_Const
#ifdef RAD_DEBUG_BLOCKKEEP
#include <complex>	// must be here when SACADO_NAMESPACE is #defined
#endif
#endif

namespace Sacado {
namespace Rad2 {

// -DRAD_NO_USING_STDCC is needed, e.g., with Sun CC 5.7
#ifndef RAD_NO_USING_STDCC
  // Bring math functions into scope
  using std::exp;
  using std::log;
  using std::log10;
  using std::sqrt;
  using std::cos;
  using std::sin;
  using std::tan;
  using std::acos;
  using std::asin;
  using std::atan;
  using std::cosh;
  using std::sinh;
  using std::tanh;
  using std::abs;
  using std::fabs;
  using std::atan2;
  using std::pow;
#endif //!RAD_NO_USING_STDCC

#ifdef RAD_AUTO_AD_Const
#undef RAD_DEBUG_BLOCKKEEP
#else /*!RAD_AUTO_AD_Const*/
#ifdef RAD_DEBUG_BLOCKKEEP
#if !(RAD_DEBUG_BLOCKKEEP > 0)
#undef RAD_DEBUG_BLOCKKEEP
#else
extern "C" void _uninit_f2c(void *x, int type, long len);

template <typename T>
struct UninitType {};

template <>
struct UninitType<float> {
  static const int utype = 4;
};

template <>
struct UninitType<double> {
  static const int utype = 5;
};

template <typename T>
struct UninitType< std::complex<T> > {
  static const int utype = UninitType<T>::utype + 2;
};

#endif /*RAD_DEBUG_BLOCKKEEP > 0*/
#endif /*RAD_DEBUG_BLOCKKEEP*/
#endif /*RAD_AUTO_AD_Const*/

 class RAD_DoubleIgnore {};

 template<typename T> class
DoubleAvoid {
 public:
	typedef double	dtype;
	typedef T	ttype;
	};
 template<> class
DoubleAvoid<double> {
 public:
	typedef RAD_DoubleIgnore &dtype;
	typedef RAD_DoubleIgnore &ttype;
	};

#define Dtype typename DoubleAvoid<Double>::dtype
#define Ttype typename DoubleAvoid<Double>::ttype

 template<typename Double> class IndepADvar;
 template<typename Double> class ConstADvar;
 template<typename Double> class ConstADvari;
 template<typename Double> class ADvar;
 template<typename Double> class ADvar1;
 template<typename Double> class ADvar1g;
 template<typename Double> class ADvar1s;
 template<typename Double> class ADvar2;
 template<typename Double> class ADvar2g;
 template<typename Double> class ADvar2q;
 template<typename Double> class ADvari;
 template<typename Double> class ADvari_block;
 template<typename Double> class ADvarn;
 template<typename Double> class Derp;

 template<typename Double> struct
ADmemblock {	// We get memory in ADmemblock chunks and never give it back,
		// but reuse it once computations start anew after call(s) on
		// ADcontext::Gradcomp() or ADcontext::Weighted_Gradcomp().
	ADmemblock *next;
	Double memblk[2000];
	};

 template<typename Double> class
ADvari_block {
 public:
	typedef ADvari<Double> ADVari;
	enum { Gulp = 1021 };
	ADvari_block *next, *prev;
	ADVari **limit;
	ADVari *pADvari[Gulp];
	};

 template<typename Double> class
ADcontext {	// A singleton class: one instance in radops.c
	typedef ADmemblock<Double> ADMemblock;
	typedef ADvari <Double> ADVari;
	typedef ADvar1 <Double> ADVar1;
	typedef ADvar2 <Double> ADVar2;
	typedef ADvarn <Double> ADVarn;
	typedef ADvar1g<Double> ADVar1g;
	typedef ADvar1s<Double> ADVar1s;
	typedef ADvar2g<Double> ADVar2g;
	typedef ADvar2q<Double> ADVar2q;
	typedef ADvari_block<Double> ADVari_block;
	ADMemblock *Busy, *Free;
	char *Mbase;
	size_t Mleft;
	ADVari **Ailimit, **Ainext;
	ADVari_block *Aibusy, *Aifree;
	ADMemblock *First;
	ADVari_block *AiFirst;
	double First0[(sizeof(ADMemblock) + sizeof(double) - 1) / sizeof(double)];
	double First1[(sizeof(ADVari_block) + sizeof(double) - 1) / sizeof(double)];
	void *new_ADmemblock(size_t);
	void new_ADvari_block();
	typedef ADvar<Double> ADVar;
	typedef Derp<Double> DErp;
	int rad_need_reinit;
	size_t rad_mleft_save;
#ifdef RAD_DEBUG_BLOCKKEEP
	int rad_busy_blocks;
	ADMemblock *rad_Oldcurmb;
#endif
 public:
	static const Double One, negOne;
	ADcontext();
	void *Memalloc(size_t len);
	static void Gradcomp(int wantgrad);
	static void Hvprod(int, ADVar**, Double*, Double*);
	static void aval_reset(void);
	static void Weighted_Gradcomp(int, ADVar**, Double*);
	static inline void Gradcomp() { Gradcomp(1); }
	inline void ADvari_record(ADVari *x) {
		if (Ainext >= Ailimit)
			new_ADvari_block();
		*Ainext++ = x;
		}
	};

 template<typename Double> class
CADcontext: public ADcontext<Double> {	// for possibly constant ADvar values
 protected:
	bool fpval_implies_const;
 public:
	friend class ADvar<Double>;
	CADcontext(): ADcontext<Double>() { fpval_implies_const = false; }
	};

 template<typename Double> class
Derp {		// one derivative-propagation operation
 public:
	friend class ADvarn<Double>;
	typedef ADvari<Double> ADVari;
	static Derp *LastDerp;
	Derp *next;
	const Double *a;
	const ADVari *b;
	mutable ADVari *c;
	Derp(){};
	Derp(const ADVari *);
	Derp(const Double *, const ADVari *);
	Derp(const Double *, const ADVari *, const ADVari *);
	/* c->aval += a * b->aval; */
	};


// Now we use #define to overcome bad design in the C++ templating system

#define Ai const ADvari<Double>&
#define AI const IndepADvar<Double>&
#define T template<typename Double>
#define D Double
#define T1(f) \
T F f (AI); \
T F f (Ai);
#define T2(r,f) \
 T r f(Ai,Ai); \
 T r f(Ai,D); \
 T r f(Ai,Dtype); \
 T r f(Ai,long); \
 T r f(Ai,int); \
 T r f(D,Ai); \
 T r f(Dtype,Ai); \
 T r f(long,Ai); \
 T r f(int,Ai); \
 T r f(AI,D); \
 T r f(AI,Dtype); \
 T r f(AI,long); \
 T r f(AI,int); \
 T r f(D,AI); \
 T r f(Dtype,AI); \
 T r f(long,AI); \
 T r f(int,AI); \
 T r f(Ai,AI);\
 T r f(AI,Ai);\
 T r f(AI,AI);

#define F ADvari<Double>&
T2(F, operator+)
T2(F, operator-)
T2(F, operator*)
T2(F, operator/)
T2(F, atan2)
T2(F, pow)
T2(F, max)
T2(F, min)
T2(int, operator<)
T2(int, operator<=)
T2(int, operator==)
T2(int, operator!=)
T2(int, operator>=)
T2(int, operator>)
T1(operator+)
T1(operator-)
T1(abs)
T1(acos)
T1(acosh)
T1(asin)
T1(asinh)
T1(atan)
T1(atanh)
T1(cos)
T1(cosh)
T1(exp)
T1(log)
T1(log10)
T1(sin)
T1(sinh)
T1(sqrt)
T1(tan)
T1(tanh)
T1(fabs)

T F copy(AI);
T F copy(Ai);

#undef F
#undef T2
#undef T1
#undef D
#undef T
#undef AI
#undef Ai

} /* namespace Rad2 */
} /* namespace Sacado */
#define SNS Sacado::Rad2
namespace std {	// Moved here from bottom for use in testing nesting of Rad with itself.
  using SNS::exp;
  using SNS::log;
  using SNS::log10;
  using SNS::sqrt;
  using SNS::cos;
  using SNS::sin;
  using SNS::tan;
  using SNS::acos;
  using SNS::asin;
  using SNS::atan;
  using SNS::cosh;
  using SNS::sinh;
  using SNS::tanh;
  using SNS::abs;
  using SNS::fabs;
  using SNS::atan2;
  using SNS::pow;
}
#undef SNS
namespace Sacado {
namespace Rad2 {

 template<typename Double>ADvari<Double>& ADf1(Double f, Double g, const IndepADvar<Double> &x);
 template<typename Double>ADvari<Double>& ADf2(Double f, Double gx, Double gy,
	const IndepADvar<Double> &x, const IndepADvar<Double> &y);
 template<typename Double>ADvari<Double>& ADfn(Double f, int n,
	const IndepADvar<Double> *x, const Double *g);

 template<typename Double> IndepADvar<Double>& ADvar_operatoreq(IndepADvar<Double>*,
	const ADvari<Double>&);
 template<typename Double> ADvar<Double>& ADvar_operatoreq(ADvar<Double>*, const ADvari<Double>&);
 template<typename Double> void AD_Const(const IndepADvar<Double>&);
 template<typename Double> void AD_Const1(Double*, const IndepADvar<Double>&);
 template<typename Double> ADvari<Double>& ADf1(Double, Double, const ADvari<Double>&);
 template<typename Double> ADvari<Double>& ADf2(Double, Double, Double,
	const ADvari<Double>&, const ADvari<Double>&);
 template<typename Double> ADvari<Double>& ADf2(Double, Double, Double,
	const IndepADvar<Double>&, const ADvari<Double>&);
 template<typename Double> ADvari<Double>& ADf2(Double, Double, Double,
	const ADvari<Double>&, const IndepADvar<Double>&);
 template<typename Double> Double val(const ADvari<Double>&);

 enum Advari_Opclass {
	Hv_const,
	Hv_copy,
	Hv_binary,
	Hv_unary,
	Hv_negate,
	Hv_plusLR,
	Hv_minusLR,
	Hv_timesL,
	Hv_timesLR,
	Hv_quotLR,
	Hv_nary
	};

 template<typename Double> ADvari<Double>&
ADf1(Double f, Double g, Double h, const ADvari<Double> &x);

 template<typename Double> ADvari<Double>&
ADf2(Double f, Double gx, Double gy, Double hxx,
			Double hxy, Double hyy, const ADvari<Double> &x, const ADvari<Double> &y);

 template<typename Double> ADvari<Double>&
ADfn(Double f, int n, const IndepADvar<Double> *x, const Double *g, const Double *h);

 template<typename Double> class
ADvari {	// implementation of an ADvar
 public:
	typedef Double value_type;
        typedef typename ScalarType<value_type>::type scalar_type;
	typedef IndepADvar<Double> IndepADVar;
#ifdef RAD_AUTO_AD_Const
	friend class IndepADvar<Double>;
#ifdef RAD_Const_WARN
	mutable const IndepADVar *padv;
#else
 protected:
	mutable const IndepADVar *padv;
#endif //RAD_Const_WARN
 private:
	ADvari *Next;
	static ADvari *First_ADvari, **Last_ADvari;
 public:
	inline void ADvari_padv(const IndepADVar *v) {
		*Last_ADvari = this;
		Last_ADvari = &this->Next;
		this->padv = v;
		}
#endif //RAD_AUTO_AD_Const
#ifdef RAD_DEBUG
	int gcgen;
	int opno;
	static int gcgen_cur, last_opno, zap_gcgen, zap_gcgen1, zap_opno;
	static FILE *debug_file;
#endif
	static ADcontext<Double> adc;
	Advari_Opclass opclass;
	Double Val;		// result of this operation
	mutable Double aval;	// adjoint -- partial of final result w.r.t. this Val
	mutable Double dO;	// deriv of op w.r.t. t in x + t*p
	mutable Double aO;	// adjoint (in Hv computation) of op
	mutable Double adO;	// adjoint (in Hv computation) of dO
	void *operator new(size_t len) {
#ifdef RAD_DEBUG
		ADvari *rv = (ADvari*)ADvari::adc.Memalloc(len);
		rv->gcgen = gcgen_cur;
		rv->opno = ++last_opno;
		if (last_opno == zap_opno && gcgen_cur == zap_gcgen)
			printf("");
		return rv;
#else
		return ADvari::adc.Memalloc(len);
#endif
		}
	void operator delete(void*) {} /*Should never be called.*/
	inline ADvari(Advari_Opclass oc, Double t):
		opclass(oc), Val(t), aval(0.), dO(0.)
		{ if (oc != Hv_const) ADvari::adc.ADvari_record(this); }
	inline ADvari(Advari_Opclass oc, Double t, Double ta):
		opclass(oc), Val(t), aval(ta), dO(0.)
		{ if (oc != Hv_const) ADvari::adc.ADvari_record(this); }
 private:
	inline ADvari(): Val(0.), aval(0.), dO(0.) {}	// prevent construction without value (?)
 public:
	friend class ConstADvari<Double>;
#ifdef RAD_AUTO_AD_Const
	friend class ADcontext<Double>;
	friend class ADvar<Double>;
	friend class ADvar1<Double>;
	friend class ADvar1s<Double>;
	friend class ADvar2<Double>;
	friend class ADvar2q<Double>;
	friend class ConstADvar<Double>;
	ADvari(const IndepADVar *, Double);
	ADvari(const IndepADVar *, Double, Double);
	ADvari(const IndepADVar *, Double, Double, int);
#endif
#define F friend
#define R ADvari&
#define Ai const ADvari&
#define T1(r,f) F r f <>(Ai);
#define T2(r,f) \
F r f <>(Ai,Ai); \
F r f <>(Ttype,Ai); \
F r f <>(Ai,Ttype); \
F r f <>(double,Ai); \
F r f <>(Ai,double); \
F r f <>(long,Ai); \
F r f <>(Ai,long); \
F r f <>(int,Ai); \
F r f <>(Ai,int);
	T1(R,operator+)
	T2(R,operator+)
	T1(R,operator-)
	T2(R,operator-)
	T2(R,operator*)
	T2(R,operator/)
	T1(R,abs)
	T1(R,acos)
	T1(R,acosh)
	T1(R,asin)
	T1(R,asinh)
	T1(R,atan)
	T1(R,atanh)
	T2(R,atan2)
	T2(R,max)
	T2(R,min)
	T1(R,cos)
	T1(R,cosh)
	T1(R,exp)
	T1(R,log)
	T1(R,log10)
	T2(R,pow)
	T1(R,sin)
	T1(R,sinh)
	T1(R,sqrt)
	T1(R,tan)
	T1(R,tanh)
	T1(R,fabs)
	T1(R,copy)
	T2(int,operator<)
	T2(int,operator<=)
	T2(int,operator==)
	T2(int,operator!=)
	T2(int,operator>=)
	T2(int,operator>)
#undef T2
#undef T1
#undef Ai
#undef R
#undef F

	friend ADvari& ADf1<>(Double f, Double g, Double h, const ADvari &x);
	friend ADvari& ADf2<>(Double f, Double gx, Double gy, Double hxx,
			Double hxy, Double hyy, const ADvari &x, const ADvari &y);
	friend ADvari& ADfn<>(Double f, int n, const IndepADVar *x,
			const Double *g, const Double *h);

	inline operator Double() { return this->Val; }
	inline operator Double() const { return this->Val; }
	};

 template<typename Double> class
ADvar1: public ADvari<Double> {	// simplest unary ops
 public:
	typedef ADvari<Double> ADVari;
	Derp<Double> d;
	ADvar1(Advari_Opclass oc, Double val1, const Double *a1, const ADVari *c1):
		ADVari(oc,val1), d(a1,this,c1) {}
#ifdef RAD_AUTO_AD_Const
	typedef typename ADVari::IndepADVar IndepADVar;
	typedef ADvar<Double> ADVar;
	ADvar1(const IndepADVar*, const IndepADVar&);
	ADvar1(const IndepADVar*, const ADVari&);
	ADvar1(Advari_Opclass oc, const Double val1, const Double *a1,
		const ADVari *c1, const ADVar *v):
			ADVari(oc,val1), d(a1,this,c1) {
		c1->padv = 0;
		this->ADvari_padv(v);
		}
#endif
	};


 template<typename Double> class
ConstADvari: public ADvari<Double> {
 private:
	ConstADvari *prevcad;
	ConstADvari() {};	// prevent construction without value (?)
	static ConstADvari *lastcad;
 public:
	typedef ADvari<Double> ADVari;
	typedef Derp<Double> DErp;
	static CADcontext<Double> cadc;
	inline void *operator new(size_t len) { return ConstADvari::cadc.Memalloc(len); }
	inline ConstADvari(Double t): ADVari(Hv_copy, t) { prevcad = lastcad; lastcad = this; }
	static void aval_reset(void);
	};


 template<typename Double> class
IndepADvar {		// an independent ADvar
 private:
	IndepADvar& operator=(IndepADvar&x) {
		/* private to prevent assignment */
#ifdef RAD_AUTO_AD_Const
		if (cv)
			cv->padv = 0;
		return new ADvar1<Double>(this,x);
#elif defined(RAD_EQ_ALIAS)
		this->cv = x.cv;
		return *this;
#else
		return ADvar_operatoreq(this,*x.cv);
#endif //RAD_AUTO_AD_Const
		}
 protected:
	static void AD_Const(const IndepADvar&);
	mutable ADvari<Double> *cv;
 public:
	typedef Double value_type;
	friend class ADvar<Double>;
	friend class ADcontext<Double>;
	friend class ADvar1<Double>;
	friend class ADvarn<Double>;
	typedef ADvari<Double> ADVari;
	typedef ADvar<Double> ADVar;
	IndepADvar(Ttype);
	IndepADvar(double);
	IndepADvar(int);
	IndepADvar(long);
	IndepADvar& operator= (Double);
#ifdef RAD_AUTO_AD_Const
	friend IndepADvar& ADvar_operatoreq<>(IndepADvar*, const ADVari&);
	inline IndepADvar(const IndepADvar &x) { cv = x.cv ? new ADvar1<Double>(this, x) : 0; };
	inline IndepADvar(ADVari *ncv) { cv = ncv; }
	inline IndepADvar() { cv = 0; }
	inline ~IndepADvar() {
					if (cv)
						cv->padv = 0;
					}
#else
	inline IndepADvar() {
#ifndef RAD_EQ_ALIAS
		cv = 0;
#endif
		}
	inline ~IndepADvar() {}
	friend IndepADvar& ADvar_operatoreq<>(IndepADvar*, const ADVari&);
#endif

#ifdef RAD_Const_WARN
	inline operator ADVari&() const {
		ADVari *tcv = this->cv;
		if (tcv->opno < 0)
			RAD_Const_Warn(tcv);
		return *tcv;
		}
	inline operator ADVari*() const {
		ADVari *tcv = this->cv;
		if (tcv->opno < 0)
			RAD_Const_Warn(tcv);
		return tcv;
		}
#else //RAD_Const_WARN
	inline operator ADVari&() const { return *this->cv; }
	inline operator ADVari*() const { return this->cv; }
#endif //RAD_Const_WARN

	Double val() const { return cv->Val; }
	Double adj() const { return cv->aval; }

	friend void AD_Const1<>(Double*, const IndepADvar&);

	friend ADVari& ADf1<>(Double, Double, const IndepADvar&);
	friend ADVari& ADf2<>(Double, Double, Double, const IndepADvar&, const IndepADvar&);
	friend ADVari& ADf2<>(Double, Double, Double, const ADVari&, const IndepADvar&);
	friend ADVari& ADf2<>(Double, Double, Double, const IndepADvar&, const ADVari&);

	static inline void Gradcomp(int wantgrad)
				{ ADcontext<Double>::Gradcomp(wantgrad); }
	static inline void Gradcomp()
				{ ADcontext<Double>::Gradcomp(1); }
	static inline void Hvprod(int n, ADVar **vp, Double *v, Double *hv)
				{ ADcontext<Double>::Hvprod(n, vp, v, hv); }
	static inline void aval_reset() { ConstADvari<Double>::aval_reset(); }
	static inline void Weighted_Gradcomp(int n, ADVar **v, Double *w)
				{ ADcontext<Double>::Weighted_Gradcomp(n, v, w); }

	/* We use #define to deal with bizarre templating rules that apparently */
	/* require us to spell the some conversion explicitly */


#define Ai const ADVari&
#define AI const IndepADvar&
#define D Double
#define T2(r,f) \
 r f <>(AI,AI);\
 r f <>(Ai,AI);\
 r f <>(AI,Ai);\
 r f <>(Ttype,AI);\
 r f <>(double,AI);\
 r f <>(long,AI);\
 r f <>(int,AI);\
 r f <>(AI,Ttype);\
 r f <>(AI,double);\
 r f <>(AI,long);\
 r f <>(AI,int);
#define T1(f) friend ADVari& f<> (AI);

#define F friend ADVari&
T2(F, operator+)
T2(F, operator-)
T2(F, operator*)
T2(F, operator/)
T2(F, atan2)
T2(F, max)
T2(F, min)
T2(F, pow)
#undef F
#define F friend int
T2(F, operator<)
T2(F, operator<=)
T2(F, operator==)
T2(F, operator!=)
T2(F, operator>=)
T2(F, operator>)

T1(operator+)
T1(operator-)
T1(abs)
T1(acos)
T1(acosh)
T1(asin)
T1(asinh)
T1(atan)
T1(atanh)
T1(cos)
T1(cosh)
T1(exp)
T1(log)
T1(log10)
T1(sin)
T1(sinh)
T1(sqrt)
T1(tan)
T1(tanh)
T1(fabs)
T1(copy)

#undef F
#undef T1
#undef T2
#undef D
#undef AI
#undef Ai

	};

 template<typename Double> class
ADvar: public IndepADvar<Double> {	// an "active" variable
 public:
        //! Turn ADvar into a meta-function class usable with mpl::apply
        template <typename U> struct apply { typedef ADvar<U> type; };

	typedef IndepADvar<Double> IndepADVar;
	typedef ADvari<Double> ADVari;
	typedef ConstADvari<Double> ConstADVari;
 private:
	void ADvar_ctr(Double d) {
		ADVari *x;
		if (ConstADVari::cadc.fpval_implies_const)
			x = new ConstADVari(d);
		else {
#ifdef RAD_AUTO_AD_Const
			x = new ADVari((IndepADVar*)this, d);
			x->ADvari_padv(this);
#else
			x = new ADVari(Hv_const, d);
#endif
			}
		this->cv = x;
		}
 public:
	friend class ADvar1<Double>;
	typedef ADvar1<Double> ADVar1;
	ADvar() { /* cv = 0; */ }
	ADvar(Ttype d)  { ADvar_ctr(d); }
	ADvar(double i) { ADvar_ctr(Double(i)); }
	ADvar(int i)	{ ADvar_ctr(Double(i)); }
	ADvar(long i)	{ ADvar_ctr(Double(i)); }
	inline ~ADvar() {}
#ifdef RAD_AUTO_AD_Const
	ADvar(IndepADVar &x) {
		this->cv = x.cv ? new ADVar1(this, x) : 0;
		}
	ADvar(ADVari &x) :IndepADVar(&x) { x.ADvari_padv((IndepADVar*)this);}
	inline ADvar& operator=(IndepADVar &x) {
		if (this->cv)
			this->cv->padv = 0;
		this->cv = new ADVar1(this,x);
		return *this;
		}
	inline ADvar& operator=(ADVari &x) {
		if (this->cv)
			this->cv->padv = 0;
		this->cv = new ADVar1(this, x);
		return *this;
		}
#else /*!RAD_AUTO_AD_Const*/
	friend ADvar& ADvar_operatoreq<>(ADvar*, const ADVari&);
#ifdef RAD_EQ_ALIAS
	/* allow aliasing v and w after "v = w;" */
	inline ADvar(const IndepADVar &x) { this->cv = (ADVari*)x.cv; }
	inline ADvar(const ADVari &x) { this->cv = (ADVari*)&x; }
	inline ADvar& operator=(IndepADVar &x) { this->cv = (ADVari*)x.cv; return *this; }
	inline ADvar& operator=(const ADVari &x) { this->cv = (ADVari*)&x; return *this; }
#else /*!RAD_EQ_ALIAS*/
	friend ADvar& ADvar_operatoreq<>(ADvar*, const ADVari&);
	ADvar(const IndepADVar &x) { this->cv = x.cv ? new ADVar1(x.cv->Val, &this->cv->adc.One, x.cv) : 0; }
	ADvar(const ADvar&x) { this->cv = x.cv ?
		new ADVar1(Hv_copy, x.cv->Val, &this->cv->adc.One, (ADVari*)x.cv) : 0; }
	ADvar(const ADVari &x) { this->cv = new ADVar1(Hv_copy, x.Val, &this->cv->adc.One, &x); }
	inline ADvar& operator=(IndepADVar &x) { return ADvar_operatoreq(this,*x.cv); };
	inline ADvar& operator=(const ADVari &x) { return ADvar_operatoreq(this,x); };
#endif /* RAD_EQ_ALIAS */
#endif /* RAD_AUTO_AD_Const */
	ADvar& operator=(Double);
	ADvar& operator+=(const ADVari&);
	ADvar& operator+=(Double);
	ADvar& operator-=(const ADVari&);
	ADvar& operator-=(Double);
	ADvar& operator*=(const ADVari&);
	ADvar& operator*=(Double);
	ADvar& operator/=(const ADVari&);
	ADvar& operator/=(Double);
	inline static bool get_fpval_implies_const(void)
		{ return ConstADVari::cadc.fpval_implies_const; }
	inline static void set_fpval_implies_const(bool newval)
		{ ConstADVari::cadc.fpval_implies_const = newval; }
	inline static bool setget_fpval_implies_const(bool newval) {
		bool oldval = ConstADVari::cadc.fpval_implies_const;
		ConstADVari::cadc.fpval_implies_const = newval;
		return oldval;
		}
	static inline void Gradcomp(int wantgrad)
				{ ADcontext<Double>::Gradcomp(wantgrad); }
	static inline void Gradcomp()
				{ ADcontext<Double>::Gradcomp(1); }
	static inline void aval_reset() { ConstADVari::aval_reset(); }
	static inline void Weighted_Gradcomp(int n, ADvar **v, Double *w)
				{ ADcontext<Double>::Weighted_Gradcomp(n, v, w); }
	};

template<typename Double>
 inline void AD_Const1(Double *notused, const IndepADvar<Double>&v)
{ IndepADvar<Double>::AD_Const(v); }

template<typename Double>
 inline void AD_Const(const IndepADvar<Double>&v) { AD_Const1((Double*)0, v); }

 template<typename Double> class
ConstADvar: public ADvar<Double> {
 public:
	typedef ADvar<Double> ADVar;
	typedef typename ADVar::ADVari ADVari;
	typedef typename ADVar::ConstADVari ConstADVari;
	typedef Derp<Double> DErp;
	typedef typename ADVar::IndepADVar IndepADVar;
 private: // disable op=
	ConstADvar& operator+=(ADVari&);
	ConstADvar& operator+=(Double);
	ConstADvar& operator-=(ADVari&);
	ConstADvar& operator-=(Double);
	ConstADvar& operator*=(ADVari&);
	ConstADvar& operator*=(Double);
	ConstADvar& operator/=(ADVari&);
	ConstADvar& operator/=(Double);
	void ConstADvar_ctr(Double);
 public:
	ConstADvar(Ttype d)	{ ConstADvar_ctr(d); }
	ConstADvar(double i)	{ ConstADvar_ctr(Double(i)); }
	ConstADvar(int i)	{ ConstADvar_ctr(Double(i)); }
	ConstADvar(long i)	{ ConstADvar_ctr(Double(i)); }
	ConstADvar(const IndepADVar&);
	ConstADvar(const ConstADvar&);
	ConstADvar(const ADVari&);
	inline ~ConstADvar() {}
#ifdef RAD_NO_CONST_UPDATE
 private:
#endif
	ConstADvar();
	inline ConstADvar& operator=(Double d) { this->cv->Val = d; return *this; }
	inline ConstADvar& operator=(ADVari& d) { this->cv->Val = d.Val; return *this; }
 };

 template<typename Double> class
ADvar1s: public ADvar1<Double> { // unary ops with partials
 public:
	typedef ADvar1<Double> ADVar1;
	typedef typename ADVar1::ADVari ADVari;
	Double pL;	// deriv of op w.r.t. left operand L
	ADvar1s(Double val1, Double a1, const ADVari *c1):
		ADVar1(Hv_timesL,val1,&pL,c1), pL(a1) {}
#ifdef RAD_AUTO_AD_Const
	Double a;
	typedef typename ADVar1::ADVar ADVar;
	ADvar1s(Double val1, Double a1, const ADVari *c1, const ADVar *v):
		ADVar1(Hv_timesL,val1,&a,c1,v), a(a1) {}
#endif
	};

 template<typename Double> class
ADvar1g: public ADvar1<Double> { // unary ops with partials
 public:
	typedef ADvar1<Double> ADVar1;
	typedef typename ADVar1::ADVari ADVari;
	Double pL;	// deriv of op w.r.t. left operand L
	Double pL2;	// partial of op w.r.t. L,L
	ADvar1g(Double val1, Double d1, Double d2, const ADVari *c1):
		ADVar1(Hv_unary,val1,&pL,c1), pL(d1), pL2(d2) {}
	};

 template<typename Double> class
ADvar2: public ADvari<Double> {	// basic binary ops
 public:
	typedef ADvari<Double> ADVari;
	typedef Derp<Double> DErp;
	DErp dL, dR;
	ADvar2(Advari_Opclass oc, Double val1, const ADVari *Lcv, const Double *Lc,
		const ADVari *Rcv, const Double *Rc):
			ADVari(oc,val1) {
		dR.next = DErp::LastDerp;
		dL.next = &dR;
		DErp::LastDerp = &dL;
		dL.a = Lc;
		dL.c = (ADVari*)Lcv;
		dR.a = Rc;
		dR.c = (ADVari*)Rcv;
		dL.b = dR.b = this;
		}
#ifdef RAD_AUTO_AD_Const
	typedef ADvar<Double> ADVar;
	ADvar2(Advari_Opclass oc, Double val1, const ADVari *Lcv, const Double *Lc,
		const ADVari *Rcv, const Double *Rc, ADVar *v):
			ADVari(oc,val1) {
		dR.next = DErp::LastDerp;
		dL.next = &dR;
		DErp::LastDerp = &dL;
		dL.a = Lc;
		dL.c = Lcv;
		dR.a = Rc;
		dR.c = Rcv;
		dL.b = dR.b = this;
		Lcv->padv = 0;
		this->ADvari_padv(v);
		}
#endif
	};

 template<typename Double> class
ADvar2q: public ADvar2<Double> { // binary ops with partials
 public:
	typedef ADvar2<Double> ADVar2;
	typedef typename ADVar2::ADVari ADVari;
	typedef typename ADVar2::DErp DErp;
	Double pL;	// deriv of op w.r.t. left operand L
	Double pR;	// deriv of op w.r.t. right operand R
	Double pLR;	// second partial w.r.t. L,R
	Double pR2;	// second partial w.r.t. R,R
	ADvar2q(Double val1, Double Lp, Double Rp, Double LR, Double R2,
		const ADVari *Lcv, const ADVari *Rcv):
			ADVar2(Hv_quotLR,val1,Lcv,&pL,Rcv,&pR),
			pL(Lp), pR(Rp), pLR(LR), pR2(R2) {}
#ifdef RAD_AUTO_AD_Const
	typedef typename ADVar2::ADVar ADVar;
	ADvar2q(Double val1, Double Lp, Double Rp, Double LR, Double R2,
		const ADVari *Lcv, const ADVari *Rcv, const ADVar *v):
			ADVar2(Hv_quotLR,val1,Lcv,&pL,Rcv,&pR,v),
			pL(Lp), pR(Rp), pLR(LR), pR2(R2) {}
#endif
	};

 template<typename Double> class
ADvar2g: public ADvar2<Double> { // general binary ops with partials
 public:
	typedef ADvar2<Double> ADVar2;
	typedef typename ADVar2::ADVari ADVari;
	Double pL;	// deriv of op w.r.t. left operand L
	Double pR;	// deriv of op w.r.t. right operand R
	Double pL2;	// second partial w.r.t. L,L
	Double pLR;	// second partial w.r.t. L,R
	Double pR2;	// second partial w.r.t. R,R
	ADvar2g(Double val1, Double Lp, Double Rp, Double L2, Double LR, Double R2,
		const ADVari *Lcv, const ADVari *Rcv):
			ADVar2(Hv_binary,val1,Lcv,&pL,Rcv,&pR),
			pL(Lp), pR(Rp), pL2(L2), pLR(LR), pR2(R2) { }
	};

 template<typename Double> class
ADvarn: public ADvari<Double> { // n-ary ops with partials g and
				// 2nd partials h (lower triangle, rowwise)
 public:
	typedef ADvari<Double> ADVari;
	typedef ADvar<Double> ADVar;
	typedef typename ADVari::IndepADVar IndepADVar;
	typedef Derp<Double> DErp;
	int n;
	Double *G, *H;
	DErp *D;
	ADvarn(Double val1, int n1, const IndepADVar *x, const Double *g, const Double *h):
			ADVari(Hv_nary,val1), n(n1) {
		DErp *d1, *dlast;
		Double *a1;
		int i, nh;

		a1 = G = (Double*)ADVari::adc.Memalloc(n1*sizeof(*g));
		d1 = D = (DErp*)ADVari::adc.Memalloc(n1*sizeof(DErp));
		dlast = DErp::LastDerp;
		for(i = 0; i < n1; i++, d1++) {
			d1->next = dlast;
			dlast = d1;
			a1[i] = g[i];
			d1->a = &a1[i];
			d1->b = this;
			d1->c = x[i].cv;
			}
		DErp::LastDerp = dlast;
		nh = (n1*(n1+1)) >> 1;
		a1 = H = (double*)ADVari::adc.Memalloc(nh * sizeof(*h));
		for(i = 0; i < nh; i++)
			a1[i] = h[i];
		}
	};

template<typename Double>
 inline ADvari<Double>& operator+(const ADvari<Double> &T) { return *(ADvari<Double>*)&T; }

template<typename Double>
 inline int operator<(const ADvari<Double> &L, const ADvari<Double> &R) { return L.Val < R.Val; }
template<typename Double>
 inline int operator<(const ADvari<Double> &L, Double R) { return L.Val < R; }
template<typename Double>
 inline int operator<(Double L, const ADvari<Double> &R) { return L < R.Val; }

template<typename Double>
 inline int operator<=(const ADvari<Double> &L, const ADvari<Double> &R) { return L.Val <= R.Val; }
template<typename Double>
 inline int operator<=(const ADvari<Double> &L, Double R) { return L.Val <= R; }
template<typename Double>
 inline int operator<=(Double L, const ADvari<Double> &R) { return L <= R.Val; }

template<typename Double>
 inline int operator==(const ADvari<Double> &L, const ADvari<Double> &R) { return L.Val == R.Val; }
template<typename Double>
 inline int operator==(const ADvari<Double> &L, Double R) { return L.Val == R; }
template<typename Double>
 inline int operator==(Double L, const ADvari<Double> &R) { return L == R.Val; }

template<typename Double>
 inline int operator!=(const ADvari<Double> &L, const ADvari<Double> &R) { return L.Val != R.Val; }
template<typename Double>
 inline int operator!=(const ADvari<Double> &L, Double R) { return L.Val != R; }
template<typename Double>
 inline int operator!=(Double L, const ADvari<Double> &R) { return L != R.Val; }

template<typename Double>
 inline int operator>=(const ADvari<Double> &L, const ADvari<Double> &R) { return L.Val >= R.Val; }
template<typename Double>
 inline int operator>=(const ADvari<Double> &L, Double R) { return L.Val >= R; }
template<typename Double>
 inline int operator>=(Double L, const ADvari<Double> &R) { return L >= R.Val; }

template<typename Double>
 inline int operator>(const ADvari<Double> &L, const ADvari<Double> &R) { return L.Val > R.Val; }
template<typename Double>
 inline int operator>(const ADvari<Double> &L, Double R) { return L.Val > R; }
template<typename Double>
 inline int operator>(Double L, const ADvari<Double> &R) { return L > R.Val; }

template<typename Double>
 inline void *ADcontext<Double>::Memalloc(size_t len) {
		if (Mleft >= len)
			return Mbase + (Mleft -= len);
		return new_ADmemblock(len);
		}

template<typename Double>
 inline Derp<Double>::Derp(const ADVari *c1): c((ADVari*)c1) {
		next = LastDerp;
		LastDerp = this;
		}

template<typename Double>
 inline Derp<Double>::Derp(const Double *a1, const ADVari *c1): a(a1), c((ADVari*)c1) {
		next = LastDerp;
		LastDerp = this;
		}

template<typename Double>
 inline Derp<Double>::Derp(const Double *a1, const ADVari *b1, const ADVari *c1):
	a(a1), b(b1), c((ADVari*)c1) {
		next = LastDerp;
		LastDerp = this;
		}

/**** radops ****/

template<typename Double> Derp<Double> *Derp<Double>::LastDerp = 0;
template<typename Double> ADcontext<Double> ADvari<Double>::adc;
template<typename Double> const Double ADcontext<Double>::One = 1.;
template<typename Double> const Double ADcontext<Double>::negOne = -1.;
template<typename Double> CADcontext<Double> ConstADvari<Double>::cadc;
template<typename Double> ConstADvari<Double> *ConstADvari<Double>::lastcad;

#ifdef RAD_AUTO_AD_Const
template<typename Double> ADvari<Double>*   ADvari<Double>::First_ADvari;
template<typename Double> ADvari<Double>**  ADvari<Double>::Last_ADvari = &ADvari<Double>::First_ADvari;
#endif

#ifdef RAD_DEBUG
#ifndef RAD_DEBUG_gcgen1
#define RAD_DEBUG_gcgen1 -1
#endif
template<typename Double> int ADvari<Double>::gcgen_cur;
template<typename Double> int ADvari<Double>::last_opno;
template<typename Double> int ADvari<Double>::zap_gcgen;
template<typename Double> int ADvari<Double>::zap_gcgen1 = RAD_DEBUG_gcgen1;
template<typename Double> int ADvari<Double>::zap_opno;
template<typename Double> FILE *ADvari<Double>::debug_file;
#endif


 template<typename Double>
ADcontext<Double>::ADcontext()
{
	ADVari_block *fb;

	First = (ADMemblock*)First0;
	First->next = 0;
	Busy = First;
	Free = 0;
	Mbase = (char*)First->memblk;
	Mleft = sizeof(First->memblk);
	AiFirst = Aibusy = fb = (ADVari_block*)First1;
	Aifree = 0;
	Ainext = fb->pADvari;
	fb->next = fb->prev = 0;
	fb->limit = Ailimit = fb->pADvari + ADVari_block::Gulp;
	rad_need_reinit = 0;
#ifdef RAD_DEBUG_BLOCKKEEP
	rad_busy_blocks = 0;
	rad_mleft_save = 0;
	rad_Oldcurmb = 0;
#endif
	}

template<typename Double> void*
ADcontext<Double>::new_ADmemblock(size_t len)
{
	ADMemblock *mb, *mb0, *mb1, *mbf, *x;
	ADVari_block *b;
#ifdef RAD_AUTO_AD_Const
	ADVari *a, *anext;
	IndepADvar<Double> *v;
#ifdef RAD_Const_WARN
	ADVari *cv;
	int i, j;
#endif
#endif /*RAD_AUTO_AD_Const*/

	if ((rad_need_reinit & 1) && this == &ADVari::adc) {
		rad_need_reinit &= ~1;
		DErp::LastDerp = 0;
		Aibusy = b = AiFirst;
		Aifree = b->next;
		b->next = b->prev = 0;
		Ailimit = b->limit = (Ainext = b->pADvari) + ADVari_block::Gulp;
#ifdef RAD_DEBUG_BLOCKKEEP
		Mleft = rad_mleft_save;
		if (Mleft < sizeof(First->memblk))
			_uninit_f2c(Mbase + Mleft,
				UninitType<Double>::utype,
			 	(sizeof(First->memblk) - Mleft)
				/sizeof(typename Sacado::ValueType<Double>::type));
		if ((mb = Busy->next)) {
			if (!(mb0 = rad_Oldcurmb))
				mb0 = (ADMemblock*)First0;
			for(;; mb = mb->next) {
				_uninit_f2c(mb->memblk,
					UninitType<Double>::utype,
					sizeof(First->memblk)
					/sizeof(typename Sacado::ValueType<Double>::type));
				if (mb == mb0)
					break;
				}
			}
		rad_Oldcurmb = Busy;
		if (rad_busy_blocks >= RAD_DEBUG_BLOCKKEEP) {
			rad_busy_blocks = 0;
			rad_Oldcurmb = 0;
			mb0 = (ADMemblock*)First0;
			mbf =  Free;
			for(mb = Busy; mb != mb0; mb = mb1) {
				mb1 = mb->next;
				mb->next = mbf;
				mbf = mb;
				}
			Free = mbf;
			Busy = mb;
			Mbase = (char*)First->memblk;
			Mleft = sizeof(First->memblk);
			}

#else /* !RAD_DEBUG_BLOCKKEEP */

		mb0 = First;
		mbf =  Free;
		for(mb = Busy; mb != mb0; mb = mb1) {
			mb1 = mb->next;
			mb->next = mbf;
			mbf = mb;
			}
		Free = mbf;
		Busy = mb;
		Mbase = (char*)First->memblk;
		Mleft = sizeof(First->memblk);
#ifdef RAD_AUTO_AD_Const
		if (ADVari::adc.rad_need_reinit & 2) {
			ADVari::adc.rad_need_reinit &= ~2;
			*ADVari::Last_ADvari = 0;
			ADVari::Last_ADvari = &ADVari::First_ADvari;
                        anext = ADVari::First_ADvari;
			if (anext) {
				while((a = anext)) {
					anext = a->Next;
					if ((v = (IndepADvar<Double> *)a->padv)) {
#ifdef RAD_Const_WARN
						if ((i = a->opno) > 0)
							i = -i;
						j = a->gcgen;
						v->cv = cv = new ADVari(v, a->Val, a->aval);
						cv->opno = i;
						cv->gcgen = j;
#else
						v->cv = new ADVari(v, a->Val, a->aval);
#endif
						}
					}
				}
			}
#endif /*RAD_AUTO_AD_Const*/
#endif /*RAD_DEBUG_BLOCKKEEP*/
		if (Mleft >= len)
			return Mbase + (Mleft -= len);
		}

	if ((x = Free))
		Free = x->next;
	else
		x = new ADMemblock;
#ifdef RAD_DEBUG_BLOCKKEEP
	rad_busy_blocks++;
#endif
	x->next = Busy;
	Busy = x;
	return (Mbase = (char*)x->memblk) +
		(Mleft = sizeof(First->memblk) - len);
	}

template<typename Double> void
ADcontext<Double>::new_ADvari_block()
{
	ADVari_block *ob, *nb;
	ob = Aibusy;
	ob->limit = Ailimit;	// should be unnecessary, but harmless
	if ((nb = Aifree))
		Aifree = nb->next;
	else
		nb = new ADVari_block;
	Aibusy = Aibusy->next = nb;
	nb->limit = Ailimit = (Ainext = nb->pADvari) + ADVari_block::Gulp;
	ob->next = nb;
	nb->prev = ob;
	nb->next = 0;
	}

template<typename Double> void
ADcontext<Double>::Gradcomp(int wantgrad)
{
	DErp *d;

	if (ADVari::adc.rad_need_reinit) {
		for(d = DErp::LastDerp; d; d = d->next)
			d->c->aval = 0;
		}
	if (!(ADVari::adc.rad_need_reinit & 1)) {
		ADVari::adc.rad_need_reinit = 1;
		ADVari::adc.rad_mleft_save = ADVari::adc.Mleft;
		ADVari::adc.Mleft = 0;
		}
#ifdef RAD_DEBUG
	if (ADVari::gcgen_cur == ADVari::zap_gcgen1) {
		const char *fname;
		if (!(fname = getenv("RAD_DEBUG_FILE")))
			fname = "rad_debug.out";
		else if (!*fname)
			fname = 0;
		if (fname)
			ADVari::debug_file = fopen(fname, "w");
		ADVari::zap_gcgen1 = -1;
		}
#endif
	if ((d = DErp::LastDerp) != 0) {
#ifdef RAD_AUTO_AD_Const
		 ADVari::adc.rad_need_reinit |= 2;
#endif /*RAD_AUTO_AD_Const*/
		if (wantgrad) {
			d->b->aval = 1;
#ifdef RAD_DEBUG
			if (ADVari::debug_file)
				do {
					fprintf(ADVari::debug_file, "%d\t%d\t%g + %g * %g",
						d->c->opno, d->b->opno, d->c->aval,
						*d->a, d->b->aval);
					d->c->aval += *d->a * d->b->aval;
					fprintf(ADVari::debug_file, " = %g\n", d->c->aval);
					fflush(ADVari::debug_file);
					} while((d = d->next));
			else
#endif
			do d->c->aval += *d->a * d->b->aval;
			while((d = d->next));
			DErp::LastDerp = 0;
			}
		}
#ifdef RAD_DEBUG
	if (ADVari::debug_file) {
		fclose(ADVari::debug_file);
		ADVari::debug_file = 0;
		}
#endif //RAD_DEBUG
#ifdef RAD_DEBUG
	ADVari::gcgen_cur++;
	ADVari::last_opno = 0;
#endif
	}

template<typename Double> void
ADcontext<Double>::Weighted_Gradcomp(int n, ADVar **V, Double *w)
{
	DErp *d;
	int i;
#ifdef RAD_Const_WARN
	ADVari *cv;
	int j;
#endif
#ifdef RAD_AUTO_AD_Const
	ADVari *a, *anext;
	IndepADvar<Double> *v;
#endif /*RAD_AUTO_AD_Const*/

	if (ADVari::adc.rad_need_reinit) {
		for(d = DErp::LastDerp; d; d = d->next)
			d->c->aval = 0;
		}
	if (!(ADVari::adc.rad_need_reinit & 1)) {
		ADVari::adc.rad_need_reinit = 1;
		ADVari::adc.rad_mleft_save = ADVari::adc.Mleft;
		ADVari::adc.Mleft = 0;
		}
#ifdef RAD_DEBUG
	if (ADVari::gcgen_cur == ADVari::zap_gcgen1) {
		const char *fname;
		if (!(fname = getenv("RAD_DEBUG_FILE")))
			fname = "rad_debug.out";
		else if (!*fname)
			fname = 0;
		if (fname)
			ADVari::debug_file = fopen(fname, "w");
		ADVari::zap_gcgen1 = -1;
		}
#endif
	if ((d = DErp::LastDerp) != 0) {
		for(i = 0; i < n; i++)
			V[i]->cv->aval = w[i];
#ifdef RAD_DEBUG
		if (ADVari::debug_file)
			do {
				fprintf(ADVari::debug_file, "%d\t%d\t%g + %g * %g",
					d->c->opno, d->b->opno, d->c->aval, *d->a, d->b->aval);
				d->c->aval += *d->a * d->b->aval;
				fprintf(ADVari::debug_file, " = %g\n", d->c->aval);
				fflush(ADVari::debug_file);
				} while((d = d->next));
		else
#endif
		do d->c->aval += *d->a * d->b->aval;
		while((d = d->next));
		}
#ifdef RAD_DEBUG
	if (ADVari::debug_file) {
		fclose(ADVari::debug_file);
		ADVari::debug_file = 0;
		}
#endif //RAD_DEBUG
#ifdef RAD_AUTO_AD_Const
	*ADVari::Last_ADvari = 0;
	ADVari::Last_ADvari = &ADVari::First_ADvari;
	if ((anext = ADVari::First_ADvari) && !(ADVari::adc.rad_need_reinit & 2)) {
		ADVari::adc.rad_need_reinit = 3;
		while((a = anext)) {
			anext = a->Next;
			if ((v = (IndepADvar<Double> *)a->padv)) {
#ifdef RAD_Const_WARN
				if ((i = a->opno) > 0)
					i = -i;
				j = a->gcgen;
				v->cv = cv = new ADVari(v, a->Val, a->aval);
				cv->opno = i;
				cv->gcgen = j;
#else
				v->cv = new ADVari(v, a->Val, a->aval);
#endif
				}
			}
		DErp::LastDerp = 0;
		}
#endif /*RAD_AUTO_AD_Const*/
#ifdef RAD_DEBUG
	ADVari::gcgen_cur++;
	ADVari::last_opno = 0;
#endif
		}

 template<typename Double>
IndepADvar<Double>::IndepADvar(Ttype d)
{

	ADVari *x = new ADVari(Hv_const, d);
	cv = x;
	}

 template<typename Double>
IndepADvar<Double>::IndepADvar(double i)
{

	ADVari *x = new ADVari(Hv_const, Double(i));
	cv = x;
	}

 template<typename Double>
IndepADvar<Double>::IndepADvar(int i)
{

	ADVari *x = new ADVari(Hv_const, Double(i));
	cv = x;
	}

 template<typename Double>
IndepADvar<Double>::IndepADvar(long i)
{

	ADVari *x = new ADVari(Hv_const, Double(i));
	cv = x;
	}

 template<typename Double>
ConstADvar<Double>::ConstADvar()
{
	ConstADVari *x = new ConstADVari(0.);
	this->cv = x;
	}

 template<typename Double> void
ConstADvar<Double>::ConstADvar_ctr(Double d)
{
	ConstADVari *x = new ConstADVari(d);
	this->cv = x;
	}

 template<typename Double>
ConstADvar<Double>::ConstADvar(const IndepADVar &x)
{
	ConstADVari *y = new ConstADVari(x.cv->Val);
	DErp *d = new DErp(&x.adc.One, y, x.cv);
	this->cv = y;
	}

 template<typename Double>
ConstADvar<Double>::ConstADvar(const ConstADvar &x)
{
	ConstADVari *y = new ConstADVari(x.cv->Val);
	DErp *d = new DErp(&x.cv->adc.One, y, (ADVari*)x.cv);
	this->cv = y;
	}

 template<typename Double>
ConstADvar<Double>::ConstADvar(const ADVari &x)
{
	ConstADVari *y = new ConstADVari(x.Val);
	DErp *d = new DErp(&x.adc.One, y, &x);
	this->cv = y;
	}

 template<typename Double>
 void
IndepADvar<Double>::AD_Const(const IndepADvar &v)
{
	typedef ConstADvari<Double> ConstADVari;

	ConstADVari *ncv = new ConstADVari(v.val());
#ifdef RAD_AUTO_AD_Const
	if (v.cv)
		v.cv->padv = 0;
#endif
	v.cv = ncv;
	}

 template<typename Double>
 void
ConstADvari<Double>::aval_reset()
{
	ConstADvari *x = ConstADvari::lastcad;
	while(x) {
		x->aval = 0;
		x = x->prevcad;
		}
	}

#ifdef RAD_AUTO_AD_Const

 template<typename Double>
ADvari<Double>::ADvari(const IndepADVar *x, Double d): Val(d), aval(0.)
{
	opclass = Hv_const;
	this->ADvari_padv(x);
	}

 template<typename Double>
ADvari<Double>::ADvari(const IndepADVar *x, Double d, Double g): Val(d), aval(g)
{
	opclass = Hv_const;
	this->ADvari_padv(x);
	}

 template<typename Double>
ADvar1<Double>::ADvar1(const IndepADVar *x, const IndepADVar &y):
	ADVari(Hv_copy, y.cv->Val), d((const Double*)&ADcontext<Double>::One, (ADVari*)this, y.cv)
{
	this->ADvari_padv(x);
	}

 template<typename Double>
ADvar1<Double>::ADvar1(const IndepADVar *x, const ADVari &y):
	ADVari(Hv_copy, y.Val), d((const Double*)&ADcontext<Double>::One, this, &y)
{
	this->ADvari_padv(x);
	}

#endif /* RAD_AUTO_AD_Const */

 template<typename Double>
 IndepADvar<Double>&
ADvar_operatoreq(IndepADvar<Double> *This, const ADvari<Double> &x)
{ This->cv = new ADvar1<Double>(Hv_copy, x.Val, &x.adc.One, &x);
  return *(IndepADvar<Double>*) This; }

 template<typename Double>
 ADvar<Double>&
ADvar_operatoreq(ADvar<Double> *This, const ADvari<Double> &x)
{ This->cv = new ADvar1<Double>(Hv_copy, x.Val, &x.adc.One, &x); return *(ADvar<Double>*) This; }

 template<typename Double>
 IndepADvar<Double>&
IndepADvar<Double>::operator=(Double d)
{
#ifdef RAD_AUTO_AD_Const
	ADVari *ncv = new ADVari(this, d);
	if (this->cv)
		this->cv->padv = 0;
	this->cv = ncv;
#else
	this->cv = new ADVari(Hv_const, d);
#endif
	return *this;
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator=(Double d)
{
#ifdef RAD_AUTO_AD_Const
	ADVari *nv = new ADVari(this, d);
	if (this->cv)
		this->cv->padv = 0;
	this->cv = nv;
#else
	this->cv = ConstADVari::cadc.fpval_implies_const
		? new ConstADVari(d)
		: new ADVari(Hv_const, d);
#endif
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator-(const ADvari<Double> &T) {
	return *(new ADvar1<Double>(Hv_negate, -T.Val, &T.adc.negOne, &T));
	}

 template<typename Double>
 ADvari<Double>&
operator+(const ADvari<Double> &L, const ADvari<Double> &R) {
	return *(new ADvar2<Double>(Hv_plusLR, L.Val + R.Val, &L, &L.adc.One, &R, &L.adc.One));
	}

#ifdef RAD_AUTO_AD_Const
#define RAD_ACA ,this
#else
#define RAD_ACA /*nothing*/
#endif

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator+=(const ADVari &R) {
	ADVari *Lcv = this->cv;
	this->cv = new ADvar2<Double>(Hv_plusLR, Lcv->Val + R.Val, Lcv,
					&R.adc.One, &R, &R.adc.One RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator+(const ADvari<Double> &L, Double R) {
	return *(new ADvar1<Double>(Hv_copy, L.Val + R, &L.adc.One, &L));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator+=(Double R) {
	ADVari *tcv = this->cv;
	this->cv = new ADVar1(Hv_copy, tcv->Val + R, &tcv->adc.One, tcv RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator+(Double L, const ADvari<Double> &R) {
	return *(new ADvar1<Double>(Hv_copy, L + R.Val, &R.adc.One, &R));
	}

 template<typename Double>
 ADvari<Double>&
operator-(const ADvari<Double> &L, const ADvari<Double> &R) {
	return *(new ADvar2<Double>(Hv_minusLR, L.Val - R.Val, &L, &L.adc.One, &R, &L.adc.negOne));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator-=(const ADVari &R) {
	ADVari *Lcv = this->cv;
	this->cv = new ADvar2<Double>(Hv_minusLR,Lcv->Val - R.Val, Lcv,
					&R.adc.One, &R, &R.adc.negOne RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator-(const ADvari<Double> &L, Double R) {
	return *(new ADvar1<Double>(Hv_copy, L.Val - R, &L.adc.One, &L));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator-=(Double R) {
	ADVari *tcv = this->cv;
	this->cv = new ADVar1(Hv_copy, tcv->Val - R, &tcv->adc.One, tcv RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator-(Double L, const ADvari<Double> &R) {
	return *(new ADvar1<Double>(Hv_negate, L - R.Val, &R.adc.negOne, &R));
	}

 template<typename Double>
 ADvari<Double>&
operator*(const ADvari<Double> &L, const ADvari<Double> &R) {
	return *(new ADvar2<Double>(Hv_timesLR, L.Val * R.Val, &L, &R.Val, &R, &L.Val));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator*=(const ADVari &R) {
	ADVari *Lcv = this->cv;
	this->cv = new ADvar2<Double>(Hv_timesLR, Lcv->Val * R.Val, Lcv,
					&R.Val, &R, &Lcv->Val RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator*(const ADvari<Double> &L, Double R) {
	return *(new ADvar1s<Double>(L.Val * R, R, &L));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator*=(Double R) {
	ADVari *Lcv = this->cv;
	this->cv = new ADvar1s<Double>(Lcv->Val * R, R, Lcv RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator*(Double L, const ADvari<Double> &R) {
	return *(new ADvar1s<Double>(L * R.Val, L, &R));
	}

 template<typename Double>
 ADvari<Double>&
operator/(const ADvari<Double> &L, const ADvari<Double> &R) {
	Double Lv = L.Val, Rv = R.Val, pL = 1. / Rv, q = Lv/Rv, qpL = q*pL;
	return *(new ADvar2q<Double>(q, pL, -qpL, -pL*pL, 2.*pL*qpL, &L, &R));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator/=(const ADVari &R) {
	ADVari *Lcv = this->cv;
	Double Lv = Lcv->Val, Rv = R.Val, pL = 1. / Rv, q = Lv/Rv, qpL = q*pL;
	this->cv = new ADvar2q<Double>(q, pL, -qpL, -pL*pL, 2.*pL*qpL, Lcv, &R RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
operator/(const ADvari<Double> &L, Double R) {
	return *(new ADvar1s<Double>(L.Val / R, 1./R, &L));
	}

 template<typename Double>
 ADvari<Double>&
operator/(Double L, const ADvari<Double> &R) {
	Double recip = 1. / R.Val;
	Double q = L * recip;
	Double d1 = -q*recip;
	return *(new ADvar1g<Double>(q, d1, -q*d1, &R));
	}

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator/=(Double R) {
	ADVari *Lcv = this->cv;
	this->cv = new ADvar1s<Double>(Lcv->Val / R, 1./R, Lcv RAD_ACA);
	return *this;
	}

 template<typename Double>
 ADvari<Double>&
acos(const ADvari<Double> &v) {
	Double t = v.Val, t1 = 1. - t*t, d1 = -1./std::sqrt(t1);
	return *(new ADvar1g<Double>(std::acos(t), d1, t*d1/t1, &v));
	}

 template<typename Double>
 ADvari<Double>&
acosh(const ADvari<Double> &v) {
	Double d1, t, t1, t2;
	t = v.Val;
	t1 = std::sqrt(t2 = t*t - 1.);
	d1 = 1. / t1;
	return *(new ADvar1g<Double>(std::log(t + t1), d1, -t*d1/t2, &v));
	}

 template<typename Double>
 ADvari<Double>&
asin(const ADvari<Double> &v) {
	Double d1, t, t1;
	t = v.Val;
	d1 = 1. / std::sqrt(t1 = 1. - t*t);
	return *(new ADvar1g<Double>(std::asin(t), d1, t*d1/t1, &v));
	}

 template<typename Double>
 ADvari<Double>&
asinh(const ADvari<Double> &v) {
	Double d1, t, t1, t2, td;
	t = v.Val;
	t1 = std::sqrt(t2 = t*t + 1.);
	d1 = 1. / t1;
	td = 1.;
	if (t < 0.)
		td = -1.;
	return *(new ADvar1g<Double>(td*std::log(t*td + t1), d1, -(t/t2)*d1, &v));
	}

 template<typename Double>
 ADvari<Double>&
atan(const ADvari<Double> &v) {
	Double t = v.Val, d1 = 1./(1. + t*t);
	return *(new ADvar1g<Double>(std::atan(t), d1, -(t+t)*d1*d1, &v));
	}

 template<typename Double>
 ADvari<Double>&
atanh(const ADvari<Double> &v) {
	Double t = v.Val, d1 = 1./(1. - t*t);
	return *(new ADvar1g<Double>(0.5*std::log((1.+t)/(1.-t)), d1, (t+t)*d1*d1, &v));
	}

 template<typename Double>
 ADvari<Double>&
atan2(const ADvari<Double> &L, const ADvari<Double> &R) {
	Double R2, t, t2, x, x2, y, y2;
	x = L.Val;
	y = R.Val;
	x2 = x*x;
	y2 = y*y;
	t = 1./(x2 + y2);
	t2 = t*t;
	R2 = 2.*t2*x*y;
	return *(new ADvar2g<Double>(std::atan2(x,y), y*t, -x*t, -R2, t2*(x2 - y2), R2, &L, &R));
	}

 template<typename Double>
 ADvari<Double>&
atan2(Double x, const ADvari<Double> &R) {
	Double t, x2, y, y2;
	y = R.Val;
	x2 = x*x;
	y2 = y*y;
	t = 1./(x2 + y2);
	return *(new ADvar1g<Double>(std::atan2(x,y), -x*t, 2.*t*t*x*y, &R));
	}

 template<typename Double>
 ADvari<Double>&
atan2(const ADvari<Double> &L, Double y) {
	Double t, x, x2, y2;
	x = L.Val;
	x2 = x*x;
	y2 = y*y;
	t = 1./(x2 + y2);
	return *(new ADvar1g<Double>(std::atan2(x,y), y*t, -2.*t*t*x*y, &L));
	}

 template<typename Double>
 ADvari<Double>&
max(const ADvari<Double> &L, const ADvari<Double> &R) {
	const ADvari<Double> &x = L.Val >= R.Val ? L : R;
	return *(new ADvar1<Double>(Hv_copy, x.Val, &x.adc.One, &x));
	}

 template<typename Double>
 ADvari<Double>&
max(Double L, const ADvari<Double> &R) {
	if (L >= R.Val)
		return *(new ADvari<Double>(Hv_const, L));
	return *(new ADvar1<Double>(Hv_copy, R.Val, &R.adc.One, &R));
	}

 template<typename Double>
 ADvari<Double>&
max(const ADvari<Double> &L, Double R) {
	if (L.Val >= R)
		return *(new ADvar1<Double>(Hv_copy, L.Val, &L.adc.One, &L));
	return *(new ADvari<Double>(Hv_const, R));
	}

 template<typename Double>
 ADvari<Double>&
min(const ADvari<Double> &L, const ADvari<Double> &R) {
	const ADvari<Double> &x = L.Val <= R.Val ? L : R;
	return *(new ADvar1<Double>(Hv_copy, x.Val, &x.adc.One, &x));
	}

 template<typename Double>
 ADvari<Double>&
min(Double L, const ADvari<Double> &R) {
	if (L <= R.Val)
		return *(new ADvari<Double>(Hv_const, L));
	return *(new ADvar1<Double>(Hv_copy, R.Val, &R.adc.One, &R));
	}

 template<typename Double>
 ADvari<Double>&
min(const ADvari<Double> &L, Double R) {
	if (L.Val <= R)
		return *(new ADvar1<Double>(Hv_copy, L.Val, &L.adc.One, &L));
	return *(new ADvari<Double>(Hv_const, R));
	}

 template<typename Double>
 ADvari<Double>&
cos(const ADvari<Double> &v) {
	Double t = std::cos(v.Val);
	return *(new ADvar1g<Double>(t, -std::sin(v.Val), -t, &v));
	}

 template<typename Double>
 ADvari<Double>&
cosh(const ADvari<Double> &v) {
	Double t = std::cosh(v.Val);
	return *(new ADvar1g<Double>(t, std::sinh(v.Val), t, &v));
	}

 template<typename Double>
 ADvari<Double>&
exp(const ADvari<Double> &v) {
	Double t = std::exp(v.Val);
	return *(new ADvar1g<Double>(t, t, t, &v));
	}

 template<typename Double>
 ADvari<Double>&
log(const ADvari<Double> &v) {
	Double x = v.Val, d1 = 1. / x;
	return *(new ADvar1g<Double>(std::log(x), d1, -d1*d1, &v));
	}

 template<typename Double>
 ADvari<Double>&
log10(const ADvari<Double> &v) {
	static double num = 1. / std::log(10.);
	Double d1, t, x;
	x = v.Val;
	t = 1. / x;
	d1 = num * t;
	return *(new ADvar1g<Double>(std::log10(x), d1, -d1*t, &v));
	}

 template<typename Double>
 ADvari<Double>&
pow(const ADvari<Double> &L, const ADvari<Double> &R) {
	Double dx, dy, t, x, xlog, xym1, y;
	x = L.Val;
	y = R.Val;
	t = std::pow(x,y);
	xym1 = t / x;
	xlog = std::log(x);
	dx = y*xym1;
	dy = t * xlog;
	return *(new ADvar2g<Double>(t, dx, dy, (y-1.)*dx/x, xym1*(1. + y*xlog), dy*xlog, &L, &R));
	}

 template<typename Double>
 ADvari<Double>&
pow(Double x, const ADvari<Double> &R) {
	Double dy, t, xlog, y;
	y = R.Val;
	t = std::pow(x,y);
	xlog = std::log(x);
	dy = t * xlog;
	return *(new ADvar1g<Double>(t, dy, dy*xlog, &R));
	}

 template<typename Double>
 ADvari<Double>&
pow(const ADvari<Double> &L, Double y) {
	Double dx, t, x;
	x = L.Val;
	t = std::pow(x,y);
	dx = y*t/x;
	return *(new ADvar1g<Double>(t, dx, (y-1.)*dx/x, &L));
	}

 template<typename Double>
 ADvari<Double>&
sin(const ADvari<Double> &v) {
	Double t = std::sin(v.Val);
	return *(new ADvar1g<Double>(t, std::cos(v.Val), -t, &v));
	}

 template<typename Double>
 ADvari<Double>&
sinh(const ADvari<Double> &v) {
	Double t = std::sinh(v.Val);
	return *(new ADvar1g<Double>(t, std::cosh(v.Val), t, &v));
	}

 template<typename Double>
 ADvari<Double>&
sqrt(const ADvari<Double> &v) {
	Double t = std::sqrt(v.Val);
	Double d1 = 0.5 / t;
	return *(new ADvar1g<Double>(t, d1, -0.5*d1/v.Val, &v));
	}

 template<typename Double>
 ADvari<Double>&
tan(const ADvari<Double> &v) {
	Double d1, rv, t;
	rv = std::tan(v.Val);
	t = 1. / std::cos(v.Val);
	d1 = t*t;
	return *(new ADvar1g<Double>(rv, d1, (rv+rv)*d1, &v));
	}

 template<typename Double>
 ADvari<Double>&
tanh(const ADvari<Double> &v) {
	Double d1, rv, t;
	rv = std::tanh(v.Val);
	t = 1. / std::cosh(v.Val);
	d1 = t*t;
	return *(new ADvar1g<Double>(rv, d1, -(rv+rv)*d1, &v));
	}

 template<typename Double>
 ADvari<Double>&
abs(const ADvari<Double> &v) {
	Double t, p;
	p = 1.;
	if ((t = v.Val) < 0) {
		t = -t;
		p = -p;
		}
	return *(new ADvar1g<Double>(t, p, 0., &v));
	}

 template<typename Double>
 ADvari<Double>&
fabs(const ADvari<Double> &v) {	// Synonym for "abs"
				// "fabs" is not the best choice of name,
				// but this name is used at Sandia.
	Double t, p;
	p = 1.;
	if ((t = v.Val) < 0) {
		t = -t;
		p = -p;
		}
	return *(new ADvar1g<Double>(t, p, 0., &v));
	}

 template<typename Double>
 ADvari<Double>&
ADf1(Double f, Double g, Double h, const ADvari<Double> &x) {
	return *(new ADvar1g<Double>(f, g, h, &x));
	}

 template<typename Double>
 inline ADvari<Double>&
ADf1(Double f, Double g, Double h, const IndepADvar<Double> &x) {
	return *(new ADvar1g<Double>(f, g, h, x.cv));
	}

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, Double hxx, Double hxy, Double hyy,
		const ADvari<Double> &x, const ADvari<Double> &y) {
	return *(new ADvar2g<Double>(f, gx, gy, hxx, hxy, hyy, &x, &y));
	}

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, Double hxx, Double hxy, Double hyy,
		const ADvari<Double> &x, const IndepADvar<Double> &y) {
	return *(new ADvar2g<Double>(f, gx, gy, hxx, hxy, hyy, &x, y.cv));
	}

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, Double hxx, Double hxy, Double hyy,
		const IndepADvar<Double> &x, const ADvari<Double> &y) {
	return *(new ADvar2g<Double>(f, gx, gy, hxx, hxy, hyy, x.cv, &y));
	}

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, Double hxx, Double hxy, Double hyy,
		const IndepADvar<Double> &x, const IndepADvar<Double> &y) {
	return *(new ADvar2g<Double>(f, gx, gy, hxx, hxy, hyy, x.cv, y.cv));
	}

 template<typename Double>
 ADvari<Double>&
ADfn(Double f, int n, const IndepADvar<Double> *x, const Double *g, const Double *h) {
	return *(new ADvarn<Double>(f, n, x, g, h));
	}

 template<typename Double>
 inline ADvari<Double>&
ADfn(Double f, int n, const ADvar<Double> *x, const Double *g, const Double *h) {
	return ADfn<Double>(f, n, (IndepADvar<Double>*)x, g, h);
	}

 template<typename Double>
 void
ADcontext<Double>::Hvprod(int n, ADvar<Double> **x, Double *v, Double *hv)
{
	ADVari *a, *aL, *aR, **ap, **ape;
	ADVari_block *b, *b0;
	DErp *d;
	Double aO, adO, *g, *h, *h0, t, tL, tR;
	int i, j, k, m;
	for(i = 0; i < n; i++) {
		a = x[i]->cv;
		a->dO = v[i];
		a->aO = a->adO = 0.;
		}
	ADVari::adc.Aibusy->limit = ADVari::adc.Ainext;
	a = 0;
	for(b0 = 0, b = ADVari::adc.AiFirst; b; b0 = b, b = b->next) {
		ap = b->pADvari;
		ape = b->limit;
		while(ap < ape) {
			a = *ap++;
			a->aO = a->adO = 0.;
			switch(a->opclass) {
			 case Hv_copy:
				a->dO = ((ADVar1*)a)->d.c->dO;
				break;
			 case Hv_binary:
				a->dO =   ((ADVar2g*)a)->pL * ((ADVar2g*)a)->dL.c->dO
					+ ((ADVar2g*)a)->pR * ((ADVar2g*)a)->dR.c->dO;
				break;
			 case Hv_unary:
				a->dO = ((ADVar1g*)a)->pL * ((ADVar1g*)a)->d.c->dO;
				break;
			 case Hv_negate:
				a->dO = -((ADVar1*)a)->d.c->dO;
				break;
			 case Hv_plusLR:
				a->dO = ((ADVar2*)a)->dL.c->dO + ((ADVar2*)a)->dR.c->dO;
				break;
			 case Hv_minusLR:
				a->dO = ((ADVar2*)a)->dL.c->dO - ((ADVar2*)a)->dR.c->dO;
				break;
			 case Hv_timesL:
				a->dO = ((ADVar1s*)a)->pL * ((ADVar1s*)a)->d.c->dO;
				break;
			 case Hv_timesLR:
				a->dO =   ((ADVar2*)a)->dR.c->Val * ((ADVar2*)a)->dL.c->dO
					+ ((ADVar2*)a)->dL.c->Val * ((ADVar2*)a)->dR.c->dO;
				break;
			 case Hv_quotLR:
				a->dO =   ((ADVar2q*)a)->pL * ((ADVar2q*)a)->dL.c->dO
					+ ((ADVar2q*)a)->pR * ((ADVar2q*)a)->dR.c->dO;
				break;
			 case Hv_nary:
				d = ((ADVarn*)a)->D;
				m = ((ADVarn*)a)->n;
				g = ((ADVarn*)a)->G;
				t = 0.;
				for(i = 0; i < m; i++)
					t += g[i] * d[i].c->dO;
				a->dO = t;
                                break;
                        case Hv_const:
                          ;
			 }
			}
		}
	if (a)
		a->adO = 1.;
	for(b = b0; b; b = b->prev) {
		ape = b->pADvari;
		ap = b->limit;
		while(ap > ape) {
			a = *--ap;
			aO = a->aO;
			adO = a->adO;
			switch(a->opclass) {
			 case Hv_copy:
				aL = ((ADVar1*)a)->d.c;
				aL->aO += aO;
				aL->adO += adO;
				break;
			 case Hv_binary:
				aL = ((ADVar2g*)a)->dL.c;
				aR = ((ADVar2g*)a)->dR.c;
				tL = adO*aL->dO;
				tR = adO*aR->dO;
				aL->aO += aO*((ADVar2g*)a)->pL
					+ tL*((ADVar2g*)a)->pL2
					+ tR*((ADVar2g*)a)->pLR;
				aR->aO += aO*((ADVar2g*)a)->pR
					+ tL*((ADVar2g*)a)->pLR
					+ tR*((ADVar2g*)a)->pR2;
				aL->adO += adO * ((ADVar2g*)a)->pL;
				aR->adO += adO * ((ADVar2g*)a)->pR;
				break;
			 case Hv_unary:
				aL = ((ADVar1g*)a)->d.c;
				aL->aO += aO * ((ADVar1g*)a)->pL
					+ adO * aL->dO * ((ADVar1g*)a)->pL2;
				aL->adO += adO * ((ADVar1g*)a)->pL;
				break;
			 case Hv_negate:
				aL = ((ADVar1*)a)->d.c;
				aL->aO -= aO;
				aL->adO -= adO;
				break;
			 case Hv_plusLR:
				aL = ((ADVar2*)a)->dL.c;
				aR = ((ADVar2*)a)->dR.c;
				aL->aO += aO;
				aL->adO += adO;
				aR->aO += aO;
				aR->adO += adO;
				break;
			 case Hv_minusLR:
				aL = ((ADVar2*)a)->dL.c;
				aR = ((ADVar2*)a)->dR.c;
				aL->aO += aO;
				aL->adO += adO;
				aR->aO -= aO;
				aR->adO -= adO;
				break;
			 case Hv_timesL:
				aL = ((ADVar1s*)a)->d.c;
				aL->aO += aO * (tL = ((ADVar1s*)a)->pL);
				aL->adO += adO * tL;
				break;
			 case Hv_timesLR:
				aL = ((ADVar2*)a)->dL.c;
				aR = ((ADVar2*)a)->dR.c;
				aL->aO += aO * (tL = aR->Val) + adO*aR->dO;
				aR->aO += aO * (tR = aL->Val) + adO*aL->dO;
				aL->adO += adO * tL;
				aR->adO += adO * tR;
				break;
			 case Hv_quotLR:
				aL = ((ADVar2q*)a)->dL.c;
				aR = ((ADVar2q*)a)->dR.c;
				tL = adO*aL->dO;
				tR = adO*aR->dO;
				aL->aO += aO*((ADVar2q*)a)->pL
					+ tR*((ADVar2q*)a)->pLR;
				aR->aO += aO*((ADVar2q*)a)->pR
					+ tL*((ADVar2q*)a)->pLR
					+ tR*((ADVar2q*)a)->pR2;
				aL->adO += adO * ((ADVar2q*)a)->pL;
				aR->adO += adO * ((ADVar2q*)a)->pR;
				break;
			 case Hv_nary:
				d  = ((ADVarn*)a)->D;
				m  = ((ADVarn*)a)->n;
				g  = ((ADVarn*)a)->G;
				h0 = ((ADVarn*)a)->H;
				for(i = 0; i < m; i++) {
					aL = d[i].c;
					aL->adO += adO * (t = g[i]);
					aL->aO += t*aO;
					t = adO * aL->dO;
					for(h = h0, j = 0; j <= i; j++)
						d[j].c->aO += t * *h++;
					h0 = h--;
					for(k = j; j < m; j++)
						d[j].c->aO += t * *(h += k++);
					}
                        case Hv_const:
                          ;
			 }
			}
		}
	for(i = 0; i < n; i++) {
		a = x[i]->cv;
		a->dO = 0.;
		hv[i] = a->aO;
		}
	}

 template<typename Double>
 inline Double
val(const ADvari<Double> &x) {
	return x.Val;
	}

#undef RAD_ACA
#define A (ADvari<Double>*)
#ifdef RAD_Const_WARN
#define C(x) (((x)->opno < 0) ? RAD_Const_Warn(x) : 0, *A x)
#else
#define C(x) *A x
#endif
#define T template<typename Double> inline
#define F ADvari<Double>&
#define Ai const ADvari<Double>&
#define AI const IndepADvar<Double>&
#define D Double
#define T2(r,f) \
 T r f(Ai L, AI R) { return f(L, C(R.cv)); }\
 T r f(AI L, Ai R) { return f(C(L.cv), R); }\
 T r f(AI L, AI R) { return f(C(L.cv), C(R.cv)); }\
 T r f(AI L, D R) { return f(C(L.cv), R); }\
 T r f(Ai L, Dtype R) { return f(L, (D)R); }\
 T r f(AI L, Dtype R) { return f(C(L.cv), (D)R); }\
 T r f(Ai L, long R) { return f(L, (D)R); }\
 T r f(AI L, long R) { return f(C(L.cv), (D)R); }\
 T r f(Ai L, int R) { return f(L, (D)R); }\
 T r f(AI L, int R) { return f(C(L.cv), (D)R); }\
 T r f(D L, AI R) { return f(L, C(R.cv)); }\
 T r f(Dtype L, Ai R) { return f((D)L, R); }\
 T r f(Dtype L, AI R) { return f((D)L, C(R.cv)); }\
 T r f(long L, Ai R) { return f((D)L, R); }\
 T r f(long L, AI R) { return f((D)L, C(R.cv)); }\
 T r f(int L, Ai R) { return f((D)L, R); }\
 T r f(int L, AI R) { return f((D)L, C(R.cv)); }

T2(F, operator+)
T2(F, operator-)
T2(F, operator*)
T2(F, operator/)
T2(F, atan2)
T2(F, pow)
T2(F, max)
T2(F, min)
T2(int, operator<)
T2(int, operator<=)
T2(int, operator==)
T2(int, operator!=)
T2(int, operator>=)
T2(int, operator>)

#undef T2
#undef D

#define T1(f)\
 T F f(AI x) { return f(C(x.cv)); }

T1(operator+)
T1(operator-)
T1(abs)
T1(acos)
T1(acosh)
T1(asin)
T1(asinh)
T1(atan)
T1(atanh)
T1(cos)
T1(cosh)
T1(exp)
T1(log)
T1(log10)
T1(sin)
T1(sinh)
T1(sqrt)
T1(tan)
T1(tanh)
T1(fabs)

T F copy(AI x)
{ return *(new ADvar1<Double>(Hv_copy, x.cv->Val, &ADcontext<Double>::One, (ADvari<Double>*)x.cv)); }

T F copy(Ai x)
{ return *(new ADvar1<Double>(Hv_copy, x.Val, &ADcontext<Double>::One, (ADvari<Double>*)&x)); }

#undef T1
#undef AI
#undef Ai
#undef F
#undef T
#undef A
#undef C
#undef Ttype
#undef Dtype

} /* namespace Rad2 */
} /* namespace Sacado */

#endif /* SACADO_TRAD2_H */
