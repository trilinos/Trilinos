// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// TRAD package (Templated Reverse Automatic Differentiation) --
// a package specialized for function and gradient evaluations.
// Written in 2004 and 2005 by David M. Gay at Sandia National Labs,
// Albuquerque, NM.

#ifndef SACADO_TRAD_H
#define SACADO_TRAD_H

#include "Sacado_ConfigDefs.h"
#include "Sacado_trad_Traits.hpp"
#include "Sacado_Base.hpp"

#if defined(RAD_DEBUG_BLOCKKEEP) && !defined(HAVE_SACADO_UNINIT)
#undef RAD_DEBUG_BLOCKKEEP
#endif

#ifndef RAD_REINIT
#define RAD_REINIT 0
#endif //RAD_REINIT

// RAD_ALLOW_WANTDERIV adds little overhead, so for simplicity
// we make it the default, which can be cancelled by compiling
// with -DRAD_DISALLOW_WANTDERIV:

#undef RAD_ALLOW_WANTDERIV
#ifndef RAD_DISALLOW_WANTDERIV
#define RAD_ALLOW_WANTDERIV
#endif

// Adjust names to avoid confusion when different source files
// are compiled with different RAD_REINIT settings.

#define RAD_REINIT_0(x) /*nothing*/
#define RAD_REINIT_1(x) /*nothing*/
#define RAD_REINIT_2(x) /*nothing*/
#define RAD_cvchk(x)    /*nothing*/

#if RAD_REINIT == 1 //{{
#undef  RAD_REINIT_1
#define RAD_REINIT_1(x) x
#ifdef RAD_ALLOW_WANTDERIV
#define ADvar ADvar_1n
#define IndepADvar IndepADvar_1n
#else
#define ADvar ADvar_1
#define IndepADvar IndepADvar_1
#endif //RAD_ALLOW_WANTDERIV
#elif RAD_REINIT == 2 //}{
#undef  RAD_REINIT_2
#define RAD_REINIT_2(x) x
#undef  RAD_cvchk
#define RAD_cvchk(x) if (x.gen != x.IndepADvar_root.gen) { x.cv = new ADvari<Double>(x.Val);\
        x.gen = x.IndepADvar_root.gen; }
#ifdef RAD_ALLOW_WANTDERIV
#define IndepADvar IndepADvar_2n
#define ADvar ADvar_2n
#else
#define IndepADvar IndepADvar_2
#define ADvar ADvar_2
#endif //RAD_ALLOW_WANTDERIV
#elif RAD_REINIT != 0 //}{
Botch!  Expected "RAD_REINIT" to be 0, 1, or 2.
Got "RAD_REINIT = " RAD_REINIT .
#else //}{
#undef RAD_ALLOW_WANTDERIV
#undef RAD_REINIT_0
#define RAD_REINIT_0(x) x
#endif //}}

#ifdef RAD_ALLOW_WANTDERIV
#define Allow_noderiv(x) x
#else
#define Allow_noderiv(x) /*nothing*/
#endif

#if RAD_REINIT > 0
#undef RAD_Const_WARN
#undef RAD_AUTO_AD_Const
#undef RAD_DEBUG_BLOCKKEEP
#endif

#include <stddef.h>
#include <Sacado_cmath.hpp>

#ifdef RAD_Const_WARN   // ==> RAD_AUTO_AD_Const and RAD_DEBUG
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
#include <complex>      // must come before namespace Sacado
#endif
#endif

namespace Sacado {
namespace Rad {

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
        typedef double  dtype;
        typedef long    ltype;
        typedef int     itype;
        typedef T       ttype;
        };
 template<> class
DoubleAvoid<double> {
 public:
        typedef RAD_DoubleIgnore &dtype;
        typedef long              ltype;
        typedef int               itype;
        typedef RAD_DoubleIgnore &ttype;
        };
template<> class
DoubleAvoid<int> {
 public:
        typedef double            dtype;
        typedef long              ltype;
        typedef RAD_DoubleIgnore &itype;
        typedef RAD_DoubleIgnore &ttype;
        };
template<> class
DoubleAvoid<long> {
 public:
        typedef double            dtype;
        typedef RAD_DoubleIgnore &ltype;
        typedef int               itype;
        typedef RAD_DoubleIgnore &ttype;
        };

#define Dtype typename DoubleAvoid<Double>::dtype
#define Ltype typename DoubleAvoid<Double>::ltype
#define Itype typename DoubleAvoid<Double>::itype
#define Ttype typename DoubleAvoid<Double>::ttype

 template<typename Double> class IndepADvar;
 template<typename Double> class ConstADvar;
 template<typename Double> class ConstADvari;
 template<typename Double> class ADvar;
 template<typename Double> class ADvari;
 template<typename Double> class ADvar1;
 template<typename Double> class ADvar1s;
 template<typename Double> class ADvar2;
 template<typename Double> class ADvar2q;
 template<typename Double> class ADvarn;
 template<typename Double> class Derp;

 template<typename Double> struct
ADmemblock {    // We get memory in ADmemblock chunks and never give it back,
                // but reuse it once computations start anew after call(s) on
                // ADcontext::Gradcomp() or ADcontext::Weighted_Gradcomp().
        ADmemblock *next;
   char memblk[1000*sizeof(double)];
        };

 template<typename Double> class
ADcontext {
        typedef ADmemblock<Double> ADMemblock;
        typedef ADvar<Double> ADVar;
        typedef ADvari<Double> ADVari;
        typedef Derp<Double> DErp;

        ADMemblock *Busy, *First, *Free;
        char *Mbase;
        size_t Mleft, rad_mleft_save;
        int rad_need_reinit;
#if RAD_REINIT > 0
        ADMemblock *DBusy, *DFree;
        size_t DMleft, nderps;
#endif
#ifdef RAD_DEBUG_BLOCKKEEP
        int rad_busy_blocks;
        ADMemblock *rad_Oldcurmb;
#endif
        void *new_ADmemblock(size_t);
        void do_init();
 public:
        static const Double One, negOne;
        inline ADcontext() { do_init(); }
        void *Memalloc(size_t len);
        static void Gradcomp(int wantgrad);
        static inline void Gradcomp() { Gradcomp(1); }
        static void aval_reset();
        static void free_all();
        static void re_init();
        static void zero_out();
        static void Weighted_Gradcomp(size_t, ADVar**, Double*);
        static void Outvar_Gradcomp(ADVar&);
#if RAD_REINIT > 0
        DErp *new_Derp(const Double *, const ADVari*, const ADVari*);
        RAD_REINIT_1(friend class ConstADvari<Double>;)
#endif
        };

#if RAD_REINIT == 0
 template<typename Double> class
CADcontext: public ADcontext<Double> {  // for possibly constant ADvar values
 protected:
        bool fpval_implies_const;
 public:
        friend class ADvar<Double>;
        CADcontext(): ADcontext<Double>() { fpval_implies_const = false; }
        };
#endif

 template<typename Double> class
Derp {          // one derivative-propagation operation
 public:
        friend class ADvarn<Double>;
        typedef ADvari<Double> ADVari;
#if RAD_REINIT > 0
        const Double a;
        inline void *operator new(size_t, Derp *x) { return x; }
        inline void operator delete(void*, Derp *) {}
#else
        static Derp *LastDerp;
        Derp *next;
        const Double *a;
#endif
        const ADVari *b;
        const ADVari *c;
        Derp(){};
        Derp(const ADVari *);
        Derp(const Double *, const ADVari *);
        Derp(const Double *, const ADVari *, const ADVari *);
        inline void *operator new(size_t len) { return (Derp*)ADVari::adc.Memalloc(len); }
        /* c->aval += a * b->aval; */
        };

#if RAD_REINIT > 0 //{
 template<typename Double> Derp<Double>*
ADcontext<Double>::new_Derp(const Double *a, const ADvari<Double> *b, const ADvari<Double> *c)
{
        ADMemblock *x;
        DErp *rv;
        Allow_noderiv(if (!b->wantderiv) return 0;)
        if (this->DMleft == 0) {
                if ((x = DFree))
                        DFree = x->next;
                else
                        x = new ADMemblock;
                x->next = DBusy;
                DBusy = x;
                this->DMleft = nderps;
                }
        rv = &((DErp*)DBusy->memblk)[--this->DMleft];
        new(rv) DErp(a,b,c);
        return rv;
        }
#endif //}

// Now we use #define to overcome bad design in the C++ templating system

#define Ai const Base< ADvari<Double> >&
#define AI const Base< IndepADvar<Double> >&
#define T template<typename Double>
#define D Double
#define T1(f) \
T F f (AI); \
T F f (Ai);
#define T2(r,f) \
 T r f(Ai,Ai); \
 T r f(Ai,D); \
 T r f(Ai,Dtype); \
 T r f(Ai,Ltype); \
 T r f(Ai,Itype); \
 T r f(D,Ai); \
 T r f(Dtype,Ai); \
 T r f(Ltype,Ai); \
 T r f(Itype,Ai); \
 T r f(AI,D); \
 T r f(AI,Dtype); \
 T r f(AI,Ltype); \
 T r f(AI,Itype); \
 T r f(D,AI); \
 T r f(Dtype,AI);   \
 T r f(Ltype,AI); \
 T r f(Itype,AI); \
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
T1(cbrt)

T F copy(AI);
T F copy(Ai);

#undef F
#undef T2
#undef T1
#undef D
#undef T
#undef AI
#undef Ai

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

 template<typename Double> class
ADvari : public Base< ADvari<Double> > {        // implementation of an ADvar
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
        mutable Double Val;     // result of this operation
        mutable Double aval;    // adjoint -- partial of final result w.r.t. this Val
        Allow_noderiv(mutable int wantderiv;)
        void *operator new(size_t len) {
#ifdef RAD_DEBUG
                ADvari *rv = (ADvari*)ADvari::adc.Memalloc(len);
                rv->gcgen = gcgen_cur;
                rv->opno = ++last_opno;
                if (last_opno == zap_opno && gcgen_cur == zap_gcgen)
                        printf("Got to opno %d\n", last_opno);
                return rv;
#else
                return ADvari::adc.Memalloc(len);
#endif
                }
        void operator delete(void*) {} /*Should never be called.*/
#ifdef RAD_ALLOW_WANTDERIV
        inline ADvari(Double t): Val(t), aval(0.), wantderiv(1) {}
        inline ADvari(): Val(0.), aval(0.), wantderiv(1) {}
        inline void no_deriv() { wantderiv = 0; }
#else
        inline ADvari(Double t): Val(t), aval(0.) {}
        inline ADvari(): Val(0.), aval(0.) {}
#endif
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
#endif
#define F friend
#define R ADvari&
#define Ai const Base< ADvari >&
#define D Double
#define T1(r,f) F r f <>(Ai);
#define T2(r,f) \
F r f <>(Ai,Ai); \
F r f <>(Ai,D);       \
F r f <>(D,Ai);
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
        T1(R,cbrt)
        T2(int,operator<)
        T2(int,operator<=)
        T2(int,operator==)
        T2(int,operator!=)
        T2(int,operator>=)
        T2(int,operator>)
#undef D
#undef T2
#undef T1
#undef Ai
#undef R
#undef F

        friend ADvari& ADf1<>(Double f, Double g, const ADvari &x);
        friend ADvari& ADf2<>(Double f, Double gx, Double gy, const ADvari &x, const ADvari &y);
        friend ADvari& ADfn<>(Double f, int n, const IndepADVar *x, const Double *g);

        inline operator Double() { return this->Val; }
        inline operator Double() const { return this->Val; }

        static ADcontext<Double> adc;
        };

 template<typename Double> class
ADvar1: public ADvari<Double> { // simplest unary ops
 public:
        typedef ADvari<Double> ADVari;
        ADvar1(Double val1): ADVari(val1) {}
#if RAD_REINIT > 0
        ADvar1(Double val1, const Double *a1, const ADVari *c1): ADVari(val1) {
#ifdef RAD_ALLOW_WANTDERIV
                if (!ADVari::adc.new_Derp(a1,this,c1))
                        this->no_deriv();
#else
                ADVari::adc.new_Derp(a1,this,c1);
#endif
                }
#else // RAD_REINIT == 0
        Derp<Double> d;
        ADvar1(Double val1, const ADVari *c1): ADVari(val1), d(c1) {}
        ADvar1(Double val1, const Double *a1, const ADVari *c1):
                ADVari(val1), d(a1,this,c1) {}
#ifdef RAD_AUTO_AD_Const
        typedef typename ADVari::IndepADVar IndepADVar;
        typedef ADvar<Double> ADVar;
        ADvar1(const IndepADVar*, const IndepADVar&);
        ADvar1(const IndepADVar*, const ADVari&);
        ADvar1(const Double val1, const Double *a1, const ADVari *c1, const ADVar *v):
                        ADVari(val1), d(a1,this,c1) {
                c1->padv = 0;
                this->ADvari_padv(v);
                }
#endif
#endif // RAD_REININT > 0
        };


 template<typename Double> class
ConstADvari: public ADvari<Double> {
 private:
        RAD_REINIT_0(ConstADvari *prevcad;)
        ConstADvari() {};       // prevent construction without value (?)
        RAD_REINIT_0(static ConstADvari *lastcad;)
        friend class ADcontext<Double>;
 public:
        typedef ADvari<Double> ADVari;
        typedef Derp<Double> DErp;
#if RAD_REINIT == 0
        static CADcontext<Double> cadc;
        inline void *operator new(size_t len) { return ConstADvari::cadc.Memalloc(len); }
        inline ConstADvari(Double t): ADVari(t) { prevcad = lastcad; lastcad = this; }
#else
        inline void *operator new(size_t len) { return ADVari::adc.Memalloc(len); }
        inline ConstADvari(Double t): ADVari(t) {}
#endif
        static void aval_reset(void);
        };

 template<typename Double> class
IndepADvar_base0 {
 public:
#if RAD_REINIT == 1
        IndepADvar_base0 *ADvnext, *ADvprev;
        inline IndepADvar_base0(double) { ADvnext = ADvprev = this; }
#elif RAD_REINIT == 2
        typedef unsigned long ADGenType;
        mutable ADGenType gen;
        inline IndepADvar_base0(double) { gen = 1; }
#endif
        inline IndepADvar_base0() {}
        };

 template<typename Double> class
IndepADvar_base: public IndepADvar_base0<Double> {
 public:
#if RAD_REINIT > 0
        Allow_noderiv(mutable int wantderiv;)
        static IndepADvar_base0<Double> IndepADvar_root;
        RAD_REINIT_2(static typename IndepADvar_base0<Double>::ADGenType gen0;)
#endif
        inline IndepADvar_base(Allow_noderiv(int wd)) {
#if RAD_REINIT == 1
                IndepADvar_root.ADvprev = (this->ADvprev = IndepADvar_root.ADvprev)->ADvnext = this;
                this->ADvnext = &IndepADvar_root;
////#elif RAD_REINIT == 2
////            **** Do NOT assign gen here -- it's too soon.
////            **** Must wait til after cv is assigned, as gen may be incremented
////            **** in the course of allocating cv.
////            //// this->gen = IndepADvar_root.gen;
#endif
                Allow_noderiv(this->wantderiv = wd;)
                }
        inline ~IndepADvar_base() {
#if RAD_REINIT == 1
                (this->ADvnext->ADvprev = this->ADvprev)->ADvnext = this->ADvnext;
#endif
                }
#if RAD_REINIT > 0
 private:
        inline IndepADvar_base(double d): IndepADvar_base0<Double>(d) {}
#endif
#if RAD_REINIT == 1
 protected:
        void Move_to_end();
#endif
        };

#if RAD_REINIT > 0 //{
template<typename Double> IndepADvar_base0<Double>
        IndepADvar_base<Double>::IndepADvar_root(0.);
#if RAD_REINIT == 1
 template<typename Double> void
IndepADvar_base<Double>::Move_to_end() {
        if (this != IndepADvar_root.ADvprev) {
                (this->ADvnext->ADvprev = this->ADvprev)->ADvnext = this->ADvnext;
                IndepADvar_root.ADvprev = (this->ADvprev = IndepADvar_root.ADvprev)->ADvnext = this;
                this->ADvnext = &IndepADvar_root;
                }
        }
#elif RAD_REINIT == 2
template<typename Double> typename IndepADvar_base0<Double>::ADGenType
        IndepADvar_base<Double>::gen0(1);
#endif
#endif //}

 template<typename Double> class
IndepADvar: protected IndepADvar_base<Double>, public Base< IndepADvar<Double> > {         // an independent ADvar
 protected:
        static void AD_Const(const IndepADvar&);
        mutable ADvari<Double> *cv;
 private:
        IndepADvar& operator=(IndepADvar&x) {
                /* private to prevent assignment */
                RAD_cvchk(x)
#ifdef RAD_AUTO_AD_Const
                if (cv)
                        cv->padv = 0;
                cv = new ADvar1<Double>(this,x);
                return *this;
#elif defined(RAD_EQ_ALIAS)
                this->cv = x.cv;
                return *this;
#else
                return ADvar_operatoreq(this,*x.cv);
#endif //RAD_AUTO_AD_Const
                }
 public:
        int Wantderiv(int);
        RAD_REINIT_2(Double Val;)
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
#ifdef RAD_ALLOW_WANTDERIV
        inline int Wantderiv() { return this->wantderiv; }
 protected:
        inline IndepADvar(void*, int wd): IndepADvar_base<Double>(wd){RAD_REINIT_1(cv = 0;)}
        IndepADvar(Ttype, int);
        IndepADvar(Double, int);
 public:
#else
        inline int Wantderiv() { return 1; }
#endif
#ifdef RAD_AUTO_AD_Const
        inline IndepADvar(const IndepADvar &x) { cv = x.cv ? new ADvar1<Double>(this, x) : 0; };
        inline IndepADvar() { cv = 0; }
        inline ~IndepADvar() {
                                        if (cv)
                                                cv->padv = 0;
                                        }
#else
        inline IndepADvar() Allow_noderiv(: IndepADvar_base<Double>(1)) {
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
        inline operator ADVari&() const { RAD_cvchk((*this)) return *this->cv; }
        inline operator ADVari*() const { RAD_cvchk((*this)) return this->cv; }
#endif //RAD_Const_WARN

        inline Double val() const {
#if RAD_REINIT == 2
                return Val;
#else
                return cv->Val;
#endif
                }
        inline Double adj() const {
                return
                        RAD_REINIT_2(this->gen != this->gen0 ? 0. :)
                        cv->aval;
        }

        friend void AD_Const1<>(Double*, const IndepADvar&);

        friend ADVari& ADf1<>(Double, Double, const IndepADvar&);
        friend ADVari& ADf2<>(Double, Double, Double, const IndepADvar&, const IndepADvar&);
        friend ADVari& ADf2<>(Double, Double, Double, const ADVari&, const IndepADvar&);
        friend ADVari& ADf2<>(Double, Double, Double, const IndepADvar&, const ADVari&);

        static inline void Gradcomp(int wantgrad)
                                { ADcontext<Double>::Gradcomp(wantgrad); }
        static inline void Gradcomp()
                                { ADcontext<Double>::Gradcomp(1); }
        static inline void aval_reset() { ConstADvari<Double>::aval_reset(); }
        static inline void Weighted_Gradcomp(size_t n, ADVar **v, Double *w)
                                { ADcontext<Double>::Weighted_Gradcomp(n, v, w); }
        static inline void Outvar_Gradcomp(ADVar &v)
                                { ADcontext<Double>::Outvar_Gradcomp(v); }

        /* We use #define to deal with bizarre templating rules that apparently */
        /* require us to spell the some conversion explicitly */


#define Ai const Base< ADVari >&
#define AI const Base< IndepADvar >&
#define D Double
#define T2(r,f) \
 r f <>(AI,AI);\
 r f <>(Ai,AI);\
 r f <>(AI,Ai);\
 r f <>(D,AI);\
 r f <>(AI,D);
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
T1(cbrt)

#undef D
#undef F
#undef T1
#undef T2
#undef AI
#undef Ai

        };

 template<typename Double> class
ADvar: public IndepADvar<Double> {      // an "active" variable
 public:
        //! Turn ADvar into a meta-function class usable with mpl::apply
        template <typename U> struct apply { typedef ADvar<U> type; };

        typedef IndepADvar<Double> IndepADVar;
        typedef typename IndepADVar::ADVari ADVari;
        typedef ConstADvari<Double> ConstADVari;
 private:
        void ADvar_ctr(Double d) {
#if RAD_REINIT == 0 //{
                ADVari *x;
                if (ConstADVari::cadc.fpval_implies_const)
                        x = new ConstADVari(d);
                else {
#ifdef RAD_AUTO_AD_Const //{
                        x = new ADVari((IndepADVar*)this, d);
                        x->ADvari_padv(this);
#else
                        x = new ADVari(d);
#endif //}
                        }
#else
                RAD_REINIT_1(this->cv = 0;)
                ADVari *x = new ADVari(d);
                RAD_REINIT_2(this->Val = d;)
#endif //}
                this->cv = x;
                RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                }
 public:
        friend class ADvar1<Double>;
        typedef ADvar1<Double> ADVar1;
        inline ADvar() { /* cv = 0; */ }
        inline ADvar(Ttype d)  { ADvar_ctr(d); }
        inline ADvar(double i) { ADvar_ctr(Double(i)); }
        inline ADvar(int i)     { ADvar_ctr(Double(i)); }
        inline ADvar(long i)    { ADvar_ctr(Double(i)); }
        inline ~ADvar() {}
        Allow_noderiv(inline ADvar(void *v, int wd): IndepADVar(v, wd) {})
#ifdef RAD_AUTO_AD_Const
        inline ADvar(IndepADVar &x) {
                this->cv = x.cv ? new ADVar1(this, x) : 0;
                RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                }
        inline ADvar(ADVari &x) { this->cv = &x; x.ADvari_padv(this); }
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
        inline ADvar(const IndepADVar &x) {
                RAD_cvchk(x)
                this->cv = (ADVari*)x.cv;
                }
        inline ADvar(const ADVari &x) { this->cv = (ADVari*)&x; }
        inline ADvar& operator=(IndepADVar &x) {
                RAD_cvchk(x)
                this->cv = (ADVari*)x.cv; return *this;
                }
        inline ADvar& operator=(const ADVari &x) { this->cv = (ADVari*)&x; return *this; }
#else /*!RAD_EQ_ALIAS*/
        ADvar(const IndepADVar &x) {
                if (x.cv) {
                        RAD_cvchk(x)
                        RAD_REINIT_1(this->cv = 0;)
                        this->cv = new ADVar1(x.cv->Val, &this->cv->adc.One, x.cv);
                        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                        }
                else
                        this->cv = 0;
                }
        ADvar(const ADvar&x) {
                if (x.cv) {
                        RAD_cvchk(x)
                        RAD_REINIT_1(this->cv = 0;)
                        this->cv = new ADVar1(x.cv->Val, &this->cv->adc.One, (ADVari*)x.cv);
                        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                        }
                else
                        this->cv = 0;
                }
        ADvar(const  ADVari &x) {
                RAD_REINIT_1(this->cv = 0;)
                this->cv = new ADVar1(x.Val, &this->cv->adc.One, &x);
                RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                }
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
#if RAD_REINIT == 0
        inline static bool get_fpval_implies_const(void)
                { return ConstADVari::cadc.fpval_implies_const; }
        inline static void set_fpval_implies_const(bool newval)
                { ConstADVari::cadc.fpval_implies_const = newval; }
        inline static bool setget_fpval_implies_const(bool newval) {
                bool oldval = ConstADVari::cadc.fpval_implies_const;
                ConstADVari::cadc.fpval_implies_const = newval;
                return oldval;
                }
#else
        inline static bool get_fpval_implies_const(void) { return true; }
        inline static void set_fpval_implies_const(bool newval) {}
        inline static bool setget_fpval_implies_const(bool newval) { return newval; }
#endif
        static inline void Gradcomp(int wantgrad)
                                { ADcontext<Double>::Gradcomp(wantgrad); }
        static inline void Gradcomp()
                                { ADcontext<Double>::Gradcomp(1); }
        static inline void aval_reset() { ConstADVari::aval_reset(); }
        static inline void Weighted_Gradcomp(size_t n, ADvar **v, Double *w)
                                { ADcontext<Double>::Weighted_Gradcomp(n, v, w); }
        static inline void Outvar_Gradcomp(ADvar &v)
                                { ADcontext<Double>::Outvar_Gradcomp(v); }
        };

#if RAD_REINIT == 0
template<typename Double>
 inline void AD_Const1(Double *notused, const IndepADvar<Double>&v)
{ IndepADvar<Double>::AD_Const(v); }

template<typename Double>
 inline void AD_Const(const IndepADvar<Double>&v) { AD_Const1((Double*)0, v); }
#else
template<typename Double>
 inline void AD_Const(const IndepADvar<Double>&v) {}
#endif //RAD_REINIT == 0

 template<typename Double> class
ConstADvar: public ADvar<Double> {
 public:
        typedef ADvar<Double> ADVar;
        typedef typename ADVar::ADVar1 ADVar1;
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
        ConstADvar(Ttype d)     { ConstADvar_ctr(d); }
        ConstADvar(double i)    { ConstADvar_ctr(Double(i)); }
        ConstADvar(int i)       { ConstADvar_ctr(Double(i)); }
        ConstADvar(long i)      { ConstADvar_ctr(Double(i)); }
        ConstADvar(const IndepADVar&);
        ConstADvar(const ConstADvar&);
        ConstADvar(const ADVari&);
        inline ~ConstADvar() {}
#ifdef RAD_NO_CONST_UPDATE
 private:
#endif
        ConstADvar();
        inline ConstADvar& operator=(Double d) {
#if RAD_REINIT > 0
                this->cv = new ADVari(d);
                RAD_REINIT_2(this->Val = d;)
#else
                this->cv->Val = d;
#endif
                return *this;
                }
        inline ConstADvar& operator=(ADVari& d) {
#if RAD_REINIT > 0
                this->cv = new ADVar1(this,d);
                RAD_REINIT_2(this->Val = d;)
#else
                this->cv->Val = d.Val;
#endif
                return *this;
                }
 };

#ifdef RAD_ALLOW_WANTDERIV
 template<typename Double> class
ADvar_nd: public ADvar<Double>
{
 public:
        typedef ADvar<Double> ADVar;
        typedef IndepADvar<Double> IndepADVar;
        typedef typename IndepADVar::ADVari ADVari;
        typedef ADvar1<Double> ADVar1;
 private:
        void ADvar_ndctr(Double d) {
                ADVari *x = new ADVari(d);
                this->cv = x;
                RAD_REINIT_2(this->Val = d;)
                RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                }
 public:
        inline ADvar_nd():         ADVar((void*)0,0) {}
        inline ADvar_nd(Ttype d):  ADVar((void*)0,0) { ADvar_ndctr(d); }
        inline ADvar_nd(double i): ADVar((void*)0,0) { ADvar_ndctr(Double(i)); }
        inline ADvar_nd(int i):    ADVar((void*)0,0) { ADvar_ndctr(Double(i)); }
        inline ADvar_nd(long i):   ADVar((void*)0,0) { ADvar_ndctr(Double(i)); }
        inline ADvar_nd(Double d):  ADVar((void*)0,0) { ADvar_ndctr(d); }
        inline ~ADvar_nd() {}
        ADvar_nd(const IndepADVar &x):   ADVar((void*)0,0) {
                if (x.cv) {
                        this->cv = new ADVar1(x.cv->Val, &this->cv->adc.One, x.cv);
                        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                        }
                else
                        this->cv = 0;
                }
        ADvar_nd(const ADVari &x):   ADVar((void*)0,0) {
                this->cv = new ADVar1(x.Val, &this->cv->adc.One, &x);
                RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
                }
        inline ADvar_nd& operator=(const ADVari &x) { return (ADvar_nd&)ADvar_operatoreq(this,x); };
        };
#else
#define ADvar_nd ADvar
#endif //RAD_ALLOW_WANTDERIV

 template<typename Double> class
ADvar1s: public ADvar1<Double> { // unary ops with partial "a"
 public:
        typedef ADvar1<Double> ADVar1;
        typedef typename ADVar1::ADVari ADVari;
#if RAD_REINIT > 0 //{{
        inline ADvar1s(Double val1, Double a1, const ADVari *c1): ADVar1(val1,&a1,c1) {}
#else //}{
        Double a;
        ADvar1s(Double val1, Double a1, const ADVari *c1): ADVar1(val1,&a,c1), a(a1) {}
#ifdef RAD_AUTO_AD_Const
        typedef typename ADVar1::ADVar ADVar;
        ADvar1s(Double val1, Double a1, const ADVari *c1, const ADVar *v):
                ADVar1(val1,&a,c1,v), a(a1) {}
#endif
#endif // }} RAD_REINIT > 0
        };

 template<typename Double> class
ADvar2: public ADvari<Double> { // basic binary ops
 public:
        typedef ADvari<Double> ADVari;
        typedef Derp<Double> DErp;
        ADvar2(Double val1): ADVari(val1) {}
#if RAD_REINIT > 0 //{{
        ADvar2(Double val1, const ADVari *Lcv, const Double *Lc,
                            const ADVari *Rcv, const Double *Rc): ADVari(val1) {
#ifdef RAD_ALLOW_WANTDERIV
                DErp *Ld, *Rd;
                Ld = ADVari::adc.new_Derp(Lc, this, Lcv);
                Rd = ADVari::adc.new_Derp(Rc, this, Rcv);
                if (!Ld && !Rd)
                        this->no_deriv();
#else //!RAD_ALLOW_WANTDERIV
                ADVari::adc.new_Derp(Lc, this, Lcv);
                ADVari::adc.new_Derp(Rc, this, Rcv);
#endif //RAD_ALLOW_WANTDERIV
                }
#else //}{ RAD_REINIT == 0
        DErp dL, dR;
        ADvar2(Double val1, const ADVari *Lcv, const Double *Lc,
                const ADVari *Rcv, const Double *Rc):
                        ADVari(val1) {
                dR.next = DErp::LastDerp;
                dL.next = &dR;
                DErp::LastDerp = &dL;
                dL.a = Lc;
                dL.c = Lcv;
                dR.a = Rc;
                dR.c = Rcv;
                dL.b = dR.b = this;
                }
#ifdef RAD_AUTO_AD_Const
        typedef ADvar<Double> ADVar;
        ADvar2(Double val1, const ADVari *Lcv, const Double *Lc,
                const ADVari *Rcv, const Double *Rc, ADVar *v):
                        ADVari(val1) {
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
#endif // }} RAD_REINIT
        };

 template<typename Double> class
ADvar2q: public ADvar2<Double> { // binary ops with partials "a", "b"
 public:
        typedef ADvar2<Double> ADVar2;
        typedef typename ADVar2::ADVari ADVari;
        typedef typename ADVar2::DErp DErp;
#if RAD_REINIT > 0 //{{
        inline ADvar2q(Double val1, Double Lp, Double Rp, const ADVari *Lcv, const ADVari *Rcv):
                        ADVar2(val1) {
#ifdef RAD_ALLOW_WANTDERIV
                DErp *Ld, *Rd;
                Ld = ADVari::adc.new_Derp(&Lp, this, Lcv);
                Rd = ADVari::adc.new_Derp(&Rp, this, Rcv);
                if (!Ld && !Rd)
                        this->no_deriv();
#else //!RAD_ALLOW_WANTDERIV
                ADVari::adc.new_Derp(&Lp, this, Lcv);
                ADVari::adc.new_Derp(&Rp, this, Rcv);
#endif //RAD_ALLOW_WANTDERIV
                }
#else //}{ RADINIT == 0
        Double a, b;
        ADvar2q(Double val1, Double Lp, Double Rp, const ADVari *Lcv, const ADVari *Rcv):
                        ADVar2(val1), a(Lp), b(Rp) {
                this->dR.next = DErp::LastDerp;
                this->dL.next = &this->dR;
                DErp::LastDerp = &this->dL;
                this->dL.a = &a;
                this->dL.c = Lcv;
                this->dR.a = &b;
                this->dR.c = Rcv;
                this->dL.b = this->dR.b = this;
                }
#ifdef RAD_AUTO_AD_Const
        typedef typename ADVar2::ADVar ADVar;
        ADvar2q(Double val1, Double Lp, Double Rp, const ADVari *Lcv,
                const ADVari *Rcv, const ADVar *v):
                        ADVar2(val1), a(Lp), b(Rp) {
                this->dR.next = DErp::LastDerp;
                this->dL.next = &this->dR;
                DErp::LastDerp = &this->dL;
                this->dL.a = &a;
                this->dL.c = Lcv;
                this->dR.a = &b;
                this->dR.c = Rcv;
                this->dL.b = this->dR.b = this;
                Lcv->padv = 0;
                this->ADvari_padv(v);
                }
#endif
#endif //}} RAD_REINIT > 0
        };

 template<typename Double> class
ADvarn: public ADvari<Double> { // n-ary ops with partials "a"
 public:
        typedef ADvari<Double> ADVari;
        typedef typename ADVari::IndepADVar IndepADVar;
        typedef Derp<Double> DErp;
#if RAD_REINIT > 0 // {{
        ADvarn(Double val1, int n1, const IndepADVar *x, const Double *g): ADVari(val1) {
#ifdef RAD_ALLOW_WANTDERIV
                int i, nd;
                for(i = nd = 0; i < n1; i++)
                        if (ADVari::adc.new_Derp(g+i, this, x[i].cv))
                                ++nd;
                if (!nd)
                        this->no_deriv();
#else
                for(int i = 0; i < n1; i++)
                        ADVari::adc.new_Derp(g+i, this, x[i].cv);
#endif // RAD_ALLOW_WANTDERIV
                }
#else //}{
        int n;
        Double *a;
        DErp *Da;
        ADvarn(Double val1, int n1, const IndepADVar *x, const Double *g):
                        ADVari(val1), n(n1) {
                DErp *d1, *dlast;
                Double *a1;
                int i;

                a1 = a = (Double*)ADVari::adc.Memalloc(n*sizeof(*a));
                d1 = Da =  (DErp*)ADVari::adc.Memalloc(n*sizeof(DErp));
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
                }
#endif //}} RAD_REINIT > 0
        };

template<typename Double>
inline ADvari<Double>& operator+(const Base< ADvari<Double> > &TT) {
  const ADvari<Double>& T = TT.derived();
  return *(ADvari<Double>*)&T; }

template<typename Double>
 inline int operator<(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& L = LL.derived();
  const ADvari<Double>& R = RR.derived();
  return L.Val < R.Val;
}
template<typename Double>
inline int operator<(const Base< ADvari<Double> > &LL, Double R) {
  const ADvari<Double>& L = LL.derived();
  return L.Val < R;
}
template<typename Double>
 inline int operator<(Double L, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& R = RR.derived();
  return L < R.Val;
}

template<typename Double>
 inline int operator<=(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& L = LL.derived();
  const ADvari<Double>& R = RR.derived();
  return L.Val <= R.Val;
}
template<typename Double>
 inline int operator<=(const Base< ADvari<Double> > &LL, Double R) {
  const ADvari<Double>& L = LL.derived();
  return L.Val <= R;
}
template<typename Double>
 inline int operator<=(Double L, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& R = RR.derived();
  return L <= R.Val;
}

template<typename Double>
 inline int operator==(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& L = LL.derived();
  const ADvari<Double>& R = RR.derived();
  return L.Val == R.Val;
}
template<typename Double>
 inline int operator==(const Base< ADvari<Double> > &LL, Double R) {
  const ADvari<Double>& L = LL.derived();
  return L.Val == R;
}
template<typename Double>
 inline int operator==(Double L, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& R = RR.derived();
  return L == R.Val;
}

template<typename Double>
 inline int operator!=(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& L = LL.derived();
  const ADvari<Double>& R = RR.derived();
  return L.Val != R.Val;
}
template<typename Double>
 inline int operator!=(const Base< ADvari<Double> > &LL, Double R) {
  const ADvari<Double>& L = LL.derived();
  return L.Val != R;
}
template<typename Double>
 inline int operator!=(Double L, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& R = RR.derived();
  return L != R.Val;
}

template<typename Double>
 inline int operator>=(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& L = LL.derived();
  const ADvari<Double>& R = RR.derived();
  return L.Val >= R.Val;
}
template<typename Double>
 inline int operator>=(const Base< ADvari<Double> > &LL, Double R) {
  const ADvari<Double>& L = LL.derived();
  return L.Val >= R;
}
template<typename Double>
 inline int operator>=(Double L, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& R = RR.derived();
  return L >= R.Val;
}

template<typename Double>
 inline int operator>(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& L = LL.derived();
  const ADvari<Double>& R = RR.derived();
  return L.Val > R.Val;
}
template<typename Double>
 inline int operator>(const Base< ADvari<Double> > &LL, Double R) {
  const ADvari<Double>& L = LL.derived();
  return L.Val > R;
}
template<typename Double>
 inline int operator>(Double L, const Base< ADvari<Double> > &RR) {
  const ADvari<Double>& R = RR.derived();
  return L > R.Val;
}

template<typename Double>
 inline void *ADcontext<Double>::Memalloc(size_t len) {
                if (Mleft >= len)
                        return Mbase + (Mleft -= len);
                return new_ADmemblock(len);
                }
#if RAD_REINIT > 0 //{{

template<typename Double>
 inline Derp<Double>::Derp(const ADVari *c1): c(c1) {}

template<typename Double>
 inline Derp<Double>::Derp(const Double *a1, const ADVari *c1): a(*a1), c(c1) {}


template<typename Double>
 inline Derp<Double>::Derp(const Double *a1, const ADVari *b1, const ADVari *c1):
        a(*a1), b(b1), c(c1) {}
#else //}{ RAD_REINIT == 0

template<typename Double>
 inline Derp<Double>::Derp(const ADVari *c1): c(c1) {
                next = LastDerp;
                LastDerp = this;
                }

template<typename Double>
 inline Derp<Double>::Derp(const Double *a1, const ADVari *c1): a(a1), c(c1) {
                next = LastDerp;
                LastDerp = this;
                }


template<typename Double>
 inline Derp<Double>::Derp(const Double *a1, const ADVari *b1, const ADVari *c1): a(a1), b(b1), c(c1) {
                next = LastDerp;
                LastDerp = this;
                }
#endif //}} RAD_REINIT

/**** radops ****/

#if RAD_REINIT > 0
#else
template<typename Double> Derp<Double> *Derp<Double>::LastDerp = 0;
#endif //RAD_REINIT
template<typename Double> ADcontext<Double> ADvari<Double>::adc;
template<typename Double> const Double ADcontext<Double>::One = 1.;
template<typename Double> const Double ADcontext<Double>::negOne = -1.;
RAD_REINIT_0(template<typename Double> CADcontext<Double> ConstADvari<Double>::cadc;)
RAD_REINIT_0(template<typename Double> ConstADvari<Double> *ConstADvari<Double>::lastcad;)

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

template<typename Double> void ADcontext<Double>::do_init()
{
        First = new ADMemblock;
        First->next = 0;
        Busy = First;
        Free = 0;
        Mbase = (char*)First->memblk;
        Mleft = sizeof(First->memblk);
        rad_need_reinit = 0;
#ifdef RAD_DEBUG_BLOCKKEEP
        rad_busy_blocks = 0;
        rad_mleft_save = 0;
        rad_Oldcurmb = 0;
#endif
#if RAD_REINIT > 0
        DBusy = new ADMemblock;
        DBusy->next = 0;
        DFree = 0;
        DMleft = nderps = sizeof(DBusy->memblk)/sizeof(DErp);
#endif
        }

template<typename Double> void ADcontext<Double>::free_all()
{
        typedef ConstADvari<Double> ConstADVari;
        ADMemblock *mb, *mb1;

        for(mb = ADVari::adc.Busy; mb; mb = mb1) {
                mb1 = mb->next;
                delete mb;
                }
        for(mb = ADVari::adc.Free; mb; mb = mb1) {
                mb1 = mb->next;
                delete mb;
                }
        for(mb = ConstADVari::cadc.Busy; mb; mb = mb1) {
                mb1 = mb->next;
                delete mb;
                }
        ConstADVari::cadc.Busy  = ADVari::adc.Busy  = ADVari::adc.Free = 0;
        ConstADVari::cadc.Mleft = ADVari::adc.Mleft = 0;
        ConstADVari::cadc.Mbase = ADVari::adc.Mbase = 0;
#if RAD_REINIT > 0
        for(mb = ADVari::adc.DBusy; mb; mb = mb1) {
                mb1 = mb->next;
                delete mb;
                }
        for(mb = ADVari::adc.DFree; mb; mb = mb1) {
                mb1 = mb->next;
                delete mb;
                }
        ADVari::adc.DBusy  = ADVari::adc.DFree = 0;
        ADVari::adc.DMleft = 0;
        ConstADVari::cadc.Mbase = ADVari::adc.Mbase = 0;
#else
        ConstADVari::lastcad = 0;
        Derp<Double>::LastDerp = 0;
#endif
        }

template<typename Double> void ADcontext<Double>::re_init()
{
        typedef ConstADvari<Double> ConstADVari;

        if (ConstADVari::cadc.Busy || ADVari::adc.Busy || ADVari::adc.Free
#if RAD_REINIT > 0
                || ADVari::adc.DBusy || ADVari::adc.DFree
#endif
                ) free_all();
        ADVari::adc.do_init();
        ConstADVari::cadc.do_init();
        }

template<typename Double> void*
ADcontext<Double>::new_ADmemblock(size_t len)
{
        ADMemblock *mb, *mb0, *mb1, *mbf, *x;
#ifdef RAD_AUTO_AD_Const
        ADVari *a, *anext;
        IndepADvar<Double> *v;
#ifdef RAD_Const_WARN
        ADVari *cv;
        int i, j;
#endif
#endif /*RAD_AUTO_AD_Const*/
#if RAD_REINIT == 1
        typedef IndepADvar_base0<Double> ADvb;
        typedef IndepADvar<Double> IADv;
        ADVari *tcv;
        ADvb *vb, *vb0;
#endif

        if ((rad_need_reinit & 1) && this == &ADVari::adc) {
                rad_need_reinit &= ~1;
                RAD_REINIT_0(DErp::LastDerp = 0;)
#ifdef RAD_DEBUG_BLOCKKEEP
                Mleft = rad_mleft_save;
                if (Mleft < sizeof(First->memblk))
                        _uninit_f2c(Mbase + Mleft,
                                UninitType<Double>::utype,
                                (sizeof(First->memblk) - Mleft)
                                /sizeof(typename Sacado::ValueType<Double>::type));
                if ((mb = Busy->next)) {
                        mb0 = rad_Oldcurmb;
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
                        mb0 = 0;
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
#endif /*RAD_DEBUG_BLOCKKEEP*/
#ifdef RAD_AUTO_AD_Const // {
                *ADVari::Last_ADvari = 0;
                ADVari::Last_ADvari = &ADVari::First_ADvari;
                a = ADVari::First_ADvari;
                if (a) {
                        do {
                                anext = a->Next;
                                if ((v = (IndepADvar<Double> *)a->padv)) {
#ifdef RAD_Const_WARN
                                        if ((i = a->opno) > 0)
                                                i = -i;
                                        j = a->gcgen;
                                        v->cv = cv = new ADVari(v, a->Val);
                                        cv->opno = i;
                                        cv->gcgen = j;
#else
                                        v->cv = new ADVari(v, a->Val);
#endif
                                        }
                                }
                                while((a = anext));
                        }
#endif // } RAD_AUTO_AD_Const
#if RAD_REINIT > 0 //{
                mb = mb0 = DBusy;
                while((mb1 = mb->next)) {
                        mb->next = mb0;
                        mb0 = mb;
                        mb = mb1;
                        }
                DBusy = mb;
                DFree = mb->next;
                mb->next = 0;
                DMleft = nderps;
#if RAD_REINIT == 1
                vb = vb0 = &IndepADvar_base<Double>::IndepADvar_root;
                while((vb = vb->ADvnext) != vb0)
                        if ((tcv = ((IADv*)vb)->cv))
                                ((IADv*)vb)->cv = new ADVari(tcv->Val);
#elif RAD_REINIT == 2
                ++IndepADvar<Double>::gen0;
#endif
#endif //}
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
ADcontext<Double>::Gradcomp(int wantgrad)
{
#if RAD_REINIT > 0 //{{
        ADMemblock *mb;
        DErp *d, *de;

        if (ADVari::adc.rad_need_reinit && wantgrad) {
                mb = ADVari::adc.DBusy;
                d = ((DErp*)mb->memblk) + ADVari::adc.DMleft;
                de = ((DErp*)mb->memblk) + ADVari::adc.nderps;
                for(;;) {
                        for(; d < de; d++)
                                d->c->aval = 0;
                        if (!(mb = mb->next))
                                break;
                        d = (DErp*)mb->memblk;
                        de = d + ADVari::adc.nderps;
                        }
                }
#else //}{ RAD_REINIT == 0
        DErp *d;

        if (ADVari::adc.rad_need_reinit && wantgrad) {
                for(d = DErp::LastDerp; d; d = d->next)
                        d->c->aval = 0;
                }
#endif //}} RAD_REINIT

        if (!(ADVari::adc.rad_need_reinit & 1)) {
                ADVari::adc.rad_need_reinit = 1;
                ADVari::adc.rad_mleft_save = ADVari::adc.Mleft;
                ADVari::adc.Mleft = 0;
                RAD_REINIT_2(++IndepADvar_base<Double>::IndepADvar_root.gen;)
                }
#ifdef RAD_DEBUG
        if (ADVari::gcgen_cur == ADVari::zap_gcgen1 && wantgrad) {
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
#if RAD_REINIT > 0 //{{
        if (ADVari::adc.DMleft < ADVari::adc.nderps && wantgrad) {
                mb = ADVari::adc.DBusy;
                d = ((DErp*)mb->memblk) + ADVari::adc.DMleft;
                de = ((DErp*)mb->memblk) + ADVari::adc.nderps;
                d->b->aval = 1;
                for(;;) {
#ifdef RAD_DEBUG
                        if (ADVari::debug_file) {
                            for(; d < de; d++) {
                                fprintf(ADVari::debug_file, "%d\t%d\t%g + %g * %g",
                                        d->c->opno, d->b->opno, d->c->aval, d->a, d->b->aval);
                                d->c->aval += d->a * d->b->aval;
                                fprintf(ADVari::debug_file, " = %g\n", d->c->aval);
                                }
                            }
                        else
#endif
                            for(; d < de; d++)
                                d->c->aval += d->a * d->b->aval;
                        if (!(mb = mb->next))
                                break;
                        d = (DErp*)mb->memblk;
                        de = d + ADVari::adc.nderps;
                    }
                }
#else //}{ RAD_REINIT == 0
        if ((d = DErp::LastDerp) && wantgrad) {
                d->b->aval = 1;
#ifdef RAD_DEBUG
                if (ADVari::debug_file)
                        do {
                                fprintf(ADVari::debug_file, "%d\t%d\t%g + %g * %g",
                                        d->c->opno, d->b->opno, d->c->aval, *d->a, d->b->aval);
                                d->c->aval += *d->a * d->b->aval;
                                fprintf(ADVari::debug_file, " = %g\n", d->c->aval);
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
#endif // }} RAD_REINIT
#ifdef RAD_DEBUG
        if (ADVari::debug_file)
                fflush(ADVari::debug_file);
        ADVari::gcgen_cur++;
        ADVari::last_opno = 0;
#endif
        }

 template<typename Double> void
ADcontext<Double>::Weighted_Gradcomp(size_t n, ADVar **V, Double *w)
{
        size_t i;
#if RAD_REINIT > 0 //{{
        ADMemblock *mb;
        DErp *d, *de;

        if (ADVari::adc.rad_need_reinit) {
                mb = ADVari::adc.DBusy;
                d = ((DErp*)mb->memblk) + ADVari::adc.DMleft;
                de = ((DErp*)mb->memblk) + ADVari::adc.nderps;
                for(;;) {
                        for(; d < de; d++)
                                d->c->aval = 0;
                        if (!(mb = mb->next))
                                break;
                        d = (DErp*)mb->memblk;
                        de = d + ADVari::adc.nderps;
                        }
                }
#else //}{ RAD_REINIT == 0
        DErp *d;

        if (ADVari::adc.rad_need_reinit) {
                for(d = DErp::LastDerp; d; d = d->next)
                        d->c->aval = 0;
                }
#endif //}} RAD_REINIT

        if (!(ADVari::adc.rad_need_reinit & 1)) {
                ADVari::adc.rad_need_reinit = 1;
                ADVari::adc.rad_mleft_save = ADVari::adc.Mleft;
                ADVari::adc.Mleft = 0;
                RAD_REINIT_2(++IndepADvar_base<Double>::IndepADvar_root.gen;)
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
#if RAD_REINIT > 0 //{{
        if (ADVari::adc.DMleft < ADVari::adc.nderps) {
                for(i = 0; i < n; i++)
                        V[i]->cv->aval = w[i];
                mb = ADVari::adc.DBusy;
                d = ((DErp*)mb->memblk) + ADVari::adc.DMleft;
                de = ((DErp*)mb->memblk) + ADVari::adc.nderps;
                d->b->aval = 1;
                for(;;) {
#ifdef RAD_DEBUG
                    if (ADVari::debug_file) {
                        for(; d < de; d++) {
                                fprintf(ADVari::debug_file, "%d\t%d\t%g + %g * %g",
                                        d->c->opno, d->b->opno, d->c->aval, d->a, d->b->aval);
                                d->c->aval += d->a * d->b->aval;
                                fprintf(ADVari::debug_file, " = %g\n", d->c->aval);
                                }
                        }
                    else
#endif
                        for(; d < de; d++)
                                d->c->aval += d->a * d->b->aval;
                        if (!(mb = mb->next))
                                break;
                        d = (DErp*)mb->memblk;
                        de = d + ADVari::adc.nderps;
                        }
                }
#else //}{ RAD_REINIT == 0
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
#endif // }} RAD_REINIT
#ifdef RAD_DEBUG
        if (ADVari::debug_file)
                fflush(ADVari::debug_file);
        ADVari::gcgen_cur++;
        ADVari::last_opno = 0;
#endif
        }

 template<typename Double> void
ADcontext<Double>::Outvar_Gradcomp(ADVar &V)
{
        Double w = 1;
        ADVar *v = &V;
        ADcontext<Double>::Weighted_Gradcomp(1, &v, &w);
        }

 template<typename Double> void
ADcontext<Double>::zero_out(void)
{
  for(DErp *d = DErp::LastDerp; d; d = d->next) {
    if (d->a)
      *(const_cast<Double*>(d->a)) = Double(0.0);
    if (d->b) {
      d->b->aval = Double(0.0);
      d->b->Val = Double(0.0);
    }
    if (d->c) {
      d->c->aval = Double(0.0);
      d->c->Val = Double(0.0);
    }
  }
}

 template<typename Double>
IndepADvar<Double>::IndepADvar(Ttype d) Allow_noderiv(: IndepADvar_base<Double>(1))
{
        RAD_REINIT_1(cv = 0;)
        cv = new ADVari(d);
        RAD_REINIT_2(Val = d; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
IndepADvar<Double>::IndepADvar(double i) Allow_noderiv(: IndepADvar_base<Double>(1))
{
        RAD_REINIT_1(cv = 0;)
        cv = new ADVari(Double(i));
        RAD_REINIT_2(Val = i; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
IndepADvar<Double>::IndepADvar(int i) Allow_noderiv(: IndepADvar_base<Double>(1))
{
        RAD_REINIT_1(cv = 0;)
        cv = new ADVari(Double(i));
        RAD_REINIT_2(Val = i; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
IndepADvar<Double>::IndepADvar(long i) Allow_noderiv(: IndepADvar_base<Double>(1))
{
        RAD_REINIT_1(cv = 0;)
        cv = new ADVari(Double(i));
        RAD_REINIT_2(Val = i; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
ConstADvar<Double>::ConstADvar()
{
        RAD_REINIT_1(this->cv = 0;)
        this->cv = new ConstADVari(0.);
        RAD_REINIT_2(this->Val = 0.; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double> void
ConstADvar<Double>::ConstADvar_ctr(Double d)
{
        RAD_REINIT_1(this->cv = 0;)
        this->cv = new ConstADVari(d);
        RAD_REINIT_2(this->Val = d; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
ConstADvar<Double>::ConstADvar(const IndepADVar &x)
{
        RAD_cvchk(x)
        RAD_REINIT_1(this->cv = 0;)
        ConstADVari *y = new ConstADVari(x.cv->Val);
#if RAD_REINIT > 0
        x.cv->adc.new_Derp(&x.adc.One, y, x.cv);
#else
        DErp *d = new DErp(&x.adc.One, y, x.cv);
#endif
        this->cv = y;
        RAD_REINIT_2(this->Val = y->Val; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
ConstADvar<Double>::ConstADvar(const ConstADvar &x)
{
        RAD_cvchk(x)
        RAD_REINIT_1(this->cv = 0;)
        ConstADVari *y = new ConstADVari(x.cv->Val);
#if RAD_REINIT > 0
        x.cv->adc.new_Derp(&x.cv->adc.One, y, x.cv);
#else
        DErp *d = new DErp(&x.cv->adc.One, y, (ADVari*)x.cv);
#endif
        this->cv = y;
        RAD_REINIT_2(this->Val = y->Val; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
ConstADvar<Double>::ConstADvar(const ADVari &x)
{
        RAD_REINIT_1(this->cv = 0;)
        ConstADVari *y = new ConstADVari(x.Val);
#if RAD_REINIT > 0
        x.adc.new_Derp(&x.adc.One, y, &x);
#else
        DErp *d = new DErp(&x.adc.One, y, &x);
#endif
        this->cv = y;
        RAD_REINIT_2(this->Val = y->Val; this->gen = this->IndepADvar_root.gen;)
        }

 template<typename Double>
 void
IndepADvar<Double>::AD_Const(const IndepADvar &v)
{
        typedef ConstADvari<Double> ConstADVari;

        ConstADVari *ncv = new ConstADVari(v.val());
#ifdef RAD_AUTO_AD_Const
        v.cv->padv = 0;
#endif
        v.cv = ncv;
        RAD_REINIT_2(v.gen = v.IndepADvar_root.gen;)
        }

 template<typename Double>
 int
IndepADvar<Double>::Wantderiv(int n)
{
#ifdef RAD_ALLOW_WANTDERIV
#if RAD_REINIT == 2
        if (this->gen != this->IndepADvar_root.gen) {
                cv = new ADVari(Val);
                this->gen = this->IndepADvar_root.gen;
                }
#endif
        return this->wantderiv = cv->wantderiv = n;
#else
        return 1;
#endif // RAD_ALLOW_WANTDERIV
        }

 template<typename Double>
 void
ConstADvari<Double>::aval_reset()
{
#if RAD_REINIT == 0
        ConstADvari *x = ConstADvari::lastcad;
        while(x) {
                x->aval = 0;
                x = x->prevcad;
                }
#elif RAD_REINIT == 1
        ADvari<Double>::adc.new_ADmemblock(0);
#endif
        }

#ifdef RAD_AUTO_AD_Const

 template<typename Double>
ADvari<Double>::ADvari(const IndepADVar *x, Double d): Val(d), aval(0.)
{
        this->ADvari_padv(x);
        Allow_noderiv(wantderiv = 1;)
        }

 template<typename Double>
ADvari<Double>::ADvari(const IndepADVar *x, Double d, Double g): Val(d), aval(g)
{
        this->ADvari_padv(x);
        Allow_noderiv(wantderiv = 1;)
        }

 template<typename Double>
ADvar1<Double>::ADvar1(const IndepADVar *x, const IndepADVar &y):
        ADVari(y.cv->Val), d((const Double*)&ADcontext<Double>::One, (ADVari*)this, y.cv)
{
        this->ADvari_padv(x);
        }

 template<typename Double>
ADvar1<Double>::ADvar1(const IndepADVar *x, const ADVari &y):
        ADVari(y.Val), d((const Double*)&ADcontext<Double>::One, this, &y)
{
        this->ADvari_padv(x);
        }

#else /* !RAD_AUTO_AD_Const */

 template<typename Double>
 IndepADvar<Double>&
ADvar_operatoreq(IndepADvar<Double> *This, const ADvari<Double> &x)
{
        This->cv = new ADvar1<Double>(x.Val, &x.adc.One, &x);
        RAD_REINIT_1(This->Move_to_end();)
        RAD_REINIT_2(This->Val = x.Val; This->gen = This->IndepADvar_root.gen;)
        return *(IndepADvar<Double>*) This;
        }

 template<typename Double>
 ADvar<Double>&
ADvar_operatoreq(ADvar<Double> *This, const ADvari<Double> &x)
{
        This->cv = new ADvar1<Double>(x.Val, &x.adc.One, &x);
        RAD_REINIT_1(This->Move_to_end();)
        RAD_REINIT_2(This->Val = x.Val; This->gen = This->IndepADvar_root.gen;)
        return *(ADvar<Double>*) This;
        }

#endif /* RAD_AUTO_AD_Const */


 template<typename Double>
 IndepADvar<Double>&
IndepADvar<Double>::operator=(Double d)
{
#ifdef RAD_AUTO_AD_Const
        if (this->cv)
                this->cv->padv = 0;
        this->cv = new ADVari(this,d);
#else
        this->cv = new ADVari(d);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
#endif
        return *this;
        }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator=(Double d)
{
#ifdef RAD_AUTO_AD_Const
        if (this->cv)
                this->cv->padv = 0;
        this->cv = new ADVari(this,d);
#else
        this->cv = RAD_REINIT_0(ConstADVari::cadc.fpval_implies_const
                ? new ConstADVari(d)
                : ) new ADVari(d);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = d; this->gen = this->IndepADvar_root.gen;)
#endif
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator-(const Base< ADvari<Double> > &TT) {
   const ADvari<Double>& T = TT.derived();
   return *(new ADvar1<Double>(-T.Val, &T.adc.negOne, &T));
 }

 template<typename Double>
 ADvari<Double>&
operator+(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   return *(new ADvar2<Double>(L.Val + R.Val, &L, &L.adc.One, &R, &L.adc.One));
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
        this->cv = new ADvar2<Double>(Lcv->Val + R.Val, Lcv, &R.adc.One, &R, &R.adc.One RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator+(const Base< ADvari<Double> > &LL, Double R) {
   const ADvari<Double>& L = LL.derived();
   return *(new ADvar1<Double>(L.Val + R, &L.adc.One, &L));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator+=(Double R) {
        ADVari *tcv = this->cv;
        this->cv = new ADVar1(tcv->Val + R, &tcv->adc.One, tcv RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator+(Double L, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   return *(new ADvar1<Double>(L + R.Val, &R.adc.One, &R));
 }

 template<typename Double>
 ADvari<Double>&
operator-(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   return *(new ADvar2<Double>(L.Val - R.Val, &L, &L.adc.One, &R, &L.adc.negOne));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator-=(const ADVari &R) {
        ADVari *Lcv = this->cv;
        this->cv = new ADvar2<Double>(Lcv->Val - R.Val, Lcv, &R.adc.One, &R, &R.adc.negOne RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator-(const Base< ADvari<Double> > &LL, Double R) {
   const ADvari<Double>& L = LL.derived();
   return *(new ADvar1<Double>(L.Val - R, &L.adc.One, &L));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator-=(Double R) {
        ADVari *tcv = this->cv;
        this->cv = new ADVar1(tcv->Val - R, &tcv->adc.One, tcv RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator-(Double L, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   return *(new ADvar1<Double>(L - R.Val, &R.adc.negOne, &R));
 }

 template<typename Double>
 ADvari<Double>&
operator*(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   return *(new ADvar2<Double>(L.Val * R.Val, &L, &R.Val, &R, &L.Val));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator*=(const ADVari &R) {
        ADVari *Lcv = this->cv;
        this->cv = new ADvar2<Double>(Lcv->Val * R.Val, Lcv, &R.Val, &R, &Lcv->Val RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator*(const Base< ADvari<Double> > &LL, Double R) {
   const ADvari<Double>& L = LL.derived();
   return *(new ADvar1s<Double>(L.Val * R, R, &L));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator*=(Double R) {
        ADVari *Lcv = this->cv;
        this->cv = new ADvar1s<Double>(Lcv->Val * R, R, Lcv RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator*(Double L, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   return *(new ADvar1s<Double>(L * R.Val, L, &R));
 }

 template<typename Double>
 ADvari<Double>&
operator/(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   Double Lv = L.Val, Rv = R.Val, pL = 1. / Rv, q = Lv/Rv;
   return *(new ADvar2q<Double>(q, pL, -q*pL, &L, &R));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator/=(const ADVari &R) {
        ADVari *Lcv = this->cv;
        Double Lv = Lcv->Val, Rv = R.Val, pL = 1. / Rv, q = Lv/Rv;
        this->cv = new ADvar2q<Double>(q, pL, -q*pL, Lcv, &R RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
operator/(const Base< ADvari<Double> > &LL, Double R) {
   const ADvari<Double>& L = LL.derived();
   return *(new ADvar1s<Double>(L.Val / R, 1./R, &L));
 }

 template<typename Double>
 ADvari<Double>&
operator/(Double L, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   Double recip = 1. / R.Val;
   Double q = L * recip;
   return *(new ADvar1s<Double>(q, -q*recip, &R));
 }

 template<typename Double>
 ADvar<Double>&
ADvar<Double>::operator/=(Double R) {
        ADVari *Lcv = this->cv;
        this->cv = new ADvar1s<Double>(Lcv->Val / R, 1./R, Lcv RAD_ACA);
        RAD_REINIT_1(this->Move_to_end();)
        RAD_REINIT_2(this->Val = this->cv->Val;)
        RAD_REINIT_2(this->gen = this->IndepADvar_root.gen;)
        return *this;
        }

 template<typename Double>
 ADvari<Double>&
acos(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = v.Val;
   return *(new ADvar1s<Double>(std::acos(t), -1./std::sqrt(1. - t*t), &v));
 }

 template<typename Double>
 ADvari<Double>&
acosh(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = v.Val, t1 = std::sqrt(t*t - 1.);
   return *(new ADvar1s<Double>(std::log(t + t1), 1./t1, &v));
 }

 template<typename Double>
 ADvari<Double>&
asin(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = v.Val;
   return *(new ADvar1s<Double>(std::asin(t), 1./std::sqrt(1. - t*t), &v));
 }

 template<typename Double>
 ADvari<Double>&
asinh(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = v.Val, td = 1., t1 = std::sqrt(t*t + 1.);
   if (t < 0.) {
     t = -t;
     td = -1.;
   }
   return *(new ADvar1s<Double>(td*std::log(t + t1), 1./t1, &v));
 }

 template<typename Double>
 ADvari<Double>&
atan(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = v.Val;
   return *(new ADvar1s<Double>(std::atan(t), 1./(1. + t*t), &v));
 }

 template<typename Double>
 ADvari<Double>&
atanh(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = v.Val;
   return *(new ADvar1s<Double>(0.5*std::log((1.+t)/(1.-t)), 1./(1. - t*t), &v));
 }

 template<typename Double>
 ADvari<Double>&
atan2(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   Double x = L.Val, y = R.Val, t = x*x + y*y;
   return *(new ADvar2q<Double>(std::atan2(x,y), y/t, -x/t, &L, &R));
 }

 template<typename Double>
 ADvari<Double>&
atan2(Double x, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   Double y = R.Val, t = x*x + y*y;
   return *(new ADvar1s<Double>(std::atan2(x,y), -x/t, &R));
 }

 template<typename Double>
 ADvari<Double>&
atan2(const Base< ADvari<Double> > &LL, Double y) {
   const ADvari<Double>& L = LL.derived();
   Double x = L.Val, t = x*x + y*y;
   return *(new ADvar1s<Double>(std::atan2(x,y), y/t, &L));
 }

 template<typename Double>
 ADvari<Double>&
max(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   const ADvari<Double> &x = L.Val >= R.Val ? L : R;
   return *(new ADvar1<Double>(x.Val, &x.adc.One, &x));
 }

 template<typename Double>
 ADvari<Double>&
max(Double L, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   if (L >= R.Val)
     return *(new ADvari<Double>(L));
   return *(new ADvar1<Double>(R.Val, &R.adc.One, &R));
 }

 template<typename Double>
 ADvari<Double>&
max(const Base< ADvari<Double> > &LL, Double R) {
   const ADvari<Double>& L = LL.derived();
   if (L.Val >= R)
     return *(new ADvar1<Double>(L.Val, &L.adc.One, &L));
   return *(new ADvari<Double>(R));
 }

 template<typename Double>
 ADvari<Double>&
min(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   const ADvari<Double> &x = L.Val <= R.Val ? L : R;
   return *(new ADvar1<Double>(x.Val, &x.adc.One, &x));
 }

 template<typename Double>
 ADvari<Double>&
min(Double L, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   if (L <= R.Val)
     return *(new ADvari<Double>(L));
   return *(new ADvar1<Double>(R.Val, &R.adc.One, &R));
 }

 template<typename Double>
 ADvari<Double>&
min(const Base< ADvari<Double> > &LL, Double R) {
   const ADvari<Double>& L = LL.derived();
   if (L.Val <= R)
     return *(new ADvar1<Double>(L.Val, &L.adc.One, &L));
   return *(new ADvari<Double>(R));
 }

 template<typename Double>
 ADvari<Double>&
cos(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   return *(new ADvar1s<Double>(std::cos(v.Val), -std::sin(v.Val), &v));
 }

 template<typename Double>
 ADvari<Double>&
cosh(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   return *(new ADvar1s<Double>(std::cosh(v.Val), std::sinh(v.Val), &v));
 }

 template<typename Double>
 ADvari<Double>&
exp(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
#if RAD_REINIT > 0
   Double t = std::exp(v.Val);
   return *(new ADvar1s<Double>(t, t, &v));
#else
   ADvar1<Double>* rcv = new ADvar1<Double>(std::exp(v.Val), &v);
   rcv->d.a = &rcv->Val;
   rcv->d.b = rcv;
   return *rcv;
#endif
 }

 template<typename Double>
 ADvari<Double>&
log(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double x = v.Val;
   return *(new ADvar1s<Double>(std::log(x), 1. / x, &v));
 }

 template<typename Double>
 ADvari<Double>&
log10(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   static double num = 1. / std::log(10.);
   Double x = v.Val;
   return *(new ADvar1s<Double>(std::log10(x), num / x, &v));
 }

 template<typename Double>
 ADvari<Double>&
pow(const Base< ADvari<Double> > &LL, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& L = LL.derived();
   const ADvari<Double>& R = RR.derived();
   Double x = L.Val, y = R.Val, t = std::pow(x,y);
   return *(new ADvar2q<Double>(t, y*t/x, t*std::log(x), &L, &R));
 }

 template<typename Double>
 ADvari<Double>&
pow(Double x, const Base< ADvari<Double> > &RR) {
   const ADvari<Double>& R = RR.derived();
   Double t = std::pow(x,R.Val);
   return *(new ADvar1s<Double>(t, t*std::log(x), &R));
 }

 template<typename Double>
 ADvari<Double>&
pow(const Base< ADvari<Double> > &LL, Double y) {
   const ADvari<Double>& L = LL.derived();
   Double x = L.Val, t = std::pow(x,y);
   return *(new ADvar1s<Double>(t, y*t/x, &L));
 }

 template<typename Double>
 ADvari<Double>&
sin(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   return *(new ADvar1s<Double>(std::sin(v.Val), std::cos(v.Val), &v));
 }

 template<typename Double>
 ADvari<Double>&
sinh(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   return *(new ADvar1s<Double>(std::sinh(v.Val), std::cosh(v.Val), &v));
 }

 template<typename Double>
 ADvari<Double>&
sqrt(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = std::sqrt(v.Val);
   return *(new ADvar1s<Double>(t, 0.5/t, &v));
 }

 template<typename Double>
 ADvari<Double>&
tan(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = std::cos(v.Val);
   return *(new ADvar1s<Double>(std::tan(v.Val), 1./(t*t), &v));
 }

 template<typename Double>
 ADvari<Double>&
tanh(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = 1. / std::cosh(v.Val);
   return *(new ADvar1s<Double>(std::tanh(v.Val), t*t, &v));
 }

 template<typename Double>
 ADvari<Double>&
abs(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t, p;
   p = 1;
   if ((t = v.Val) < 0) {
     t = -t;
     p = -p;
   }
   return *(new ADvar1s<Double>(t, p, &v));
 }

 template<typename Double>
 ADvari<Double>&
fabs(const Base< ADvari<Double> > &vv) {
   // Synonym for "abs"
   // "fabs" is not the best choice of name,
   // but this name is used at Sandia.
   const ADvari<Double>& v = vv.derived();
   Double t, p;
   p = 1;
   if ((t = v.Val) < 0) {
     t = -t;
     p = -p;
   }
   return *(new ADvar1s<Double>(t, p, &v));
 }

template<typename Double>
 ADvari<Double>&
cbrt(const Base< ADvari<Double> > &vv) {
   const ADvari<Double>& v = vv.derived();
   Double t = std::cbrt(v.Val);
   return *(new ADvar1s<Double>(t, 1.0/(3.0*t*t), &v));
 }

 template<typename Double>
 ADvari<Double>&
ADf1(Double f, Double g, const ADvari<Double> &x) {
        return *(new ADvar1s<Double>(f, g, &x));
        }

 template<typename Double>
 inline ADvari<Double>&
ADf1(Double f, Double g, const IndepADvar<Double> &x) {
        return *(new ADvar1s<Double>(f, g, x.cv));
        }

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, const ADvari<Double> &x, const ADvari<Double> &y) {
        return *(new ADvar2q<Double>(f, gx, gy, &x, &y));
        }

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, const ADvari<Double> &x, const IndepADvar<Double> &y) {
        return *(new ADvar2q<Double>(f, gx, gy, &x, y.cv));
        }

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, const IndepADvar<Double> &x, const ADvari<Double> &y) {
        return *(new ADvar2q<Double>(f, gx, gy, x.cv, &y));
        }

 template<typename Double>
 ADvari<Double>&
ADf2(Double f, Double gx, Double gy, const IndepADvar<Double> &x, const IndepADvar<Double> &y) {
        return *(new ADvar2q<Double>(f, gx, gy, x.cv, y.cv));
        }

 template<typename Double>
 ADvari<Double>&
ADfn(Double f, int n, const IndepADvar<Double> *x, const Double *g) {
        return *(new ADvarn<Double>(f, n, x, g));
        }

 template<typename Double>
 inline ADvari<Double>&
ADfn(Double f, int n, const ADvar<Double> *x, const Double *g) {
        return ADfn<Double>(f, n, (IndepADvar<Double>*)x, g);
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
#define Ai const Base< ADvari<Double> >&
#define AI const Base< IndepADvar<Double> >&
#define D Double
#define CAI(x,y) const IndepADvar<Double> & x = y.derived()
#define CAi(x,y) const ADvari<Double> & x = y.derived()
#define T2(r,f) \
 T r f(Ai LL, AI RR) { CAi(L,LL); CAI(R,RR); RAD_cvchk(R) return f(L, C(R.cv)); } \
 T r f(AI LL, Ai RR) { CAI(L,LL); CAi(R,RR); RAD_cvchk(L) return f(C(L.cv), R); }\
 T r f(AI LL, AI RR) { CAI(L,LL); CAI(R,RR); RAD_cvchk(L) RAD_cvchk(R) return f(C(L.cv), C(R.cv)); }\
 T r f(AI LL, D R) { CAI(L,LL); RAD_cvchk(L) return f(C(L.cv), R); } \
 T r f(D  L,  AI RR) { CAI(R,RR); RAD_cvchk(R) return f(L, C(R.cv)); } \
 T r f(Ai L,  Dtype R) { return f(L, (D)R); }\
 T r f(AI L,  Dtype R) { return f(L, (D)R); }\
 T r f(Ai L,  Ltype R) { return f(L, (D)R); }\
 T r f(AI L,  Ltype R) { return f(L, (D)R); }\
 T r f(Ai L,  Itype R) { return f(L, (D)R); }\
 T r f(AI L,  Itype R) { return f(L, (D)R); }\
 T r f(Dtype L, Ai R) { return f((D)L, R); }\
 T r f(Dtype L, AI R) {return f((D)L, R); }\
 T r f(Ltype L, Ai R) { return f((D)L, R); }\
 T r f(Ltype L, AI R) { return f((D)L, R); }\
 T r f(Itype L, Ai R) { return f((D)L, R); }\
 T r f(Itype L, AI R) { return f((D)L, R); }

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
  T F f(AI xx) { CAI(x,xx); RAD_cvchk(x) return f(C(x.cv)); }

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
T1(cbrt)

T F copy(AI xx)
{
  CAI(x,xx);
  RAD_cvchk(x)
    return *(new ADvar1<Double>(x.cv->Val, &ADcontext<Double>::One, (ADvari<Double>*)x.cv));
}

T F copy(Ai xx)
{ CAi(x,xx); return *(new ADvar1<Double>(x.Val, &ADcontext<Double>::One, (ADvari<Double>*)&x)); }

#undef RAD_cvchk
#undef T1
#undef AI
#undef Ai
#undef F
#undef T
#undef A
#undef C
#undef Ttype
#undef Dtype
#undef Ltype
#undef Itype
#undef CAI
#undef CAi
#undef D

} /* namespace Rad */
} /* namespace Sacado */

#undef SNS
#undef RAD_REINIT_2
#undef Allow_noderiv

#endif /* SACADO_TRAD_H */
