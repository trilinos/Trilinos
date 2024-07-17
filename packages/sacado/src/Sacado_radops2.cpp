// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// Support routines for Hessian-vector products via the RAD package
// (Reverse Automatic Differentiation) -- a package specialized for
// function and gradient evaluations, and herewith extended for Hessian-
// vector products.  This specialization is for several Hessian-vector
// products at one point, where the vectors are determined on the fly,
// e.g., by the conjugate-gradient algorithm.  Where the vectors are known
// in advance, nestings of RAD and FAD might be more efficient.
// RAD Hessian-vector product extension written in 2007 by David M. Gay at
// Sandia National Labs, Albuquerque, NM.

#include "Sacado_rad2.hpp"

#ifndef SACADO_NO_NAMESPACE
namespace Sacado {
namespace Rad2d { // "2" for 2nd derivatives, "d" for "double"
#endif

#ifdef RAD_AUTO_AD_Const
ADvari *ADvari::First_ADvari, **ADvari::Last_ADvari = &ADvari::First_ADvari;
#undef RAD_DEBUG_BLOCKKEEP
#else
#ifdef RAD_DEBUG_BLOCKKEEP
#if !(RAD_DEBUG_BLOCKKEEP > 0)
#undef RAD_DEBUG_BLOCKKEEP
#else
extern "C" void _uninit_f2c(void *x, int type, long len);
#define TYDREAL 5
static ADmemblock *rad_Oldcurmb;
static int rad_busy_blocks;
#endif /*RAD_DEBUG_BLOCKKEEP > 0*/
#endif /*RAD_DEBUG_BLOCKKEEP*/
#endif /*RAD_AUTO_AD_Const*/

Derp *Derp::LastDerp = 0;
ADcontext ADvari::adc;
CADcontext ConstADvari::cadc;
ConstADvari *ConstADvari::lastcad;
static int rad_need_reinit;
#ifdef RAD_DEBUG_BLOCKKEEP
static size_t rad_mleft_save;
#endif

const double CADcontext::One = 1., CADcontext::negOne = -1.;

ADcontext::ADcontext()
{
	First.next = 0;
	Busy = &First;
	Free = 0;
	Mbase = (char*)First.memblk;
	Mleft = sizeof(First.memblk);
	Aibusy = &AiFirst;
	Aifree = 0;
	Ainext = AiFirst.pADvari;
	AiFirst.next = AiFirst.prev = 0;
	AiFirst.limit = Ailimit = AiFirst.pADvari + ADvari_block::Gulp;
	}

 void*
ADcontext::new_ADmemblock(size_t len)
{
	ADmemblock *mb, *mb0, *mb1, *mbf, *x;
	ADvari_block *b;
#ifdef RAD_AUTO_AD_Const
	ADvari *a, *anext;
	IndepADvar *v;
#endif /* RAD_AUTO_AD_Const */

	if (rad_need_reinit && this == &ADvari::adc) {
		rad_need_reinit = 0;
		Derp::LastDerp = 0;
		Aibusy = b = &AiFirst;
		Aifree = b->next;
		b->next = b->prev = 0;
		Ailimit = b->limit = (Ainext = b->pADvari) + ADvari_block::Gulp;
#ifdef RAD_DEBUG_BLOCKKEEP
		Mleft = rad_mleft_save;
		if (Mleft < sizeof(First.memblk))
			_uninit_f2c(Mbase + Mleft, TYDREAL,
			 (sizeof(First.memblk) - Mleft)/sizeof(double));
		if ((mb = Busy->next)) {
			if (!(mb0 = rad_Oldcurmb))
				mb0 = &First;
			for(;; mb = mb->next) {
				_uninit_f2c(mb->memblk, TYDREAL,
					sizeof(First.memblk)/sizeof(double));
				if (mb == mb0)
					break;
				}
			}
		rad_Oldcurmb = Busy;
		if (rad_busy_blocks >= RAD_DEBUG_BLOCKKEEP) {
			rad_busy_blocks = 0;
			rad_Oldcurmb = 0;
			mb0 = &First;
			mbf =  Free;
			for(mb = Busy; mb != mb0; mb = mb1) {
				mb1 = mb->next;
				mb->next = mbf;
				mbf = mb;
				}
			Free = mbf;
			Busy = mb;
			Mbase = (char*)First.memblk;
			Mleft = sizeof(First.memblk);
			}

#else /* !RAD_DEBUG_BLOCKKEEP */

		mb0 = &ADvari::adc.First;
		mbf =  ADvari::adc.Free;
		for(mb = ADvari::adc.Busy; mb != mb0; mb = mb1) {
			mb1 = mb->next;
			mb->next = mbf;
			mbf = mb;
			}
		ADvari::adc.Free = mbf;
		ADvari::adc.Busy = mb;
		ADvari::adc.Mbase = (char*)ADvari::adc.First.memblk;
		ADvari::adc.Mleft = sizeof(ADvari::adc.First.memblk);
#ifdef RAD_AUTO_AD_Const
		*ADvari::Last_ADvari = 0;
		ADvari::Last_ADvari = &ADvari::First_ADvari;
		if ((anext = ADvari::First_ADvari)) {
			while((a = anext)) {
				anext = a->Next;
				if ((v = a->padv))
					v->cv = new ADvari(v, a->Val);
				}
			}
#endif /*RAD_AUTO_AD_Const*/
#endif /* RAD_DEBUG_BLOCKKEEP */
		if (Mleft >= len)
			return Mbase + (Mleft -= len);
		}

	if ((x = Free))
		Free = x->next;
	else
		x = new ADmemblock;
#ifdef RAD_DEBUG_BLOCKKEEP
	rad_busy_blocks++;
#endif
	x->next = Busy;
	Busy = x;
	return (Mbase = (char*)x->memblk) +
		(Mleft = sizeof(First.memblk) - len);
	}

 void
ADcontext::new_ADvari_block()
{
	ADvari_block *ob, *nb;
	ob = Aibusy;
	ob->limit = Ailimit;	// should be unnecessary, but harmless
	if ((nb = Aifree))
		Aifree = nb->next;
	else
		nb = new ADvari_block;
	Aibusy = Aibusy->next = nb;
	nb->limit = Ailimit = (Ainext = nb->pADvari) + ADvari_block::Gulp;
	ob->next = nb;
	nb->prev = ob;
	nb->next = 0;
	}

 void
ADcontext::Gradcomp(int wantgrad)
{
	Derp *d;

	if (rad_need_reinit && wantgrad) {
		for(d = Derp::LastDerp; d; d = d->next)
			d->c->aval = 0;
		}
	else {
		rad_need_reinit = 1;
#ifdef RAD_DEBUG_BLOCKKEEP
		rad_mleft_save = ADvari::adc.Mleft;
#endif
		ADvari::adc.Mleft = 0;
		}

	if ((d = Derp::LastDerp) && wantgrad) {
		d->b->aval = 1;
		do d->c->aval += *d->a * d->b->aval;
		while((d = d->next));
		}
	}

 void
ADcontext::Weighted_Gradcomp(int n, ADvar **v, double *w)
{
	Derp *d;
	int i;

	if (rad_need_reinit) {
		for(d = Derp::LastDerp; d; d = d->next)
			d->c->aval = 0;
		}
	else {
		rad_need_reinit = 1;
#ifdef RAD_DEBUG_BLOCKKEEP
		rad_mleft_save = ADvari::adc.Mleft;
#endif
		ADvari::adc.Mleft = 0;
		}
	if ((d = Derp::LastDerp)) {
		for(i = 0; i < n; i++)
			v[i]->cv->aval = w[i];
		do d->c->aval += *d->a * d->b->aval;
		while((d = d->next));
		}
	}

#ifdef RAD_AUTO_AD_Const

IndepADvar::IndepADvar(double d)
{
	cv = new ADvari(this, d);
	}

IndepADvar::IndepADvar(int d)
{
	cv = new ADvari(this, (double)d);
	}

IndepADvar::IndepADvar(long d)
{
	cv = new ADvari(this, (double)d);
	}

#else /*!RAD_AUTO_AD_Const*/

IndepADvar::IndepADvar(double d)
{
	ADvari *x = new ADvari(Hv_const, d);
	cv = x;
	}

IndepADvar::IndepADvar(int d)
{
	ADvari *x = new ADvari(Hv_const, (double)d);
	cv = x;
	}

IndepADvar::IndepADvar(long d)
{
	ADvari *x = new ADvari(Hv_const, (double)d);
	cv = x;
	}

#endif /*RAD_AUTO_AD_Const*/

 void
ADvar::ADvar_ctr(double d)
{
#ifdef RAD_AUTO_AD_Const
	ADvari *x = new ADvari(this, d);
#else
	ADvari *x = ConstADvari::cadc.fpval_implies_const
		? new ConstADvari(d)
		: new ADvari(Hv_const, d);
#endif
	cv = x;
	}

ConstADvar::ConstADvar()
{
	ConstADvari *x = new ConstADvari(0.);
	cv = x;
	}

 void
ConstADvar::ConstADvar_ctr(double d)
{
	ConstADvari *x = new ConstADvari(d);
	cv = x;
	}
ConstADvar::ConstADvar(const ADvari &x)
{
	ConstADvari *y = new ConstADvari(x.Val);
	new Derp(&CADcontext::One, y, &x); /* for side effect; value not used */
	cv = y;
	}

 void
IndepADvar::AD_Const(const IndepADvar &v)
{
	ConstADvari *ncv = new ConstADvari(v.val());
	((IndepADvar*)&v)->cv = ncv;
	}

 void
ConstADvari::aval_reset()
{
	ConstADvari *x = ConstADvari::lastcad;
	while(x) {
		x->aval = 0;
		x = x->prevcad;
		}
	}

#ifdef RAD_AUTO_AD_Const

ADvari::ADvari(const IndepADvar *x, double d): Val(d), aval(0.)
{
	opclass = Hv_const;
	*Last_ADvari = this;
	Last_ADvari = &Next;
	padv = (IndepADvar*)x;
	}

ADvar1::ADvar1(const IndepADvar *x, const IndepADvar &y):
	ADvari(Hv_copy, y.cv->Val), d(&CADcontext::One, this, y.cv)
{
	*ADvari::Last_ADvari = this;
	ADvari::Last_ADvari = &Next;
	padv = (IndepADvar*)x;
	}

ADvar1::ADvar1(const IndepADvar *x, const ADvari &y):
	ADvari(Hv_copy, y.Val), d(&CADcontext::One, this, &y)
{
	*ADvari::Last_ADvari = this;
	ADvari::Last_ADvari = &Next;
	padv = (IndepADvar*)x;
	}

#else

 ADvar&
ADvar_operatoreq(ADvar *This, const ADvari &x)
{ This->cv = new ADvar1(Hv_copy, x.Val, &CADcontext::One, (ADvari*)&x); return *This; }

 IndepADvar&
ADvar_operatoreq(IndepADvar *This, const ADvari &x)
{ This->cv = new ADvar1(Hv_copy, x.Val, &CADcontext::One, (ADvari*)&x); return *This; }

#endif /*RAD_AUTO_AD_Const*/

ADvar2q::ADvar2q(double val1, double Lp, double Rp,
		double LR, double R2, const ADvari *Lcv, const ADvari *Rcv):
			ADvar2(Hv_quotLR,val1,Lcv,&pL,Rcv,&pR),
			pL(Lp), pR(Rp), pLR(LR), pR2(R2) { }

ADvar2g::ADvar2g(double val1, double Lp, double Rp,
		double L2, double LR, double R2, const ADvari *Lcv, const ADvari *Rcv):
			ADvar2(Hv_binary,val1,Lcv,&pL,Rcv,&pR),
			pL(Lp), pR(Rp), pL2(L2), pLR(LR), pR2(R2) { }

 IndepADvar&
IndepADvar::operator=(double d)
{
#ifdef RAD_AUTO_AD_Const
	if (cv)
		cv->padv = 0;
	cv = new ADvari(this,d);
#else
	cv = new ADvari(Hv_const, d);
#endif
	return *this;
	}

 ADvar&
ADvar::operator=(double d)
{
#ifdef RAD_AUTO_AD_Const
	if (cv)
		cv->padv = 0;
	cv = new ADvari(this, d);
#else
	cv = ConstADvari::cadc.fpval_implies_const
		? new ConstADvari(d)
		: new ADvari(Hv_const, d);
#endif
	return *this;
	}

 ADvari&
operator-(const ADvari &T) {
	return *(new ADvar1(Hv_negate, -T.Val, &CADcontext::negOne, &T));
	}

 ADvari&
operator+(const ADvari &L, const ADvari &R) {
	return *(new ADvar2(Hv_plusLR,L.Val + R.Val, &L, &CADcontext::One, &R, &CADcontext::One));
	}

 ADvar&
ADvar::operator+=(const ADvari &R) {
	ADvari *Lcv = cv;
#ifdef RAD_AUTO_AD_Const
	Lcv->padv = 0;
#endif
	cv = new ADvar2(Hv_plusLR, Lcv->Val + R.Val, Lcv, &CADcontext::One, &R, &CADcontext::One);
	return *this;
	}

 ADvari&
operator+(const ADvari &L, double R) {
	return *(new ADvar1(Hv_copy, L.Val + R, &CADcontext::One, &L));
	}

 ADvar&
ADvar::operator+=(double R) {
	ADvari *tcv = cv;
#ifdef RAD_AUTO_AD_Const
	tcv->padv = 0;
#endif
	cv = new ADvar1(Hv_copy, tcv->Val + R, &CADcontext::One, tcv);
	return *this;
	}

 ADvari&
operator+(double L, const ADvari &R) {
	return *(new ADvar1(Hv_copy, L + R.Val, &CADcontext::One, &R));
	}

 ADvari&
operator-(const ADvari &L, const ADvari &R) {
	return *(new ADvar2(Hv_minusLR,L.Val - R.Val, &L, &CADcontext::One, &R, &CADcontext::negOne));
	}

 ADvar&
ADvar::operator-=(const ADvari &R) {
	ADvari *Lcv = cv;
#ifdef RAD_AUTO_AD_Const
	Lcv->padv = 0;
#endif
	cv = new ADvar2(Hv_minusLR, Lcv->Val - R.Val, Lcv, &CADcontext::One, &R, &CADcontext::negOne);
	return *this;
	}

 ADvari&
operator-(const ADvari &L, double R) {
	return *(new ADvar1(Hv_copy, L.Val - R, &CADcontext::One, &L));
	}

 ADvar&
ADvar::operator-=(double R) {
	ADvari *tcv = cv;
#ifdef RAD_AUTO_AD_Const
	tcv->padv = 0;
#endif
	cv = new ADvar1(Hv_copy, tcv->Val - R, &CADcontext::One, tcv);
	return *this;
	}

 ADvari&
operator-(double L, const ADvari &R) {
	return *(new ADvar1(Hv_negate, L - R.Val, &CADcontext::negOne, &R));
	}

 ADvari&
operator*(const ADvari &L, const ADvari &R) {
	return *(new ADvar2(Hv_timesLR, L.Val * R.Val, &L, &R.Val, &R, &L.Val));
	}

 ADvar&
ADvar::operator*=(const ADvari &R) {
	ADvari *Lcv = cv;
#ifdef RAD_AUTO_AD_Const
	Lcv->padv = 0;
#endif
	cv = new ADvar2(Hv_timesLR, Lcv->Val * R.Val, Lcv, &R.Val, &R, &Lcv->Val);
	return *this;
	}

 ADvari&
operator*(const ADvari &L, double R) {
	return *(new ADvar1s(L.Val * R, R, &L));
	}

 ADvar&
ADvar::operator*=(double R) {
	ADvari *Lcv = cv;
#ifdef RAD_AUTO_AD_Const
	Lcv->padv = 0;
#endif
	cv = new ADvar1s(Lcv->Val * R, R, Lcv);
	return *this;
	}

 ADvari&
operator*(double L, const ADvari &R) {
	return *(new ADvar1s(L * R.Val, L, &R));
	}

 ADvari&
operator/(const ADvari &L, const ADvari &R) {
	double Lv = L.Val, Rv = R.Val, pL = 1. / Rv, q = Lv/Rv, qpL = q*pL;
	return *(new ADvar2q(q, pL, -qpL, -pL*pL, 2.*pL*qpL, &L, &R));
	}

 ADvar&
ADvar::operator/=(const ADvari &R) {
	ADvari *Lcv = cv;
#ifdef RAD_AUTO_AD_Const
	Lcv->padv = 0;
#endif
	double Lv = Lcv->Val, Rv = R.Val, pL = 1. / Rv, q = Lv/Rv, qpL = q*pL;
	cv = new ADvar2q(q, pL, -qpL, -pL*pL, 2.*pL*qpL, Lcv, &R);
	return *this;
	}

 ADvari&
operator/(const ADvari &L, double R) {
	return *(new ADvar1s(L.Val / R, 1./R, &L));
	}

 ADvari&
operator/(double L, const ADvari &R) {
	double recip = 1. / R.Val;
	double q = L * recip;
	double d1 = -q*recip;
	return *(new ADvar1g(q, d1, -q*d1, &R));
	}

 ADvar&
ADvar::operator/=(double R) {
	ADvari *Lcv = cv;
#ifdef RAD_AUTO_AD_Const
	Lcv->padv = 0;
#endif
	cv = new ADvar1s(Lcv->Val / R, 1./R, Lcv);
	return *this;
	}

 ADvari&
acos(const ADvari &v) {
	double t = v.Val;
	double t1 = 1. - t*t;
	double d1 = -1. / sqrt(t1);
	return *(new ADvar1g(acos(t), d1, t*d1/t1, &v));
	}

 ADvari&
acosh(const ADvari &v) {
	double d1, t, t1, t2;
	t = v.Val;
	t1 = sqrt(t2 = t*t - 1.);
	d1 = 1./t1;
	return *(new ADvar1g(log(t + t1), d1, -t*d1/t2, &v));
	}

 ADvari&
asin(const ADvari &v) {
	double d1, t, t1;
	t = v.Val;
	d1 = 1. / sqrt(t1 = 1. - t*t);
	return *(new ADvar1g(asin(t), d1, t*d1/t1, &v));
	}

 ADvari&
asinh(const ADvari &v) {
	double d1, t, t1, t2, td;
	t = v.Val;
	t1 = sqrt(t2 = t*t + 1.);
	d1 = 1. / t1;
	td = 1.;
	if (t < 0.)
		td = -1.;
	return *(new ADvar1g(td*log(t*td + t1), d1, -(t/t2)*d1, &v));
	}

 ADvari&
atan(const ADvari &v) {
	double d1, t;
	t = v.Val;
	d1 = 1./(1. + t*t);
	return *(new ADvar1g(atan(t), d1, -(t+t)*d1*d1, &v));
	}

 ADvari&
atanh(const ADvari &v) {
	double t = v.Val;
	double d1 = 1./(1. - t*t);
	return *(new ADvar1g(0.5*log((1.+t)/(1.-t)), d1, (t+t)*d1*d1, &v));
	}

 ADvari&
max(const ADvari &L, const ADvari &R) {
	const ADvari &x = L.Val >= R.Val ? L : R;
	return *(new ADvar1(Hv_copy, x.Val, &CADcontext::One, &x));
	}

 ADvari&
max(double L, const ADvari &R) {
	if (L >= R.Val)
		return *(new ADvari(Hv_const, L));
	return *(new ADvar1(Hv_copy, R.Val, &CADcontext::One, &R));
	}

 ADvari&
max(const ADvari &L, double R) {
	if (L.Val >= R)
		return *(new ADvar1(Hv_copy, L.Val, &CADcontext::One, &L));
	return *(new ADvari(Hv_const, R));
	}

 ADvari&
min(const ADvari &L, const ADvari &R) {
	const ADvari &x = L.Val <= R.Val ? L : R;
	return *(new ADvar1(Hv_copy, x.Val, &CADcontext::One, &x));
	}

 ADvari&
min(double L, const ADvari &R) {
	if (L <= R.Val)
		return *(new ADvari(Hv_const, L));
	return *(new ADvar1(Hv_copy, R.Val, &CADcontext::One, &R));
	}

 ADvari&
min(const ADvari &L, double R) {
	if (L.Val <= R)
		return *(new ADvar1(Hv_copy, L.Val, &CADcontext::One, &L));
	return *(new ADvari(Hv_const, R));
	}

 ADvari&
atan2(const ADvari &L, const ADvari &R) {
	double x = L.Val, y = R.Val;
	double x2 = x*x, y2 = y*y;
	double t = 1./(x2 + y2);
	double t2 = t*t;
	double R2 = 2.*t2*x*y;
	return *(new ADvar2g(atan2(x,y), y*t, -x*t, -R2, t2*(x2 - y2), R2, &L, &R));
	}

 ADvari&
atan2(double x, const ADvari &R) {
	double y = R.Val;
	double x2 = x*x, y2 = y*y;
	double t = 1./(x2 + y2);
	return *(new ADvar1g(atan2(x,y), -x*t, 2.*t*t*x*y, &R));
	}

 ADvari&
atan2(const ADvari &L, double y) {
	double x = L.Val;
	double x2 = x*x, y2 = y*y;
	double t = 1./(x2 + y2);
	return *(new ADvar1g(atan2(x,y), y*t, -2.*t*t*x*y, &L));
	}

 ADvari&
cos(const ADvari &v) {
	double t = cos(v.Val);
	return *(new ADvar1g(t, -sin(v.Val), -t, &v));
	}

 ADvari&
cosh(const ADvari &v) {
	double t = cosh(v.Val);
	return *(new ADvar1g(t, sinh(v.Val), t, &v));
	}

 ADvari&
exp(const ADvari &v) {
	double t = exp(v.Val);
	return *(new ADvar1g(t, t, t, &v));
	}

 ADvari&
log(const ADvari &v) {
	double x = v.Val;
	double d1 = 1. / x;
	return *(new ADvar1g(log(x), d1, -d1*d1, &v));
	}

 ADvari&
log10(const ADvari &v) {
	static double num = 1. / log(10.);
	double x = v.Val, t = 1. / x;
	double d1 = num * t;
	return *(new ADvar1g(log10(x), d1, -d1*t, &v));
	}

 ADvari&
pow(const ADvari &L, const ADvari &R) {
	double x = L.Val, y = R.Val, t = pow(x,y);
	double xym1 = t/x;
	double xlog = log(x);
	double dx = y*xym1;
	double dy = t * xlog;
	return *(new ADvar2g(t, dx, dy, (y-1.)*dx/x, xym1*(1. + y*xlog), dy*xlog, &L, &R));
	}

 ADvari&
pow(double x, const ADvari &R) {
	double y = R.Val, t = pow(x,y);
	double xlog = log(x);
	double dy = t * xlog;
	return *(new ADvar1g(t, dy, dy*xlog, &R));
	}

 ADvari&
pow(const ADvari &L, double y) {
	double x = L.Val, t = pow(x,y);
	double dx = y*t/x;
	return *(new ADvar1g(t, dx, (y-1.)*dx/x, &L));
	}

 ADvari&
sin(const ADvari &v) {
	double t = sin(v.Val);
	return *(new ADvar1g(t, cos(v.Val), -t, &v));
	}

 ADvari&
sinh(const ADvari &v) {
	double t = sinh(v.Val);
	return *(new ADvar1g(t, cosh(v.Val), t, &v));
	}

 ADvari&
sqrt(const ADvari &v) {
	double t = sqrt(v.Val);
	double d1 = 0.5/t;
	return *(new ADvar1g(t, d1, -0.5*d1/v.Val, &v));
	}

 ADvari&
tan(const ADvari &v) {
	double d1, rv, t;
	rv = tan(v.Val);
	t = 1. / cos(v.Val);
	d1 = t*t;
	return *(new ADvar1g(rv, d1, (rv+rv)*d1, &v));
	}

 ADvari&
tanh(const ADvari &v) {
	double d1, rv, t;
	rv = tanh(v.Val);
	t = 1. / cosh(v.Val);
	d1 = t*t;
	return *(new ADvar1g(rv, d1, -(rv+rv)*d1, &v));
	}

 ADvari&
fabs(const ADvari &v) {	// "fabs" is not the best choice of name,
			// but this name is used at Sandia.
	double t, p;
	p = 1;
	if ((t = v.Val) < 0) {
		t = -t;
		p = -p;
		}
	return *(new ADvar1g(t, p, 0., &v));
	}

 ADvari&
ADf1(double f, double g, double h, const ADvari &x) {
	return *(new ADvar1g(f, g, h, &x));
	}

 ADvari&
ADf2(double f, double gx, double gy, double hxx, double hxy, double hyy,
		const ADvari &x, const ADvari &y) {
	return *(new ADvar2g(f, gx, gy, hxx, hxy, hyy, &x, &y));
	}

 ADvarn::ADvarn(double val1, int n1, const ADvar *x, const double *g, const double *h):
		ADvari(Hv_nary,val1), n(n1)
{
	Derp *d1, *dlast;
	double *a1;
	int i, nh;

	a1 = G = (double*)ADvari::adc.Memalloc(n1*sizeof(*g));
	d1 = D =  (Derp*)ADvari::adc.Memalloc(n1*sizeof(Derp));
	dlast = Derp::LastDerp;
	for(i = 0; i < n1; i++, d1++) {
		d1->next = dlast;
		dlast = d1;
		a1[i] = g[i];
		d1->a = &a1[i];
		d1->b = this;
		d1->c = x[i].cv;
		}
	Derp::LastDerp = dlast;
	nh = (n1*(n1+1)) >> 1;
	a1 = H = (double*)ADvari::adc.Memalloc(nh * sizeof(*h));
	for(i = 0; i < nh; i++)
		a1[i] = h[i];
	}

 ADvari&
ADfn(double f, int n, const ADvar *x, const double *g, const double *h) {
	return *(new ADvarn(f, n, x, g, h));
	}

 void
ADcontext::Hvprod(int n, ADvar **x, double *v, double *hv)
{
	ADvari *a, *aL, *aR, **ap, **ape;
	ADvari_block *b, *b0;
	Derp *d;
	double aO, adO, *g, *h, *h0, t, tL, tR;
	int i, j, k, m;
	for(i = 0; i < n; i++) {
		a = x[i]->cv;
		a->dO = v[i];
		a->aO = a->adO = 0.;
		}
	ADvari::adc.Aibusy->limit = ADvari::adc.Ainext;
	a = 0;
	for(b0 = 0, b = &ADvari::adc.AiFirst; b; b0 = b, b = b->next) {
		ap = b->pADvari;
		ape = b->limit;
		while(ap < ape) {
			a = *ap++;
			a->aO = a->adO = 0.;
			switch(a->opclass) {
			 case Hv_copy:
				a->dO = ((ADvar1*)a)->d.c->dO;
				break;
			 case Hv_binary:
				a->dO =   ((ADvar2g*)a)->pL * ((ADvar2g*)a)->dL.c->dO
					+ ((ADvar2g*)a)->pR * ((ADvar2g*)a)->dR.c->dO;
				break;
			 case Hv_unary:
				a->dO = ((ADvar1g*)a)->pL * ((ADvar1g*)a)->d.c->dO;
				break;
			 case Hv_negate:
				a->dO = -((ADvar1*)a)->d.c->dO;
				break;
			 case Hv_plusLR:
				a->dO = ((ADvar2*)a)->dL.c->dO + ((ADvar2*)a)->dR.c->dO;
				break;
			 case Hv_minusLR:
				a->dO = ((ADvar2*)a)->dL.c->dO - ((ADvar2*)a)->dR.c->dO;
				break;
			 case Hv_timesL:
				a->dO = ((ADvar1s*)a)->pL * ((ADvar1s*)a)->d.c->dO;
				break;
			 case Hv_timesLR:
				a->dO =   ((ADvar2*)a)->dR.c->Val * ((ADvar2*)a)->dL.c->dO
					+ ((ADvar2*)a)->dL.c->Val * ((ADvar2*)a)->dR.c->dO;
				break;
			 case Hv_quotLR:
				a->dO =   ((ADvar2q*)a)->pL * ((ADvar2q*)a)->dL.c->dO
					+ ((ADvar2q*)a)->pR * ((ADvar2q*)a)->dR.c->dO;
			 case Hv_const:	// This case does not arise; the intent here
					// is to eliminate a red-herring compiler warning.
				break;
			 case Hv_nary:
				d = ((ADvarn*)a)->D;
				m = ((ADvarn*)a)->n;
				g = ((ADvarn*)a)->G;
				t = 0.;
				for(i = 0; i < m; i++)
					t += g[i] * d[i].c->dO;
				a->dO = t;
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
				aL = ((ADvar1*)a)->d.c;
				aL->aO += aO;
				aL->adO += adO;
				break;
			 case Hv_binary:
				aL = ((ADvar2g*)a)->dL.c;
				aR = ((ADvar2g*)a)->dR.c;
				tL = adO*aL->dO;
				tR = adO*aR->dO;
				aL->aO += aO*((ADvar2g*)a)->pL
					+ tL*((ADvar2g*)a)->pL2
					+ tR*((ADvar2g*)a)->pLR;
				aR->aO += aO*((ADvar2g*)a)->pR
					+ tL*((ADvar2g*)a)->pLR
					+ tR*((ADvar2g*)a)->pR2;
				aL->adO += adO * ((ADvar2g*)a)->pL;
				aR->adO += adO * ((ADvar2g*)a)->pR;
				break;
			 case Hv_unary:
				aL = ((ADvar1g*)a)->d.c;
				aL->aO += aO * ((ADvar1g*)a)->pL
					+ adO * aL->dO * ((ADvar1g*)a)->pL2;
				aL->adO += adO * ((ADvar1g*)a)->pL;
				break;
			 case Hv_negate:
				aL = ((ADvar1*)a)->d.c;
				aL->aO -= aO;
				aL->adO -= adO;
				break;
			 case Hv_plusLR:
				aL = ((ADvar2*)a)->dL.c;
				aR = ((ADvar2*)a)->dR.c;
				aL->aO += aO;
				aL->adO += adO;
				aR->aO += aO;
				aR->adO += adO;
				break;
			 case Hv_minusLR:
				aL = ((ADvar2*)a)->dL.c;
				aR = ((ADvar2*)a)->dR.c;
				aL->aO += aO;
				aL->adO += adO;
				aR->aO -= aO;
				aR->adO -= adO;
				break;
			 case Hv_timesL:
				aL = ((ADvar1s*)a)->d.c;
				aL->aO += aO * (tL = ((ADvar1s*)a)->pL);
				aL->adO += adO * tL;
				break;
			 case Hv_timesLR:
				aL = ((ADvar2*)a)->dL.c;
				aR = ((ADvar2*)a)->dR.c;
				aL->aO += aO * (tL = aR->Val) + adO*aR->dO;
				aR->aO += aO * (tR = aL->Val) + adO*aL->dO;
				aL->adO += adO * tL;
				aR->adO += adO * tR;
				break;
			 case Hv_quotLR:
				aL = ((ADvar2q*)a)->dL.c;
				aR = ((ADvar2q*)a)->dR.c;
				tL = adO*aL->dO;
				tR = adO*aR->dO;
				aL->aO += aO*((ADvar2q*)a)->pL
					+ tR*((ADvar2q*)a)->pLR;
				aR->aO += aO*((ADvar2q*)a)->pR
					+ tL*((ADvar2q*)a)->pLR
					+ tR*((ADvar2q*)a)->pR2;
				aL->adO += adO * ((ADvar2q*)a)->pL;
				aR->adO += adO * ((ADvar2q*)a)->pR;
			 case Hv_const:	// This case does not arise; the intent here
					// is to eliminate a red-herring compiler warning.
				break;
			 case Hv_nary:
				d  = ((ADvarn*)a)->D;
				m  = ((ADvarn*)a)->n;
				g  = ((ADvarn*)a)->G;
				h0 = ((ADvarn*)a)->H;
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
			 }
			}
		}
	for(i = 0; i < n; i++) {
		a = x[i]->cv;
		a->dO = 0.;
		hv[i] = a->aO;
		}
	}

#ifndef SACADO_NO_NAMESPACE
} // namespace Rad2d
} // namespace Sacado
#endif
