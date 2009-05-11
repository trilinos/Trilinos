// ----------------------------------------------------------------------------
//	Organization:	Natural Resource Ecology Laboratory
//			Colorado State University, Fort Collins, CO 80523 USA
//			www.nrel.colostate.edu
//	Project:  Century Soil Organic Matter Model
//	File:	  Random.cpp
//	Class:	  nrel::stats::Random
//
//	Description:
//	Random number generator using the code from zufall.h/.c.
//	See header file for more details.
// ----------------------------------------------------------------------------
//	Author:	Tom Hilinski, tom.hilinski@colostate.edu, Mar07
//	History: See header file.
// ----------------------------------------------------------------------------

#include "Random.h"
#include <cmath>
#include <stdexcept>
#include <string>

using std::exp;
using std::cos;
using std::sin;
using std::sqrt;
using std::log;


nrel::stats::Random::Random (
	int const useSeed /* = 1802 */ )
{
	InitializeSeedBuffer (useSeed);
}

// ----------------------------------------------------------------------------
//	InitializeSeedBuffer (formerly zufalli)
//	Generates initial seed buffer by linear congruential
//	method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
//	Variable seed should be 0 < seed < 31328 ( numeric_limits<int> ? )
// ----------------------------------------------------------------------------
void nrel::stats::Random::InitializeSeedBuffer (
	int const seed)
{
    int const ij = ( seed != 0 ? seed : 1802 );
    int i = ij / 177 % 177 + 2;
    int j = ij % 177 + 2;
    int k = 9373 / 169 % 178 + 1;
    int l = 9373 % 169;
    for (int ii = 0; ii < klotz0_1.size; ++ii)
    {
	double s = 0.0;
	double t = 0.5;
	for (int jj = 1; jj <= 24; ++jj)
	{
	    int m = i * j % 179 * k % 179;
	    i = j;
	    j = k;
	    k = m;
	    l = (l * 53 + 1) % 169;
	    if (l * m % 64 >= 32)
		s += t;
	    t *= 0.5;
	}
	klotz0_1.buff[ii] = s;
    }
}

// ----------------------------------------------------------------------------
//	Uniform (formerly zufall)
//	portable lagged Fibonacci series uniform random number
//	generator with "lags" -273 und -607:
//	W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
// ----------------------------------------------------------------------------
void nrel::stats::Random::Uniform (
	int const n,
	std::vector<double> & a)
{
    if (n <= 0)
	return;

    a.resize (n);
    a.assign (n, 0.0);

    int left, bptr, aptr0;
    double t;
    int vl, k273, k607, kptr;
    int aptr = 0;
    int nn = n;

  L1:

    if (nn <= 0)
	return;

    /* factor nn = q*607 + r */

    int q = (nn - 1) / klotz0_1.size;
    left = klotz0_1.size - klotz0_1.ptr;

    if (q <= 1)
    {

	/* only one or fewer full segments */

	if (nn < left)
	{
            kptr = klotz0_1.ptr;
	    for (int i = 0; i < nn; ++i)
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    klotz0_1.ptr += nn;
	    return;
	}
	else
	{
            kptr = klotz0_1.ptr;
#pragma _CRI ivdep
	    for (int i = 0; i < left; ++i)
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    klotz0_1.ptr = 0;
	    aptr += left;
	    nn -= left;
	    /*  buff -> buff case */
	    vl = 273;
	    k273 = 334;
	    k607 = 0;
	    for (int k = 0; k < 3; ++k)
	    {
#pragma _CRI ivdep
		for (int i = 0; i < vl; ++i)
		{
		   t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		   klotz0_1.buff[k607+i] = t - (double) ((int) t);
		}
		k607 += vl;
		k273 += vl;
		vl = 167;
		if (k == 0)
		{
		    k273 = 0;
		}
	    }
	    goto L1;
	}
    }
    else
    {

	/* more than 1 full segment */

        kptr = klotz0_1.ptr;
#pragma _CRI ivdep
	for (int i = 0; i < left; ++i)
	    a[i + aptr] = klotz0_1.buff[kptr + i];
	nn -= left;
	klotz0_1.ptr = 0;
	aptr += left;

	/* buff -> a(aptr0) */

	vl = 273;
	k273 = 334;
	k607 = 0;
	for (int k = 0; k < 3; ++k)
	{
	    if (k == 0)
	    {
#pragma _CRI ivdep
		for (int i = 0; i < vl; ++i)
		{
		    t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		    a[aptr + i] = t - (double) ((int) t);
		}
		k273 = aptr;
		k607 += vl;
		aptr += vl;
		vl = 167;
	    }
	    else
	    {
#pragma _CRI ivdep
		for (int i = 0; i < vl; ++i)
		{
		    t = a[k273 + i] + klotz0_1.buff[k607 + i];
		    a[aptr + i] = t - (double) ((int) t);
		}
		k607 += vl;
		k273 += vl;
		aptr += vl;
	    }
	}
	nn += -klotz0_1.size;

	/* a(aptr-607) -> a(aptr) for last of the q-1 segments */

	aptr0 = aptr - klotz0_1.size;
	vl = klotz0_1.size;

	for (int qq = 0; qq < q-2; ++qq)
	{
	    k273 = aptr0 + 334;
#pragma _CRI ivdep
	    for (int i = 0; i < vl; ++i)
	    {
		t = a[k273 + i] + a[aptr0 + i];
		a[aptr + i] = t - (double) ((int) t);
	    }
	    nn += -klotz0_1.size;
	    aptr += vl;
	    aptr0 += vl;
	}

	/* a(aptr0) -> buff, last segment before residual */

	vl = 273;
	k273 = aptr0 + 334;
	k607 = aptr0;
	bptr = 0;
	for (int k = 0; k < 3; ++k)
	{
	    if (k == 0)
	    {
#pragma _CRI ivdep
		for (int i = 0; i < vl; ++i)
		{
		    t = a[k273 + i] + a[k607 + i];
		    klotz0_1.buff[bptr + i] = t - (double) ((int) t);
		}
		k273 = 0;
		k607 += vl;
		bptr += vl;
		vl = 167;
	    }
	    else
	    {
#pragma _CRI ivdep
		for (int i = 0; i < vl; ++i)
		{
		    t = klotz0_1.buff[k273 + i] + a[k607 + i];
		    klotz0_1.buff[bptr + i] = t - (double) ((int) t);
		}
		k607 += vl;
		k273 += vl;
		bptr += vl;
	    }
	}
	goto L1;
    }
}

// ----------------------------------------------------------------------------
//	Normal (formerly normalen)
//	Box-Muller method for Gaussian random numbers.
// ----------------------------------------------------------------------------
void nrel::stats::Random::Normal (
	int const n,
	std::vector<double> & x)
{
    if (n <= 0)
	return;

    x.resize (n);
    x.assign (n, 0.0);

    if (klotz1_1.first == 0)
    {
	normal00();
	klotz1_1.first = 1;
    }
    int ptr = 0;
    int nn = n;

  L1:
    int const left = klotz1_1.size - klotz1_1.xptr;
    int const kptr = klotz1_1.xptr;
    if (nn < left)
    {
#pragma _CRI ivdep
	for (int i = 0; i < nn; ++i)
	    x[i + ptr] = klotz1_1.xbuff[kptr + i];
	klotz1_1.xptr += nn;
	return;
    }
    else
    {
#pragma _CRI ivdep
	for (int i = 0; i < left; ++i)
	    x[i + ptr] = klotz1_1.xbuff[kptr + i];
	klotz1_1.xptr = 0;
	ptr += left;
	nn -= left;
	normal00();
	goto L1;
    }
}

// ----------------------------------------------------------------------------
//	Poisson (formerly fische)
//	Poisson generator for distribution function of p's:
//	q(mu,p) = exp(-mu) mu**p/p!
// ----------------------------------------------------------------------------
void nrel::stats::Random::Poisson (
	int const n,
	double const mu,
	std::vector<int> & p)
{
    if (n <= 0)
	return;

    p.resize (n);
    p.assign (n, 0);

    double const pmu = exp(-mu);
    int p0 = 0;
    int nsegs = (n - 1) / klotz1_1.size;
    int left = n - (nsegs << 10);
    ++nsegs;
    int nl0 = left;

    std::vector<int> indx (klotz1_1.size);
    std::vector<double> q (klotz1_1.size);
    std::vector<double> u (klotz1_1.size);

    for (int k = 0; k < nsegs; ++k)
    {
	for (int i = 0; i < left; ++i)
	{
	    indx[i] = i;
	    p[p0 + i] = 0;
	    q[i] = 1.;
	}

	/* Begin iterative loop on segment of p's */
	do
	{
	    Uniform (left, u);	// Get the needed uniforms
	    int jj = 0;
	    for (int i = 0; i < left; ++i)
	    {
		int const ii = indx[i];
		double const q0 = q[ii] * u[i];
		q[ii] = q0;
		if (q0 > pmu)
		{
		    indx[jj++] = ii;
		    ++p[p0 + ii];
		}
	    }
	    left = jj;	// any left in this segment?
	}
	while (left > 0);

	p0  += nl0;
	nl0  = klotz1_1.size;
	left = klotz1_1.size;
    }
}

// ----------------------------------------------------------------------------
//	zufallsv
//	saves common blocks klotz0, containing seeds and
//	pointer to position in seed block. IMPORTANT: svblk must be
//	dimensioned at least 608 in driver. The entire contents
//	of klotz0 (pointer in buff, and buff) must be saved.
// ----------------------------------------------------------------------------
void nrel::stats::Random::zufallsv (
	std::vector<double> & saveBuffer)	//  size = klotz0_1_.size + 1
{
    saveBuffer.resize (klotz0_1_::size + 1);
    saveBuffer[0] = (double) klotz0_1.ptr;
#pragma _CRI ivdep
    for (short i = 0; i < klotz0_1.size; ++i)
	saveBuffer[i + 1] = klotz0_1.buff[i];
}

// ----------------------------------------------------------------------------
//	zufallrs
//	restores common block klotz0, containing seeds and pointer
//	to position in seed block. IMPORTANT: saveBuffer must be
//	dimensioned at least 608 in driver. The entire contents
//	of klotz0 must be restored.
// ----------------------------------------------------------------------------
void nrel::stats::Random::zufallrs (
	std::vector<double> const & saveBuffer)	//  size = klotz0_1_.size + 1
{
    if ( saveBuffer.size() != klotz0_1_::size + 1 )
    {
	std::string msg;
	msg = "nrel::stats::Random::zufallrs, restore of uninitialized block.";
	throw std::runtime_error (msg);
    }
    klotz0_1.ptr = (int) saveBuffer[0];
#pragma _CRI ivdep
    for (short i = 0; i < klotz0_1.size; ++i)
	klotz0_1.buff[i] = saveBuffer[i + 1];
}

// ----------------------------------------------------------------------------
//	normalsv
//	save zufall block klotz0
//	IMPORTANT: svbox must be dimensioned at
//	least 1634 in driver. The entire contents of blocks
//	klotz0 (via zufallsv) and klotz1 must be saved.
// ----------------------------------------------------------------------------
void nrel::stats::Random::normalsv (
	std::vector<double> & saveBuffer)		// size = 1634
{
    if (klotz1_1.first == 0)
    {
	std::string msg;
	msg = "nrel::stats::Random::normalsv, save of uninitialized block.";
	throw std::runtime_error (msg);
    }

    saveBuffer.resize (1634);
    zufallsv (saveBuffer);

    saveBuffer[klotz0_1.size + 1] = (double) klotz1_1.first;	// [608]
    saveBuffer[klotz0_1.size + 2] = (double) klotz1_1.xptr;	// [609]
    int const k = klotz0_1.size + 3;				// 610
#pragma _CRI ivdep
    for (short i = 0; i < klotz1_1.size; ++i)
	saveBuffer[i + k] = klotz1_1.xbuff[i];
}

// ----------------------------------------------------------------------------
//	normalrs
//	restore zufall blocks klotz0 and klotz1
//	IMPORTANT: saveBuffer must be dimensioned at
//	least 1634 in driver. The entire contents
//	of klotz0 and klotz1 must be restored.
// ----------------------------------------------------------------------------
void nrel::stats::Random::normalrs (
	std::vector<double> const & saveBuffer)		// size = 1634
{
    zufallrs (saveBuffer);

    klotz1_1.first = (int) saveBuffer[klotz0_1.size + 1];	// [608]
    if (klotz1_1.first == 0)
    {
	std::string msg;
	msg = "nrel::stats::Random::normalrs, restore of uninitialized block.";
	throw std::runtime_error (msg);
    }

    klotz1_1.xptr = (int) saveBuffer[klotz0_1.size + 2];	// [609]
    int const k = klotz0_1.size + 3;				// 610
#pragma _CRI ivdep
    for (short i = 0; i < klotz1_1.size; ++i)
	klotz1_1.xbuff[i] = saveBuffer[i + k];
}

// ----------------------------------------------------------------------------
//	normal00
// ----------------------------------------------------------------------------
void nrel::stats::Random::normal00 ()
{
    /* Builtin functions */
    /* double cos(), sin(), log(), sqrt(); */

    double const twopi = 6.2831853071795862;

    Uniform (klotz1_1.size, klotz1_1.xbuff);

#pragma _CRI ivdep
    for (short i = 0; i < klotz1_1.size - 1; i += 2)
    {
	double const r1 = twopi * klotz1_1.xbuff[i];
	double const t1 = cos(r1);
	double const t2 = sin(r1);
	double const r2 = sqrt(-2.*(log(1. - klotz1_1.xbuff[i+1])));
	klotz1_1.xbuff[i]   = t1 * r2;
	klotz1_1.xbuff[i+1] = t2 * r2;
    }
}

//--- end of file ---
