// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IQR_GMRES_TOOLS_H
#define IQR_GMRES_TOOLS_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

#include <iostream>
#include <shylu_internal_gmres.h>

#include <Epetra_Operator.h>

/*
    This typename provides some tools to handle with gmres.
    In particular it stores informations of a computed gmres which are useful
    for the followings gmres, such as initial guess or preconditioner.
    It is even possible to have nested right preconditioners using
    the pointer to shylu_internal_gmres_tools *P.
    Data stored are:
        v : array of basis vectors (for "AP^-1 y" and  "Ax")
        w : for future use: same as v but basis of "x" and "P^-1y".
           for the moment v and w are he same.
        Q,R : QR-factorization of the hessenberg matrix
*/

namespace IQR
{

template < typename Map, typename MultiVector, typename LocalMatrix, typename LocalVector >
class GMRESManager
{
public:
    // Constructor and destructor

    // flex != 0 means that we will use flexible gmres and
    // we have to use w
    GMRESManager(const Map& map, const int &rest, const int flex = 0, const bool scaling = true);
    ~GMRESManager();

    // Public methods
    int initialGuess(const  MultiVector &b,
		     MultiVector& X) const;
    // solve min_{x\inV} || Ax -b||
    int ApplyInverse(const MultiVector &b, MultiVector& X) const;
    // delete P and restart the gmres.
    int start();

    // Public data
    int restart;
    int  m;
    int isFlex;
    bool doScaling;
    bool isFirst;

    // Q
    LocalVector cs;
    LocalVector sn;
    // R , with abuse of notation
    LocalMatrix H;
    // v basis of AP^-1 y = Ax and of y
    MultiVector *v;
    // w basis of x = P^-1 y
    MultiVector *w;
    GMRESManager *P; // preconditioner(s)
    Epetra_Operator* P2;

    // Related to parallel operation
    const Map& map_;
};

template < typename Map, typename MultiVector, typename LocalMatrix, typename LocalVector >
GMRESManager< Map, MultiVector, LocalMatrix, LocalVector >::GMRESManager(const Map& map, const int &rest, const int flex, const bool scaling)
    : restart(rest),
      m(-1),
      isFlex(flex),
      doScaling(scaling),
      cs(rest + 1, 0.0),
      sn(rest + 1, 0.0),
      H(rest + 1, LocalVector(rest, 0.0)),
      P(0),
      map_(map)
{
    if (rest < 0) {
        std::cout << "shylu_internal_gmres_tools: restart must be positive" << std::endl;
        exit(1);
    }

    // After you get it working, try removing these...
    v = new MultiVector(map_, rest + 1, true);

    if (isFlex == 0) {
        w = v;
    } else {
        w = new MultiVector(map_, rest + 1, true);
    }
}

template < typename Map, typename MultiVector, typename LocalMatrix, typename LocalVector >
GMRESManager< Map, MultiVector, LocalMatrix, LocalVector >::~GMRESManager()
{
    delete v;
    if (isFlex != 0) {
        delete w;
    }
    //delete P;
}

template < typename Map, typename MultiVector, typename LocalMatrix, typename LocalVector >
int GMRESManager< Map, MultiVector, LocalMatrix, LocalVector >::initialGuess(const  MultiVector &b,
										      MultiVector& X) const
{
    int mm(m);

    if (m > restart) {
        std::cout << "GMRES_TOOLS.initialGuess: mm>m, not expected using  m"
                  << std::endl;
        mm = restart;
    }

    if (P != 0 && isFlex == 0) {// if there is a (right) preconditioner,  use it
        this->ApplyInverse(b, X);
        // x= P->initialGuess(x);
        return 0;
    }

    if (mm <= -1 ) {
        X = b;
        return 0;
    }

    //std::cout << "GMRES_TOOLS.initialGuess: guess with m = " << mm
    //          << std::endl;
    int i;

    LocalVector bp(mm + 1, 0.0);
    // bp = Q_{n+1}^T b and x = x - bp Q_{n+1} 
    for (i = 0; i < mm + 1; i++) {
        b.Dot(*(*v)(i), &bp[i]);
    }

    // applying Q: bp = \Omega_n bp. At this time, bp = \Omega_n Q_{n+1}^T b
    for (i = 0; i < mm; i++) {
        ApplyPlaneRotation(bp[i], bp[i + 1], cs[i], sn[i]);
    }

    // solving R_{n} y_{n} = bp and setting x = x + Q_{n} y_{n} 
    Update(X, mm-1, H, bp, *w);
    //Update(y, i, G.H, , *(G.v));

    return 0;
}


template < typename Map, typename MultiVector, typename LocalMatrix, typename LocalVector >
int GMRESManager< Map, MultiVector, LocalMatrix, LocalVector >::start()
{
    m = -1;

    if (P != 0 ) {
        delete P;
        P = 0;
    }

    return 0;
}


template < typename Map, typename MultiVector, typename LocalMatrix, typename LocalVector >
int GMRESManager< Map, MultiVector,
                  LocalMatrix, LocalVector >::ApplyInverse(const  MultiVector &b,
                                                           MultiVector& X) const
{
    int mm(m);

    if (m > restart) {
        std::cout << "GMRES_TOOLS.solve: mm>m, not expected using  m"
                  << std::endl;
        mm = restart;
    }

    if (mm <= -1 ) {
        std::cout << "GMRES_TOOLS.solve: no preconditioner" << std::endl;
        X = b;
        return 1;
    }

    //std::cout << "GMRES_TOOLS.solve: prec with m = " << mm << std::endl;
//    int myPID = map_.Comm().MyPID();
//    if (! myPID) std::cout << "Apply inverse " << id << "\n";

    int i;
    LocalVector bp(mm + 1, 0.0);
    MultiVector x(map_, 1, false);
    x = b;

    // bp = Q_{n+1}^T b and x = x - bp Q_{n+1} 
    for (i = 0; i < mm + 1; i++) { // using " V is an orthogonal basis
        MultiVector* vi = (*v)(i);
        b.Dot(*vi, &bp[i]);
        //std::cout << "bp: " << bp[i] << std::endl;
        //    x-= bp(i) * w[i]; // computing the orthogonal component of b
        x.Update(-bp[i], *vi, 1.0);
    }

    // applying Q: bp = \Omega_n bp. At this time, bp = \Omega_n Q_{n+1}^T b
    for (i = 0; i < mm; i++) {
        ApplyPlaneRotation(bp[i], bp[i + 1], cs[i], sn[i]);
    }

    typedef typename  LocalVector::value_type Scalar;
    Scalar scaling(0);
    for (i = 0; i < mm; i++)
    {
        scaling += H[i][i];
    }
    scaling = static_cast<Scalar>(mm) / scaling;

    LocalVector bq(mm + 1, 0.0);
    bq[mm] = bp[mm];

    // bq = \omega_n^T (0...0 bp(n+1))^T
    for (i = mm-1; i >= 0; i--) {
    	ApplyPlaneRotation(bq[i], bq[i + 1], cs[i], -sn[i]);
//        bq[i] = - sn[i] * bq[i + 1];
//        bq[i + 1] = cs[i] * bq[i + 1];
    }

    // x = x + Q_n bq
    for (i = 0; i <= mm; i++) {
        MultiVector* vi = (*v)(i);
        //    x+= bq(i)*w[i];
        x.Update(bq[i], *vi, 1.0);
    }

    if (doScaling) {
        x.Scale(scaling);
    }

    // solving R_{n} y_{n} = bp and setting x = x + Q_{n} y_{n} 
    Update(x, mm-1, H, bp, *w); // summing to x the solution of Rx=bp

    if (isFirst) {
    	if (P2 !=0 && isFlex == 0) {
    		P2->ApplyInverse(x, X);
    	} else {
    		X = x;
    	}
    } else {
		if (P != 0 && isFlex == 0) {
			P->ApplyInverse(x, X);
		} else {
			X = x;
		}
    }
    return 0;
}

} // namespace IQR

#endif // IQR_GMRES_TOOLS_H
