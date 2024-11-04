//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

#ifndef IQR_GMRES_H
#define IQR_GMRES_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

#include <cmath>
#include <iostream>

#include <shylu_internal_gmres_tools.h>

namespace IQR
{

struct IdPreconditioner
{
    void ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
    {
        Y = X;
    }
};

//! Generate i-th Given rotation Gi. Note Qn= G0 * ... * Gn
template <typename Scalar>
void GeneratePlaneRotation(const Scalar &dx, const Scalar &dy, Scalar &cs,
                           Scalar &sn)
{
    if (dy == 0.0) {
        cs = 1.0;
        sn = 0.0;
    } else if (std::abs(dy) > std::abs(dx)) {
        Scalar temp = dx / dy;
        sn = 1.0 / std::sqrt( 1.0 + temp * temp );
        cs = temp * sn;
    } else {
        Scalar temp = dy / dx;
        cs = 1.0 / std::sqrt( 1.0 + temp * temp );
        sn = temp * cs;
    }
}

//! Apply i-th Given rotation Gi. Note Qn= G0 * ... * Gn
//! Gi = \left( \begin{array}[ccc] I 0 0
template <typename Scalar>
void ApplyPlaneRotation(Scalar &dx, Scalar &dy, const Scalar &cs,
                        const Scalar &sn)
{
    Scalar temp  =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}


//! solving R_{k+1} y_{k+1} = bp and setting x = x + Q_{k+1} y_{k+1}
template < typename LocalMatrix, typename LocalVector, typename MultiVector >
void Update(MultiVector &x, const int k, const LocalMatrix &h,
            const LocalVector &s, const MultiVector &v)
{
    LocalVector y(s);

    // Backsolve:
    for (int i = k; i >= 0; i--) {
        y[i] /= h[i][i];
        for (int j = i - 1; j >= 0; j--) {
            y[j] -= h[j][i] * y[i];
        }
    }

    for (int j = 0; j <= k; j++) {
        x.Update(y[j], *v(j), 1.0);
    }
}

template < typename Operator, typename MultiVector, typename LeftPrec,
           typename RightPrec, typename GMRESManager, typename LocalVector,
           typename Scalar>
int GMRES(const Operator &A, MultiVector &x, const MultiVector &b,
          LeftPrec *L, RightPrec *M, GMRESManager &G, int &max_iter,
          Scalar &tol)
{
    // Storing a reference to the parallel map
  //auto& b.Map() = b.Map();
  //int myPID = b.Map().Comm().MyPID();

    Scalar resid;
    int i(0), j(1), k(0);
    // The following initial guess was wrong! Indeed it was altering the whole QR factorization
    // initial guess from previous solves : compute x
    //M->ApplyInverse(b, x);

    LocalVector s(G.restart + 1);
    MultiVector w(b.Map(), 1, true);

    Scalar normb;
    L->ApplyInverse(b, w);
    w.Norm2(&normb);

    MultiVector t(b.Map(), 1, true);
    A.Apply(x, t);
    w.Update(1.0, b, -1.0, t, 0.0);

    MultiVector r(b.Map(), 1, true);
    L->ApplyInverse(w, r);
    Scalar beta;
    r.Norm2(&beta);

    if (normb == 0.0) {
        normb = 1;
    }

    if ((resid = beta / normb) <= tol) { // qui ho migliorato
        tol = resid;
        max_iter = 0;
        return 0;
    }

    while (j <= max_iter) {
        MultiVector* v0 = (*G.v)(0);
        v0->Update(1.0 / beta, r, 0.0);
        s.assign(G.restart + 1, 0.0);
        s[0] = beta;

        for (i = 0; i < G.restart && j <= max_iter; i++, j++) {
            M->ApplyInverse(*((*G.v)(i)), t);
            A.Apply(t, r);
            L->ApplyInverse(r, w);

            for (k = 0; k <= i; k++) {
                MultiVector* vk = (*G.v)(k);
                w.Dot(*vk, &(G.H[k][i]));
                w.Update(-G.H[k][i], *vk, 1.0);
            }
            w.Norm2(&(G.H[i + 1][i]));
            MultiVector* vi1 = (*G.v)(i + 1);
            // Set (*G.v)(i + 1) to w/||w||
            vi1->Scale(1.0 / G.H[i + 1][i], w);

            for (k = 0; k < i; k++) {
                ApplyPlaneRotation(G.H[k][i], G.H[k + 1][i], G.cs[k], G.sn[k]);
            }

      // Generate i-th Given rotation Gi. Note Qn= G0 * ... * Gn
            GeneratePlaneRotation(G.H[i][i], G.H[i + 1][i], G.cs[i], G.sn[i]);
      // Apply Gi
            ApplyPlaneRotation(G.H[i][i], G.H[i + 1][i], G.cs[i], G.sn[i]);
            ApplyPlaneRotation(s[i], s[i + 1], G.cs[i], G.sn[i]);

//            if (! myPID) std::cout << "iter: " << j << ", residual: " << resid << std::endl;

            if ((resid = abs(s[i + 1]) / normb) < tol) {
                MultiVector y(b.Map(), 1, true);
                Update(y, i, G.H, s, *(G.v));
                M->ApplyInverse(y, t);
                x.Update(1.0, t, 1.0);
                tol = resid;
                max_iter = j;
                G.m = i;
                return 2;
            }
        }

        MultiVector y(b.Map(), 1, true);
        Update(y, i - 1, G.H, s, *(G.v));
        M->ApplyInverse(y, t);
        x.Update(1.0, t, 1.0);
        A.Apply(x, t);
        w.Update(1.0, b, -1.0, t, 0.0);
        L->ApplyInverse(w, r);
        r.Norm2(&beta);

        if ((resid = beta / normb) < tol) {
            tol = resid;
            max_iter = j;
            G.m = i;
            return 3;
        }
    }

    tol = resid;
    G.m = i;
    return 1;
}

} // namespace IQR

#endif // IQR_GMRES_H
