// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KLU2_ORDINALTRAITS_H
#define KLU2_ORDINALTRAITS_H

template <typename T>
struct KLU_OrdinalTraits
{
    /* None of this methods should technically return anything.
     * If the OrdinalType is no supported and there is no
     * specialization available, we will error out. These
     * returns are just to keep the compiler happy and stop
     * spewing warnings.
     * */
    static inline T btf_order (T n, T *Ap, T *Ai, double maxwork,
        double *work, T *P, T *Q, T *R, T *nmatch, T *Work)
    {
        return n;
    }

    static inline T btf_strongcomp (T n, T *Ap, T *Ai, T *Q, T *P,
        T *R, T *Work)
    {
        return n;
    }

    static inline T amd_order (T n, T *Ap, T *Ai, T *P, double *Control,
        double *Info)
    {
        return n;
    }

    static inline T colamd (T n_row, T n_col, T Alen, T *A, T *p,
        double *knobs, T *stats)
    {
        return n_row;
    }

    static inline T colamd_recommended (T nnz, T n_row, T n_col)
    {
        return nnz;
    }
};

template<>
struct KLU_OrdinalTraits<int>
{
    static inline int btf_order (int n, int *Ap, int *Ai,
        double maxwork, double *work, int *P, int *Q, int *R, int *nmatch,
        int *Work)
    {
        return (trilinos_btf_order (n, Ap, Ai, maxwork, work, P, Q, R, nmatch,
                    Work)) ;
    }

    static inline int btf_strongcomp (int n, int *Ap, int *Ai, int *Q, int *P,
        int *R, int *Work)
    {
        return (trilinos_btf_strongcomp (n, Ap, Ai, Q, P, R, Work)) ;
    }

    static inline int amd_order (int n, int *Ap, int *Ai, int *P,
        double *Control, double *Info)
    {
        return (trilinos_amd_order(n, Ap, Ai, P, Control, Info)) ;
    }

    static inline int colamd (int n_row, int n_col, int Alen, int *A, int *p,
        double *knobs, int *stats)
    {
        return(trilinos_colamd (n_row, n_col, Alen, A, p, knobs, stats));
    }

    // TODO : return size_t
    static inline int colamd_recommended (int nnz, int n_row, int n_col)
    {
        return(trilinos_colamd_recommended(nnz, n_row, n_col));
    }
};

template<>
struct KLU_OrdinalTraits<ptrdiff_t>
{
// These should all be UF_long, which I presume is resolving to ptrdiff_t
// ptrdiff_t is ptrdiff_t on Linux64, but just to be safe
    static inline ptrdiff_t btf_order (ptrdiff_t n, ptrdiff_t *Ap, ptrdiff_t *Ai,
        double maxwork, double *work, ptrdiff_t *P, ptrdiff_t *Q, ptrdiff_t *R, ptrdiff_t *nmatch,
        ptrdiff_t *Work)
    {
        return (trilinos_btf_l_order (n, Ap, Ai, maxwork, work, P, Q, R, nmatch,
                    Work));
    }

    static inline ptrdiff_t btf_strongcomp (ptrdiff_t n, ptrdiff_t *Ap, ptrdiff_t *Ai, ptrdiff_t *Q,
        ptrdiff_t *P, ptrdiff_t *R, ptrdiff_t *Work)
    {
        return(trilinos_btf_l_strongcomp (n, Ap, Ai, Q, P, R, Work)) ;
    }

    static inline ptrdiff_t amd_order (ptrdiff_t n, ptrdiff_t *Ap, ptrdiff_t *Ai, ptrdiff_t *P,
        double *Control, double *Info)
    {
        return (trilinos_amd_l_order(n, Ap, Ai, P, Control, Info)) ;
    }

    static inline ptrdiff_t colamd (ptrdiff_t n_row, ptrdiff_t n_col, ptrdiff_t Alen, ptrdiff_t *A,
        ptrdiff_t *p, double *knobs, ptrdiff_t *stats)
    {
        return(trilinos_colamd_l (n_row, n_col, Alen, A, p, knobs, stats));
    }

    static inline ptrdiff_t colamd_recommended (ptrdiff_t nnz, ptrdiff_t n_row, ptrdiff_t n_col)
    {
        return(trilinos_colamd_l_recommended(nnz, n_row, n_col));
    }
};

#endif
