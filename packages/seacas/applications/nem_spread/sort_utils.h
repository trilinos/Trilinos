/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SORT_UTILS_H
#define SORT_UTILS_H
template <typename INT> void gds_iqsort(INT v[], INT iv[], size_t N);

template <typename INT> void gds_qsort(INT v[], size_t N);

template <typename INT> void indexed_sort(INT v[], INT iv[], size_t N);
#endif
