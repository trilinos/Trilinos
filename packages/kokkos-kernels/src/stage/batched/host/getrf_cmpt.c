#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mkl.h"
#include "sys/time.h"

#define idx2(i, j, ldi) ((j * ldi) + i)

#define min(X,Y) (((X) < (Y)) ? (X) : (Y))
#define max(X,Y) (((X) > (Y)) ? (X) : (Y))

#define GETRF dgetrf
#define GETRF_COMPUTE_BATCH LAPACKE_dgetrf_compute_batch
#define fabst fabs
#define FPTYPE double
#define VECTOR_LENGTH 8

#define REPEAT 1000
#define MEM_ALIGNMENT 64

#define GRP_COUNT 1

int main(int argc, char *argv[])
{

  int grp, i, j, ii, jj, nthr;
  int g, i_matrix, num_matrices;
#pragma omp parallel
  {
    nthr = omp_get_num_threads();
  }
  unsigned long op_count = 0;
  double startt, stopt, mint, mints;
  mint = 1e5;
  mints = 1e5;
  int group_count = GRP_COUNT;
  int pack_length = VECTOR_LENGTH;

  int m_init, grp_init;
  if (argc > 1) m_init = atoi(argv[1]);
  else m_init = 5;
  if (argc > 2) grp_init = atoi(argv[2]);
  else grp_init = 512;

  MKL_INT m[GRP_COUNT] = {m_init};
  MKL_INT n[GRP_COUNT] = {m_init};

  MKL_INT lda[GRP_COUNT] = {m_init};

  MKL_INT size_per_grp[GRP_COUNT] = {grp_init};

  int num_pointers = 0;
  for (i = 0; i < GRP_COUNT; i++) num_pointers += size_per_grp[i];

  int a_total = 0;
  for (i = 0; i < GRP_COUNT; i++) a_total += m[i] * n[i] * size_per_grp[i];

  FPTYPE *a, *d;
  a = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  d = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  MKL_INT *ipiv, *info;
  ipiv = (MKL_INT *)_mm_malloc( a_total * sizeof(MKL_INT), MEM_ALIGNMENT );
  info = (MKL_INT *)_mm_malloc( num_pointers * sizeof(MKL_INT), MEM_ALIGNMENT );

  for (i = 0; i < a_total; i++) a[i] = rand() / (FPTYPE) RAND_MAX + .5;
  for (i = 0; i < grp_init; i++)
      for (j = 0; j < m_init * m_init; j+= m_init+1)
          a[i*(m_init*m_init) + j] += 10.0;
  for (i = 0; i < a_total; i++) d[i] = a[i];
  for (i = 0; i < a_total; i++) ipiv[i] = 0;

  FPTYPE *a_array[num_pointers], *d_array[num_pointers];
  MKL_INT *ipiv_array[num_pointers];
  int a_idx = 0;
  int p_num = 0;
  for (i = 0; i < GRP_COUNT; i++) {
    double emn = (double) min(m[i], n[i]);
    for (j = 0; j < size_per_grp[i]; j++) {
      a_array[p_num] = &a[ a_idx ];
      d_array[p_num] = &d[ a_idx ];
      ipiv_array[p_num] = &ipiv[ a_idx ];
      p_num++;
      a_idx += m[i] * n[i];
      op_count += (2.0 * m[i] * m[i] * m[i]) / 3.0;
    }
  }

  // setup packed arrays
  int grp_idx = 0;
  int A_p_idx = 0;
  int i_grp, i_format, format_tail;
  FPTYPE *a_packed;
  a_packed = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );

  // setup compact_t structs
  compact_t A_p;
  A_p.layout = CblasColMajor;
  A_p.rows = m;
  A_p.cols = n;
  A_p.stride = lda;
  A_p.group_count = GRP_COUNT;
  A_p.size_per_group = size_per_grp;
  A_p.format = VECTOR_LENGTH;
  A_p.mat = a_packed;
  for (i_grp=0; i_grp<GRP_COUNT; i_grp++) {
    for (i_format=0; i_format<(size_per_grp[i_grp]/VECTOR_LENGTH)*VECTOR_LENGTH; i_format+=VECTOR_LENGTH) {
      for (j=0; j<A_p.cols[i_grp]; j++) {
        for (i=0; i<A_p.rows[i_grp]; i++) {
          for (ii=0; ii<VECTOR_LENGTH; ii++) {
            a_packed[ii + A_p_idx] = a_array[ii + grp_idx + i_format][A_p.rows[i_grp]*j + i];
          }
          A_p_idx+=VECTOR_LENGTH;
        }
      }
    }
    // tail handling
    format_tail = size_per_grp[i_grp] - i_format;
    i_format = (size_per_grp[i_grp]/VECTOR_LENGTH)*VECTOR_LENGTH;

    if (format_tail > 0) {
      for (j=0; j<A_p.cols[i_grp]; j++) {
        for (i=0; i<A_p.rows[i_grp]; i++) {
          for (ii=0; ii<format_tail; ii++) {
            a_packed[ii + A_p_idx] = a_array[ii + grp_idx + i_format][A_p.rows[i_grp]*j + i];
          }
          A_p_idx+=format_tail;
        }
      }
    }
    grp_idx += size_per_grp[i_grp];
  }

  printf("\n THREADS --- m = k = n --- GETRF_BATCH --- GETRF_BATCH_COMPUTE\n");
  for (i = 0; i < REPEAT; i++) {
    num_matrices = 0;
    startt = omp_get_wtime();
    for (g=0; g < GRP_COUNT; g++) {
#pragma omp parallel for shared(m, n, d_array, lda, ipiv_array, info, g) private(i_matrix) schedule(static,1)
        for (i_matrix=num_matrices; i_matrix<size_per_grp[g] + num_matrices; i_matrix++) {
            GETRF( &m[g], &n[g], d_array[i_matrix], &lda[g], ipiv_array[i_matrix], &info[i_matrix] );
        }
        num_matrices += size_per_grp[g];
    }
    stopt = omp_get_wtime();
    mint = min(mint, (stopt - startt));
  }

  int one=1;
  for (i = 0; i < REPEAT; i++) {
    startt = omp_get_wtime();
    GETRF_COMPUTE_BATCH( &A_p );
    stopt = omp_get_wtime();
    mints = min(mints, (stopt - startt));
  }

  printf(" %d      --- %d         --- %5.2f      --- %5.2f\n", nthr, m_init, (double)op_count/(mint*1e9), (double)op_count/(mints*1e9));

  _mm_free(a_packed);
  _mm_free(a);
  _mm_free(d);
  _mm_free(ipiv);
  _mm_free(info);

  return 0;
}
