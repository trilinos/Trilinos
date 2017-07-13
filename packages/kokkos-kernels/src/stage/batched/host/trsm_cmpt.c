#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mkl.h"
#include "sys/time.h"
#include "math.h"

#define min(X,Y) (((X) < (Y)) ? (X) : (Y))

#define TRSM_BATCH cblas_dtrsm_batch
#define TRSM_COMPUTE_BATCH cblas_dtrsm_compute_batch
#define fabst fabs
#define FPTYPE double
#define VECTOR_LENGTH 8

#define REPEAT 1000
#define MEM_ALIGNMENT 64

#define GRP_COUNT 1

int main(int argc, char *argv[])
{

  int i, j, ii, jj, nthr, grp;
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

  MKL_INT lda[GRP_COUNT] = {m_init};
  MKL_INT ldb[GRP_COUNT] = {m_init};
  FPTYPE alpha[GRP_COUNT] = {1.0};

  MKL_INT size_per_grp[GRP_COUNT] = {grp_init};

  CBLAS_SIDE SIDE[GRP_COUNT] = {CblasRight};
  CBLAS_UPLO UPLO[GRP_COUNT] = {CblasUpper};
  CBLAS_TRANSPOSE TRANSA[GRP_COUNT] = {CblasNoTrans};
  CBLAS_DIAG DIAG[GRP_COUNT] = {CblasNonUnit};

  int num_pointers = 0;
  for (i = 0; i < GRP_COUNT; i++) num_pointers += size_per_grp[i];

  int a_total = 0;
  int b_total = 0;
  for (i = 0; i < GRP_COUNT; i++) a_total += m[i] * m[i] * size_per_grp[i];
  for (i = 0; i < GRP_COUNT; i++) b_total += m[i] * m[i] * size_per_grp[i];

  FPTYPE *a, *b, *c;
  a = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  b = (FPTYPE *)_mm_malloc( b_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  c = (FPTYPE *)_mm_malloc( b_total * sizeof(FPTYPE), MEM_ALIGNMENT );

  for (i = 0; i < a_total; i++) a[i] = (FPTYPE) rand() / (FPTYPE) RAND_MAX + .5;
  for (i = 0; i < b_total; i++) b[i] = (FPTYPE) rand() / (FPTYPE) RAND_MAX + .5;
  for (i = 0; i < b_total; i++) c[i] = b[i];

  FPTYPE *a_array[num_pointers], *b_array[num_pointers];
  int a_idx = 0;
  int b_idx = 0;
  int p_num = 0;
  for (i = 0; i < GRP_COUNT; i++) {
    for (j = 0; j < size_per_grp[i]; j++) {
      a_array[p_num] = &a[ a_idx ];
      b_array[p_num] = &b[ b_idx ];
      p_num++;
      a_idx += m[i] * m[i];
      b_idx += m[i] * m[i];
      op_count += m[i] * m[i] * m[i];
    }
  }
  for (i = 0; i < GRP_COUNT; i++) {
    int idx_matrix = 0;
    for (ii = 0; ii < i; ii++) idx_matrix += size_per_grp[ii];
    for (j = 0; j < size_per_grp[i]; j++) {
      for (ii = 0; ii < m[i]; ii++) {
        a_array[idx_matrix][ii*lda[i] + ii] = 10.0;
      }
      idx_matrix++;
    }
  }

  // setup packed arrays
  int grp_idx = 0;
  int A_p_idx = 0;
  int B_p_idx = 0;
  int i_grp, i_format, format_tail;
  FPTYPE *a_packed, *b_packed;
  a_packed = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  b_packed = (FPTYPE *)_mm_malloc( b_total * sizeof(FPTYPE), MEM_ALIGNMENT );

  // setup compact_t structs
  compact_t A_p, B_p, C_p;
  A_p.layout = CblasColMajor;
  A_p.rows = m;
  A_p.cols = m;
  A_p.stride = lda;
  A_p.group_count = GRP_COUNT;
  A_p.size_per_group = size_per_grp;
  A_p.format = VECTOR_LENGTH;
  A_p.mat = a_packed;
  B_p.layout = CblasColMajor;
  B_p.rows = m;
  B_p.cols = m;
  B_p.stride = ldb;
  B_p.group_count = GRP_COUNT;
  B_p.size_per_group = size_per_grp;
  B_p.format = VECTOR_LENGTH;
  B_p.mat = b_packed;
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
      for (j=0; j<B_p.cols[i_grp]; j++) {
        for (i=0; i<B_p.rows[i_grp]; i++) {
          for (ii=0; ii<VECTOR_LENGTH; ii++) {
            b_packed[ii + B_p_idx] = b_array[ii + grp_idx + i_format][B_p.rows[i_grp]*j + i];
          }
          B_p_idx+=VECTOR_LENGTH;
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
      for (j=0; j<B_p.cols[i_grp]; j++) {
        for (i=0; i<B_p.rows[i_grp]; i++) {
          for (ii=0; ii<format_tail; ii++) {
            b_packed[ii + B_p_idx] = b_array[ii + grp_idx + i_format][B_p.rows[i_grp]*j + i];
          }
          B_p_idx+=format_tail;
        }
      }
    }
    grp_idx += size_per_grp[i_grp];
  }

  printf("\n THREADS --- m = n --- TRSM_BATCH --- TRSM_BATCH_COMPUTE\n");
  for (i = 0; i < REPEAT; i++) {
    startt = omp_get_wtime();
    TRSM_BATCH( CblasColMajor, SIDE, UPLO, TRANSA, DIAG, m, m, alpha, (const FPTYPE**)a_array, lda, b_array, ldb, GRP_COUNT, size_per_grp );
    stopt = omp_get_wtime();
    mint = min(mint, (stopt - startt));
  }

  for (i = 0; i < REPEAT; i++) {
    startt = omp_get_wtime();
    TRSM_COMPUTE_BATCH( SIDE, UPLO, TRANSA, DIAG, alpha, &A_p, &B_p );
    stopt = omp_get_wtime();
    mints = min(mints, (stopt - startt));
  }

  printf(" %d      --- %d         --- %5.2f      --- %5.2f\n", nthr, m_init, (double)op_count/(mint*1e9), (double)op_count/(mints*1e9));

  _mm_free(a_packed);
  _mm_free(b_packed);
  _mm_free(a);
  _mm_free(b);
  _mm_free(c);

  return 0;
}
