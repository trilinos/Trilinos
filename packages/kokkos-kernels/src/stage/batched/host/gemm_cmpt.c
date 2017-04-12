#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mkl.h"
#include "sys/time.h"

#define min(X,Y) (((X) < (Y)) ? (X) : (Y))

#define GEMM_BATCH cblas_dgemm_batch
#define GEMM_COMPUTE_BATCH cblas_dgemm_compute_batch
#define fabst fabs
#define FPTYPE double
#define VECTOR_LENGTH 8

#define REPEAT 1000
#define MEM_ALIGNMENT 64

#define GRP_COUNT 1

int main(int argc, char *argv[])
{

  int grp, i, j, ii, jj, nthr;
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
  MKL_INT k[GRP_COUNT] = {m_init};
  MKL_INT n[GRP_COUNT] = {m_init};

  MKL_INT lda[GRP_COUNT] = {m_init};
  MKL_INT ldb[GRP_COUNT] = {m_init};
  MKL_INT ldc[GRP_COUNT] = {m_init};
  FPTYPE alpha[GRP_COUNT] = {1.0};
  FPTYPE beta[GRP_COUNT] = {1.0};

  MKL_INT size_per_grp[GRP_COUNT] = {grp_init};
  MKL_INT format_tail;

  CBLAS_TRANSPOSE TRANSA[GRP_COUNT] = {CblasNoTrans};
  CBLAS_TRANSPOSE TRANSB[GRP_COUNT] = {CblasNoTrans};

  int num_pointers = 0;
  for (i = 0; i < GRP_COUNT; i++) num_pointers += size_per_grp[i];

  int a_total = 0;
  int b_total = 0;
  int c_total = 0;
  for (i = 0; i < GRP_COUNT; i++) a_total += m[i] * k[i] * size_per_grp[i];
  for (i = 0; i < GRP_COUNT; i++) b_total += k[i] * n[i] * size_per_grp[i];
  for (i = 0; i < GRP_COUNT; i++) c_total += m[i] * n[i] * size_per_grp[i];

  FPTYPE *a, *b, *c, *d;
  a = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  b = (FPTYPE *)_mm_malloc( b_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  c = (FPTYPE *)_mm_malloc( c_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  d = (FPTYPE *)_mm_malloc( c_total * sizeof(FPTYPE), MEM_ALIGNMENT );

  for (i = 0; i < a_total; i++) a[i] = rand() / (FPTYPE) RAND_MAX + .5;
  for (i = 0; i < b_total; i++) b[i] = rand() / (FPTYPE) RAND_MAX + .5;
  int ci = 0;
  for (i = 0; i < c_total; i++) {
    c[i] = ci;
    d[i] = c[i];
    ci++;
  }

  FPTYPE *a_array[num_pointers], *b_array[num_pointers], *c_array[num_pointers], *d_array[num_pointers];
  int a_idx = 0;
  int b_idx = 0;
  int c_idx = 0;
  int p_num = 0;
  for (i = 0; i < GRP_COUNT; i++) {
    for (j = 0; j < size_per_grp[i]; j++) {
      a_array[p_num] = &a[ a_idx ];
      b_array[p_num] = &b[ b_idx ];
      c_array[p_num] = &c[ c_idx ];
      d_array[p_num] = &d[ c_idx ];
      p_num++;
      a_idx += m[i] * k[i];
      b_idx += k[i] * n[i];
      c_idx += m[i] * n[i];
      op_count += 2 * m[i] * n[i] * k[i];
    }
  }

  // setup packed arrays
  int grp_idx = 0;
  int A_p_idx = 0;
  int B_p_idx = 0;
  int C_p_idx = 0;
  int i_grp, i_format;
  FPTYPE *a_packed, *b_packed, *c_packed;
  a_packed = (FPTYPE *)_mm_malloc( a_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  b_packed = (FPTYPE *)_mm_malloc( b_total * sizeof(FPTYPE), MEM_ALIGNMENT );
  c_packed = (FPTYPE *)_mm_malloc( c_total * sizeof(FPTYPE), MEM_ALIGNMENT );

  // setup compact_t structs
  compact_t A_p, B_p, C_p;
  A_p.layout = CblasColMajor;
  A_p.rows = m;
  A_p.cols = k;
  A_p.stride = lda;
  A_p.group_count = GRP_COUNT;
  A_p.size_per_group = size_per_grp;
  A_p.format = VECTOR_LENGTH;
  A_p.mat = a_packed;
  B_p.layout = CblasColMajor;
  B_p.rows = k;
  B_p.cols = n;
  B_p.stride = ldb;
  B_p.group_count = GRP_COUNT;
  B_p.size_per_group = size_per_grp;
  B_p.format = VECTOR_LENGTH;
  B_p.mat = b_packed;
  C_p.layout = CblasColMajor;
  C_p.rows = m;
  C_p.cols = n;
  C_p.stride = ldc;
  C_p.group_count = GRP_COUNT;
  C_p.size_per_group = size_per_grp;
  C_p.format = VECTOR_LENGTH;
  C_p.mat = c_packed;

  // pack the matrices in VECTOR_LENGTH chunks
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
      for (j=0; j<C_p.cols[i_grp]; j++) {
        for (i=0; i<C_p.rows[i_grp]; i++) {
          for (ii=0; ii<VECTOR_LENGTH; ii++) {
            c_packed[ii + C_p_idx] = c_array[ii + grp_idx + i_format][C_p.rows[i_grp]*j + i];
          }
          C_p_idx+=VECTOR_LENGTH;
        }
      }
    }
    // tail handling
    format_tail = size_per_grp[i_grp] - i_format;
    i_format = (size_per_grp[i_grp]/VECTOR_LENGTH)*VECTOR_LENGTH;

    if(format_tail > 0) {
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
      for (j=0; j<C_p.cols[i_grp]; j++) {
        for (i=0; i<C_p.rows[i_grp]; i++) {
          for (ii=0; ii<format_tail; ii++) {
            c_packed[ii + C_p_idx] = c_array[ii + grp_idx + i_format][C_p.rows[i_grp]*j + i];
          }
          C_p_idx+=format_tail;
        }
      }
    }
    grp_idx += size_per_grp[i_grp];
  }

  printf("\n THREADS --- m = k = n --- GEMM_BATCH --- GEMM_BATCH_COMPUTE\n");
  for (i = 0; i < REPEAT; i++) {
    startt = omp_get_wtime();
    GEMM_BATCH( CblasColMajor, TRANSA, TRANSB, m, n, k, alpha, (const FPTYPE**)a_array, lda, (const FPTYPE**)b_array, ldb, beta, d_array, ldc, GRP_COUNT, size_per_grp );
    stopt = omp_get_wtime();
    mint = min(mint, (stopt - startt));
  }
  int one=1;
  for (i = 0; i < REPEAT; i++) {
    startt = omp_get_wtime();
    GEMM_COMPUTE_BATCH( TRANSA, TRANSB, alpha, &A_p, &B_p, beta, &C_p );
    stopt = omp_get_wtime();
    mints = min(mints, (stopt - startt));
  }

  printf(" %d      --- %d         --- %5.2f      --- %5.2f\n", nthr, m_init, (double)op_count/(mint*1e9), (double)op_count/(mints*1e9));

  _mm_free(a_packed);
  _mm_free(b_packed);
  _mm_free(c_packed);
  _mm_free(a);
  _mm_free(b);
  _mm_free(c);
  _mm_free(d);

  return 0;
}
