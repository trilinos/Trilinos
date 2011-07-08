/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/* Redistribution and use in source and binary forms, with or without     */
/* modification, are permitted provided that the following conditions are */
/* met:                                                                   */
/*                                                                        */
/* 1. Redistributions of source code must retain the above copyright      */
/* notice, this list of conditions and the following disclaimer.          */
/*                                                                        */
/* 2. Redistributions in binary form must reproduce the above copyright   */
/* notice, this list of conditions and the following disclaimer in the    */
/* documentation and/or other materials provided with the distribution.   */
/*                                                                        */
/* 3. Neither the name of the Corporation nor the names of the            */
/* contributors may be used to endorse or promote products derived from   */
/* this software without specific prior written permission.               */
/*                                                                        */
/* THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY        */
/* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR     */
/* PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE    */
/* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,    */
/* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR     */
/* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF */
/* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING   */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS     */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 *
 *  Multi-array 'axpby'
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>

#if defined( HAVE_MPI )
#include <mpi.h>
#endif

int test_c_tpi_dnax( int , int );

int main( int argc , char ** argv )
{
  int num_thread[] = { 1 , 2 , 4 , 6 , 8 , 12 , 16 };
  int num_test = sizeof(num_thread) / sizeof(int);
 
  const int ntrial = 1 < argc ? atoi( argv[1] ) : 2 ;
  int i ;

#if defined( HAVE_MPI )
  int rank ;
 
  MPI_Init( & argc , & argv );
  MPI_Comm_rank( MPI_COMM_WORLD , & rank );
  if ( 0 == rank ) {
#endif


  fprintf( stdout , "\"TESTING Multiarray 'axpby' with: %s\"\n" ,
           TPI_Version() );
 
  for ( i = 0 ; i < num_test ; ++i ) {
    test_c_tpi_dnax( num_thread[i] , ntrial );
  }

#if defined( HAVE_MPI )
  }
  MPI_Finalize();
#endif
 
  return 0 ;
}

/*------------------------------------------------------------------------*/

typedef double SCALAR ;

/*------------------------------------------------------------------------*/

struct TestTPI_DNAX {
  SCALAR * coef ;
  SCALAR * array ;
  unsigned number ;
  unsigned length ;
  unsigned stride ;
  unsigned chunk_length ;
};

/*------------------------------------------------------------------------*/

static
void test_dnax_column( const unsigned num_array , 
                       const unsigned stride ,
                       const unsigned length , 
                       const SCALAR * const coef ,
                       SCALAR * const array )
{
  unsigned i = 0 ;
  for ( ; i < length ; ++i ) {
    SCALAR * const a = array + i ;
    SCALAR tmp = 0 ;
    unsigned j = 0 ;
    for ( ; j < num_array ; ++j ) { tmp += coef[j] * a[ j * stride ] ; }
    a[0] = tmp ;
  }
}

static
void test_dnax_row( const unsigned num_array , 
                    const unsigned stride ,
                    const unsigned length , 
                    const SCALAR * const coef ,
                    SCALAR * const array )
{
  unsigned i = 0 ;
  for ( ; i < length ; ++i ) {
    SCALAR * const a = array + i * stride ;
    SCALAR tmp = 0 ;
    unsigned j = 0 ;
    for ( ; j < num_array ; ++j ) { tmp += coef[j] * a[j] ; }
    a[0] = tmp ;
  }
}

/*------------------------------------------------------------------------*/
/*  The multi-array storage is flat: every array is fully contiguous.
 *  Work corresponds to a span of the array.
 */
static
void test_dnax_flat_work( TPI_Work * work )
{
  const struct TestTPI_DNAX * const info =
    (struct TestTPI_DNAX *) work->info ;

  const unsigned which_chunk = work->rank ;
  const unsigned beg_local   = info->chunk_length * which_chunk ;
  const unsigned max_local   = info->length - beg_local ;
  const unsigned len_local   = info->chunk_length < max_local ?
                               info->chunk_length : max_local ;

  test_dnax_column( info->number ,
                    info->stride ,
                    len_local ,
                    info->coef ,
                    info->array + beg_local );

  return ;
}

/*  The multi-array storage is chunked: each array has a contiguous chunk;
 *  but chunk-subarrays are contiguously grouped.
 */
static
void test_dnax_column_work( TPI_Work * work )
{
  const struct TestTPI_DNAX * const info =
    (struct TestTPI_DNAX *) work->info ;

  const unsigned which_chunk = work->rank ;
  const unsigned beg_local   = info->chunk_length * which_chunk ;
  const unsigned max_local   = info->length - beg_local ;
  const unsigned len_local   = info->chunk_length < max_local ?
                               info->chunk_length : max_local ;

  const unsigned chunk_size = info->chunk_length * info->number ;

  test_dnax_column( info->number ,
                    info->chunk_length ,
                    len_local ,
                    info->coef ,
                    info->array + which_chunk * chunk_size );

  return ;
}

static
void test_dnax_row_work( TPI_Work * work )
{
  const struct TestTPI_DNAX * const info =
    (struct TestTPI_DNAX *) work->info ;

  const unsigned which_chunk = work->rank ;
  const unsigned beg_local   = info->chunk_length * which_chunk ;
  const unsigned max_local   = info->length - beg_local ;
  const unsigned len_local   = info->chunk_length < max_local ?
                               info->chunk_length : max_local ;

  const unsigned chunk_size = info->chunk_length * info->number ;

  test_dnax_row( info->number ,
                 info->number ,
                 len_local ,
                 info->coef ,
                 info->array + which_chunk * chunk_size );

  return ;
}

/*------------------------------------------------------------------------*/
/* Process identical block of allocated memory as a
 * as a flat array, chunked-column, and chunked-row.
 */

static
void test_tpi_dnax_driver( const int nthread ,
                           const unsigned Mflop_target ,
                           const unsigned num_trials ,
                           const unsigned num_test ,
                           const unsigned num_test_array[] ,
                           const unsigned length_array ,
                           const unsigned length_chunk )
{
  const unsigned max_array = num_test_array[ num_test - 1 ];

  const unsigned num_chunk =
    ( length_array + length_chunk - 1 ) / length_chunk ;

  const unsigned stride_array = num_chunk * length_chunk ;
  const unsigned size_alloc   = max_array * stride_array ;

  SCALAR * const coef  = (SCALAR *) malloc( max_array * sizeof(SCALAR) );
  SCALAR * const array = (SCALAR *) malloc( size_alloc * sizeof(SCALAR) );

  struct TestTPI_DNAX data = { NULL , NULL , 0 , 0 , 0 , 0 };

  unsigned i_test , i , j ;

  data.coef = coef ;

  if ( NULL == array ) {
    fprintf(stderr,"allocation failure for %u\n",size_alloc);
    abort();
  }

  for ( i = 0 ; i < max_array ; ++i ) { coef[i] = 0 ; }

  printf("\n\"test_tpi_dnax[%d]( length_array = %u , stride_array = %u )\"\n",
         nthread , length_array , stride_array );
  printf("\"NUMBER OF THREADS\" , %d\n" , nthread );
  printf("\"NUMBER OF CHUNKS\" , %u\n" , num_chunk );
  printf("\"NUMBER OF TRIALS\" , %u \n", num_trials );

  printf("\"TEST\" , \"#ARRAY\" \"DT-MEAN\" , \"DT-STDDEV\" , \"MFLOP-MEAN\" , \"MFLOP-STDDEV\"\n");

  /*----------------------------------------------------------------------*/

  for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
    const unsigned num_array = num_test_array[ i_test ];
    const unsigned num_sets  = max_array / num_array ;

    const double mflop_cycle =
      ((double)( 2 * num_array * length_array )) / 1.0e6 ;

    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

    double dt_sum = 0 ;
    double dt_sum_2 = 0 ;

    data.length       = length_array ;
    data.number       = num_array ;
    data.stride       = stride_array ;
    data.chunk_length = length_chunk ;

    for ( i = 0 ; i < size_alloc ; ++i ) { array[i] = 0 ; }

    for ( j = 0 ; j < num_trials ; ++j ) {

      double dt_tmp = TPI_Walltime();
      for ( i = 0 ; i < ncycle ; ++i ) {
        data.array = array + stride_array * num_array * ( i % num_sets );
        TPI_Run( & test_dnax_flat_work , & data , num_chunk , 0 );
      }
      dt_tmp = TPI_Walltime() - dt_tmp ;

      dt_sum += dt_tmp ;
      dt_sum_2 += dt_tmp * dt_tmp ;
    }

    {
      const double dt_mean = dt_sum / num_trials ;
      const double dt_sdev = sqrt( ( num_trials * dt_sum_2 - dt_sum * dt_sum ) / ( num_trials * ( num_trials - 1 ) ) );
      const double mflop_mean = mflop_cycle * ncycle / dt_mean ;
      const double mflop_sdev = mflop_mean * dt_sdev / ( dt_mean + dt_sdev );

      printf("\"FLAT  ARRAY\"  , %6u , %9.5g , %9.3g , %9.5g , %9.3g\n",
             num_array, dt_mean, dt_sdev, mflop_mean, mflop_sdev );
    }
  }

  /*----------------------------------------------------------------------*/

  for ( i_test = 0 ; i_test < num_test ; ++i_test ) {

    const unsigned num_array = num_test_array[ i_test ];
    const unsigned num_sets  = max_array / num_array ;

    const double mflop_cycle =
      ((double)( 2 * num_array * length_array )) / 1.0e6 ;

    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

    double dt_sum = 0 ;
    double dt_sum_2 = 0 ;

    data.length       = length_array ;
    data.number       = num_array ;
    data.stride       = stride_array ;
    data.chunk_length = length_chunk ;

    for ( i = 0 ; i < size_alloc ; ++i ) { array[i] = 0 ; }

    for ( j = 0 ; j < num_trials ; ++j ) {

      double dt_tmp = TPI_Walltime();
      for ( i = 0 ; i < ncycle ; ++i ) {
        data.array = array + stride_array * num_array * ( i % num_sets );
        TPI_Run( & test_dnax_column_work , & data , num_chunk , 0 );
      }
      dt_tmp = TPI_Walltime() - dt_tmp ;

      dt_sum += dt_tmp ;
      dt_sum_2 += dt_tmp * dt_tmp ;
    }

    {
      const double dt_mean = dt_sum / num_trials ;
      const double dt_sdev = sqrt( ( num_trials * dt_sum_2 - dt_sum * dt_sum ) / ( num_trials * ( num_trials - 1 ) ) );
      const double mflop_mean = mflop_cycle * ncycle / dt_mean ;
      const double mflop_sdev = mflop_mean * dt_sdev / ( dt_mean + dt_sdev );

      printf("\"CHUNK COLUMN\" , %6u , %9.5g , %9.3g , %9.5g , %9.3g\n",
             num_array, dt_mean, dt_sdev, mflop_mean, mflop_sdev );
    }
  }

  /*----------------------------------------------------------------------*/

  for ( i_test = 0 ; i_test < num_test ; ++i_test ) {

    const unsigned num_array = num_test_array[ i_test ];
    const unsigned num_sets  = max_array / num_array ;

    const double mflop_cycle =
      ((double)( 2 * num_array * length_array )) / 1.0e6 ;

    const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

    double dt_sum = 0 ;
    double dt_sum_2 = 0 ;

    data.length       = length_array ;
    data.number       = num_array ;
    data.stride       = stride_array ;
    data.chunk_length = length_chunk ;

    for ( i = 0 ; i < size_alloc ; ++i ) { array[i] = 0 ; }

    for ( j = 0 ; j < num_trials ; ++j ) {

      double dt_tmp = TPI_Walltime();

      for ( i = 0 ; i < ncycle ; ++i ) {
        data.array = array + stride_array * num_array * ( i % num_sets );
        TPI_Run( & test_dnax_row_work , & data , num_chunk , 0 );
      }
      dt_tmp = TPI_Walltime() - dt_tmp ;

      dt_sum += dt_tmp ;
      dt_sum_2 += dt_tmp * dt_tmp ;
    }

    {
      const double dt_mean = dt_sum / num_trials ;
      const double dt_sdev = sqrt( ( num_trials * dt_sum_2 - dt_sum * dt_sum ) / ( num_trials * ( num_trials - 1 ) ) );
      const double mflop_mean = mflop_cycle * ncycle / dt_mean ;
      const double mflop_sdev = mflop_mean * dt_sdev / ( dt_mean + dt_sdev );

      printf("\"CHUNK ROW\"    , %6u , %9.5g , %9.3g , %9.5g , %9.3g\n",
             num_array, dt_mean, dt_sdev, mflop_mean, mflop_sdev );
    }
  }

  /*----------------------------------------------------------------------*/

  free( array );
  free( coef );
}

/*------------------------------------------------------------------------*/

int test_c_tpi_dnax( int nthread , int ntrial )
{
  const unsigned Mflop_target = 10 ;
  const unsigned num_array[6] = { 2 , 5 , 10 , 20 , 50 , 100 };
  const unsigned ntest = sizeof(num_array) / sizeof(unsigned);

  if ( ntrial <= 0 ) { ntrial = 7 ; }

  TPI_Init( nthread );

  test_tpi_dnax_driver( nthread ,
                        Mflop_target * nthread ,
                        ntrial    /* number trials */ ,
                        ntest     /* number of tests */ ,
                        num_array /* number of arrays for each test */ ,
                        1e6       /* array computation length */ ,
                        1000      /* chunk length */ );

  TPI_Finalize();

  return 0 ;
}



