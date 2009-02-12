/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
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
  const struct TestTPI_DNAX * const data =
    (struct TestTPI_DNAX *) work->shared ;

  const unsigned which_chunk = work->work_rank ;
  const unsigned beg_local   = data->chunk_length * which_chunk ;
  const unsigned max_local   = data->length - beg_local ;
  const unsigned len_local   = data->chunk_length < max_local ?
                               data->chunk_length : max_local ;

  test_dnax_column( data->number ,
                    data->stride ,
                    len_local ,
                    data->coef ,
                    data->array + beg_local );

  return ;
}

/*  The multi-array storage is chunked: each array has a contiguous chunk;
 *  but chunk-subarrays are contiguously grouped.
 */
static
void test_dnax_column_work( TPI_Work * work )
{
  const struct TestTPI_DNAX * const data =
    (struct TestTPI_DNAX *) work->shared ;

  const unsigned which_chunk = work->work_rank ;
  const unsigned beg_local   = data->chunk_length * which_chunk ;
  const unsigned max_local   = data->length - beg_local ;
  const unsigned len_local   = data->chunk_length < max_local ?
                               data->chunk_length : max_local ;

  const unsigned chunk_size = data->chunk_length * data->number ;

  test_dnax_column( data->number ,
                    data->chunk_length ,
                    len_local ,
                    data->coef ,
                    data->array + which_chunk * chunk_size );

  return ;
}

static
void test_dnax_row_work( TPI_Work * work )
{
  const struct TestTPI_DNAX * const data =
    (struct TestTPI_DNAX *) work->shared ;

  const unsigned which_chunk = work->work_rank ;
  const unsigned beg_local   = data->chunk_length * which_chunk ;
  const unsigned max_local   = data->length - beg_local ;
  const unsigned len_local   = data->chunk_length < max_local ?
                               data->chunk_length : max_local ;

  const unsigned chunk_size = data->chunk_length * data->number ;

  test_dnax_row( data->number ,
                 data->number ,
                 len_local ,
                 data->coef ,
                 data->array + which_chunk * chunk_size );

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



