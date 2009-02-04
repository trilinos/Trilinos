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

  SCALAR coef[ max_array ];

  SCALAR * const array = (SCALAR *) malloc( size_alloc * sizeof(SCALAR) );

  struct TestTPI_DNAX data = { coef , NULL , 0 , 0 , 0 , 0 };

  if ( NULL == array ) {
    fprintf(stderr,"allocation failure for %u\n",size_alloc);
    abort();
  }

  {
    unsigned i = 0 ;
    for ( ; i < max_array ; ++i ) { coef[i] = 0 ; }
  }

  printf("\n\"test_tpi_dnax[%d]( length_array = %u , stride_array = %u )\"\n",
         nthread , length_array , stride_array );
  printf("\"NUMBER OF THREADS\" , %d\n" , nthread );
  printf("\"NUMBER OF CHUNKS\" , %u\n" , num_chunk );
  printf("\"NUMBER OF TRIALS\" , %u \n", num_trials );

  printf("\"NUMBER OF ARRAYS\"");
  {
    unsigned i_test = 0 ;
    for ( ; i_test < num_test ; ++i_test ) {
      printf(" , %u", num_test_array[ i_test ] );
    }
  }
  printf("\n");

  /*----------------------------------------------------------------------*/
  {
    double dt_min[ num_test ];
    double dt_max[ num_test ];
    double dt_mean[ num_test ];
    double mflops_min[ num_test ];
    double mflops_max[ num_test ];
    double mflops_mean[ num_test ];

    unsigned i_test = 0 ;

    for ( ; i_test < num_test ; ++i_test ) {
      const unsigned num_array = num_test_array[ i_test ];
      const unsigned num_sets  = max_array / num_array ;

      const double mflop_cycle =
        ((double)( 2 * num_array * length_array )) / 1.0e6 ;

      const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

      data.length       = length_array ;
      data.number       = num_array ;
      data.stride       = stride_array ;
      data.chunk_length = length_chunk ;

      { unsigned i = 0 ; for ( ; i < size_alloc ; ++i ) { array[i] = 0 ; } }

      {
        double dt_tmp ;
        unsigned repeat = 0 ;
        unsigned i ;
        for ( ; repeat < num_trials ; ++repeat ) {

          dt_tmp = TPI_Walltime();
          for ( i = 0 ; i < ncycle ; ++i ) {
            data.array = array + stride_array * num_array * ( i % num_sets );
            TPI_Run( & test_dnax_flat_work , & data , num_chunk , 0 );
          }
          dt_tmp = TPI_Walltime() - dt_tmp ;

          if ( 0 == repeat ) {
            dt_min[ i_test ] = dt_tmp ;
            dt_max[ i_test ] = dt_tmp ;
            dt_mean[ i_test ] = dt_tmp / num_trials ;
          }
          else {
            dt_mean[ i_test ] += dt_tmp / num_trials ;
          }
          if ( dt_tmp < dt_min[ i_test ] ) { dt_min[ i_test ] = dt_tmp ; }
          if ( dt_tmp > dt_max[ i_test ] ) { dt_max[ i_test ] = dt_tmp ; }
        }
      }

      mflops_max[ i_test ] = mflop_cycle * ncycle / dt_min[ i_test ];
      mflops_min[ i_test ] = mflop_cycle * ncycle / dt_max[ i_test ];
      mflops_mean[ i_test ] = mflop_cycle * ncycle / dt_mean[ i_test ];
    }

    printf("\"FLAT ARRAY max  time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_max[ i_test ] );
    }
    printf("\n");

    printf("\"FLAT ARRAY mean time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_mean[ i_test ] );
    }
    printf("\n");

    printf("\"FLAT ARRAY min  time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_min[ i_test ] );
    }
    printf("\n");

    printf("\"FLAT ARRAY max  Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_max[ i_test ] );
    }
    printf("\n");

    printf("\"FLAT ARRAY mean Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_mean[ i_test ] );
    }
    printf("\n");

    printf("\"FLAT ARRAY min  Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_min[ i_test ] );
    }
    printf("\n");
  }

  /*----------------------------------------------------------------------*/
  {
    double dt_min[ num_test ];
    double dt_max[ num_test ];
    double dt_mean[ num_test ];
    double mflops_min[ num_test ];
    double mflops_max[ num_test ];
    double mflops_mean[ num_test ];

    unsigned i_test = 0 ;

    for ( ; i_test < num_test ; ++i_test ) {
      const unsigned num_array = num_test_array[ i_test ];
      const unsigned num_sets  = max_array / num_array ;

      const double mflop_cycle =
        ((double)( 2 * num_array * length_array )) / 1.0e6 ;

      const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

      data.length       = length_array ;
      data.number       = num_array ;
      data.stride       = stride_array ;
      data.chunk_length = length_chunk ;

      { unsigned i = 0 ; for ( ; i < size_alloc ; ++i ) { array[i] = 0 ; } }

      {
        double dt_tmp ;
        unsigned i ;
        unsigned repeat = 0 ;
        for ( ; repeat < num_trials ; ++repeat ) {

          dt_tmp = TPI_Walltime();
          for ( i = 0 ; i < ncycle ; ++i ) {
            data.array = array + stride_array * num_array * ( i % num_sets );
            TPI_Run( & test_dnax_column_work , & data , num_chunk , 0 );
          }
          dt_tmp = TPI_Walltime() - dt_tmp ;

          if ( 0 == repeat ) {
            dt_min[ i_test ] = dt_tmp ;
            dt_max[ i_test ] = dt_tmp ;
            dt_mean[ i_test ] = dt_tmp / num_trials ;
          }
          else {
            dt_mean[ i_test ] += dt_tmp / num_trials ;
          }
          if ( dt_tmp < dt_min[ i_test ] ) { dt_min[ i_test ] = dt_tmp ; }
          if ( dt_tmp > dt_max[ i_test ] ) { dt_max[ i_test ] = dt_tmp ; }
        }
      }

      mflops_max[ i_test ] = mflop_cycle * ncycle / dt_min[ i_test ];
      mflops_min[ i_test ] = mflop_cycle * ncycle / dt_max[ i_test ];
      mflops_mean[ i_test ] = mflop_cycle * ncycle / dt_mean[ i_test ];
    }

    printf("\"CHUNK COLUMN max  time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_max[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK COLUMN mean time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_mean[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK COLUMN min  time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_min[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK COLUMN max  Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_max[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK COLUMN mean Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_mean[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK COLUMN min  Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_min[ i_test ] );
    }
    printf("\n");
  }

  /*----------------------------------------------------------------------*/
  {
    double dt_min[ num_test ];
    double dt_max[ num_test ];
    double dt_mean[ num_test ];
    double mflops_min[ num_test ];
    double mflops_max[ num_test ];
    double mflops_mean[ num_test ];

    unsigned i_test = 0 ;

    for ( ; i_test < num_test ; ++i_test ) {
      const unsigned num_array = num_test_array[ i_test ];
      const unsigned num_sets  = max_array / num_array ;

      const double mflop_cycle =
        ((double)( 2 * num_array * length_array )) / 1.0e6 ;

      const unsigned ncycle = 1 + (unsigned)( Mflop_target / mflop_cycle );

      data.length       = length_array ;
      data.number       = num_array ;
      data.stride       = stride_array ;
      data.chunk_length = length_chunk ;

      { unsigned i = 0 ; for ( ; i < size_alloc ; ++i ) { array[i] = 0 ; } }

      {
        double dt_tmp ;
        unsigned i ;
        unsigned repeat = 0 ;
        for ( ; repeat < num_trials ; ++repeat ) {

          dt_tmp = TPI_Walltime();
          for ( i = 0 ; i < ncycle ; ++i ) {
            data.array = array + stride_array * num_array * ( i % num_sets );
            TPI_Run( & test_dnax_row_work , & data , num_chunk , 0 );
          }
          dt_tmp = TPI_Walltime() - dt_tmp ;

          if ( 0 == repeat ) {
            dt_min[ i_test ] = dt_tmp ;
            dt_max[ i_test ] = dt_tmp ;
            dt_mean[ i_test ] = dt_tmp / num_trials ;
          }
          else {
            dt_mean[ i_test ] += dt_tmp / num_trials ;
          }
          if ( dt_tmp < dt_min[ i_test ] ) { dt_min[ i_test ] = dt_tmp ; }
          if ( dt_tmp > dt_max[ i_test ] ) { dt_max[ i_test ] = dt_tmp ; }
        }
      }

      mflops_max[ i_test ] = mflop_cycle * ncycle / dt_min[ i_test ];
      mflops_min[ i_test ] = mflop_cycle * ncycle / dt_max[ i_test ];
      mflops_mean[ i_test ] = mflop_cycle * ncycle / dt_mean[ i_test ];
    }

    printf("\"CHUNK ROW max  time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_max[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK ROW mean time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_mean[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK ROW min  time (sec)\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", dt_min[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK ROW max  Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_max[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK ROW mean Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_mean[ i_test ] );
    }
    printf("\n");

    printf("\"CHUNK ROW min  Mflops\"");
    for ( i_test = 0 ; i_test < num_test ; ++i_test ) {
      printf(" , %lf", mflops_min[ i_test ] );
    }
    printf("\n");
  }

  /*----------------------------------------------------------------------*/

  free( array );
}

/*------------------------------------------------------------------------*/

int test_c_tpi_dnax( int nthread )
{
  const unsigned Mflop_target = 10 ;
  const unsigned num_array[6] = { 2 , 5 , 10 , 20 , 50 , 100 };

  TPI_Init( nthread );

  test_tpi_dnax_driver( nthread ,
                        Mflop_target * nthread ,
                        5         /* number trials */ ,
                        6         /* number of tests */ ,
                        num_array /* number of arrays for each test */ ,
                        1e6       /* array computation length */ ,
                        1000      /* chunk length */ );

  TPI_Finalize();

  return 0 ;
}



