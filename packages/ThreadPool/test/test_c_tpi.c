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
 */

#include <stdio.h>
#include <TPI.h>


/*------------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void test_tpi_noop( void * , TPI_ThreadPool );

int test_c_tpi_noop( int num_test , int * num_thread )
{
  const unsigned n = 1e5 ;
  const unsigned n_trial = 5 ;
  int itest ;

  fprintf(stdout,"\n\"TPI_Run(noop) test\"\n");
  fprintf(stdout,
          "\"# Thread\" , \"Min microsec\" , \"Mean microsec\" , \"Max microsec\" , \"Min Use #\" , \"Mean Use #\" , \"Max Use #\"\n");

  for ( itest = 0 ; itest < num_test ; ++itest ) {
    const unsigned num = num_thread[ itest ];

    int count[ num ];

    unsigned i , j ;
    double dt ;
    double dt_min = 0 , dt_max = 0 ;
    double use_min = 0 , use_max = 0 ;
    double dt_mean = 0 , use_mean = 0 ;

    TPI_Init( num );

    /* Time many tries and trials, get the min and max time */

    TPI_Run( & test_tpi_noop , NULL );

    for ( j = 0 ; j < num ; ++j ) { count[j] = 0 ; }

    for ( j = 0 ; j < n_trial ; ++j ) {
  
      TPI_Run_count( num , NULL );

      dt = TPI_Walltime();
      for ( i = 0 ; i < n ; ++i ) { TPI_Run( & test_tpi_noop , NULL ); }
      dt = TPI_Walltime() - dt ;

      TPI_Run_count( num , count );

      dt_mean += dt ;

      if ( ! j ) {
        dt_min = dt_max = dt ;
        use_min = use_max = count[0] ;
      }

      if ( dt < dt_min ) { dt_min = dt ; }
      if ( dt > dt_max ) { dt_max = dt ; }
      for ( i = 0 ; i < num ; ++i ) {
        use_mean += count[i] ;
        if ( use_min > count[i] ) { use_min = count[i] ; }
        if ( use_max < count[i] ) { use_max = count[i] ; }
      }
    }

    dt_min  *= 1.0e6 / n ;
    dt_mean *= 1.0e6 / ( n * n_trial );
    dt_max  *= 1.0e6 / n ;

    use_min  /= n ;
    use_mean /= n * n_trial * num ;
    use_max  /= n ;

    fprintf(stdout,
            "%u , %g , %g , %g , %g , %g , %g\n",
            num, dt_min, dt_mean, dt_max, use_min, use_mean, use_max );

    fflush(stdout);

    TPI_Finalize();
  }

  return 0 ;
}

static void test_tpi_noop( void * arg , TPI_ThreadPool pool )
{ while ( arg && pool ) { arg = *((void**)arg) ; } return ; }

/*------------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

struct TestTPI {
  int total ;
  int count ;
};

static void test_tpi_loop( void * arg , TPI_ThreadPool pool )
{
  static const char name[] = "test_tpi_loop" ;

  struct TestTPI * const data = (struct TestTPI *) arg ;

  int size  = -1 ;
  int rank  = -1 ;
  int result = 0 ;
  int begin  = 0 ;
  int number = 0 ;

  if ( ( result = TPI_Rank( pool , & rank , & size ) ) ) {
    fprintf(stderr,"\n%s: TPI_Pool_rank = %d\n",name,result);
  }

  if ( ( result = TPI_Partition(rank,size,data->total, & begin, & number) ) ) {
    fprintf(stderr,"\n%s: TPI_Partition = %d\n",name,result);
  }
  else {
    int count = 0 ;
    int i , j ;
    for ( j = 0 ; j < 101 ; ++j ) {
      if ( j % 2 ) {
        for ( i = 0 ; i < number ; ++i ) { --count ; }
      }
      else {
        for ( i = 0 ; i < number ; ++i ) { ++count ; }
      }
    }

    TPI_Lock( pool , 0 );
    data->count += count ;
    TPI_Unlock( pool , 0 );
  }

  return ;
}

int test_c_tpi_single( int size )
{
  static const char name[] = "test_c_tpi_single" ;
  int result = 0 ;

  struct TestTPI data = { 10000 /* 1000000 */ , 0 };

  TPI_Init( size );

  fprintf(stdout,"\"%s[%d] starting...",name,size);
  fflush(stdout);

  /*--------------------------------*/

  {
    int n ;
    for ( n = 1 ; n < 64 ; ++n ) {
      data.count = 0 ;
      if ( ( result = TPI_Set_lock_size( 1 ) ) ) {
        fprintf(stderr,"\n%s: TPI_Set_lock_size = %d\n",name,result);
      }
      if ( ( result = TPI_Run( & test_tpi_loop , & data ) ) ) {
        fprintf(stderr,"\n%s: TPI_Run = %d\n",name,result);
      }
      else {
        if ( ( result = data.count != data.total ) ) {
          fprintf(stderr,"\n%s: test_tpi_loop : %d != %d\n",name,
                  data.count , data.total );
        }
      }
    }
  }

  /*--------------------------------*/

  {
    TPI_parallel_subprogram func = & test_tpi_loop ;
    void * ptr = & data ;
    int n ;
    for ( n = 1 ; n < 64 ; ++n ) {
      data.count = 0 ;

      if ( ( result = TPI_Set_lock_size( 1 ) ) ) {
        fprintf(stderr,"\n%s: TPI_Set_lock_size = %d\n",name,result);
      }
      if ( ( result = TPI_Run_many( 1 , & func , & ptr , & n ) ) ) {
        fprintf(stderr,"\n%s: TPI_Run_many[%d] = %d\n",name,n,result);
      }
      else {
        if ( ( result = data.count != data.total ) ) {
          fprintf(stderr,"\n%s: test_tpi_loop[%d] : %d != %d\n",name,n,
                  data.count , data.total );
        }
      }
    }
  }

  /*--------------------------------*/

  fprintf(stdout,"completed successfully\"\n");
  fflush(stdout);

  TPI_Finalize();

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

struct TestTPIMany {
  int * flag ;
  int size ;
};

static void test_tpi_many_one( void * arg , TPI_ThreadPool pool )
{
  struct TestTPIMany * const data = (struct TestTPIMany *) arg ;

  int rank ;
  int size ;
  int result ;

  if ( ( result = TPI_Rank( pool , & rank , & size ) ) ) {
    fprintf(stderr,"\ntest_tpi_many_one failed TPI_Rank = %d\n",result);
  }
  else if ( rank < 0 || size <= rank || size != data->size ) {
    fprintf(stderr,"\ntest_tpi_many_one failed rank = %d, size = %d, data->size = %d\n",
           rank , size , data->size );
  }
  else {
    data->flag[rank] += 1 ;
  }
}

static void test_tpi_many_two( void * arg , TPI_ThreadPool pool )
{
  struct TestTPIMany * const data = (struct TestTPIMany *) arg ;

  int rank ;
  int size ;
  int result ;

  if ( ( result = TPI_Rank( pool , & rank , & size ) ) ) {
    fprintf(stderr,"\ntest_tpi_many_two failed TPI_Rank = %d\n",result);
  }
  else if ( rank < 0 || size <= rank || size != data->size ) {
    fprintf(stderr,"\ntest_tpi_many_two failed rank = %d, size = %d, data->size = %d\n",
           rank , size , data->size );
  }
  else {
    data->flag[rank] += 2 ;
  }
}

static void test_tpi_many_three( void * arg , TPI_ThreadPool pool )
{
  struct TestTPIMany * const data = (struct TestTPIMany *) arg ;

  int rank ;
  int size ;
  int result ;

  if ( ( result = TPI_Rank( pool , & rank , & size ) ) ) {
    fprintf(stderr,"\ntest_tpi_many_three failed TPI_Rank = %d\n",result);
  }
  else if ( rank < 0 || size <= rank || size != data->size ) {
    fprintf(stderr,"\ntest_tpi_many_three failed rank = %d, size = %d, data->size = %d\n",
           rank , size , data->size );
  }
  else {
    data->flag[rank] += 3 ;
  }
}

int test_c_tpi_many( int size )
{
  enum { N = 100 };
  static const char name[] = "test_c_tpi_many" ;

  struct TestTPIMany data[ 3 ] ;

  void * pointers[3] = { & data[0] , & data[1] , & data[2] };

  TPI_parallel_subprogram func[3] =
    { & test_tpi_many_one , & test_tpi_many_two , & test_tpi_many_three };

  int length[ 3 ] = { 0 , 0 , 0 };

  int flags[ N ];

  int result = 0 ;

  int i , j ;

  fprintf(stdout,"\"%s[%d] starting...",name,size);
  fflush(stdout);

  TPI_Init( size );

  for ( i = 1 ; i < N - 1 ; ++i ) {

    const int two_begin   = i ;
    const int three_begin = i + ( N - i ) / 2 ;

    const int one_length   = two_begin ;
    const int two_length   = three_begin - two_begin ;
    const int three_length = N - three_begin ;

    for ( j = 0 ; j < N ; ++j ) { flags[j] = 0 ; }

    data[0].size = length[0] = one_length ;
    data[1].size = length[1] = two_length ;
    data[2].size = length[2] = three_length ;

    data[0].flag = flags ;
    data[1].flag = flags + two_begin ;
    data[2].flag = flags + three_begin ;

    if ( ( result = TPI_Run_many( 3 , func , pointers , length ) ) ) {
      fprintf(stderr,"\n%s: TPI_Run_many(%d,%d,%d) = %d\n",
              name,length[0],length[1],length[2],result);
    }
    else {
      for ( j = 0 ; j < two_begin ; ++j ) {
        if ( flags[j] != 1 ) {
          printf("test_tpi_many failed at = %d\n",j);
        }
      }
      for ( ; j < three_begin ; ++j ) {
        if ( flags[j] != 2 ) {
          printf("test_tpi_many failed at = %d\n",j);
        }
      }
      for ( ; j < N ; ++j ) {
        if ( flags[j] != 3 ) {
          printf("test_tpi_many failed at = %d\n",j);
        }
      }
    }
  }

  TPI_Finalize();

  fprintf(stdout,"completed successfully\"\n");
  fflush(stdout);

  return result ;
}

/*--------------------------------------------------------------------*/

