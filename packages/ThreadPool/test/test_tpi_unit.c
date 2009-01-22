

#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>

static void test_task_unit( TPI_Work * );
static void test_task_lock( TPI_Work * );

int test_c_tpi_unit( const int nthread , const int nwork )
{
  const int ntrial  = 100 ;
  int result ;

  result = TPI_Init( nthread );

  printf("  %d = TPI_Init( %d )\n",result,nthread);

  /* Unit test */
  {
    int flags[ nwork ];
    int i , k ;
    int np = 0 ;
    for ( i = 0 ; i < ntrial ; ++i ) {
      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }
      np += TPI_Run( test_task_unit , flags , nwork , 0 );
      for ( k = 0 ; k < nwork ; ++k ) {
        if ( flags[k] != 1 ) {
          printf("  test_task_unit failed at trial = %d\n",i);
          abort();
        }
      }
    }
    np /= ntrial ;
    printf("  test_task_unit passed, mean parallelism = %d\n",np);
  }

  /* Test locking */
/*
  {
    const int nlock = 2 ;
    int i , j , k ;

    for ( i = 0 ; i < ntrial ; ++i ) {
      for ( j = 1 ; j <= nwork ; ++j ) {
        for ( k = 1 ; k <= nlock ; ++k ) {
          int ncount = 0 ;
          TPI_Run( test_task_lock , & ncount , nwork , nlock );
        }
      }
    }
  }
*/

  TPI_Finalize();

  return 0 ;
}

static void test_task_unit( TPI_Work * task )
{
  const int nloop = 1000 ;
  int * const flag = (int *)( task->shared );
  int i ;
  for ( i = 0 ; i < nloop ; ++i ) {
    flag[ task->work_rank ] += 1 ;
  }
  flag[ task->work_rank ] /= nloop ;
  
}

static void test_task_lock( TPI_Work * work )
{}

