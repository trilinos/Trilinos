

#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>

static void test_task_unit( TPI_Work * );
static void test_task_lock( TPI_Work * );

int test_c_tpi_unit( const int nthread , const int nwork )
{
  const int ntrial = 10 ;
  int result ;

  result = TPI_Init( nthread );

  printf("  %d = TPI_Init( %d )\n",result,nthread);

  /* Unit test */
  {
    int flags[ nwork ];
    int i , k ;
    double dt = 0 ;
    for ( i = 0 ; i < ntrial ; ++i ) {
      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }
      {
        const double t = TPI_Walltime();
        TPI_Run( test_task_unit , flags , nwork , 0 );
        dt += TPI_Walltime() - t ;
      }
      for ( k = 0 ; k < nwork ; ++k ) {
        if ( flags[k] != 1 ) {
          printf("  test_task_unit[flag] failed at trial = %d\n",i);
          abort();
        }
      }
    }
    dt /= ntrial ;
    printf("  test_task_unit[flag] passed, work = %d , mean time = %f\n",nwork,dt);
  }

  /* Test locking */
  {
    const int nlock = 2 ;
    int i ;
    double dt = 0 ;

    for ( i = 0 ; i < ntrial ; ++i ) {
      int ncount = 0 ;
      const double t = TPI_Walltime();
      TPI_Run( test_task_lock , & ncount , nwork , nlock );
      dt += TPI_Walltime() - t ;
      if ( ncount != nwork * nlock ) {
        printf("  test_task_unit[lock] failed at trial = %d\n",i);
        abort();
      }
    }
    dt /= ntrial ;
    printf("  test_task_unit[lock] passed, work = %d , mean time = %f\n",nwork,dt);
  }

  TPI_Finalize();

  return 0 ;
}

static void test_task_unit( TPI_Work * task )
{
  const int nloop = 10000 ;
  int * const flag = (int *)( task->shared );
  int i ;
  for ( i = 0 ; i < nloop ; ++i ) {
    flag[ task->work_rank ] += 1 ;
  }
  flag[ task->work_rank ] /= nloop ;
  
}

static void test_task_lock( TPI_Work * work )
{
  int * const ncount = (int*)( work->shared );
  int i ;
  for ( i = 0 ; i < work->lock_count ; ++i ) {
    TPI_Lock( 0 );
    ++*ncount ;
    TPI_Unlock( 0 );
  }
  return ;
}

