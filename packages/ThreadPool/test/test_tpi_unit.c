

#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>

static void test_flag( TPI_Work * );
static void test_lock( TPI_Work * );
static void test_reduce_work( TPI_Work * );
static void test_reduce_init( TPI_Work * );
static void test_reduce_join( TPI_Work * , void * );

int test_c_tpi_unit( const int nthread , const int nwork )
{
  int * const flags = (int *) malloc( nwork * sizeof(int) );
  const int ntrial = 100 ;

  const int result = TPI_Init( nthread );

  if ( result != nthread ) {
    printf("%d != TPI_Init(%d) : FAILED\n",result,nthread);
    abort();
  }

  /* Unit test */
  {
    int i , k ;
    double dt = 0 ;
    for ( i = 0 ; i < ntrial ; ++i ) {
      for ( k = 0 ; k < nthread ; ++k ) { flags[k] = 0 ; }
      {
        const double t = TPI_Walltime();
        TPI_Run_threads( test_flag , flags , 0 );
        dt += TPI_Walltime() - t ;
      }
      for ( k = 0 ; k < nthread ; ++k ) {
        if ( flags[k] != 1 ) {
          printf("  TPI_Run_threads(test_flag) failed at trial = %d\n",i);
          abort();
        }
      }
    }
    dt /= ntrial ;
    printf("\"  TPI_Run_threads(test_flag,*,0) passed: N = %d , mean time = %g\"\n",nthread,dt);
  }

  /* Unit test */
  {
    int i , k ;
    double dt = 0 ;
    for ( i = 0 ; i < ntrial ; ++i ) {
      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }
      {
        const double t = TPI_Walltime();
        TPI_Run( test_flag , flags , nwork , 0 );
        dt += TPI_Walltime() - t ;
      }
      for ( k = 0 ; k < nwork ; ++k ) {
        if ( flags[k] != 1 ) {
          printf("  TPI_Run(test_flag) failed at trial = %d\n",i);
          abort();
        }
      }
    }
    dt /= ntrial ;
    printf("\"  TPI_Run(test_flag,*,%d,0) passed: mean time = %g\"\n",nwork,dt);
  }

  /* Test locking */
  {
    const int nlock = 2 ;
    int i ;
    double dt = 0 ;

    for ( i = 0 ; i < ntrial ; ++i ) {
      int ncount = 0 ;
      const double t = TPI_Walltime();
      TPI_Run( test_lock , & ncount , nwork , nlock );
      dt += TPI_Walltime() - t ;
      if ( ncount != nwork * nlock ) {
        printf("  TPI_Run(test_lock) failed at trial = %d\n",i);
        abort();
      }
    }
    dt /= ntrial ;
    printf("\"  TPI_Run(test_lock,*,%d,%d) passed: mean time = %g\"\n",nwork,nlock,dt);
  }

  /* Test reduction */
  {
    int i ;
    double dt = 0 ;
    for ( i = 0 ; i < ntrial ; ++i ) {
      int work_count = 0 ;
      const double t = TPI_Walltime();
      TPI_Run_reduce( test_reduce_work , NULL , nwork ,
                      test_reduce_join , test_reduce_init ,
                      sizeof(int) , & work_count );
      dt += TPI_Walltime() - t ;
      if ( work_count != nwork ) {
        printf("  TPI_Run_reduce(test_reduce) failed at trial = %d\n",i);
        abort();
      }
    }
    dt /= ntrial ;
    printf("\"  TPI_Run_reduce(test_reduce,NULL,%d,...) passed: mean time = %g\"\n",nwork,dt);
  }

  /* Test blocking */
  {
    int i ;
    double t = 0 , dt = 0 ;
    for ( i = 0 ; i < ntrial ; ++i ) {
      TPI_Block();
      t = TPI_Walltime();
      TPI_Unblock();
      dt += TPI_Walltime() - t ;
    }
    dt /= ntrial ;
    printf("\"  TPI_Block() & TPI_Unblock() passed: mean time = %g\"\n",dt);
  }

  TPI_Finalize();

  free( flags );

  return 0 ;
}

static void test_reduce_work( TPI_Work * work )
{
  int * const count = (int *) work->reduce ;
  ++*count ;
}

static void test_reduce_init( TPI_Work * work )
{
  int * const d = (int *) ( work->reduce );
  *d = 0 ;
}

static void test_reduce_join( TPI_Work * work , void * src )
{
        int * const d = (int *) ( work->reduce );
  const int * const s = (const int *) src ;

  *d += *s ;
}

static void test_flag( TPI_Work * task )
{
  int * const ncount = (int *)( task->info );
  ncount[ task->rank ] += 1 ;
}

static void test_lock( TPI_Work * work )
{
  int * const ncount = (int*)( work->info );
  int i ;
  for ( i = 0 ; i < work->lock_count ; ++i ) {
    TPI_Lock( 0 );
    ++*ncount ;
    TPI_Unlock( 0 );
  }
  return ;
}

