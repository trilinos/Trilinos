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

static void test_tpi_noop( TPI_Work * );

int test_c_tpi_noop( int num_test , int * num_thread )
{
  void * const p = NULL ;
  const unsigned n_work = 10000 ;
  const unsigned n_trial = 100 ;
  int itest ;

  fprintf(stdout,"\n\"TPI_Run(noop) test\"\n");
  fprintf(stdout,"\"NUMBER OF TRIALS\" , %u\n", n_trial );
  fprintf(stdout,"\"WORK COUNT\" , %u\n", n_work );
  fprintf(stdout,
          "\"# Thread\" , \"Min microsec\" , \"Mean microsec\" , \"Max microsec\"\n");

  for ( itest = 0 ; itest < num_test ; ++itest ) {
    const unsigned num = num_thread[ itest ];

    unsigned j ;
    double dt_min = 0 , dt_max = 0 , dt_mean = 0 ;
    double dt ;

    /* Time many tries and trials, get the min and max time */

    for ( j = 0 ; j < n_trial ; ++j ) {
  
      TPI_Init( num );

      dt = TPI_Walltime();
      TPI_Run( & test_tpi_noop, p, n_work, 0 );
      dt = TPI_Walltime() - dt ;

      dt_mean += dt ;

      if ( ! j ) {
        dt_min = dt_max = dt ;
      }

      if ( dt < dt_min ) { dt_min = dt ; }
      if ( dt > dt_max ) { dt_max = dt ; }

      TPI_Finalize();
    }

    dt_min  *= 1.0e6 ;
    dt_mean *= 1.0e6 / n_trial ;
    dt_max  *= 1.0e6 ;

    fprintf(stdout, "%u , %g , %g , %g\n", num, dt_min, dt_mean, dt_max );

    fflush(stdout);
  }

  return 0 ;
}

static void test_tpi_noop( TPI_Work * work )
{
  void * p = work->shared ;
  while ( p ) { p = *((void**)p) ; } return ;
}

/*------------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

struct TestTPI {
  int total ;
  int count ;
};

static void test_tpi_loop( TPI_Work * work )
{
  struct TestTPI * const data = (struct TestTPI *) work->shared ;

  TPI_Lock( 0 );
  data->count += 1 ;
  TPI_Unlock( 0 );

  return ;
}

int test_c_tpi_single( int nthread )
{
  static const char name[] = "test_c_tpi_single" ;
  int result = 0 ;

  struct TestTPI data = { 10000 /* 1000000 */ , 0 };

  TPI_Init( nthread );

  fprintf(stdout,"\"%s[%d] starting...",name,nthread);
  fflush(stdout);

  /*--------------------------------*/

  {
    int n ;
    for ( n = 1 ; 0 <= result && n < 64 ; ++n ) {
      data.count = 0 ;
      result = TPI_Run( & test_tpi_loop , & data , data.total , 1 );
      if ( result < 0 ) {
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

  fprintf(stdout,"completed successfully\"\n");
  fflush(stdout);

  TPI_Finalize();

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

