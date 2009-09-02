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

#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>

int test_pthreads_performance( int , int * );
int test_c_tpi_dnax( int , int );
int test_c_tpi_unit( int nthread , int nwork );

int main( int argc , char ** argv )
{
  int num_thread[] = { 1 , 2 , 4 , 6 , 8 , 12 , 16 };
  int num_test = sizeof(num_thread) / sizeof(int);

  const int num_trial = 1 < argc ? atoi( argv[1] ) : 0 ;
  const int nwork     = 2 < argc ? atoi( argv[2] ) : 1000 ;

  int i ;

/*  test_pthreads_performance( num_test , num_thread ); */

  for ( i = 0 ; i < num_test ; ++i ) {
    test_c_tpi_unit( num_thread[i], nwork );
  }

  if ( num_trial ) {
    for ( i = 0 ; i < num_test ; ++i ) {
      test_c_tpi_dnax( num_thread[i] , num_trial );
    }
  }

  return 0 ;
}

