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

#include <unistd.h>
#include <sys/types.h>

#include <TPI.h>

#ifndef TPI_NO_SCHED
#define TPI_NO_SCHED 0
#endif

/*--------------------------------------------------------------------*/

int TPI_Concurrency()
{
#if TPI_NO_SCHED
  return 0 ;
#else
  enum { NTMP = 8 };

  extern int sched_getaffinity( pid_t, unsigned int , unsigned long * );

  int count = 0 ;
  unsigned long tmp[ NTMP ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  if ( ! sched_getaffinity( 0 , sizeof(tmp) , tmp ) ) {

    int i ;

    for ( i = 0 ; i < NTMP ; ++i ) {
      unsigned long t = tmp[i] ;
      for ( ; t ; t >>= 1 ) { if ( t & 01 ) { ++count ; } }
    }
  }
  return count ;
#endif
}

