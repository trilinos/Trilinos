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



#include <stdio.h>
#include <stdlib.h>

#include <BoxPartitionIB.h>

/*--------------------------------------------------------------------*/
/* Recursively split a box into into (up-ip) sub-boxes */

typedef const int RangeInput[2] ;
typedef       int RangeOutput[2] ;
typedef RangeInput  * const BoxInput ;
typedef RangeOutput * const BoxOutput ;

static 
void box_partition( int ip , int up , int axis ,
                    BoxInput box ,
                    int (* const p_box)[3][2] )
{
  const int np = up - ip ;
  if ( 1 == np ) {
    p_box[ip][0][0] = box[0][0] ; p_box[ip][0][1] = box[0][1] ;
    p_box[ip][1][0] = box[1][0] ; p_box[ip][1][1] = box[1][1] ;
    p_box[ip][2][0] = box[2][0] ; p_box[ip][2][1] = box[2][1] ;
  }
  else {
    const int n = box[ axis ][1] - box[ axis ][0] ;
    const int np_low = np / 2 ;  /* Rounded down */
    const int np_upp = np - np_low ;

    const int n_upp = (int) (((double) n) * ( ((double)np_upp) / ((double)np)));
    const int n_low = n - n_upp ;
    const int next_axis = ( axis + 2 ) % 3 ;

    if ( np_low ) { /* P = [ip,ip+np_low) */
      int dbox[3][2] ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      dbox[ axis ][1] = dbox[ axis ][0] + n_low ;

      box_partition( ip, ip + np_low, next_axis,
                     (const int (*)[2]) dbox, p_box );
    }

    if ( np_upp ) { /* P = [ip+np_low,ip+np_low+np_upp) */
      int dbox[3][2] ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      ip += np_low ;
      dbox[ axis ][0] += n_low ;
      dbox[ axis ][1]  = dbox[ axis ][0] + n_upp ;

      box_partition( ip, ip + np_upp, next_axis,
                     (const int (*)[2]) dbox, p_box );
    }
  }
}

void box_partition_rcb( const int np ,
                        const int root_box[3][2] ,
                        int    pbox[][3][2] )
{
  box_partition( 0 , np , 2 , root_box , pbox );
}

/*--------------------------------------------------------------------*/

static int box_intersect( BoxInput a , BoxInput b , BoxOutput c )
{
  int i ;
  for ( i = 0 ; i < 3 ; ++i ) {
    c[i][0] = a[i][0] < b[i][0] ? b[i][0] : a[i][0] ;
    c[i][1] = a[i][1] < b[i][1] ? a[i][1] : b[i][1] ;
  }

  return c[0][0] < c[0][1] && c[1][0] < c[1][1] && c[2][0] < c[2][1] ;
}
 

/*--------------------------------------------------------------------*/

static void global_to_use_box( BoxInput gbox ,
                               BoxInput pbox ,
                               const int ghost ,
                                     BoxOutput interiorBox ,
                                     BoxOutput useBox )
{
  int i = 0 ;

  for ( i = 0 ; i < 3 ; ++i ) {
    const int n = pbox[i][1] - pbox[i][0] ;

    if ( n < 0 ) {
      abort();
    }

    interiorBox[i][0] = gbox[i][0] == pbox[i][0]
                      ? gbox[i][0] :  pbox[i][0] + ghost ;

    interiorBox[i][1] = gbox[i][1] == pbox[i][1]
                      ? gbox[i][1] :  pbox[i][1] - ghost ;

    if ( interiorBox[i][1] < pbox[i][0] ) {
      interiorBox[i][1] = pbox[i][0] ;
    }

    if ( interiorBox[i][0] > pbox[i][1] ) {
      interiorBox[i][0] = pbox[i][1] ;
    }

    if ( interiorBox[i][1] < interiorBox[i][0] ) {
      interiorBox[i][1] = interiorBox[i][0] ;
    }

    useBox[i][0] = pbox[i][0] - ghost ;
    useBox[i][1] = pbox[i][1] + ghost ;

    if ( useBox[i][0] < gbox[i][0] ) { useBox[i][0] = gbox[i][0] ; }
    if ( useBox[i][1] > gbox[i][1] ) { useBox[i][1] = gbox[i][1] ; }
  }
}


/*  A use-box is the owned box plus the ghost layers.
 *  Map a global (x,y,z) to a local integer ordinate.
 */
static int map_global_to_use_box( BoxInput useBox ,
                                  const int global_x ,
                                  const int global_y ,
                                  const int global_z )
{
  const int nx = useBox[0][1] - useBox[0][0] ;
  const int ny = useBox[1][1] - useBox[1][0] ;
  const int nz = useBox[2][1] - useBox[2][0] ;
  const int ix = global_x     - useBox[0][0] ;
  const int iy = global_y     - useBox[1][0] ;
  const int iz = global_z     - useBox[2][0] ;

  const int good = 0 <= ix && ix < nx &&
                   0 <= iy && iy < ny &&
                   0 <= iz && iz < nz ;

  if ( nx < 0 || ny < 0 || nz < 0 ) {
    abort();
  }
  if ( ! good ) {
    abort();
  }

  return good ? ix + iy * nx + iz * nx * ny : -1 ;
}

int box_map_local( const int local_uses[3][2] ,
                   const int map_local_id[] ,
                   const int global_x ,
                   const int global_y ,
                   const int global_z )
{
  int i = map_global_to_use_box( local_uses , global_x , global_y , global_z );

  if ( 0 <= i ) { i = map_local_id[i] ; }

  return i ;
}


/*--------------------------------------------------------------------*/

static void resize_int( int ** a , int * allocLen , int newLen )
{
  int k = 32;
  while ( k < newLen ) { k <<= 1 ; }
  if ( NULL == *a )
    { *a = malloc( sizeof(int)*(*allocLen = k) ); }
  else if ( *allocLen < k ) 
    { *a = realloc(*a , sizeof(int)*(*allocLen = k)); }
}

void box_partition_map( 
  const int np ,
  const int my_p ,
  const int gbox[3][2] ,
  const int pbox[][3][2] ,
  const int ghost ,

  int    map_use_box[3][2] ,
  int    map_local_id[] ,
  int *  map_count_interior ,
  int *  map_count_owns ,
  int *  map_count_uses ,
  int ** map_recv_pc ,
  int ** map_send_pc ,
  int ** map_send_id )
{
  int * recv_pc = (int *) malloc( ( np + 1 ) * sizeof(int) );
  int * send_pc = (int *) malloc( ( np + 1 ) * sizeof(int) );

  int   id_length = 0 ;

  int * send_id  = NULL ;
  int   send_id_size = 0 ;

  int own_length , use_length , int_length ;
  int count_interior , count_parallel ;
  int iSend ;
  int g_ix , g_iy , g_iz ;
  int i ;

  int my_int_box[3][2] ;

  global_to_use_box( gbox , pbox[my_p] , ghost , my_int_box , map_use_box );

  own_length = ( pbox[my_p][0][1] - pbox[my_p][0][0] ) *
               ( pbox[my_p][1][1] - pbox[my_p][1][0] ) *
               ( pbox[my_p][2][1] - pbox[my_p][2][0] );

  use_length = ( map_use_box[0][1] - map_use_box[0][0] ) *
               ( map_use_box[1][1] - map_use_box[1][0] ) *
               ( map_use_box[2][1] - map_use_box[2][0] );

  int_length = ( my_int_box[0][1] - my_int_box[0][0] ) *
               ( my_int_box[1][1] - my_int_box[1][0] ) *
               ( my_int_box[2][1] - my_int_box[2][0] );

  for ( i = 0 ; i < id_length ; ++i ) { map_local_id[i] = -1 ; }

  /* Fill in locally owned portion: { interior , parallel } */

  count_interior = 0 ;
  count_parallel = int_length ;

  for ( g_iz = pbox[my_p][2][0] ; g_iz < pbox[my_p][2][1] ; ++g_iz ) {
  for ( g_iy = pbox[my_p][1][0] ; g_iy < pbox[my_p][1][1] ; ++g_iy ) {
  for ( g_ix = pbox[my_p][0][0] ; g_ix < pbox[my_p][0][1] ; ++g_ix ) {

    const int local =
      map_global_to_use_box( (BoxInput) map_use_box, g_ix, g_iy, g_iz );

    if ( local < 0 ) { 
      abort();
    }

    if ( my_int_box[2][0] <= g_iz && g_iz < my_int_box[2][1] &&
         my_int_box[1][0] <= g_iy && g_iy < my_int_box[1][1] &&
         my_int_box[0][0] <= g_ix && g_ix < my_int_box[0][1] ) {
      /* Interior */
      map_local_id[ local ] = count_interior++ ;
    }
    else {
      /* Parallel */
      map_local_id[ local ] = count_parallel++ ;
    }
  }
  }
  }

  if ( count_interior != int_length ) { abort(); }
  if ( count_parallel != own_length ) { abort(); }

  /* Fill in off-process received portion: { ( i + my_p ) % np } */

  recv_pc[0] = count_parallel ;
  recv_pc[1] = count_parallel ;
  send_pc[0] = 0 ;
  send_pc[1] = 0 ;
  iSend = 0 ;

  for ( i = 1 ; i < np ; ++i ) {
    const int ip = ( i + my_p ) % np ;
    int recv_box[3][2] ;
    int send_box[3][2] ;
    int other_int_box[3][2] ;
    int other_use_box[3][2] ;

    /* Received portions */

    if ( box_intersect( (BoxInput) map_use_box , (BoxInput) pbox[ip] , recv_box ) ) {

      for ( g_iz = recv_box[2][0] ; g_iz < recv_box[2][1] ; ++g_iz ) {
      for ( g_iy = recv_box[1][0] ; g_iy < recv_box[1][1] ; ++g_iy ) {
      for ( g_ix = recv_box[0][0] ; g_ix < recv_box[0][1] ; ++g_ix ) {

        const int local = map_global_to_use_box( (BoxInput) map_use_box, g_ix, g_iy, g_iz );

        map_local_id[ local ] = count_parallel++ ;
      }
      }
      }
    }
    recv_pc[i+1] = count_parallel ;

    /* Sent items */

    global_to_use_box( gbox, pbox[ip], ghost, other_int_box, other_use_box );

    if ( box_intersect( (BoxInput) other_use_box , (BoxInput) pbox[my_p] , send_box ) ) {

      int nSend = ( send_box[0][1] - send_box[0][0] ) *
                  ( send_box[1][1] - send_box[1][0] ) *
                  ( send_box[2][1] - send_box[2][0] );

      resize_int( & send_id , & send_id_size , (iSend + nSend ) );

      for ( g_iz = send_box[2][0] ; g_iz < send_box[2][1] ; ++g_iz ) {
      for ( g_iy = send_box[1][0] ; g_iy < send_box[1][1] ; ++g_iy ) {
      for ( g_ix = send_box[0][0] ; g_ix < send_box[0][1] ; ++g_ix ) {

        const int local = map_global_to_use_box( (BoxInput) map_use_box, g_ix, g_iy, g_iz );

        if ( map_local_id[ local ] < count_interior ) { abort(); }

        send_id[ iSend ] = map_local_id[ local ] ;
        ++iSend ;
      }
      }
      }
    }
    send_pc[i+1] = iSend ;
  }

  if ( count_parallel != use_length ) { abort(); }

  *map_count_interior = int_length ;
  *map_count_owns     = own_length ;
  *map_count_uses     = use_length ;
  *map_recv_pc        = recv_pc ;
  *map_send_pc        = send_pc ;
  *map_send_id        = send_id ;
}

/*--------------------------------------------------------------------*/

#ifdef UNIT_TEST

static int box_contain( const int a[3][2] , const int b[3][2] )
{
  return a[0][0] <= b[0][0] && b[0][1] <= a[0][1] &&
         a[1][0] <= b[1][0] && b[1][1] <= a[1][1] &&
         a[2][0] <= b[2][0] && b[2][1] <= a[2][1] ;
}

static void box_print( FILE * fp , const int a[][2] )
{
  fprintf(fp,"{ [ %d , %d ) , [ %d , %d ) , [ %d , %d ) }",
                a[0][0] , a[0][1] ,  
                a[1][0] , a[1][1] ,  
                a[2][0] , a[2][1] );
}

static int box_disjoint( BoxInput a , BoxInput b )
{
  return a[0][1] <= b[0][0] || b[0][1] <= a[0][0] ||
         a[1][1] <= b[1][0] || b[1][1] <= a[1][0] ||
         a[2][1] <= b[2][0] || b[2][1] <= a[2][0] ;
}


static void test_box( const int box[3][2] , const int np )
{
  const int ncell_box = box[0][1] * box[1][1] * box[2][1] ;
  int ncell_total = 0 ;
  int ncell_min = ncell_box ;
  int ncell_max = 0 ;
  int (*pbox)[3][2] ;
  int i , j ;

  pbox = (int (*)[3][2]) malloc( sizeof(int) * np * 3 * 2 );

  box_partition( 0 , np , 2 , box , pbox );

  for ( i = 0 ; i < np ; ++i ) {
    const int ncell = ( pbox[i][0][1] - pbox[i][0][0] ) *
                      ( pbox[i][1][1] - pbox[i][1][0] ) *
                      ( pbox[i][2][1] - pbox[i][2][0] );

    if ( ! box_contain( box , (const int (*)[2]) pbox[i] ) ) {
      fprintf(stdout,"  OUT OF BOUNDS pbox[%d/%d] = ",i,np);
      box_print(stdout,(const int (*)[2]) pbox[i]);
      fprintf(stdout,"\n");
      abort();
    }

    for ( j = i + 1 ; j < np ; ++j ) {
      if ( ! box_disjoint( (const int (*)[2]) pbox[i] ,
                           (const int (*)[2]) pbox[j] ) ) {
        fprintf(stdout,"  NOT DISJOINT pbox[%d/%d] = ",i,np);
        box_print(stdout, (const int (*)[2]) pbox[i]);
        fprintf(stdout,"\n");
        fprintf(stdout,"               pbox[%d/%d] = ",j,np);
        box_print(stdout, (const int (*)[2]) pbox[j]);
        fprintf(stdout,"\n");
        abort();
      }
    }
    ncell_total += ncell ;

    if ( ncell_max < ncell ) { ncell_max = ncell ; }
    if ( ncell < ncell_min ) { ncell_min = ncell ; }
  }

  if ( ncell_total != ncell_box ) {
    fprintf(stdout,"  WRONG CELL COUNT NP = %d\n",np);
    abort();
  }
  fprintf(stdout,"NP = %d, total = %d, avg = %d, min = %d, max = %d\n",
          np,ncell_box,ncell_box/np,ncell_min,ncell_max);

  free( pbox );
}

/*--------------------------------------------------------------------*/

static void test_maps( const int root_box[][2] , const int np )
{
  const int ghost = 1 ;
  const int nx_global = root_box[0][1] - root_box[0][0] ;
  const int ny_global = root_box[1][1] - root_box[1][0] ;
  int map_count_interior , map_count_owns , map_count_uses ;
  int map_use_box[3][2] ;
  int ieq , i , j ;
  int (*pbox)[3][2] ;
  int **local_values ;
  int **map_local_id ;
  int **map_recv_pc ;
  int **map_send_pc ;
  int **map_send_id ;
  
  pbox = (int (*)[3][2]) malloc( sizeof(int) * np * 3 * 2 );

  box_partition( 0 , np , 2 , root_box , pbox );

  local_values = (int **) malloc( sizeof(int*) * np );
  map_local_id = (int **) malloc( sizeof(int*) * np );
  map_recv_pc  = (int **) malloc( sizeof(int*) * np );
  map_send_pc  = (int **) malloc( sizeof(int*) * np );
  map_send_id  = (int **) malloc( sizeof(int*) * np );

  /* Set each local value to the global equation number */

  for ( ieq = i = 0 ; i < np ; ++i ) {
    const int (*mybox)[2] = (const int (*)[2]) pbox[i] ;
    const int nx = mybox[0][1] - mybox[0][0] ;
    const int ny = mybox[1][1] - mybox[1][0] ;
    const int nz = mybox[2][1] - mybox[2][0] ;
    int ix , iy , iz ;

    map_local_id[i] = (int *) malloc( sizeof(int) *
                                      ( nx + 2 * ghost ) *
                                      ( ny + 2 * ghost ) *
                                      ( nz + 2 * ghost ) );

    /* Generate the partition maps for this rank */
    box_partition_map( np , i , root_box ,
                        (const int (*)[3][2])  pbox , ghost ,
                        map_use_box ,
                        map_local_id[i] ,
                        & map_count_interior ,
                        & map_count_owns ,
                        & map_count_uses ,
                        & map_recv_pc[i] , 
                        & map_send_pc[i] , & map_send_id[i] );

    if ( map_count_uses != map_recv_pc[i][np] ) { abort(); }

    local_values[i] = (int *) malloc( sizeof(int) * map_count_uses );

    for ( iz = map_use_box[2][0] ; iz < map_use_box[2][1] ; ++iz ) {
    for ( iy = map_use_box[1][0] ; iy < map_use_box[1][1] ; ++iy ) {
    for ( ix = map_use_box[0][0] ; ix < map_use_box[0][1] ; ++ix ) {

      const int igrid = map_global_to_use_box((BoxInput)map_use_box,ix,iy,iz);
      const int ieq   = map_local_id[i][ igrid ];

      if ( 0 <= ieq ) {
        local_values[i][ ieq ] =
          ix + iy * nx_global + iz * nx_global * ny_global ;
      }
    }
    }
    }
  }

  /* Pair-wise compare the local values */
  /* i  == receiving processor rank */
  /* ip == sending   processor rank */
  /* j  == receiving processor data entry for message from 'ip' */
  /* jp == sending   processor data entry for message to   'i' */

  for ( i = 0 ; i < np ; ++i ) {
    for ( j = 1 ; j < np ; ++j ) {
      const int ip = ( i + j ) % np ;
      const int jp = ( i + np - ip ) % np ;
      const int nrecv = map_recv_pc[i] [j+1]  - map_recv_pc[i] [j] ;
      const int nsend = map_send_pc[ip][jp+1] - map_send_pc[ip][jp] ;
      int k ;
      if ( nrecv != nsend ) {
        fprintf(stderr,"P%d recv %d from P%d\n",i,nrecv,ip);
        fprintf(stderr,"P%d send %d to   P%d\n",ip,nsend,i);
        abort();
      }
      for ( k = 0 ; k < nrecv ; ++k ) {
        const int irecv = map_recv_pc[i][j] + k ;
        const int isend = map_send_pc[ip][jp] + k ;
        const int val_irecv = local_values[i][irecv] ;
        const int val_isend = local_values[ip][ map_send_id[ip][isend] ] ;
        if ( val_irecv != val_isend ) {
          fprintf(stderr,"P%d recv[%d] = %d , from P%d\n",i,k,val_irecv,ip);
          fprintf(stderr,"P%d send[%d] = %d , to   P%d\n",ip,k,val_isend,i);
          abort();
        }
      }
    }
  }

  for ( i = 0 ; i < np ; ++i ) {
    free( map_local_id[i] );
    free( map_recv_pc[i] );
    free( map_send_pc[i] );
    free( map_send_id[i] );
    free( local_values[i] );
  }
  free( map_send_id );
  free( map_send_pc );
  free( map_recv_pc );
  free( map_local_id );
  free( local_values );
  free( pbox );
}

/*--------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
  int np_max = 256 ;
  int box[3][2] = { { 0 , 64 } , { 0 , 64 } , { 0 , 64 } };
  int np = 0 ;

  switch( argc ) {
  case 3:
    sscanf(argv[1],"%d",&np);
    sscanf(argv[2],"%dx%dx%d",& box[0][1] , & box[1][1] , & box[2][1] );
    if ( 0 < np ) { test_box(  (const int (*)[2]) box , np ); }
    if ( 0 < np ) { test_maps( (const int (*)[2]) box , np ); }
    break ;
  default:
    for ( np = 1 ; np <= np_max ; ++np ) {
      test_box(  (const int (*)[2]) box , np );
      test_maps( (const int (*)[2]) box , np );
    }
    break ;
  }
  return 0 ;
}

#endif


