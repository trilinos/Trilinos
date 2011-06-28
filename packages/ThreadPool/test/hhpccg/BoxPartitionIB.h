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



#ifndef BoxPartionIB_h
#define BoxPartionIB_h

/** \brief  Partition a { [ix,jx) X [iy,jy) X [iz,jz) } box.
 *
 *  Use recursive coordinate bisection to partition a box 
 *  into np disjoint sub-boxes.  Allocate (via malloc) and
 *  populate the sub-boxes, mapping the local (x,y,z) to
 *  a local ordinal, and mappings for the send-recv messages
 *  to update the ghost cells.
 *
 *  Order local ordinates as follows:
 *    {
 *      interior ,
 *      boundary ,
 *      remote[ ( my_p + i ) % np ]
 *    } 
 *      where i = 1..(np-1)
 *
 *  usage:
 *
 *  my_nx = pbox[my_p][0][1] - pbox[my_p][0][0] ;
 *  my_ny = pbox[my_p][1][1] - pbox[my_p][1][0] ;
 *  my_nz = pbox[my_p][2][1] - pbox[my_p][2][0] ;
 *
 *  for ( x = -ghost ; x < my_nx + ghost ; ++x ) {
 *  for ( y = -ghost ; y < my_ny + ghost ; ++y ) {
 *  for ( z = -ghost ; z < my_nz + ghost ; ++z ) {
 *    const int x_global = x + pbox[my_p][0][0] ;
 *    const int y_global = y + pbox[my_p][1][0] ;
 *    const int z_global = z + pbox[my_p][2][0] ;
 *
 *    const int local_ordinal =
 *      box_map_local( pbox[my_p], ghost, map_local_id, x, y, z );
 *
 *    if ( 0 <= local_ordinal ) {
 *    }
 *  }
 *  
 *  for ( i = 1 ; i < np ; ++i ) {
 *    const int recv_processor = ( my_p + i ) % np ;
 *    const int recv_ordinal_begin = map_recv_pc[i];
 *    const int recv_ordinal_end   = map_recv_pc[i+1];
 *  }
 *
 *  for ( i = 1 ; i < np ; ++i ) {
 *    const int send_processor = ( my_p + i ) % np ;
 *    const int send_map_begin = map_send_pc[i];
 *    const int send_map_end   = map_send_pc[i+1];
 *    for ( j = send_map_begin ; j < send_map_end ; ++j ) {
 *      send_ordinal = map_send_id[j] ;
 *    }
 *  }
 */


void box_partition_rcb(
  const int np              /**< [in] Number of partitions */ ,
  const int root_box[3][2]  /**< [in] Global 3D box to partition  */ ,
  int       pbox[][3][2]    /**< [out] Partition of global 3D boxes */ );

void box_partition_map(
  const int np            /**< [in] Number of partitions */ ,
  const int my_p          /**< [in] My partition */ ,
  const int gbox[3][2]    /**< [in] Global 3D box */ ,
  const int pbox[][3][2]  /**< [in] Partitions of global 3D box */ ,
  const int ghost         /**< [in] Number of grid points to ghost */ ,

  int    map_uses_box[3][2]  /**< [out] Local box expanded by ghosting */ ,
  int    map_local_id[]      /**< [out] Mapping for local points */ ,
  int *  map_count_interior  /**< [out] Number of my interior points */ ,
  int *  map_count_owns      /**< [out] Number of points I own */ ,
  int *  map_count_uses      /**< [out] Number of points I access */ ,
  int ** map_recv_pc         /**< [out] Received prefix spans per process */ ,
  int ** map_send_pc         /**< [out] Send prefix counts per process */ ,
  int ** map_send_id         /**< [out] Send grid points */ );

/* \brief  Map a global (x,y,z) to a local ordinal.  */
int box_map_local( const int local_uses[3][2] ,
                   const int map_local_id[] ,
                   const int global_x ,
                   const int global_y ,
                   const int global_z );

#endif

