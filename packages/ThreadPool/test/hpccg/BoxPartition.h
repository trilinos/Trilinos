
/** \brief  Partition a { [ix,jx) X [iy,jy) X [iz,jz) } box.
 *
 *  Use recursive coordinate bisection to partition a box 
 *  into np disjoint sub-boxes.  Allocate (via malloc) and
 *  populate the sub-boxes, mapping the local (x,y,z) to
 *  a local ordinal, and mappings for the send-recv messages
 *  to update the ghost cells.
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
  const int np            /**< [in]  Number of partitions */ ,
  const int my_p          /**< [in]  My partition rank    */ ,
  const int root_box[][2] /**< [in]  3D Box to partition  */ ,
  const int ghost         /**< [in]  Ghost cell boundary  */ ,
  int (**pbox)[3][2]      /**< [out] Partition's 3D boxes */ ,
  int ** map_local_id     /**< [out] Map local cells */ ,
  int ** map_recv_pc      /**< [out] Receive spans per processor */ ,
  int ** map_send_pc      /**< [out] Send prefix counts per processor */ ,
  int ** map_send_id      /**< [out] Send message ordinals */ );

/* \brief  Map a local (x,y,z) to a local ordinal.
 */
int box_map_local( const int box_local[][2] ,
                   const int ghost ,
                   const int map_local_id[] ,
                   const int local_x ,
                   const int local_y ,
                   const int local_z );

