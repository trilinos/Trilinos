/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/* Function Prototypes */

int MPI_Init (int *argc, char **argv[]);
int MPI_Finalize (void);

int MPI_Comm_rank (MPI_Comm comm, int* rank);
int MPI_Comm_size (MPI_Comm comm, int* size);

int MPI_Send (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Recv (void* message, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status);

int MPI_Irecv (void* message, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Isend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);

int MPI_Get_count (MPI_Status* status, MPI_Datatype datatype, int* count);

int MPI_Bsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Rsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Ssend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

int MPI_Comm_create(MPI_Comm comm, MPI_Group new_group, MPI_Comm* new_comm);
int MPI_Comm_group(MPI_Comm comm, MPI_Group* group);
int MPI_Comm_free(MPI_Comm* comm);

int MPI_Abort(MPI_Comm comm, int errorcode);

int MPI_Reduce ( void *sendbuf, void *recvbuf, int count,
   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm );

int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root,
  MPI_Comm comm );

int MPI_Op_create( MPI_User_function* function, int commute, MPI_Op* op);
int MPI_Op_free( MPI_Op *op );

int MPI_Sendrecv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
    int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
    int source, int recvtag, MPI_Comm comm, MPI_Status *status );

int MPI_Wait ( MPI_Request* request, MPI_Status* status);
int MPI_Waitany( int count, MPI_Request array_of_requests[], int* index, MPI_Status* status);
int MPI_Waitall(
        int count, 
        MPI_Request array_of_requests[], 
        MPI_Status array_of_statuses[] );


int MPI_Scatter (void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
    int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Gather ( void *sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

int MPI_Ibsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Irsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Issend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);

/* ----------- */
/* STUBBED OUT */
int MPIO_Test(MPIO_Request *request, int *flag, MPI_Status *status);
int MPIO_Wait(MPIO_Request *request, MPI_Status *status);
int MPI_Address( void *location, MPI_Aint *address);
int MPI_Allgather ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                    void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm );
int MPI_Allgatherv ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                     void *recvbuf, int *recvcounts, int *displs,
                    MPI_Datatype recvtype, MPI_Comm comm );
int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
int MPI_Alltoall( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                 MPI_Comm comm );
int MPI_Alltoallv (
        void *sendbuf,
        int *sendcnts,
        int *sdispls,
        MPI_Datatype sendtype,
        void *recvbuf,
        int *recvcnts,
        int *rdispls,
        MPI_Datatype recvtype,
        MPI_Comm comm );
int MPI_Attr_delete ( MPI_Comm comm, int keyval );
int MPI_Attr_get (
        MPI_Comm comm,
        int keyval,
        void *attr_value,
        int *flag );
int MPI_Attr_put ( MPI_Comm comm, int keyval, void *attr_value );
int MPI_Barrier ( MPI_Comm comm);
int MPI_Bsend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Buffer_attach( void *buffer, int size );
int MPI_Buffer_detach( void *bufferptr, int *size);
int MPI_Cancel( MPI_Request *request );
int MPI_Cart_coords ( MPI_Comm comm, int rank, int maxdims, int *coords );
int MPI_Cart_create ( MPI_Comm comm_old, int ndims, int *dims, int *periods,
                     int reorder, MPI_Comm *comm_cart );
int MPI_Cart_get (
        MPI_Comm comm,
        int maxdims,
        int *dims,
        int *periods,
        int *coords );
int MPI_Cart_map (
        MPI_Comm comm_old,
        int ndims,
        int *dims,
        int *periods,
        int *newrank);
int MPI_Cart_rank (
        MPI_Comm comm,
        int *coords,
        int *rank );
int MPI_Cart_shift ( MPI_Comm comm, int direction, int displ,
                    int *source, int *dest );
int MPI_Cart_sub ( MPI_Comm comm, int *remain_dims, MPI_Comm *comm_new );
int MPI_Cartdim_get ( MPI_Comm comm, int *ndims );

int MPI_Comm_compare (
        MPI_Comm  comm1,
        MPI_Comm  comm2,
        int *result);
int MPI_Comm_dup (
        MPI_Comm comm,
        MPI_Comm *comm_out );
int MPI_Comm_get_name( MPI_Comm comm, char *namep, int *reslen );
int MPI_Comm_remote_group ( MPI_Comm comm, MPI_Group *group );
int MPI_Comm_remote_size ( MPI_Comm comm, int *size );
int MPI_Comm_set_name( MPI_Comm com, char *name );
int MPI_Comm_split ( MPI_Comm comm, int color, int key, MPI_Comm *comm_out );
int MPI_Comm_test_inter ( MPI_Comm comm, int *flag );
int MPI_Dims_create(
        int nnodes,
        int ndims,
        int *dims);

int MPI_Errhandler_create(
        MPI_Handler_function *function,
        MPI_Errhandler       *errhandler);
int MPI_Errhandler_free( MPI_Errhandler *errhandler );
int MPI_Errhandler_get( MPI_Comm comm, MPI_Errhandler *errhandler );
int MPI_Errhandler_set( MPI_Comm comm, MPI_Errhandler errhandler );
int MPI_Error_class(
        int errorcode,
        int *errorclass);
int MPI_Error_string( int errorcode, char *string, int *resultlen );
int MPI_File_close(MPI_File *fh);
MPI_Fint MPI_File_c2f(MPI_File fh);
int MPI_File_delete(char *filename, MPI_Info info);
MPI_File MPI_File_f2c(MPI_Fint fh);
int MPI_File_get_amode(MPI_File fh, int *amode);
int MPI_File_get_atomicity(MPI_File fh, int *flag);
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp);
int MPI_File_get_errhandler(MPI_File fh, MPI_Errhandler *errhandler);
int MPI_File_get_group(MPI_File fh, MPI_Group *group);
int MPI_File_get_info(MPI_File fh, MPI_Info *info_used);
int MPI_File_get_position(MPI_File fh, MPI_Offset *offset);
int MPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset);
int MPI_File_get_size(MPI_File fh, MPI_Offset *size);
int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype,
                             MPI_Aint *extent);
int MPI_File_get_view(MPI_File fh, MPI_Offset *disp, MPI_Datatype *etype,
                MPI_Datatype *filetype, char *datarep);
int MPI_File_iread(MPI_File fh, void *buf, int count,
                   MPI_Datatype datatype, MPIO_Request *request);
int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype,
                      MPIO_Request *request);
int MPI_File_iread_shared(MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPIO_Request *request);
int MPI_File_iwrite(MPI_File fh, void *buf, int count,
                    MPI_Datatype datatype, MPIO_Request *request);
int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, void *buf,
                       int count, MPI_Datatype datatype,
                       MPIO_Request *request);
int MPI_File_iwrite_shared(MPI_File fh, void *buf, int count,
                       MPI_Datatype datatype, MPIO_Request *request);
int MPI_File_open(MPI_Comm comm, char *filename, int amode,
                  MPI_Info info, MPI_File *fh);
int MPI_File_preallocate(MPI_File fh, MPI_Offset size);
int MPI_File_read(MPI_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_all(MPI_File fh, void *buf, int count,
                      MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_all_begin(MPI_File fh, void *buf, int count,
                            MPI_Datatype datatype);
int MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
                   int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype,
                         MPI_Status *status);
int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype);
int MPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_read_ordered(MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count,
                             MPI_Datatype datatype);
int MPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_read_shared(MPI_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);

int MPI_File_set_atomicity(MPI_File fh, int flag);
int MPI_File_set_errhandler(MPI_File fh, MPI_Errhandler errhandler);
int MPI_File_set_info(MPI_File fh, MPI_Info info);
int MPI_File_set_size(MPI_File fh, MPI_Offset size);
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,MPI_Datatype filetype, char *datarep, MPI_Info info);

int MPI_File_sync(MPI_File fh);
int MPI_File_write(MPI_File fh, void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status);
int MPI_File_write_all(MPI_File fh, void *buf, int count,
                       MPI_Datatype datatype, MPI_Status *status);
int MPI_File_write_all_begin(MPI_File fh, void *buf, int count,
                            MPI_Datatype datatype);
int MPI_File_write_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_at(MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype,
                      MPI_Status *status);
int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                          int count, MPI_Datatype datatype,
                          MPI_Status *status);
int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype);
int MPI_File_write_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_ordered(MPI_File fh, void *buf, int count,
                         MPI_Datatype datatype, MPI_Status *status);
int MPI_File_write_ordered_begin(MPI_File fh, void *buf, int count,
                              MPI_Datatype datatype);
int MPI_File_write_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_shared(MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status);
int MPI_Finalized( int *flag );
int MPI_Gatherv ( void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                  void *recvbuf, int *recvcnts, int *displs,
                 MPI_Datatype recvtype,
                  int root, MPI_Comm comm );
int MPI_Get_count(
        MPI_Status *status,
        MPI_Datatype datatype,
        int *count );
int MPI_Get_elements ( MPI_Status *status, MPI_Datatype datatype,
                      int *elements );
int MPI_Get_processor_name(
        char *name,
        int *resultlen);
int MPI_Get_version(
        int *version,
        int *subversion );
int MPI_Graph_create ( MPI_Comm comm_old, int nnodes, int *index, int *edges,
                      int reorder, MPI_Comm *comm_graph );
int MPI_Graph_get ( MPI_Comm comm, int maxindex, int maxedges,
                   int *index, int *edges );
int MPI_Graph_map ( MPI_Comm comm_old, int nnodes, int *index, int *edges,
                   int *newrank );
int MPI_Graph_neighbors ( MPI_Comm comm, int rank, int maxneighbors,
                        int *neighbors );
int MPI_Graph_neighbors_count ( MPI_Comm comm, int rank, int *nneighbors );
int MPI_Graphdims_get ( MPI_Comm comm, int *nnodes, int *nedges );
int MPI_Group_compare ( MPI_Group group1, MPI_Group group2, int *result );
int MPI_Group_difference ( MPI_Group group1, MPI_Group group2,
                         MPI_Group *group_out );
int MPI_Group_excl ( MPI_Group group, int n, int *ranks, MPI_Group *newgroup );
int MPI_Group_free ( MPI_Group *group );
int MPI_Group_incl ( MPI_Group group, int n, int *ranks, MPI_Group *group_out );
int MPI_Group_intersection ( MPI_Group group1, MPI_Group group2,
                           MPI_Group *group_out );
int MPI_Group_range_excl ( MPI_Group group, int n, int ranges[][3],
                         MPI_Group *newgroup );
int MPI_Group_range_incl ( MPI_Group group, int n, int ranges[][3],
                         MPI_Group *newgroup );
int MPI_Group_rank ( MPI_Group group, int *rank );
int MPI_Group_size ( MPI_Group group, int *size );
int MPI_Group_translate_ranks ( MPI_Group group_a, int n, int *ranks_a,
                             MPI_Group group_b, int *ranks_b );
int MPI_Group_union ( MPI_Group group1, MPI_Group group2,
                     MPI_Group *group_out );
MPI_Fint MPI_Info_c2f(MPI_Info info);
int MPI_Info_create(MPI_Info *info);
int MPI_Info_delete(MPI_Info info, char *key);
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
MPI_Info MPI_Info_f2c(MPI_Fint info);
int MPI_Info_free(MPI_Info *info);
int MPI_Info_get(MPI_Info info, char *key, int valuelen, char *value, int *flag);
int MPI_Info_get_nkeys(MPI_Info info, int *nkeys);
int MPI_Info_get_nthkey(MPI_Info info, int n, char *key);
int MPI_Info_get_valuelen(MPI_Info info, char *key, int *valuelen, int *flag);
int MPI_Info_set(MPI_Info info, char *key, char *value);
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided );
int MPI_Initialized( int *flag );
MPI_Handle_type MPI_Int2handle( MPI_Fint f_handle, MPI_Handle_enum handle_kind );
int MPI_Intercomm_create ( MPI_Comm local_comm, int local_leader,
                         MPI_Comm peer_comm, int remote_leader, int tag,
                         MPI_Comm *comm_out );
int MPI_Intercomm_merge ( MPI_Comm comm, int high, MPI_Comm *comm_out );
int MPI_Iprobe( int source, int tag, MPI_Comm comm, int *flag,
               MPI_Status *status );
int MPI_Keyval_create (
        MPI_Copy_function *copy_fn,
        MPI_Delete_function *delete_fn,
        int *keyval,
        void *extra_state );
int MPI_Keyval_free ( int *keyval );
int MPI_Pack ( void *inbuf, int incount, MPI_Datatype datatype,
               void *outbuf, int outcount, int *position, MPI_Comm comm );
int MPI_Pack_size ( int incount, MPI_Datatype datatype, MPI_Comm comm,
                   int *size );
int MPI_Pcontrol( int level );
int MPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status );
int MPI_Recv_init( void *buf, int count, MPI_Datatype datatype, int source,
                  int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Reduce_scatter ( void *sendbuf, void *recvbuf, int *recvcnts,
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
MPI_Fint MPI_Request_c2f( MPI_Request c_request );
int MPI_Request_free( MPI_Request *request );
int MPI_Rsend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Scan ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm );
int MPI_Scatterv (
        void *sendbuf,
        int *sendcnts,
        int *displs,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcnt,
        MPI_Datatype recvtype,
        int root,
        MPI_Comm comm );
int MPI_Send_init( void *buf, int count, MPI_Datatype datatype, int dest,
                  int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Sendrecv_replace( void *buf, int count, MPI_Datatype datatype,
                        int dest, int sendtag, int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status );
int MPI_Ssend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Start(
        MPI_Request *request);
int MPI_Startall( int count, MPI_Request array_of_requests[] );
int MPI_Status_c2f( MPI_Status *c_status, MPI_Fint *f_status );
int MPI_Status_set_cancelled( MPI_Status *status, int flag );
int MPI_Status_set_elements( MPI_Status *status, MPI_Datatype datatype,
                           int count );
int MPI_Test (
        MPI_Request  *request,
        int          *flag,
        MPI_Status   *status);
int MPI_Test_cancelled(
        MPI_Status *status,
        int        *flag);
int MPI_Testall(
        int count,
        MPI_Request array_of_requests[],
        int *flag,
        MPI_Status array_of_statuses[] );
int MPI_Testany(
        int count,
        MPI_Request array_of_requests[],
        int *index, int *flag,
        MPI_Status *status );
int MPI_Testsome(
        int incount,
        MPI_Request array_of_requests[],
        int *outcount,
        int array_of_indices[],
        MPI_Status array_of_statuses[] );
int MPI_Topo_test ( MPI_Comm comm, int *top_type );
int MPI_Type_commit ( MPI_Datatype *datatype );
int MPI_Type_contiguous(
        int count,
        MPI_Datatype old_type,
        MPI_Datatype *newtype);
int MPI_Type_create_darray(int size, int rank, int ndims,
                           int *array_of_gsizes, int *array_of_distribs,
                           int *array_of_dargs, int *array_of_psizes,
                           int order, MPI_Datatype oldtype,
                           MPI_Datatype *newtype);
int MPI_Type_create_indexed_block(
        int count,
        int blocklength,
        int array_of_displacements[],
        MPI_Datatype old_type,
        MPI_Datatype *newtype );
int MPI_Type_create_subarray(
        int ndims,
        int *array_of_sizes,
        int *array_of_subsizes,
        int *array_of_starts,
        int order,
        MPI_Datatype oldtype,
        MPI_Datatype *newtype);
int MPI_Type_extent( MPI_Datatype datatype, MPI_Aint *extent );
int MPI_Type_free ( MPI_Datatype *datatype );
int MPI_Type_get_contents(
        MPI_Datatype datatype,
        int max_integers,
        int max_addresses,
        int max_datatypes,
        int *array_of_integers,
        MPI_Aint *array_of_addresses,
        MPI_Datatype *array_of_datatypes);
int MPI_Type_get_envelope(
        MPI_Datatype datatype,
        int *num_integers,
        int *num_addresses,
        int *num_datatypes,
        int *combiner);
int MPI_Type_hindexed(
        int count,
        int blocklens[],
        MPI_Aint indices[],
        MPI_Datatype old_type,
        MPI_Datatype *newtype );
int MPI_Type_hvector(
        int count,
        int blocklen,
        MPI_Aint stride,
        MPI_Datatype old_type,
        MPI_Datatype *newtype );
int MPI_Type_indexed(
        int count,
        int blocklens[],
        int indices[],
        MPI_Datatype old_type,
        MPI_Datatype *newtype );
int MPI_Type_lb ( MPI_Datatype datatype, MPI_Aint *displacement );
int MPI_Type_size ( MPI_Datatype datatype, int *size );
int MPI_Type_struct(
        int count,
        int blocklens[],
        MPI_Aint indices[],
        MPI_Datatype old_types[],
        MPI_Datatype *newtype );
int MPI_Type_ub ( MPI_Datatype datatype, MPI_Aint *displacement );
int MPI_Type_vector(
        int count,
        int blocklen,
        int stride,
        MPI_Datatype old_type,
        MPI_Datatype *newtype );
int MPI_Unpack ( void *inbuf, int insize, int *position,
                void *outbuf, int outcount, MPI_Datatype datatype,
                MPI_Comm comm );
int MPI_Waitsome(
        int incount,
        MPI_Request array_of_requests[],
        int *outcount,
        int array_of_indices[],
        MPI_Status array_of_statuses[] );
double MPI_Wtick();
double MPI_Wtime();
