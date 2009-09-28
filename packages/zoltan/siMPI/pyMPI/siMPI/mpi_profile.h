/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: rrdrake $
 *    Date: 2009/07/14 22:54:56 $
 *    Revision: 1.3 $
 ****************************************************************************/
/* Profiling */
int PMPI_Init (int *argc, char **argv[]);
int PMPI_Finalize (void);

int PMPI_Comm_rank (MPI_Comm comm, int* rank);
int PMPI_Comm_size (MPI_Comm comm, int* size);

int PMPI_Send (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int PMPI_Msend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int PMPI_Recv (void* message, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status);

int PMPI_Irecv (void* message, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Isend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);

int PMPI_Bsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int PMPI_Rsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int PMPI_Ssend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

int PMPI_Comm_create(MPI_Comm comm, MPI_Group new_group, MPI_Comm* new_comm);
int PMPI_Comm_group(MPI_Comm comm, MPI_Group* group);
int PMPI_Comm_free(MPI_Comm* comm);

int PMPI_Abort(MPI_Comm comm, int errorcode);

int PMPI_Reduce ( void *sendbuf, void *recvbuf, int count,
   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm );

int PMPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root,
  MPI_Comm comm );

int PMPI_Op_create( MPI_User_function* function, int commute, MPI_Op* op);
int PMPI_Op_free( MPI_Op *op );

int PMPI_Sendrecv( void *sendbuf, int sendcount, MPI_Datatype sendtype,
    int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
    int source, int recvtag, MPI_Comm comm, MPI_Status *status );

int PMPI_Wait ( MPI_Request* request, MPI_Status* status);
int PMPI_Waitany( int count, MPI_Request array_of_requests[], int* index, MPI_Status* status);
int PMPI_Waitall(
        int count, 
        MPI_Request array_of_requests[], 
        MPI_Status array_of_statuses[] );


int PMPI_Scatter (void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
    int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
int PMPI_Gather ( void *sendbuf, int sendcnt, MPI_Datatype sendtype,
   void *recvbuf, int recvcnt, MPI_Datatype recvtype,
   int root, MPI_Comm comm );

int PMPI_Ibsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Irsend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Issend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);

/* ----------- */
/* STUBBED OUT */
int PMPIO_Test(MPIO_Request *request, int *flag, MPI_Status *status);
int PMPIO_Wait(MPIO_Request *request, MPI_Status *status);
int PMPI_Address( void *location, MPI_Aint *address);
int PMPI_Allgather ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                    void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm );
int PMPI_Allgatherv ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                     void *recvbuf, int *recvcounts, int *displs,
                    MPI_Datatype recvtype, MPI_Comm comm );
int PMPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
int PMPI_Alltoall( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                 MPI_Comm comm );
int PMPI_Alltoallv (
        void *sendbuf,
        int *sendcnts,
        int *sdispls,
        MPI_Datatype sendtype,
        void *recvbuf,
        int *recvcnts,
        int *rdispls,
        MPI_Datatype recvtype,
        MPI_Comm comm );
int PMPI_Attr_delete ( MPI_Comm comm, int keyval );
int PMPI_Attr_get (
        MPI_Comm comm,
        int keyval,
        void *attr_value,
        int *flag );
int PMPI_Attr_put ( MPI_Comm comm, int keyval, void *attr_value );
int PMPI_Barrier ( MPI_Comm comm);
int PMPI_Bsend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request );
int PMPI_Buffer_attach( void *buffer, int size );
int PMPI_Buffer_detach( void *bufferptr, int *size);
int PMPI_Cancel( MPI_Request *request );
int PMPI_Cart_coords ( MPI_Comm comm, int rank, int maxdims, int *coords );
int PMPI_Cart_create ( MPI_Comm comm_old, int ndims, int *dims, int *periods,
                     int reorder, MPI_Comm *comm_cart );
int PMPI_Cart_get (
        MPI_Comm comm,
        int maxdims,
        int *dims,
        int *periods,
        int *coords );
int PMPI_Cart_map (
        MPI_Comm comm_old,
        int ndims,
        int *dims,
        int *periods,
        int *newrank);
int PMPI_Cart_rank (
        MPI_Comm comm,
        int *coords,
        int *rank );
int PMPI_Cart_shift ( MPI_Comm comm, int direction, int displ,
                    int *source, int *dest );
int PMPI_Cart_sub ( MPI_Comm comm, int *remain_dims, MPI_Comm *comm_new );
int PMPI_Cartdim_get ( MPI_Comm comm, int *ndims );


int PMPI_Comm_compare (
        MPI_Comm  comm1,
        MPI_Comm  comm2,
        int *result);
int PMPI_Comm_dup (
        MPI_Comm comm,
        MPI_Comm *comm_out );
int PMPI_Comm_get_name( MPI_Comm comm, char *namep, int *reslen );
int PMPI_Comm_remote_group ( MPI_Comm comm, MPI_Group *group );
int PMPI_Comm_remote_size ( MPI_Comm comm, int *size );
int PMPI_Comm_set_name( MPI_Comm com, char *name );
int PMPI_Comm_split ( MPI_Comm comm, int color, int key, MPI_Comm *comm_out );
int PMPI_Comm_test_inter ( MPI_Comm comm, int *flag );
int PMPI_Dims_create(
        int nnodes,
        int ndims,
        int *dims);

int PMPI_Errhandler_create(
        MPI_Handler_function *function,
        MPI_Errhandler       *errhandler);
int PMPI_Errhandler_free( MPI_Errhandler *errhandler );
int PMPI_Errhandler_get( MPI_Comm comm, MPI_Errhandler *errhandler );
int PMPI_Errhandler_set( MPI_Comm comm, MPI_Errhandler errhandler );
int PMPI_Error_class(
        int errorcode,
        int *errorclass);
int PMPI_Error_string( int errorcode, char *string, int *resultlen );
int PMPI_File_close(MPI_File *fh);
MPI_Fint PMPI_File_c2f(MPI_File fh);
int PMPI_File_delete(char *filename, MPI_Info info);
MPI_File PMPI_File_f2c(MPI_Fint fh);
int PMPI_File_get_amode(MPI_File fh, int *amode);
int PMPI_File_get_atomicity(MPI_File fh, int *flag);
int PMPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp);
int PMPI_File_get_errhandler(MPI_File fh, MPI_Errhandler *errhandler);
int PMPI_File_get_group(MPI_File fh, MPI_Group *group);
int PMPI_File_get_info(MPI_File fh, MPI_Info *info_used);
int PMPI_File_get_position(MPI_File fh, MPI_Offset *offset);
int PMPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset);
int PMPI_File_get_size(MPI_File fh, MPI_Offset *size);
int PMPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype,
                             MPI_Aint *extent);
int PMPI_File_get_view(MPI_File fh, MPI_Offset *disp, MPI_Datatype *etype,
                MPI_Datatype *filetype, char *datarep);
int PMPI_File_iread(MPI_File fh, void *buf, int count,
                   MPI_Datatype datatype, MPIO_Request *request);
int PMPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype,
                      MPIO_Request *request);
int PMPI_File_iread_shared(MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPIO_Request *request);
int PMPI_File_iwrite(MPI_File fh, void *buf, int count,
                    MPI_Datatype datatype, MPIO_Request *request);
int PMPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, void *buf,
                       int count, MPI_Datatype datatype,
                       MPIO_Request *request);
int PMPI_File_iwrite_shared(MPI_File fh, void *buf, int count,
                       MPI_Datatype datatype, MPIO_Request *request);
int PMPI_File_open(MPI_Comm comm, char *filename, int amode,
                  MPI_Info info, MPI_File *fh);
int PMPI_File_preallocate(MPI_File fh, MPI_Offset size);
int PMPI_File_read(MPI_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_read_all(MPI_File fh, void *buf, int count,
                      MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_read_all_begin(MPI_File fh, void *buf, int count,
                            MPI_Datatype datatype);
int PMPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
                   int count, MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype,
                         MPI_Status *status);
int PMPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype);
int PMPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_read_ordered(MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_read_ordered_begin(MPI_File fh, void *buf, int count,
                             MPI_Datatype datatype);
int PMPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_read_shared(MPI_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int PMPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);

int PMPI_File_set_atomicity(MPI_File fh, int flag);
int PMPI_File_set_errhandler(MPI_File fh, MPI_Errhandler errhandler);
int PMPI_File_set_info(MPI_File fh, MPI_Info info);
int PMPI_File_set_size(MPI_File fh, MPI_Offset size);
int PMPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,MPI_Datatype filetype, char *datarep, MPI_Info info);

int PMPI_File_sync(MPI_File fh);
int PMPI_File_write(MPI_File fh, void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_write_all(MPI_File fh, void *buf, int count,
                       MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_write_all_begin(MPI_File fh, void *buf, int count,
                            MPI_Datatype datatype);
int PMPI_File_write_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_write_at(MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype,
                      MPI_Status *status);
int PMPI_File_write_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                          int count, MPI_Datatype datatype,
                          MPI_Status *status);
int PMPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype);
int PMPI_File_write_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_write_ordered(MPI_File fh, void *buf, int count,
                         MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_write_ordered_begin(MPI_File fh, void *buf, int count,
                              MPI_Datatype datatype);
int PMPI_File_write_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_write_shared(MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status);
int PMPI_Finalized( int *flag );
int PMPI_Gatherv ( void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                  void *recvbuf, int *recvcnts, int *displs,
                 MPI_Datatype recvtype,
                  int root, MPI_Comm comm );
int PMPI_Get_count(
        MPI_Status *status,
        MPI_Datatype datatype,
        int *count );
int PMPI_Get_elements ( MPI_Status *status, MPI_Datatype datatype,
                      int *elements );
int PMPI_Get_processor_name(
        char *name,
        int *resultlen);
int PMPI_Get_version(
        int *version,
        int *subversion );
int PMPI_Graph_create ( MPI_Comm comm_old, int nnodes, int *index, int *edges,
                      int reorder, MPI_Comm *comm_graph );
int PMPI_Graph_get ( MPI_Comm comm, int maxindex, int maxedges,
                   int *index, int *edges );
int PMPI_Graph_map ( MPI_Comm comm_old, int nnodes, int *index, int *edges,
                   int *newrank );
int PMPI_Graph_neighbors ( MPI_Comm comm, int rank, int maxneighbors,
                        int *neighbors );
int PMPI_Graph_neighbors_count ( MPI_Comm comm, int rank, int *nneighbors );
int PMPI_Graphdims_get ( MPI_Comm comm, int *nnodes, int *nedges );
int PMPI_Group_compare ( MPI_Group group1, MPI_Group group2, int *result );
int PMPI_Group_difference ( MPI_Group group1, MPI_Group group2,
                         MPI_Group *group_out );
int PMPI_Group_excl ( MPI_Group group, int n, int *ranks, MPI_Group *newgroup );
int PMPI_Group_free ( MPI_Group *group );
int PMPI_Group_incl ( MPI_Group group, int n, int *ranks, MPI_Group *group_out );
int PMPI_Group_intersection ( MPI_Group group1, MPI_Group group2,
                           MPI_Group *group_out );
int PMPI_Group_range_excl ( MPI_Group group, int n, int ranges[][3],
                         MPI_Group *newgroup );
int PMPI_Group_range_incl ( MPI_Group group, int n, int ranges[][3],
                         MPI_Group *newgroup );
int PMPI_Group_rank ( MPI_Group group, int *rank );
int PMPI_Group_size ( MPI_Group group, int *size );
int PMPI_Group_translate_ranks ( MPI_Group group_a, int n, int *ranks_a,
                             MPI_Group group_b, int *ranks_b );
int PMPI_Group_union ( MPI_Group group1, MPI_Group group2,
                     MPI_Group *group_out );
MPI_Fint PMPI_Info_c2f(MPI_Info info);
int PMPI_Info_create(MPI_Info *info);
int PMPI_Info_delete(MPI_Info info, char *key);
int PMPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
MPI_Info PMPI_Info_f2c(MPI_Fint info);
int PMPI_Info_free(MPI_Info *info);
int PMPI_Info_get(MPI_Info info, char *key, int valuelen, char *value, int *flag);
int PMPI_Info_get_nkeys(MPI_Info info, int *nkeys);
int PMPI_Info_get_nthkey(MPI_Info info, int n, char *key);
int PMPI_Info_get_valuelen(MPI_Info info, char *key, int *valuelen, int *flag);
int PMPI_Info_set(MPI_Info info, char *key, char *value);
int PMPI_Init_thread(int *argc, char ***argv, int required, int *provided );
int PMPI_Initialized( int *flag );
MPI_Handle_type PMPI_Int2handle( MPI_Fint f_handle, MPI_Handle_enum handle_kind );
int PMPI_Intercomm_create ( MPI_Comm local_comm, int local_leader,
                         MPI_Comm peer_comm, int remote_leader, int tag,
                         MPI_Comm *comm_out );
int PMPI_Intercomm_merge ( MPI_Comm comm, int high, MPI_Comm *comm_out );
int PMPI_Iprobe( int source, int tag, MPI_Comm comm, int *flag,
               MPI_Status *status );
int PMPI_Keyval_create (
        MPI_Copy_function *copy_fn,
        MPI_Delete_function *delete_fn,
        int *keyval,
        void *extra_state );
int PMPI_Keyval_free ( int *keyval );
int PMPI_Pack ( void *inbuf, int incount, MPI_Datatype datatype,
               void *outbuf, int outcount, int *position, MPI_Comm comm );
int PMPI_Pack_size ( int incount, MPI_Datatype datatype, MPI_Comm comm,
                   int *size );
int PMPI_Pcontrol( int level );
int PMPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status );
int PMPI_Recv_init( void *buf, int count, MPI_Datatype datatype, int source,
                  int tag, MPI_Comm comm, MPI_Request *request );
int PMPI_Reduce_scatter ( void *sendbuf, void *recvbuf, int *recvcnts,
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
MPI_Fint PMPI_Request_c2f( MPI_Request c_request );
int PMPI_Request_free( MPI_Request *request );
int PMPI_Rsend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request );
int PMPI_Scan ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm );
int PMPI_Scatterv (
        void *sendbuf,
        int *sendcnts,
        int *displs,
        MPI_Datatype sendtype,
        void *recvbuf,
        int recvcnt,
        MPI_Datatype recvtype,
        int root,
        MPI_Comm comm );
int PMPI_Send_init( void *buf, int count, MPI_Datatype datatype, int dest,
                  int tag, MPI_Comm comm, MPI_Request *request );
int PMPI_Sendrecv_replace( void *buf, int count, MPI_Datatype datatype,
                        int dest, int sendtag, int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status );
int PMPI_Ssend_init( void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Request *request );
int PMPI_Start(
        MPI_Request *request);
int PMPI_Startall( int count, MPI_Request array_of_requests[] );
int PMPI_Status_c2f( MPI_Status *c_status, MPI_Fint *f_status );
int PMPI_Status_set_cancelled( MPI_Status *status, int flag );
int PMPI_Status_set_elements( MPI_Status *status, MPI_Datatype datatype,
                           int count );
int PMPI_Test (
        MPI_Request  *request,
        int          *flag,
        MPI_Status   *status);
int PMPI_Test_cancelled(
        MPI_Status *status,
        int        *flag);
int PMPI_Testall(
        int count,
        MPI_Request array_of_requests[],
        int *flag,
        MPI_Status array_of_statuses[] );
int PMPI_Testany(
        int count,
        MPI_Request array_of_requests[],
        int *index, int *flag,
        MPI_Status *status );
int PMPI_Testsome(
        int incount,
        MPI_Request array_of_requests[],
        int *outcount,
        int array_of_indices[],
        MPI_Status array_of_statuses[] );
int PMPI_Topo_test ( MPI_Comm comm, int *top_type );
int PMPI_Type_commit ( MPI_Datatype *datatype );
int PMPI_Type_contiguous( 
        int count,
        MPI_Datatype old_type,
        MPI_Datatype *newtype);
int PMPI_Type_create_darray(int size, int rank, int ndims, 
                           int *array_of_gsizes, int *array_of_distribs, 
                           int *array_of_dargs, int *array_of_psizes, 
                           int order, MPI_Datatype oldtype, 
                           MPI_Datatype *newtype);
int PMPI_Type_create_indexed_block( 
        int count, 
        int blocklength, 
        int array_of_displacements[], 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype );
int PMPI_Type_create_subarray(
        int ndims, 
        int *array_of_sizes, 
        int *array_of_subsizes, 
        int *array_of_starts, 
        int order, 
        MPI_Datatype oldtype, 
        MPI_Datatype *newtype);
int PMPI_Type_extent( MPI_Datatype datatype, MPI_Aint *extent );
int PMPI_Type_free ( MPI_Datatype *datatype );
int PMPI_Type_get_contents(
        MPI_Datatype datatype, 
        int max_integers, 
        int max_addresses, 
        int max_datatypes, 
        int *array_of_integers, 
        MPI_Aint *array_of_addresses, 
        MPI_Datatype *array_of_datatypes);
int PMPI_Type_get_envelope(
        MPI_Datatype datatype, 
        int *num_integers, 
        int *num_addresses, 
        int *num_datatypes, 
        int *combiner);
int PMPI_Type_hindexed( 
        int count, 
        int blocklens[], 
        MPI_Aint indices[], 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype );
int PMPI_Type_hvector( 
        int count, 
        int blocklen, 
        MPI_Aint stride, 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype );
int PMPI_Type_indexed( 
        int count, 
        int blocklens[], 
        int indices[], 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype );
int PMPI_Type_lb ( MPI_Datatype datatype, MPI_Aint *displacement );
int PMPI_Type_size ( MPI_Datatype datatype, int *size );
int PMPI_Type_struct( 
        int count, 
        int blocklens[], 
        MPI_Aint indices[], 
        MPI_Datatype old_types[], 
        MPI_Datatype *newtype );
int PMPI_Type_ub ( MPI_Datatype datatype, MPI_Aint *displacement );
int PMPI_Type_vector( 
        int count, 
        int blocklen, 
        int stride, 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype );
int PMPI_Unpack ( void *inbuf, int insize, int *position, 
                void *outbuf, int outcount, MPI_Datatype datatype, 
                MPI_Comm comm );
int PMPI_Waitsome( 
        int incount, 
        MPI_Request array_of_requests[], 
        int *outcount, 
        int array_of_indices[], 
        MPI_Status array_of_statuses[] );
double PMPI_Wtick(void);
double PMPI_Wtime(void);
