
#define MSG_EXPORT 333

void receive_import_message(
        char **buf, 
        int *num_vertices, 
        VERTEX **buf_vertices,
        ID **buf_edges
)
{
MPI_Status mpi_status;
int size;
int *buf_header;
int header_size = sizeof(int);  /* THIS MUST BE SAME AS IN SEND !!!! */

  MPI_Probe(MPI_ANY_SOURCE, MSG_EXPORT, MPI_COMM_WORLD, &mpi_status);

  MPI_Get_Count(&mpi_status, MPI_BYTE, &size);

  *buf = (char *) Array_Alloc(1, size, sizeof(char);

  error = MPI_Recv(*buf, size, MPI_BYTE, mpi_status.MPI_SOURCE, 
                   mpi_status.MPI_TAG, MPI_COMM_WORLD, &mpi_status);

  /*
   *  Set pointers for data layout.
   *  NOTE:  these pointers must match layout of SEND buffer.
   */

  buf_header = (int *) (*buf);
  *num_vertices = buf_header[0];
  *buf_vertices = (VERTEX *) (*buf + header_size);
  *buf_edges = (ID *) (*buf + header_size + *num_vertices * sizeof(VERTEX);

}

/****************************************************************************/
void update_receive(GRAPH *graph, NBORHD *nborhd)
{
int num_msgs_read;
int nborhd_size = nborhd->Num_Nbors;
char *buf;
VERTEX *buf_vertices;
ID *buf_edges;


  for (num_msgs_read = 0; num_msgs_read < nborhd_size-1; num_msgs_read++) {
    
    receive_import_message(&buf, &num_vertices, &buf_vertices, &buf_edges);

    for (i = 0, edge_cnt = 0; i < num_vertices; i++) {

      recv_vertex = buf_vertices[i];

      if (recv_vertex->New_Proc == Proc) {
        /*
         *  This vertex is a vertex to be imported!
         *  See whether the vertex is already a ghost; if so
         *  convert it to a regular vertex.  Otherwise, alloc
         *  space for a new vertex.
         */

        if ((vertex = find_vertex(graph, recv_vertex->ID)) == NULL) {
          vertex = create_new_vertex(recv_vertex->ID);
          add_vertex_to_graph(graph, vertex);
        }
        
        /*
         *  Update edge lists and other data for vertex.
         */

        update_vertex_data(vertex, recv_vertex, buf_edges[edge_cnt]);
 

        /*
         *  Add ghost cells where needed.
         */

        num_edges = get_edge_list(vertex, edge_list);
        for (edge = 0; edge < num_edges; edge++) {
          edge_vertex_ID = edge_list[edge];
          if ((edge_vertex = find_vertex(graph, edge_vertex_ID)) == NULL) {
            /*
             *  This vertex is not yet in the graph on this processor.
             *  Need to add it as a ghost node.  
             *  DON'T KNOW ITS LOCATION YET, but should get that info from
             *  wherever it is owned (as a concerned message).
             */
            edge_vertex = create_new_vertex(edge_vertex_ID);
            add_vertex_to_graph(graph, edge_vertex);
          }
        }
      }  /* end IMPORT */
      else {
        /*
         *  This vertex is a vertex that the processor is concerned about.
         */

        if ((vertex = find_vertex(graph, recv_vertex->ID)) == NULL) {
          /*
           *  This vertex is a new ghost vertex.
           */

          vertex = create_new_vertex(recv_vertex->ID);
          add_vertex_to_graph(graph, vertex);
        }
        update_vertex_data(vertex, recv_vertex, buf_edges[edge_cnt]);
      } /* end CONCERNED */
      edge_cnt += recv_vertex->num_edges;
    } /* end loop over vertices received */
    LB_Free((void **) &buf);
  } /* end loop over messages */
}

