
/*
 *  Routines to update the graph due to data migration.
 *  These are the updates due to EXPORTS.
 */

/*****************************************************************************/
static void build_send_lists(
        GRAPH *graph, 
        EXPORT_LIST *export_list, 
        LIST **proc_exp_list,
        int *proc_exp_cnt
)
{
int new_proc, edge_proc, edge_new_proc;
int edge, num_edges;
VERTEX *export_vertex, *edge_vertex;
LIST *list_ptr;

  /*
   *  Build lists of vertices to send to importing processors and
   *  concerned processors.
   */

  for (list_ptr = export_list; list_ptr != NULL; list_ptr = list_ptr->Next) {
    export_vertex = list_ptr->Data;
    new_proc = export_vertex->New_Proc;

    /*
     *  Send the vertex to the importer.
     */

    add_to_list(proc_exp_list[new_proc], export_vertex);
    proc_exp_cnt[new_proc]++;
    
    /*
     *  Loop over edges to alert concerned processors.
     */

    num_edges = get_edge_list(export_vertex);
    proc_exp_edge_cnt[new_proc] += num_edges;

    for (edge = 0; edge < num_edges; edge++) {

      edge_vertex = edge_list[edge];
      edge_proc = edge_vertex->Proc;
      edge_new_proc = edge_vertex->New_Proc;
      if (edge_proc != Proc && new_proc != edge_proc) {
        /*
         *  The edge vertex is in a concerned (but not importer) processor.
         */
 
        add_to_list(proc_exp_list[edge_proc], export_vertex);
        proc_exp_cnt[edge_proc]++;
        proc_exp_edge_cnt[edge_proc] += num_edges;
      }
      if (edge_proc != edge_new_proc && edge_new_proc != Proc && 
          new_proc != new_edge_proc) {
        /*
         *  New edge proc is concerned since a neighboring vertex
         *  is going to a new processor; new edge proc will need to
         *  set up ghost vertex.
         */

        add_to_list(proc_exp_list[new_edge_proc], export_vertex);
        proc_exp_cnt[new_edge_proc]++;
        proc_exp_edge_cnt[new_edge_proc] += num_edges;
      }
    }
  }
}

/*****************************************************************************/
static void update_graph_due_to_exports(GRAPH *graph, EXPORT_LIST *export_list)
{
int edge;
VERTEX *export_vertex, *edge_vertex;
LIST *list_ptr;
int cnt_export_edges, cnt_edge_edges;
int num_edges;

  /*
   *  Update the graph to handle changes due to exports.
   *  The export vertices may become ghost vertices, or they may
   *  be removed.  If the export vertices' edge vertices were
   *  ghost vertices, they may also be removed.
   */

  for (list_ptr = export_list; list_ptr != NULL; list_ptr = list_ptr->Next) {
    export_vertex = list_ptr->Data;

    num_edges = get_edge_list(export_vertex);
    cnt_export_edges = num_edges;;

    for (edge = 0; edge < num_edges; edge++) {
      edge_vertex = edge_list[edge];
      cnt_edge_edges = get_num_edges(edge_vertex);
      if (edge_vertex->Proc != Proc) {
        /*
         *  edge_vertex is a ghost cell; remove edge_vertex from edge list of
         *  export_vertex.
         */
        remove_edge(export_vertex, edge_vertex);
        cnt_export_edges--;

        /*
         *  Also remove export_vertex from the edge list of edge_vertex.
         */

        remove_edge(edge_vertex, export_vertex);
        cnt_edge_edges--;

        if (cnt_edge_edges == 0) {
          /*
           *  This ghost cell is no longer needed.  Remove it from the graph.
           */
          remove_vertex_from_graph(graph, edge_vertex);
        }
      }
      else if (edge_vertex->New_Proc != Proc) {
        /*
         *  This edge_vertex is also being exported; we do not need an
         *  edge from export_vertex to this edge_vertex.  Remove 
         *  edge_vertex from the edge list of export_vertex.
         */

        remove_edge(export_vertex, edge_vertex);
        cnt_export_edges--;
      }
  
    }
    if (cnt_export_edges == 0) {
      /*
       *  All of export_vertex's edge neighbors are also being
       *  exported; 
      remove_vertex_from_graph(graph, export_vertex);
    }
  }
}

/*****************************************************************************/
void send_vertices_to_processors(
        NBORHD *nborhd,
        LIST **proc_exp_list,
        int *proc_exp_cnt,
        int *proc_exp_edge_cnt
)
{
/*
 *  Send data in export lists.
 *  Note:  this routine will have to be changed to handle data typing
 *  for heterogeneous processors.
 *  Send message buffer format:
 *     number of vertices sent (int)
 *     list of vertex structures
 *     list of all edge lists for vertices sent.
 */

int i;
char *buf; 
int *buf_header;
VERTEX *buf_vertices, *vertex;
ID *buf_edges;
LIST *ptr;
int header_size = sizeof(int);   /* assume header is only num of vertices */
  
  max_size = 0;
  for (i = 0; i < nborhd->Num_Proc; i++) {
    index = nborhd->Nbor_Proc[i];
    size = proc_exp_cnt[index] * sizeof(VERTEX) +
           proc_exp_edge_cnt[index] * sizeof(ID);
    if (size > max_size) max_size = size;
  }
   
  buf = (char *) Array_Alloc(1, max_size + header_size, sizeof(char));

  for (i = 0; i < nborhd->Num_Proc; i++) {
    index = nborhd->Nbor_Proc[i];
   
    buf_header = (int *) buf;
    buf_vertices = (VERTEX *) (buf + header_size);
    buf_edges = (ID *) (buf + header_size + proc_exp_cnt[index]*sizeof(VERTEX));
    size = header_size + proc_exp_cnt[index] * sizeof(VERTEX)
         + proc_exp_edge_cnt[index] * sizeof(ID);

    buf_header[0] = proc_exp_cnt[index];

    for (ptr = proc_exp_list[index]->Data, cnt = 0, edge_cnt = 0; ptr != NULL; 
         ptr = ptr->Next, cnt++) {
      vertex = ptr->Data;
      
      copy_vertex(&(buf_vertices[cnt]), vertex);
      copy_edge_list(&(buf_edges[edge_cnt]), vertex);
       
      edge_cnt += get_num_edges
    }

    error = MPI_Send(buf, size, MPI_BYTE, index, MSG_EXPORT, MPI_COMM_WORLD);
  }

  LB_FREE(&buf);
}

/*****************************************************************************/
void update_send(GRAPH *graph, NBORHD *nborhd, EXPORT_LIST *export_list) 
{ 
LIST **nbor_list;
int *nbor_cnt;
VERTEX *export_vertex;
EXPORT_LIST *list_ptr;
LIST **proc_exp_list;
int *proc_exp_cnt, proc_exp_edge_cnt;

  proc_exp_list = (LIST **) Array_Alloc(1, Num_Proc, sizeof(LIST *));
  proc_exp_cnt = (int *) Array_Alloc(1, 2 * Num_Proc, sizeof(int));
  proc_exp_edge_cnt = proc_exp_cnt + Num_Proc;
  
  build_send_lists(graph, export_list, proc_exp_list, proc_exp_cnt, 
                   proc_exp_edge_cnt);
  
  /*
   *  Send vertices to the processors.
   */

  send_vertices_to_processors(nborhd, proc_exp_list, proc_exp_cnt, 
                              proc_exp_edge_cnt);

  update_graph_due_to_exports(graph, export_list);

  LB_FREE(&proc_exp_list);
  LB_FREE(&proc_exp_cnt);
}

