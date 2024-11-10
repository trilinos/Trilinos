// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "makeparfiles.H"
#include <assert.h>       // for assert
#include <exodusII.h>     // for ex_put_node_set_param, ex_put_node_set, ex_...
#include <mpi.h>
#include <ne_nemesisI.h>  // for ne_put_cmap_params, ne_put_eb_info_global
#include <stdio.h>        // for sprintf, printf
#include <stdlib.h>       // for exit
#include <cstring>        // for strcpy
#include <iostream>       // for operator<<, basic_ostream, char_traits, cerr
#include <sstream>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk_tools 
{

void SetFileName(char fn[], int total_subdomains, int i, const char* rootdir);
int OpenExoFile(char filename[300] );
void getMyNodeNeighbors(int my_node_neighbor[26],
int xc, int yc, int zc,
int ncuts_x,  int  ncuts_y, int  ncuts_z,
int & num_node_cmaps);
void GenerateMapIds(const int my_node_neighbor[26], int ctype[],
int node_map_node_cnts[], const int num_nodes_neighbors[26],
int node_cmap_ids[],
int num_node_cmaps, int num_elem_cmaps,
int exoid, const int my_proc_id);

void CreateNodeIdList(int num_node_cmaps,
                      const int node_map_node_cnts[],
                      const int node_cmap_ids[],
                      const int ctype[],
                      int edge_node_ids[],
                      int face_node_ids[],
                      int cids[],
                      int numnodes_per_edge,
                      int num_nodes_in_set,
                      int exoid,
                      int my_proc_id);
int WriteCubeCoords(int xc, double lenx_subcube, int yc, double leny_subcube, int zc, double lenz_subcube,
                    int num_nodes, double *xt, double *x_coord, double *yt, double *y_coord,
                    double *zt, double *z_coord, int exoid);
int* AssignIdLists(int ctype,
                   int edge_node_ids[],
                   int face_node_ids[],
                   int cids[],
                   int numnodes_per_edge,
                   int num_nodes_in_set);
int CreateConnectivity(int conn[], int nelem_per_edge, int numnodes_per_edge);
void CreateNodeLists( int numnodes_per_edge, int *nodelist1, int *nodelist2, int num_nodes_in_set);
void CreateCoordTemplate(int numnodes_per_edge, double * x_coord, double dx,
                         double * y_coord, double dy,
                         double * z_coord, double dz);
void FindNodalIdsOfFaceNeighbors(int q, int face_node_ids[], int z, int num_nodes_in_set);
void WriteNodeElemMaps(int num_nodes, int nelem_per_edge, int num_elem,
                       int xc, int yc, int zc, int exoid,
                       int ncuts_x, int ncuts_y, int ncuts_z,
                       int num_nodes_globalx, int num_nodes_globaly, int num_nodes_globalz);
void WriteNodesets(int xc, int exoid, int num_nodes_in_set,
                   int *nodelist1, double *df_ns, int ncuts_x, int *nodelist2);


void getCubeCoord(int id, int &x_c, int &y_c, int &z_c, const int& n_x, const int& n_y, const int& n_z);

int getCubeId(int &x_c, int &y_c, int &z_c, const int& n_x, const int& n_y, const int& n_z);

// This function uses the following information to create a partitioned input file
// cube.par.${num_procs}.${file_id}
// 
// my_proc_id: This processor's id
// num_procs:  Total number of processors being used to create files
// ncuts_x:    Number of cuts to make in the x-direction of the cube
// ncuts_y:    Number of cuts to make in the y-direction of the cube
// ncuts_z:    Number of cuts to make in the z-direction of the cube
// nelem_per_edge: Number of elements per subdomains cube's edge
// lenx,leny,lenz: Length of geometry in each direction. Even though these
//                 lengths can vary, they are current coded to be 1.0, hence
//                 I use the term 'cube' to describe the geometry

void MakeParFile(const int& my_proc_id, const int& num_procs, const int& ncuts_x, const int& ncuts_y,
                 const int& ncuts_z, const int& nelem_per_edge, const double& lenx, const double& leny, const double& lenz,
                 const char* rootdir)
{

  int i = 0;

  // Notice, that total_subdomains and num_procs are NOT the same

  int total_subdomains = ncuts_x * ncuts_y * ncuts_z;

  int ierr = 0;
  char dirloc[300];

  const char* title = "Scaling Study Problem";
  const char* ff = "p";
  char ftype[10];
  strcpy(ftype, ff);

  // numnodes_per_edge:   Number of nodes along each subdomain cube's edge
  // num_dim:             Dimensionality of problem. Always equal to 3.
  // num_nodes:           Total number of nodes for each subdomain cube
  // num_elem:            Total number of elements for each subdomain cube
  // num_elem_blk:        Number of element blocks ( is 1 for this case )
  // num_node_sets:       Number of node sets ( always 2 for us )
  // num_side_sets:       Number of side sets ( always 0 for us )

  int numnodes_per_edge = nelem_per_edge + 1;
  int num_dim = 3;
  int num_nodes = numnodes_per_edge * numnodes_per_edge * numnodes_per_edge;
  int num_elem = nelem_per_edge * nelem_per_edge * nelem_per_edge;
  int num_elem_blk = 1;
  int num_node_sets = 2;
  int num_side_sets = 0;

  // The following are global values, not subdomain values
  // num_elems_global:    Number of elements of the global cube geometry
  // num_nodes_globalx:   Number of nodes along the x-edge of the global cube
  // num_nodes_globaly:   Number of nodes along the y-edge of the global cube
  // num_nodes_globalz:   Number of nodes along the z-edge of the global cube
  // num_nodes_global:    Total number of nodes of the global cube

  int num_elems_global = num_elem * total_subdomains;
  int num_nodes_globalx = nelem_per_edge * ncuts_x + 1;
  int num_nodes_globaly = nelem_per_edge * ncuts_y + 1;
  int num_nodes_globalz = nelem_per_edge * ncuts_z + 1;
  int num_nodes_global = num_nodes_globalx * num_nodes_globaly * num_nodes_globalz;

  // x, y, z:             *Almost* coordinates of the subdomain cube. These are more
  //                      like template values to be used to fill in xt,yt,zt
  // xt, yt, zt:          *Actual* coordinates of the subdomain cube. The values in
  //                      x, y, z are used, and an offset is added to obtain these.

  double *x_coord = new double[num_nodes];
  double *y_coord = new double[num_nodes];
  double *z_coord = new double[num_nodes];
  double *xt = new double[num_nodes];
  double *yt = new double[num_nodes];
  double *zt = new double[num_nodes];

  // Create coordinates for nodes
  // VERY IMPORTANT: Indicies and unique integer ids are used heavily in this program.
  //                 For example, assume a 2x2x2 cube decomposition for a total of 8
  //                 subdomains. Each subdomain has a unique triplet to describe it.
  //                 Subdomain 0 is (0,0,0), and subdomain 1 is (1,0,0). We increase
  //                 ids along the x-axis first, then the y-axis, and then the z-axis.
  //                 Therefore, subdomain 'X' will have a tripled (xc,yc,zc), where
  //                 'X' and (xc,yc,zc) are BOTH unique. Hence, each subdomain cube can
  //                 be described using either a unique id or triplet.
  //                 There are two functions at the end of this file:
  //                 getCubeId(...), and getCubeCoord(...). These functions will help
  //                 generate a unique id given a triplet or vice-versa:
  //                 In addition, this same labelling of each subdomain cube is also
  //                 used in many other case, e.g., the node labelling.

  // Assume nodes increase first in x-dir, then y-dir, then z-dir

  // Calculate lengths of each subdomain cube, used later to determine nodal coordinates

  double lenx_subcube = lenx / ncuts_x;
  double leny_subcube = leny / ncuts_y;
  double lenz_subcube = lenz / ncuts_z;

  double dx = lenx_subcube / nelem_per_edge;
  double dy = leny_subcube / nelem_per_edge;
  double dz = lenz_subcube / nelem_per_edge;

  // block_id:            the element block id is 1 for everything
  // elem_type:           always using hex8s
  // num_nodes_per_elem:  since always using hex8s, this is always 8
  // num_attr:            always 0

  int block_id = 1;
  const char* elem_type = "hex8";
  int num_nodes_per_elem = 8;
  int num_attr = 0;

  int *conn = new int[num_nodes_per_elem * num_elem];
  int counter = CreateConnectivity( conn,  nelem_per_edge,  numnodes_per_edge);
  assert ( counter == num_nodes_per_elem*num_elem );
  (void)counter;

  // Create nodesets and dfs
  // The following section creates the nodeset info. Again, since the local ids are used
  // for each subdomain cube, this can be done once before we get into the subdomain loop

  int num_nodes_in_set = numnodes_per_edge * numnodes_per_edge;
  double *df_ns = new double[num_nodes_in_set];

  for(i = 0; i < num_nodes_in_set; i++)
    df_ns[i] = 1.0;

  int *nodelist1 = new int[num_nodes_in_set];
  int *nodelist2 = new int[num_nodes_in_set];

  CreateNodeLists( numnodes_per_edge, nodelist1, nodelist2, num_nodes_in_set );

  // Create Coordinates as template
  CreateCoordTemplate(numnodes_per_edge, x_coord, dx, y_coord, dy, z_coord, dz);


  // Create NEMESIS INFO here. Again, all the *stuff* here is invariant to which subdomain
  // you are are, so it is calculated before we get into the subdomain loop. Almost all.
  // The following are some of the variables that need to be defined to create a *proper*
  // decomposed par file.
  // num_internal_nodes:      I don't fully understand this, but using an example, the
  //                          value for this is all the nodes on the subdomain
  // num_border_nodes:        equal to 0
  // num_external_nodes:      equal to 0
  // num_interal_elems:       Since this is an element based decomposition, no single element
  //                          is shared. Therefore, this value is set to the total number
  //                          of elements on this subdomain cube.
  // num_node_cmaps:          Number of maps that will have to created for each subdomain cube
  //                          to communicate information about the nodes that belong to more than
  //                          one processor. (Now you know why I am confused abotu num_internal_nodes)
  // num_elem_cmaps:          Similar to num_node_cmaps. This is set to zero.
  int num_internal_nodes = num_nodes;
  int num_border_nodes = 0;
  int num_external_nodes = 0;
  int num_internal_elems = num_elem;
  int num_border_elems = 0;
  int num_node_cmaps = 0;
  int num_elem_cmaps = 0;
  // Nemesis info. needs an element map and node map. Currently these are set to 1...num_elem, and
  // 1..num_nodes.
  int *dum_elem_map = new int[num_elem];
  int *dum_node_map = new int[num_nodes];
  for(i = 0;i < num_elem;i++)
    dum_elem_map[i] = i + 1;

  for(i = 0;i < num_nodes;i++)
    dum_node_map[i] = i + 1;

  // The following is something very interesting. Each subdomain cube can have a maximum of
  // 26 neighbors. These are 8 corner neighbors, 12 edge neighbors, and 6 face neighbors.
  // Using triplet to integer id mapping, the process of finding then nodal ids for each
  // possible neighbor are easy.
  // Generate node_ids for communication between neighbors
  // Since they use local ids, they should be generated just once,
  // and reused for all cubes
  //
  // Use numnodes_per_edge for cube
  // Find nodal ids of corner neighbors
  int corners[24] = {0, 0, 0, numnodes_per_edge - 1, 0, 0, 0, numnodes_per_edge - 1, 0, numnodes_per_edge - 1, numnodes_per_edge - 1, 0, 0, 0, numnodes_per_edge - 1, numnodes_per_edge - 1, 0, numnodes_per_edge - 1, 0, numnodes_per_edge - 1, numnodes_per_edge - 1, numnodes_per_edge - 1, numnodes_per_edge - 1, numnodes_per_edge - 1};
  int cids[8];
  int a, b, c;
  for(i = 0;i < 8;i++){
    a = corners[3 * i];
    b = corners[3 * i + 1];
    c = corners[3 * i + 2];
    cids[i] = 1 + getCubeId(a, b, c, numnodes_per_edge, numnodes_per_edge, numnodes_per_edge);
  }
  // Find nodal ids of edge neighbors
  int *edge_node_ids = new int[12 * numnodes_per_edge];
  int q = numnodes_per_edge; // to make things easier to type
  int qq = q - 1;
  int z = 0;
  for(i = 0;i < q;i++){
    // only x varies
    edge_node_ids[0 * q + i] = 1 + getCubeId(i, z, z, q, q, q);
    edge_node_ids[3 * q + i] = 1 + getCubeId(i, qq, z, q, q, q);
    edge_node_ids[8 * q + i] = 1 + getCubeId(i, z, qq, q, q, q);
    edge_node_ids[11 * q + i] = 1 + getCubeId(i, qq, qq, q, q, q);
    // only y varies
    edge_node_ids[1 * q + i] = 1 + getCubeId(z, i, z, q, q, q);
    edge_node_ids[2 * q + i] = 1 + getCubeId(qq, i, z, q, q, q);
    edge_node_ids[9 * q + i] = 1 + getCubeId(z, i, qq, q, q, q);
    edge_node_ids[10 * q + i] = 1 + getCubeId(qq, i, qq, q, q, q);
    // on z varies
    edge_node_ids[6 * q + i] = 1 + getCubeId(z, qq, i, q, q, q);
    edge_node_ids[7 * q + i] = 1 + getCubeId(qq, qq, i, q, q, q);
    edge_node_ids[4 * q + i] = 1 + getCubeId(z, z, i, q, q, q);
    edge_node_ids[5 * q + i] = 1 + getCubeId(qq, z, i, q, q, q);
  }
  // Find nodal ids of face neighbors
  int *face_node_ids = new int[6 * num_nodes_in_set]; // same as numnodes^2
  FindNodalIdsOfFaceNeighbors(q, face_node_ids, z, num_nodes_in_set);

  // This is where the main subdomain loop starts. Now each processor for this code
  // will generate a cube.par.X.Y file. The loop for doing this is below, starts here.
  for(i = my_proc_id; i < total_subdomains; i += num_procs)
  {


    SetFileName(dirloc,total_subdomains, i, rootdir);


    // Now have the filename figured out and stored in *dirloc*
    // Open exodusfile named *dirloc*

    int exoid = OpenExoFile(dirloc);

    // This is described earlier also, and maybe a little better before. But I am
    // leaving the following comments in (too lazy to revise).
    // Need to come up with a cube ordering scheme:
    // Starting at the origin: cube id's increase in the x-direction first,
    // then the y-direction, and finally the z-direction.
    // Nodeset 1 is on the x=0 face, and Nodeset 2 is on the x=1 face (ASSUMING lenx=1.0)
    // Each cube can have a maximum of 26 neighbors.
    // Therefore each cube can have a unique id, and also a unique (triplet) specifying what
    // the *coordinates* of the cube are. x_cube_c = { 0 to ncuts_x -1 }, same for y_cube_c,
    // and z_cube_c. Need methods that will provide an integer id given the triplet, vice-versa.
    // Suppose cube_c = { 0, 3, 2 }, then id = z_cube_c* (ncuts_x*ncuts_y) + y_cube_c*ncuts_x + x_cube_c
    // therfore id = 2*( ) + 3*( ) + 0
    // Suppose id = 32, what are the cube coordinates cube_c?
    // z_c = (id+1)/(ncuts_x*ncuts_y)
    // id = id - z_c*(ncuts_x*ncuts_y);
    // y_c = (id+1)/ncuts_x;
    // id = id - y_c*ncuts_x;
    // z_c = id;
    int xc = 0, yc = 0, zc = 0;
    // Get coordinates of cube withing global model as index (xc, yc, zc)
    getCubeCoord(i, xc, yc, zc, ncuts_x, ncuts_y, ncuts_z);
    // Start writing data to exodus file
    ierr = ex_put_init(exoid, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets);
    assert ( ierr>=0 );
    ierr = ne_put_init_global(exoid, num_nodes_global, num_elems_global, num_elem_blk, num_node_sets, num_side_sets);
    assert ( ierr>=0 );
    ierr = ne_put_init_info(exoid, total_subdomains, 1, ftype);
    assert ( ierr>=0 );
    int gids[2] = {1, 2};
    int num_global_node_counts[2];
    num_global_node_counts[0] = num_global_node_counts[1] = (nelem_per_edge * ncuts_y + 1) * (nelem_per_edge * ncuts_z + 1);
    ierr = ne_put_ns_param_global(exoid, gids, num_global_node_counts, num_global_node_counts);
    assert ( ierr>= 0 );
    // In nemesis API, but not in ACCESS area
    ierr = ne_put_eb_info_global(exoid, &block_id, &num_elems_global);
    assert ( ierr>=0 );
    // The following section of code is where the nodal communication maps are figured out.
    // This is a little messy, but it works.
    num_node_cmaps = 0;
    num_elem_cmaps = 0;

    // The num_nodes_neighbors array determines how many nodes need to be communicated
    // for each cube's neighbor, assuming the neighbor exists
    static int num_nodes_neighbors[26] = {1, numnodes_per_edge, 1, numnodes_per_edge, num_nodes_in_set, numnodes_per_edge, 1, numnodes_per_edge, 1, numnodes_per_edge, num_nodes_in_set, numnodes_per_edge, num_nodes_in_set, num_nodes_in_set, numnodes_per_edge, num_nodes_in_set, numnodes_per_edge, 1, numnodes_per_edge, 1, numnodes_per_edge, num_nodes_in_set, numnodes_per_edge, 1, numnodes_per_edge, 1};
    int my_node_neighbor[26];
    getMyNodeNeighbors(my_node_neighbor, xc, yc, zc,
                       ncuts_x, ncuts_y, ncuts_z, num_node_cmaps);

    // For this case, num_border_nodes=num_external_nodes=0 and
    // num_border_elems = 0 also. I don't fully understand this, and
    // I am just copying from an example i saw from a decomposed par file.
    ierr = ne_put_loadbal_param(exoid, num_internal_nodes, num_border_nodes, num_external_nodes, num_internal_elems, num_border_elems, num_node_cmaps, num_elem_cmaps, my_proc_id);
    assert ( ierr>= 0 );
    // Write out dummy maps
    ierr = ne_put_elem_map(exoid, dum_elem_map, 0, my_proc_id);
    assert ( ierr>= 0 );
    ierr = ne_put_node_map(exoid, dum_node_map, 0, 0, my_proc_id);
    assert ( ierr>= 0 );
    if ( ierr < 0 ) {
      printf("failed to put node map.\n");
    }
    // Generate node_map_ids and elem_map_ids, and their corresponding counts
    int *node_cmap_ids = new int[num_node_cmaps];
    int *node_map_node_cnts = new int[num_node_cmaps];
    int *ctype = new int[num_node_cmaps];
    GenerateMapIds( my_node_neighbor, ctype, node_map_node_cnts,
                    num_nodes_neighbors, node_cmap_ids, num_node_cmaps,
                    num_elem_cmaps, exoid, my_proc_id);

    // For each node communication map, create a list of node ids, and the processors that
    // own them. In this case, the map_id is the same as the proc_id.
    //
    CreateNodeIdList(num_node_cmaps,node_map_node_cnts,
                     node_cmap_ids, ctype,
                     edge_node_ids, face_node_ids, cids,
                     numnodes_per_edge, num_nodes_in_set,
                     exoid, my_proc_id);

    delete [] ctype;
    ctype = 0;
    delete [] node_map_node_cnts;
    node_map_node_cnts = 0;
    delete [] node_cmap_ids;
    node_cmap_ids = 0;
    // Write out coordinates for cube 'i'
    // Using coordinates of cube(0,0,0)'s add offset for this cube's coordinates
    ierr = WriteCubeCoords(xc, lenx_subcube, yc, leny_subcube, zc, lenz_subcube,
                           num_nodes, xt, x_coord, yt, y_coord, zt, z_coord, exoid);
    // Write node map
    WriteNodeElemMaps(num_nodes, nelem_per_edge,  num_elem,
                      xc,  yc,  zc,  exoid,
                      ncuts_x,  ncuts_y,  ncuts_z,
                      num_nodes_globalx,  num_nodes_globaly,  num_nodes_globalz);
    // Write element block attributes. Only 1 block here.

    ierr = ex_put_block(exoid, EX_ELEM_BLOCK, block_id, elem_type, num_elem, num_nodes_per_elem, 0, 0, num_attr);

    assert ( ierr>=0 );

    // Write element block connectivity

    ierr = ex_put_conn(exoid, EX_ELEM_BLOCK, block_id, conn, nullptr, nullptr);

    assert ( ierr>=0 );

    // Write nodeset info.
    // Only if the cube is on the left or right edge (based on x-axis)
    // will it have any nodesets
    WriteNodesets(xc, exoid, num_nodes_in_set, nodelist1, df_ns, ncuts_x, nodelist2);

    // Write proc_id as element variable

    ierr = ex_put_variable_param(exoid, EX_ELEM_BLOCK, 1);
    assert ( ierr >= 0 );

    ierr = ex_put_variable_name(exoid, EX_ELEM_BLOCK, 1, "Domain #");
    assert ( ierr >= 0 );

    double time_value = 1.0;
    ierr = ex_put_time(exoid, 1, &time_value);
    assert ( ierr >= 0 );

    int total = xc+yc+zc;
    double colorVal = total%2;

    if ( zc%2 == 0 && colorVal == 1)
      colorVal = 2;

    double* values = new double[num_elem];
    for(int ii = 0; ii < num_elem; ii++)
    {
      values[ii] = colorVal;
    }

    ierr = ex_put_var(exoid, 1, EX_ELEM_BLOCK, 1, 1, num_elem, values);
    assert ( ierr >= 0 );

    delete[] values; values = 0;

    // Close exodus file

    ierr = ex_close(exoid);
    assert ( ierr >= 0 );
    if ( ierr<0 ){
      printf("failed to close file.\n");
    }

  }

  delete[] face_node_ids; face_node_ids = 0;
  delete[] edge_node_ids; edge_node_ids = 0;
  delete[] dum_elem_map;  dum_elem_map = 0;
  delete[] dum_node_map;  dum_node_map = 0;
  delete[] nodelist2;     nodelist2 = 0;
  delete[] nodelist1;     nodelist1 = 0;
  delete[] df_ns;         df_ns = 0;
  delete[] conn;          conn = 0;
  delete[] zt;            zt = 0;
  delete[] yt;            yt = 0;
  delete[] xt;            xt = 0;
  delete[] z_coord;       z_coord = 0;
  delete[] y_coord;       y_coord = 0;
  delete[] x_coord;       x_coord = 0;

  return; // Exit function
}

void getCubeCoord(int id, int &x_c, int &y_c, int &z_c, const int& n_x, const int& n_y, const int& n_z)
{
  z_c = id / (n_x * n_y);
  id -= z_c * n_x * n_y;
  y_c = id / n_x;
  x_c = id - y_c * n_x;

  return;
}

int getCubeId(int &x_c, int &y_c, int &z_c, const int& n_x, const int& n_y, const int& n_z)
{
  int id = -1;
  // id = z_c * n_z * n_y + y_c * n_x + x_c;
  id = z_c * n_x * n_y + y_c * n_x + x_c;
  return id;
}


int GetNumDigits(int n)
{
  if (n==0) {return 1;}
  if (n<0) {n = -n;}
  int d = 0;
  while (n != 0) {
    n = n / 10;
    ++d;
  }
  return d;
}


void SetFileName(char dirloc[], int total_subdomains, int i, 
                 const char* rootdir)
// The following is some ugly code to arrive at the full path and name for
// each subdomain cube. If root = /ufs/tmp and nsubdomains is 8, then
// the 'dirloc', which contains the full path and name, will be
// /ufs/tmp/cube.exo.8.0, for cube 0

// First create the filename, which depends on the
// number of processors and current file id.
{
  const char* base_fn = "cube.exo.";
  char fn[50];
  int digit = GetNumDigits(total_subdomains);
  switch ( digit )
  {
  case 1:
    sprintf(fn, "%s%d.%.1d", base_fn, total_subdomains, i);
    break;
  case 2:
    sprintf(fn, "%s%d.%.2d", base_fn, total_subdomains, i);
    break;
  case 3:
    sprintf(fn, "%s%d.%.3d", base_fn, total_subdomains, i);
    break;
  case 4:
    sprintf(fn, "%s%d.%.4d", base_fn, total_subdomains, i);
    break;
  case 5:
    sprintf(fn, "%s%d.%.5d", base_fn, total_subdomains, i);
    break;
  case 6:
    sprintf(fn, "%s%d.%.6d", base_fn, total_subdomains, i);
    break;
  default:
    throw std::runtime_error("no more than 6 digits are supported");
  }


  std::string directory(rootdir);
  if (directory.empty()) {
    sprintf(dirloc, "%s", fn);
  } else {
    sprintf(dirloc, "%s/%s", rootdir, fn);
  }
}



int OpenExoFile(char filename[300] )
{
  // Now have the filename figured out and stored in *dirloc*
  // Open exodusfile named *filename*
  int comp_ws = sizeof (double);
  int io_ws = sizeof (double);
  int exoid = ex_create( filename, EX_CLOBBER, &comp_ws, &io_ws );
  // Check for errors in creating the Exodus file
  if(exoid < 0)
  {
    std::ostringstream os;
    os<<"\nCould not create file " << filename << "\n.";
    throw std::runtime_error(os.str());
  }
  return exoid;
}


void getMyNodeNeighbors(int my_node_neighbor[26],
int xc, int yc, int zc,
int ncuts_x, int ncuts_y, int  ncuts_z,
int & num_node_cmaps)
{
  // all_neighbors: there are a maximum of 26 possible neighbors for each cube.
  // and I am using a grid system where x varies from 0 to ncuts_x-1,
  // and y varies from 0 to ncuts_y-1, and same for z. Therefore, I can
  // use the offsets below 3*26 to determine if a neighbor exists.
  static int all_neighbors[78] = {-1, -1, -1, // Corner 1
                                  0, -1, -1, // Edge 1
                                  1, -1, -1, // Corner 2
                                  -1, 0, -1, // Edge 2
                                  0, 0, -1, // Face 1 z-
                                  1, 0, -1, // Edge 3
                                  -1, 1, -1, // Corner 3
                                  0, 1, -1, // Edge 4
                                  1, 1, -1, // Corner 4
                                  //=========================
                                  -1, -1, 0, // Edge 5
                                  0, -1, 0, // Face 2 y-
                                  1, -1, 0, // Edge 6
                                  -1, 0, 0, // Face 3 x-
                                  1, 0, 0, // Face 4 x+
                                  -1, 1, 0, // Edge 7
                                  0, 1, 0, // Face 5 y+
                                  1, 1, 0, // Edge 8
                                  //=========================
                                  -1, -1, 1, // Corner 5
                                  0, -1, 1, // Edge 9
                                  1, -1, 1, // Corner 6
                                  -1, 0, 1, // Edge 10
                                  0, 0, 1, // Face 6 z+
                                  1, 0, 1, // Edge 11
                                  -1, 1, 1, // Corner 7
                                  0, 1, 1, // Edge 12
                                  1, 1, 1};

  for(int j = 0;j < 26;j++){
    my_node_neighbor[j] = -1;
  }
  for(int j = 0;j < 26;j++){
    int *cc = &all_neighbors[3 * j];
    int junk1 = xc + cc[0];
    int junk2 = yc + cc[1];
    int junk3 = zc + cc[2];
    if(junk1 < 0 || junk1 >= ncuts_x)
      continue;

    if(junk2 < 0 || junk2 >= ncuts_y)
      continue;

    if(junk3 < 0 || junk3 >= ncuts_z)
      continue;

    // At this point, the cube does have a neighbor 'j'. So add it to the list.
    num_node_cmaps++;
    my_node_neighbor[j] = getCubeId(junk1, junk2, junk3, ncuts_x, ncuts_y, ncuts_z);
  }
}



void GenerateMapIds(const int my_node_neighbor[26], int ctype[],
int node_map_node_cnts[], const int num_nodes_neighbors[26],
int node_cmap_ids[],
int num_node_cmaps, int num_elem_cmaps,
int exoid, const int my_proc_id)
{
  int counter = 0;
  for(int j = 0;j < 26;j++){
    if(my_node_neighbor[j] != -1){
      ctype[counter] = j;
      node_map_node_cnts[counter] = num_nodes_neighbors[j];
      node_cmap_ids[counter++] = my_node_neighbor[j];
    }
  }

  assert ( counter == num_node_cmaps );
  assert ( 0 == num_elem_cmaps );
  ne_put_cmap_params(exoid, node_cmap_ids, node_map_node_cnts, 0, 0, my_proc_id);
  // assert ( ierr>= 0 );
}



void CreateNodeIdList(int num_node_cmaps,
                      const int node_map_node_cnts[],
                      const int node_cmap_ids[],
                      const int ctype[],
                      int edge_node_ids[],
                      int face_node_ids[],
                      int cids[],
                      int numnodes_per_edge,
                      int num_nodes_in_set,
                      int exoid,
                      int my_proc_id)
{

  for(int j = 0; j < num_node_cmaps; j++)
  {
    int* proc_ids = new int[node_map_node_cnts[j]];

    // proc_ids is the same at the id of this communication map
    for(int k = 0; k < node_map_node_cnts[j]; k++)
    {
      proc_ids[k] = node_cmap_ids[j];
    }

    // Now we just need the node ids. We should be able to determine
    // them using edge info.
    //
    // cids, edge_node_ids, face_node_ids;
    int* node_ids =  AssignIdLists(ctype[j],
                                   edge_node_ids,
                                   face_node_ids,
                                   cids,
                                   numnodes_per_edge,
                                   num_nodes_in_set);
    int ierr = ne_put_node_cmap(exoid, node_cmap_ids[j], node_ids, proc_ids, my_proc_id);

    assert( ierr >= 0);
    (void)ierr;

    delete[] proc_ids;
    proc_ids = 0;

  }
}


int WriteCubeCoords(int xc, double lenx_subcube, int yc, double leny_subcube, int zc, double lenz_subcube,
                    int num_nodes, double *xt, double *x_coord, double *yt, double *y_coord,
                    double *zt, double *z_coord, int exoid)
{
  // Write out coordinates for cube 'i'
  // Using coordinates of cube(0,0,0)'s add offset for this cube's coordinates
  double dxx = xc * lenx_subcube;
  double dyy = yc * leny_subcube;
  double dzz = zc * lenz_subcube;
  for(int j = 0;j < num_nodes;j++){
    xt[j] = x_coord[j] + dxx;
    yt[j] = y_coord[j] + dyy;
    zt[j] = z_coord[j] + dzz;
  }
  int ierr = ex_put_coord(exoid, xt, yt, zt);
  assert ( ierr >= 0 );
  return ierr;
}


int* AssignIdLists(int ctype,
                   int edge_node_ids[],
                   int face_node_ids[],
                   int corner_node_ids[],
                   int numnodes_per_edge,
                   int num_nodes_in_set)
{
  enum geom { corner=0, edge=1, face=2};
  static int mapCtyp[][2] = {
    { corner, 1},
    { edge, 1},
    { corner, 2},
    { edge, 2},
    { face, 1},
    { edge, 3},
    { corner, 3},
    { edge, 4},
    { corner, 4},
    { edge, 5},
    { face, 2},
    { edge, 6},
    { face, 3},
    { face, 4},
    { edge, 7},
    { face, 5},
    { edge, 8},
    { corner, 5},
    { edge, 9},
    { corner, 6},
    { edge, 10},
    { face, 6},
    { edge, 11},
    { corner, 7},
    { edge, 12},
    { corner, 8}
  };

  assert(ctype>=0);
  assert(ctype<26); // 8 corners + 12 edges + 6 faces

  int* node_ids = 0;

  enum geom mygeom = (enum geom)mapCtyp[ctype][0];
  int myindex = mapCtyp[ctype][1]-1;

  switch ( mygeom )
  {
  case corner:
    node_ids = &corner_node_ids[myindex];
    break;
  case edge:
    node_ids = &edge_node_ids[myindex * numnodes_per_edge];
    break;
  case face:
    node_ids = &face_node_ids[myindex * num_nodes_in_set];
    break;
  }

  return node_ids;
}



int CreateConnectivity(int conn[], int nelem_per_edge, int numnodes_per_edge)
{
  // Create Element Block Connectivity
  // The following section creates the connectivity of the subdomain cube. Interestingly enough,
  // since the local ids are used, the connectivity is the same for each subdomain cube. This is
  // why it is calculated before we get into the subdomain loop.

  // Again, vary elements from 1 to .. , x-dir first, y-dir, then z-dir
  int counter = 0;

  for(int k = 0; k < nelem_per_edge; k++)
  {
    for(int j = 0; j < nelem_per_edge; j++)
    {
      for(int i = 0; i < nelem_per_edge; i++)
      {
        // Given a triplet, calculate the unique id. Hey, the id happens to be the
        // nodal id needed to define the connectivity. Again, the ids here are local ids.

        int nid = getCubeId(i, j, k, numnodes_per_edge, numnodes_per_edge, numnodes_per_edge);
        nid += 1;
        int nid1 = nid + numnodes_per_edge;
        conn[counter + 0] = nid;
        conn[counter + 1] = nid + 1;
        conn[counter + 2] = nid1 + 1;
        conn[counter + 3] = nid1;
        conn[counter + 5] = nid + 1 + (numnodes_per_edge * numnodes_per_edge);
        conn[counter + 4] = conn[counter + 5] - 1;
        conn[counter + 6] = nid1 + 1 + (numnodes_per_edge * numnodes_per_edge);
        conn[counter + 7] = conn[counter + 6] - 1;
        counter += 8;
      }
    }
  }
  return counter;
}



void CreateNodeLists( int numnodes_per_edge, int *nodelist1, int *nodelist2, int num_nodes_in_set)
{
  int counter = 0;
  int lface = 0;
  int rface = numnodes_per_edge - 1;
  for(int k = 0;k < numnodes_per_edge;k++)
    for(int j = 0;j < numnodes_per_edge;j++){
      // Once again, given a triplet, find the unique id. Use to determine nodal ids for each
      // nodeset. These are in local id numbering.
      nodelist1[counter] = getCubeId(lface, j, k, numnodes_per_edge, numnodes_per_edge, numnodes_per_edge) + 1;
      nodelist2[counter] = getCubeId(rface, j, k, numnodes_per_edge, numnodes_per_edge, numnodes_per_edge) + 1;
      counter++;
    }

  assert ( counter == num_nodes_in_set );
}



void CreateCoordTemplate(int numnodes_per_edge, double * x_coord, double dx,
                         double * y_coord, double dy,
                         double * z_coord, double dz)
{
  // Create Coordinates as template
  for(int k = 0;k < numnodes_per_edge;k++)
  {
    for(int j = 0;j < numnodes_per_edge;j++)
    {
      for(int i = 0;i < numnodes_per_edge;i++){
        // Knowing the index of the local id, get coordinates of cube (0,0,0)
        int id = getCubeId(i, j, k, numnodes_per_edge, numnodes_per_edge, numnodes_per_edge);
        x_coord[id] = dx * i;
        y_coord[id] = dy * j;
        z_coord[id] = dz * k;
      }
    }
  }
}




void FindNodalIdsOfFaceNeighbors(int q, int face_node_ids[], int z, int num_nodes_in_set)
{
  int counter = 0;
  int qq = q - 1;
  //face 1
  for(int i = 0;i < q;i++)
    for(int j = 0;j < q;j++)
      face_node_ids[counter++] = 1 + getCubeId(j, i, z, q, q, q);


  assert ( counter == num_nodes_in_set );
  //face 2
  for(int i = 0;i < q;i++)
    for(int j = 0;j < q;j++)
      face_node_ids[counter++] = 1 + getCubeId(j, z, i, q, q, q);


  assert ( counter == 2*num_nodes_in_set );
  //face 3
  for(int i = 0;i < q;i++)
    for(int j = 0;j < q;j++)
      face_node_ids[counter++] = 1 + getCubeId(z, j, i, q, q, q);


  assert ( counter == 3*num_nodes_in_set );
  //face 4
  for(int i = 0;i < q;i++)
    for(int j = 0;j < q;j++)
      face_node_ids[counter++] = 1 + getCubeId(qq, j, i, q, q, q);


  assert ( counter == 4*num_nodes_in_set );
  //face 5
  for(int i = 0;i < q;i++)
    for(int j = 0;j < q;j++)
      face_node_ids[counter++] = 1 + getCubeId(j, qq, i, q, q, q);


  assert ( counter == 5*num_nodes_in_set );
  //face 6
  for(int i = 0;i < q;i++)
    for(int j = 0;j < q;j++)
      face_node_ids[counter++] = 1 + getCubeId(j, i, qq, q, q, q);


  assert ( counter == 6*num_nodes_in_set );
}


void WriteNodeElemMaps(int num_nodes, int nelem_per_edge, int num_elem,
                       int xc, int yc, int zc, int exoid,
                       int ncuts_x, int ncuts_y, int ncuts_z,
                       int num_nodes_globalx, int num_nodes_globaly, int num_nodes_globalz)
{
  int *nodemap = new int[num_nodes];
  // We have (xc,yc,zc) coordiante index of cube.
  // Translate to nodeindex (n_xc, n_yc, n_zc )
  int n_xc = nelem_per_edge * xc;
  int n_yc = nelem_per_edge * yc;
  int n_zc = nelem_per_edge * zc;
  int i1, j1, k1;
  int counter = 0;
  for(k1 = n_zc; k1 <= (n_zc + nelem_per_edge); k1++)
    for(j1 = n_yc; j1 <= (n_yc + nelem_per_edge); j1++)
      for(i1 = n_xc; i1 <= (n_xc + nelem_per_edge); i1++)
      {
        int nid = getCubeId(i1, j1, k1, num_nodes_globalx, num_nodes_globaly, num_nodes_globalz);
        nodemap[counter++] = nid + 1;
      }

  assert( counter == num_nodes);
  int ierr = ex_put_id_map(exoid, EX_NODE_MAP, nodemap);
  assert(ierr>=0);
  (void)(ierr);
  delete[] nodemap;
  nodemap = 0;
  int *elemmap = new int[num_elem];

  n_xc = nelem_per_edge * xc;
  n_yc = nelem_per_edge * yc;
  n_zc = nelem_per_edge * zc;

  counter = 0;

  for(k1 = n_zc; k1 < (n_zc + nelem_per_edge); k1++)
    for(j1 = n_yc; j1 < (n_yc + nelem_per_edge); j1++)
      for(i1 = n_xc; i1 < (n_xc + nelem_per_edge); i1++)
      {
        int nid = getCubeId(i1, j1, k1, nelem_per_edge * ncuts_x, nelem_per_edge * ncuts_y, nelem_per_edge * ncuts_z);
        elemmap[counter++] = nid + 1;
      }

  assert( counter == num_elem);

  ierr = ex_put_id_map(exoid, EX_ELEM_MAP, elemmap);

  assert( ierr>=0);

  delete[] elemmap;
  elemmap = 0;
}




void WriteNodesets(int xc, int exoid, int num_nodes_in_set,
                   int *nodelist1, double *df_ns, int ncuts_x, int *nodelist2)
{
  // Write nodeset info.
  // Only if the cube is on the left or right edge (based on x-axis)
  // will it have any nodesets
  int ierr = 0;
  if(xc == 0) // left_face
  {
    ierr = ex_put_set_param(exoid, EX_NODE_SET, 1, num_nodes_in_set, num_nodes_in_set);
    assert( ierr >= 0);
    ierr = ex_put_set(exoid, EX_NODE_SET, 1, nodelist1, nullptr);
    assert( ierr>= 0);
    ierr = ex_put_set_dist_fact(exoid, EX_NODE_SET, 1, df_ns);
    assert( ierr>= 0);
  }
  else // if no nodeset on this subdomain
  {
    ierr = ex_put_set_param(exoid, EX_NODE_SET, 1, 0, 0);
    assert( ierr >= 0);
  }
  if(xc == ncuts_x - 1) // right_face
  {
    ierr = ex_put_set_param(exoid, EX_NODE_SET, 2, num_nodes_in_set, num_nodes_in_set);
    assert( ierr >= 0);
    ierr = ex_put_set(exoid, EX_NODE_SET, 2, nodelist2, nullptr);
    assert( ierr >= 0);
    ierr = ex_put_set_dist_fact(exoid, EX_NODE_SET, 2, df_ns);
    assert( ierr >= 0);
  }
  else // if no nodeset on this subdomain
  {
    ierr = ex_put_set_param(exoid, EX_NODE_SET, 2, 0, 0);
    if (ierr < 0) {
      assert( ierr >= 0);
    }
  }
}


} // end namespace 
