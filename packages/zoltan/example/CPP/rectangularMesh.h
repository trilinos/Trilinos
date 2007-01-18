// Class rectangularMesh 
//
//    A regular rectangular mesh, vertices are initially distributed 
//      across processes in a round robin manner.
//
//      Global IDs of vertices begin at lower left with vertex 1,
//      increase to the right, then begin at the next higher row
//      on the left, and so on.  Edges connect vertices
//      horizontally and vertically.
//
//    Y
//    ^  21----22----23----24---25
//    |  |     |     |     |    |
//    |  16----17----18----19---20
//    |  |     |     |     |    |
//    |  11----12----13----14---15
//    |  |     |     |     |    |
//    |  6-----7-----8-----9----10
//    |  |     |     |     |    |
//    |  1-----2-----3-----4----5
//    o---------------------------------> X
//
//    parameters: X and Y dimensions, X and Y stride, 
//                2D coordinates of lower left corner
//
//    If vertex weights are desired, we set the weight to the
//    number of neighbors a vertex has in the mesh.  No edge weights.
//    Only query functions for Zoltan's geometric methods are included.
//    No update (migration) is implemented.  The query functions only
//    respond with the initial partitioning.
//

#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <vector.h>
#include <zoltan_cpp.h>

//
//  MPICPP - Define this if your C++ bindings for MPI work.
//  NAMESPACES_OK - Define this if your system uses namespaces.
//

//#define MPICPP
//#define NAMESPACES_OK

#ifdef NAMESPACES_OK
using namespace std;

#ifdef MPICPP
using namespace MPI;
#endif
#endif


class rectangularMesh{
  int x_dim;       /* number of vertices in horizontal direction */
  int y_dim;       /* number of vertices in vertical direction */
  int x_stride;    /* distance between vertices in hor. direction */
  int y_stride;    /* distance between vertices in vert. direction */
  int x_1;         /* X coordinate of left most vertices */
  int y_1;         /* Y coordinate of lowest vertices */
  int nprocs;
  int rank;       
  vector<int> myIDs;

  void reset_totals() { myIDs.erase(myIDs.begin(), myIDs.end());}

  void set_par_env()
    {
#ifdef MPICPP
    this->rank = MPI::COMM_WORLD.Get_rank();
    this->nprocs = MPI::COMM_WORLD.Get_size();
#else
    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &this->nprocs);
#endif
    }

protected:

  // Return number of vertices in global mesh
  int get_num_IDs()
    {
    return this->x_dim * this->y_dim;
    }

  // Return number of vertices initially in my partition
  int get_num_my_IDs()
    {
    // Initial distribution of vertices is done in a round robin manner

    if (myIDs.size() == 0)  
      {
      int numIDs = this->x_dim * this->y_dim;
      for (int i=0; i<numIDs; i++)
        {
        if (i % this->nprocs == this->rank)
          {
          // global ID of vertex in my partition
          myIDs.push_back(i+1);
          }
        }
      }

    return myIDs.size();
    }

  // Return the global ID of my i'th vertex
  int get_my_ID(int i)
    {
    if (myIDs.size() == 0) get_num_my_IDs(); 
    return this->myIDs[i];
    }

  // Return the global ID of the vertex in the given row and column
  int get_ID(int row, int col)
    {
    // rows start at 0 at the bottom, columns start at 0 on the left
    return (row*this->x_dim + col + 1);
    }

  // Return the row and column of the given vertex global ID
  void get_row_col(int ID, int &row, int &col)
    {
    row = (ID - 1) / this->x_dim;
    col = (ID - 1) - (row * this->x_dim);
    return;
    }

  // Return the number of edges in the mesh of a vertex
  int get_num_edges(int ID)
    {
    int row=0, col=0, nedges=4;
    this->get_row_col(ID, row, col);

    if ((row == 0) || (row == this->y_dim-1)) 
      {
      nedges--;
      } 
    if ((col == 0) || (col == this->x_dim-1))
      {
      nedges--;
      } 
    return nedges; 
    }

  // Return the geometric coordinates of the vertex in the
  // given row and column.
  void get_coords(int row, int col, int &x, int &y)
    {
    x = this->x_1 + col * this->x_stride;
    y = this->y_1 + row * this->y_stride;
    return;
    }

  // Return the geometric coordinates of the vertex with the
  // given global ID.
  void get_coords(int ID, int &x, int &y)
    {
    int row=0, col=0;
    this->get_row_col(ID, row, col);
    get_coords(row, col, x, y);
    return;
    }

public:
  // Constructors
  rectangularMesh::rectangularMesh()
    {
    this->x_dim   =5;
    this->y_dim   =5;
    this->x_stride=1;
    this->y_stride=1;
    this->x_1     =0;
    this->y_1     =0;
    
    this->set_par_env(); 
    this->reset_totals();
    }
  rectangularMesh(int xdim, int ydim, 
             int xstride, int ystride,
             int x1, int y1)
    {
    this->x_dim=xdim;
    this->y_dim=ydim;
    this->x_stride=xstride;
    this->y_stride=ystride;
    this->x_1=x1;
    this->y_1=y1;

    this->set_par_env(); 
    this->reset_totals();
    }

  // Destructor
  rectangularMesh::~rectangularMesh()
    { 
    reset_totals();
    }

  // can add copy operator, copy constructor, ostream, istream

  // Set/reset parameters
  void set_x_dim(int xdim){this->x_dim = xdim; this->reset_totals();}
  void set_y_dim(int ydim){this->y_dim = ydim; this->reset_totals();}
  void set_x_stride(int xstride){this->x_stride = xstride; }
  void set_y_stride(int ystride){this->y_stride = ystride; }
  void set_x_1(int x1){this->x_1 = x1; }
  void set_y_1(int y1){this->y_1 = y1; }

  // Query functions for Zoltan library

  static int get_number_of_objects(void *data, int *ierr)
    {
    // Prototype: ZOLTAN_NUM_OBJ_FN
    // Return the number of objects I own.

    rectangularMesh *mesh = (rectangularMesh *)data;
    *ierr = ZOLTAN_OK;
    return mesh->get_num_my_IDs();
    }

  static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
    {
    // Prototype: ZOLTAN_OBJ_LIST_FN
    // Return list of my objects, with optional object weights.

    rectangularMesh *mesh = (rectangularMesh *)data;
    *ierr = ZOLTAN_OK;

    for (int i=0; i<mesh->get_num_my_IDs(); i++)
      { 
      globalID[i] = mesh->get_my_ID(i);
      localID[i] = i;
      if (wgt_dim == 1)
        {
        obj_wgts[i] = mesh->get_num_edges(globalID[i]);
        }
      }
    return;
    }

  static int get_num_geometry(void *data, int *ierr)
    {
    // Prototype: ZOLTAN_NUM_GEOM_FN
    // Return the dimension of a vertex, for geometric methods

    *ierr = ZOLTAN_OK;
    return 2;
    }

  static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                    int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
    {
    // Prototype: ZOLTAN_GEOM_MULTI_FN
    // Return the coordinates of my objects (vertices), for geometric methods.

    int x, y;
    rectangularMesh *mesh = (rectangularMesh *)data;
    *ierr = ZOLTAN_OK;

    for (int i=0; i<num_obj; i++)
      {
      mesh->get_coords(globalID[i], x, y);
      geom_vec[2*i] = (double)x;
      geom_vec[2*i + 1] = (double)y;
      }
    return;
    }

  // Draw a picture of the partitioning.  All processes send their
  // global ID list to process 0, which draws the picture.

  void draw_partitions(char *title, int len, int *gids, int want_wgts)
    {
    int root = 0;
    int iam_root = (this->rank == root);

    if ((this->x_dim > 75) || (this->y_dim > 75))
      {
      if (iam_root) cout << "too big to draw" << endl;
      return;
      }

    int *incounts = NULL, *indisp = NULL, *ingids = NULL;
    if (iam_root) incounts = new int [this->nprocs];

#ifdef MPICPP
    MPI::COMM_WORLD.Gather(&len, 1, MPI::INT, incounts, 1, MPI_INT, root);
#else
    MPI_Gather(&len, 1, MPI_INT, incounts, 1, MPI_INT, root, MPI_COMM_WORLD);
#endif

    if (iam_root){
      indisp = new int[this->nprocs + 1];
      indisp[0] = 0;
      for (int i=1; i<=this->nprocs; i++){
        indisp[i] = indisp[i-1] + incounts[i-1];
      }
      ingids = new int [indisp[this->nprocs]];
    }

#ifdef MPICPP
    MPI::COMM_WORLD.Gatherv(gids, len, MPI::INT, 
                            ingids, incounts, indisp, MPI::INT, root);
#else
    MPI_Gatherv(gids, len, MPI_INT, 
                ingids, incounts, indisp, MPI_INT, root, MPI_COMM_WORLD);
#endif

  if (!iam_root){
    return;
  }

  int ngids = indisp[this->nprocs];
 
  int *parts = new int [ngids];
  float *wgt_total = NULL;
  int nextId = 0;

  if (want_wgts){
    wgt_total = new float [this->nprocs];
    memset(wgt_total, 0, sizeof(float) * this->nprocs);
  }

  for (int i=0; i < this->nprocs; i++) {
    for (int j=0; j<incounts[i]; j++){
      parts[ingids[nextId]-1] = i;
      if (want_wgts){
        wgt_total[i] += this->get_num_edges(ingids[nextId]);
      }
      nextId++;
    }
  }

  printf("\n\n%s\n\n",title);

  int row, col;

  for (row=this->y_dim-1; row >= 0; row--){
    for (col=0; col < this->x_dim; col++){
      int local_id = get_ID(row, col) - 1;
      printf("%2d ", parts[local_id]);
      if (this->x_stride > this->y_stride) printf(" ");
    }
    printf("\n");
    if (this->y_stride > this->x_stride) printf("\n");
  }

  if (want_wgts){
    for (int i=0; i<this->nprocs; i++){
      printf("Total weight partition %d: %f\n",i, wgt_total[i]);
    }
    delete [] wgt_total;
  }

  this->get_coords(1, col, row);
  printf("Geometry: (%d, %d) X ",col, row);

  this->get_coords(ngids, col, row);
  printf("(%d, %d)\n", col, row);


  delete [] incounts;
  delete [] indisp;
  delete [] ingids;
  delete [] parts;
  }
};
