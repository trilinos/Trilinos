#include "zoltan_cpp.h"

void inferPartition(
  Zoltan &zz,               // Zoltan class populated by a geometric
                            // partition with parameter KEEP_CUTS=1
  int nCoords,              // Number of coordinates whose part 
                            // assignment should be inferred
  double *x,                // x coordinates (length nCoords)
  double *y,                // y coordinates (length nCoords)
  double *z,                // z coordinates (length nCoords)
  int *part_assignments     // Output:  part assignments inferred for 
                            // each coordinate (x_i,y_i,z_i) 
                            // (length nCoords)
)
{
  int part, proc;
  double tmp[3] = {0.,0.,0.};

  // Loop over coordinates
  for (int i = 0; i < nCoords; i++) {

    // Gather coordinates into temporary array
    tmp[0] = x[i];
    if (y) tmp[1] = y[i];
    if (z) tmp[2] = z[i];

    // Compute the coordinate's part assignment
    zz.LB_Point_PP_Assign(tmp, proc, part);

    // Store the part assignment in the return array
    part_assignments[i] = part;
  }
}

////////////////////////////////////////////////////////////////////////////


class Mesh {
public:
  int ncoords;
  double *x, *y, *z;
  Mesh(int n_, double *x_, double *y_, double *z_) : 
    ncoords(n_), x(x_), y(y_), z(z_) {}
};

int nObj(void *data, int *ierr) {
  *ierr = ZOLTAN_OK; 
  return ((Mesh *) data)->ncoords;
}

void objMulti(void *data, int ngid, int nlid, 
              ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, int wdim, float *wgt,
              int *ierr) 
{
  Mesh *mesh = (Mesh *) data;
  for (int i = 0; i < mesh->ncoords; i++) {
    lid[i] = i;
    gid[i] = i; 
  }
}

int nGeom(void *data, int *ierr) { *ierr = ZOLTAN_OK; return 3; }
void geomMulti(void *data, int ngid, int nlid, int nobj, 
               ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, int ndim,
               double *coords, int *ierr)
{
  Mesh *mesh = (Mesh *) data;

  for (int i = 0; i < nobj; i++) {
    coords[i*ndim]   = mesh->x[lid[i]];
    coords[i*ndim+1] = mesh->y[lid[i]];
    coords[i*ndim+2] = mesh->z[lid[i]];
  }
  *ierr = ZOLTAN_OK;
}


int main(int narg, char **arg) {
  MPI_Init(&narg, &arg);

  double xOne[27] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};
  double yOne[27] = {0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2};
  double zOne[27] = {0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2};
  Mesh meshOne(27, xOne, yOne, zOne);

  double xTwo[8] = {0.1, 1.1, 2.1, 1.1, 0.1, 0.1, 2.1, 2.1};
  double yTwo[8] = {0.1, 0.1, 0.1, 1.1, 1.1, 0.1, 1.1, 2.1};
  double zTwo[8] = {0.1, 0.1, 0.1, 1.1, 1.1, 1.1, 1.1, 2.1};

  Zoltan zz(MPI_COMM_WORLD);

  zz.Set_Num_Obj_Fn(nObj, &meshOne);
  zz.Set_Obj_List_Fn(objMulti, &meshOne);
  zz.Set_Num_Geom_Fn(nGeom, &meshOne);
  zz.Set_Geom_Multi_Fn(geomMulti, &meshOne);

  zz.Set_Param("LB_METHOD", "RCB");
  zz.Set_Param("KEEP_CUTS", "1");
  zz.Set_Param("RETURN_LISTS", "PART");
  zz.Set_Param("NUM_GLOBAL_PARTS", "5");

  int nChanges;
  int nGid, nLid;
  int nImp, nExp;
  ZOLTAN_ID_PTR iGid, iLid;
  ZOLTAN_ID_PTR eGid, eLid;
  int *iPart, *iProc;
  int *partOne, *eProc;

  zz.LB_Partition(nChanges, nGid, nLid, nImp, iGid, iLid, iProc, iPart,
                                        nExp, eGid, eLid, eProc, partOne);

  std::cout << "MeshOne " << std::endl;
  for (int i = 0; i < 27; i++)
    std::cout << "  (" << xOne[i] << ", "
                       << yOne[i] << ", "
                       << zOne[i] << ") in " 
                       << partOne[eLid[i]] << std::endl;

  
  int partTwo[8];
  inferPartition(zz, 8, xTwo, yTwo, zTwo, partTwo);

  std::cout << "MeshTwo " << std::endl;
  for (int i = 0; i < 8; i++)
    std::cout << "  (" << xTwo[i] << ", "
                       << yTwo[i] << ", "
                       << zTwo[i] << ") in " 
                       << partTwo[i] << std::endl;
  
  zz.LB_Free_Part(&iGid, &iLid, &iProc, &iPart);
  zz.LB_Free_Part(&eGid, &eLid, &eProc, &partOne);
}
