#ifndef FNX_MESH_H
#define FNX_MESH_H


#include <iostream>
#include "hdf5.h"
#include "mpi.h"
#include <vector>
#define HAVE_CONFIG_H
#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_LocalMap.h"
#include "parmetis.h"
#include <map>

using namespace std;

void read_data_size(string what, hid_t& file_id, int& num_vectors, int& length)
{
  herr_t	status;

  hid_t dataset_id = H5Dopen(file_id, what.c_str());
  hid_t dataspace_id = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(dataspace_id);
  vector<hsize_t> dims(rank);
  status = H5Sget_simple_extent_dims(dataspace_id, &dims[0], NULL);

  assert (dims.size() == 2);

  num_vectors = dims[0]; length = dims[1];

  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
}

void read_data(string what, hid_t& file_id, hid_t type, Epetra_Map& Map, int size,
               char* data)
{
  herr_t	status;

  hid_t dataset_id = H5Dopen(file_id, what.c_str());
  hid_t dataspace_id = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(dataspace_id);
  vector<hsize_t> dims(rank);
  status = H5Sget_simple_extent_dims(dataspace_id, &dims[0], NULL);
  
  vector<hsize_t> offset(2);
  vector<hsize_t> count(2);

  offset[0] = 0; offset[1] = Map.MinMyGID();
  count[0] = size; count[1] = Map.NumMyElements();

  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, &offset[0],
                               NULL, &count[0], NULL);

  hid_t memspace_id = H5Screate_simple (2, &count[0], NULL);  

  status = H5Dread (dataset_id, type, memspace_id, dataspace_id,
                    H5P_DEFAULT, data);

  status = H5Sclose(memspace_id);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
}

int main(int argc, char* argv[])
{
  hid_t       file_id;         /* file and dataset identifiers */
  hid_t	plist_id;        /* property list identifier( access template) */
  /*
   * MPI variables
   */
  int mpi_size, mpi_rank;
  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;
  
  /*
   * * Initialize MPI
   * */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  
#if 0
  if (false) {
    int i = 0, j = 0;
    char buf[80];
    char go = ' ';
    char hostname[80];

    if (Comm.MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
    for (i = 0; i <Comm.NumProc() ; i++) {
      if (i == Comm.MyPID() ) {
        gethostname(hostname, sizeof(hostname));
        sprintf(buf, "Host: %s\tComm().MyPID(): %d\tPID: %d", 
                hostname, Comm.MyPID(), getpid());
        printf("%s\n",buf);
        fflush(stdout);
        sleep(1);
      }
    }
    if(Comm.MyPID() == 0) {
      printf("\n");
      printf("** Pausing because environment variable ML_BREAK_FOR_DEBUGGER has been set,\n");
      puts("** or file ML_debug_now has been created");
      printf("**\n");
      printf("** You may now attach debugger to the processes listed above.\n");
      printf( "**\n");
      printf( "** Enter a character to continue > "); fflush(stdout);
      scanf("%c",&go);
    }
  }
#endif

  /* 
   * Set up file access property list with parallel I/O access
   */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  
  // Create a new file collectively.
  file_id = H5Fopen("square.h5", H5F_ACC_RDWR, H5P_DEFAULT);
  
  int num_vectors, length;

  // =============================== //
  // start reading nodal coordinates //
  // =============================== //
  
  read_data_size("/mesh/vertices", file_id, num_vectors, length);

  Epetra_Map lin_v_map(length, 0, Comm);
  int lin_v_MyElements = lin_v_map.NumMyElements();
  std::vector<double> lin_v_data(num_vectors * lin_v_MyElements);

  read_data("/mesh/vertices", file_id, H5T_IEEE_F64BE, lin_v_map, num_vectors,
            (char*)&lin_v_data[0]);

  // wrap data as Epetra_MultiVector to redistribute them
  Epetra_MultiVector lin_v(View, lin_v_map, &lin_v_data[0], lin_v_MyElements, 2);

  // ===================== //
  // reading boundary data //
  // ===================== //
  
  read_data_size("/mesh/boundaries", file_id, num_vectors, length);

  Epetra_Map lin_b_map(length, 0, Comm);
  int lin_b_NumMyElements = lin_b_map.NumMyElements();
  std::vector<int> lin_b_data(num_vectors * lin_b_NumMyElements);

  read_data("/mesh/boundaries", file_id, H5T_STD_U32BE, lin_b_map, num_vectors,
            (char*)&lin_b_data[0]);

  // wrap data as Epetra_MultiVector to redistribute them
  Epetra_MultiVector lin_b(lin_b_map, num_vectors);

  for (int i = 0 ; i < num_vectors ; ++i)
    for (int j = 0; j < lin_b_NumMyElements; ++j)
      lin_b[i][j] = (double)lin_b_data[j + i * lin_b_NumMyElements];

  // ========================= //
  // reading mesh connectivity //
  // ========================= //
  
  read_data_size("/mesh/elements", file_id, num_vectors, length);

  Epetra_Map lin_e_map(length, 0, Comm);
  int lin_e_NumMyElements = lin_e_map.NumMyElements();
  std::vector<int> lin_e_data(num_vectors * lin_e_NumMyElements);

  read_data("/mesh/elements", file_id, H5T_STD_U32BE, lin_e_map, num_vectors,
            (char*)&lin_e_data[0]);

  // wrap data as Epetra_MultiVector to redistribute them
  Epetra_MultiVector lin_e(lin_e_map, num_vectors);

  for (int i = 0 ; i < num_vectors ; ++i)
    for (int j = 0; j < lin_e_NumMyElements; ++j)
      lin_e[i][j] = (double)lin_e_data[j + i * lin_e_NumMyElements];

  // ==================================================================== //
  // Call ParMETIS to equilibrate the workload and minimize communication //
  // ==================================================================== //

  // 1) build the graph of the mesh
  vector<int> elmdist(Comm.NumProc() + 1); elmdist[0] = 0;

  Comm.GatherAll(&lin_e_NumMyElements, &(elmdist[1]), 1);

  for (int i = 2; i < Comm.NumProc() + 1; ++i)
    elmdist[i] += elmdist[i - 1];

  vector<int> eptr(lin_e_NumMyElements + 1);
  eptr[0] = 0;
  for (int i = 0; i < lin_e_NumMyElements; ++i)
    eptr[i + 1] = eptr[i] + num_vectors;
  vector<int> eind(lin_e_NumMyElements * 3);
  int icount = 0;
  for (int i = 0; i < lin_e_NumMyElements; ++i)
    for (int j = 0; j < 3; ++j)
      eind[icount++] = (int)lin_e[j][i];

  int ncommonnodes = 2;
  MPI_Comm parmetis_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &parmetis_comm);

  int numflag = 0;
  int wgtflag = 0;
  int ncon = 0;
  int nparts = Comm.NumProc();
  vector<float> tpwgts(nparts);
  vector<float> ubvec(nparts);
  int options[4]; options[0] = 0;
  vector<int> part(lin_e_NumMyElements);

  for (int i = 0; i < nparts; ++i) {
    tpwgts[i] = 1.0 / nparts;
    ubvec[i] = 1.05;
  }
  int edgecut;

  ParMETIS_V3_PartMeshKway(&elmdist[0], &eptr[0], &eind[0], 0, &wgtflag, 
                           &numflag, &ncon, &ncommonnodes, &nparts, &tpwgts[0], 
                           &ubvec[0], options, &edgecut, &part[0], &parmetis_comm);

  cout << ncommonnodes << " " << edgecut << endl;

  // MPI_Comm_destroy(&parmetis_comm);
  
  Epetra_Map proc_map(Comm.NumProc(), 0, Comm);
  Epetra_FECrsGraph proc_graph(Copy, proc_map, 0);

  for (int i = 0; i < lin_e_NumMyElements; ++i) {
    int owner = part[i];
    assert (owner >= 0 && owner < Comm.NumProc());
    int index = lin_e_map.GID(i);
    assert (index >=0 && index < lin_e_map.NumGlobalElements());
    proc_graph.InsertGlobalIndices(1, &owner, 1, &index);
  }

  Epetra_LocalMap all(lin_e_map.NumGlobalElements(), 0, Comm);
  proc_graph.GlobalAssemble(all, proc_map, false);

  // now build the final element map
  int NumMyElements; int* MyGlobalElements;
  EPETRA_CHK_ERR(proc_graph.ExtractGlobalRowView (Comm.MyPID(), NumMyElements, MyGlobalElements));

  Epetra_Map FEMap(-1, NumMyElements, MyGlobalElements, 0, Comm);
  Epetra_MultiVector FEData(FEMap, num_vectors);

  Epetra_Import FEImporter(FEMap, lin_e_map);
  FEData.Import(lin_e, FEImporter, Insert);

  Epetra_Vector tmp(FEMap); tmp.PutScalar(Comm.MyPID());
  Epetra_Vector L_FE_tmp(lin_e_map);
  L_FE_tmp.Export(tmp, FEImporter, Insert);

  cout << L_FE_tmp;

  // ============================================================ //
  // at this point we need to gather all vertices and faces that  //
  // will be used by the local elements. First we create the list //
  // of them, then we use Import statements.                      //
  // ============================================================ //

  map<int, bool> hash;
  for (int i = 0; i < FEMap.NumMyElements() ; ++i)
    for (int j = 0 ; j < FEData.NumVectors() - 1; ++j)
      hash[(int)(FEData[j][i])] = true;

  vector<int> MyGlobalVertices(hash.size());

  map<int, bool>::iterator iter;
  int count = 0;
  for (iter = hash.begin(); iter != hash.end(); ++iter)
    MyGlobalVertices[++count] = iter->first;

  // repeated entries here
  Epetra_Map VTX_Map(-1, count, &MyGlobalVertices[0], 0, Comm);
  Epetra_Import VTX_Importer(VTX_Map, lin_v_map);
  Epetra_MultiVector VTX_Data(VTX_Map, lin_v.NumVectors());
  VTX_Data.Import(lin_v, VTX_Importer, Insert);

  // still to do boundaries...





  /*
   * Close property list.
   */
  H5Pclose(plist_id);
  
  /*
   * Close the file.
   */
  H5Fclose(file_id);

  MPI_Finalize();
  return(0);
}

