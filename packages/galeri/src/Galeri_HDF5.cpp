#include "Galeri_ConfigDefs.h"
#ifdef HAVE_GALERI_HDF5
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Galeri_Exception.h"
#include "Galeri_HDF5.h"

// ==========================================================================
// data container and iterators to find a dataset with a given name
struct FindDataset_t
{
  string name;
  bool found;
};

static herr_t FindDataset(hid_t loc_id, const char *name, void *opdata)
{
  string& token = ((FindDataset_t*)opdata)->name;
  if (token == name)
    ((FindDataset_t*)opdata)->found = true;

  return(0);
}

// ==========================================================================
// This function copied from Roman Geus' FEMAXX code
static void WriteParameterListRecursive(const Teuchos::ParameterList& params, 
                                        hid_t group_id) 
{
  Teuchos::ParameterList::ConstIterator it = params.begin();
  for (; it != params.end(); ++ it) 
  {
    string key(params.name(it));
    if (params.isSublist(key)) 
    {
      // Sublist

      // Create subgroup for sublist.
      hid_t child_group_id = H5Gcreate(group_id, key.c_str(), 0);
      WriteParameterListRecursive(params.sublist(key), child_group_id);
      H5Gclose(child_group_id);
    } 
    else 
    {
      //
      // Regular parameter
      //

      // Create dataspace/dataset.
      herr_t status;
      hsize_t one = 1;
      hid_t dataspace_id, dataset_id;
      bool found = false; // to avoid a compiler error on MAC OS X GCC 4.0.0

      // Write the dataset.
      if (params.isType<string>(key)) 
      {
        string value = params.get<string>(key);
        hsize_t len = value.size() + 1;
        dataspace_id = H5Screate_simple(1, &len, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_C_S1, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, value.c_str());
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
        found = true;
      } 

      if (params.isType<bool>(key)) 
      {
        // Use H5T_NATIVE_USHORT to store a bool value
        unsigned short value = params.get<bool>(key) ? 1 : 0;
        dataspace_id = H5Screate_simple(1, &one, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, &value);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
        found = true;
      } 
      
      if (params.isType<int>(key)) 
      {
        int value = params.get<int>(key);
        dataspace_id = H5Screate_simple(1, &one, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, &value);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
        found = true;
      } 

      if (params.isType<double>(key)) 
      {
        double value = params.get<double>(key);
        dataspace_id = H5Screate_simple(1, &one, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, &value);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
        found = true;
      } 

      if (!found)
      {
        throw(Galeri::Exception(__FILE__, __LINE__, 
                                "type for parameter " + key + " not supported"));
      }
    }
  }
}

// ==========================================================================
// Recursive Operator function called by H5Giterate for each entity in group.
// This function copied from Roman Geus' FEMAXX code
static herr_t f_operator(hid_t loc_id, const char *name, void *opdata) 
{
  H5G_stat_t statbuf;
  hid_t dataset_id, space_id, type_id;
  Teuchos::ParameterList* sublist;
  Teuchos::ParameterList* params = (Teuchos::ParameterList*)opdata;
  /*
   * Get type of the object and display its name and type.
   * The name of the object is passed to this function by 
   * the Library. Some magic :-)
   */
  H5Gget_objinfo(loc_id, name, 0, &statbuf);
  switch (statbuf.type) {
  case H5G_GROUP:
    sublist = &params->sublist(name);
    H5Giterate(loc_id, name , NULL, f_operator, sublist);
    break;
  case H5G_DATASET:
    hsize_t len;
    dataset_id = H5Dopen(loc_id, name);
    space_id = H5Dget_space(dataset_id);
    if (H5Sget_simple_extent_ndims(space_id) != 1)
      throw(Galeri::Exception(__FILE__, __LINE__, 
                              "dimensionality of parameters must be 1."));
    H5Sget_simple_extent_dims(space_id, &len, NULL);
    type_id = H5Dget_type(dataset_id);
    if (H5Tequal(type_id, H5T_NATIVE_DOUBLE) > 0) {
      double value;
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      params->set(name, value);
    } else if (H5Tequal(type_id, H5T_NATIVE_INT) > 0) {
      int value;
      H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      params->set(name, value);
    } else if (H5Tequal(type_id, H5T_C_S1) > 0) {
      char* buf = new char[len];
      H5Dread(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
      params->set(name, std::string(buf));
      delete[] buf;
    } else if (H5Tequal(type_id, H5T_NATIVE_USHORT) > 0) {
      unsigned short value;
      H5Dread(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      params->set(name, value != 0 ? true : false);
    } else {
      throw(Galeri::Exception(__FILE__, __LINE__,
                              "unsupported datatype")); // FIXME
    }
    H5Tclose(type_id);
    H5Sclose(space_id);  
    H5Dclose(dataset_id);  
    break;
  default:
    throw(Galeri::Exception(__FILE__, __LINE__,
                            "unsupported datatype")); // FIXME
  }
  return 0;
}

// ==========================================================================
void Galeri::HDF5::Create(const string FileName)
{
  if (IsOpen())
    throw(Exception(__FILE__, __LINE__,
                    "an HDF5 is already open, first close the current one",
                    "using method Close(), then open/create a new one"));

  FileName_ = FileName;

#ifdef HAVE_MPI
  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);
  // Create property list for collective dataset write.

#if 0
  unsigned int boh = H5Z_FILTER_MAX;
  cout << "............." << H5Pset_filter(plist_id_, H5Z_FILTER_DEFLATE, 
                                           H5Z_FILTER_MAX, 0, &boh);
#endif


  H5Pset_fapl_mpio(plist_id_, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  cerr << "Not yet implemented" << endl;
  exit(EXIT_FAILURE);
#endif
  // create the file collectively and release property list identifier.
  file_id_ = H5Fcreate(FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, 
                      plist_id_);
  H5Pclose(plist_id_);

  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);

  IsOpen_ = true;
}

// ==========================================================================
void Galeri::HDF5::Open(const string FileName, int AccessType)
{
  if (IsOpen())
    throw(Exception(__FILE__, __LINE__,
                    "an HDF5 is already open, first close the current one",
                    "using method Close(), then open/create a new one"));

  FileName_ = FileName;

  // create the file collectively and release property list identifier.
  file_id_ = H5Fopen(FileName.c_str(), AccessType, H5P_DEFAULT);

  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);

  IsOpen_ = true;
}

// ==========================================================================
bool Galeri::HDF5::IsContained(string Name)
{
  FindDataset_t data;
  data.name = Name;
  data.found = false;

  int idx_f = H5Giterate(file_id_, "/", NULL, FindDataset, (void*)&data);

  return(data.found);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& GroupName, const string& DataSetName,
                         const int type, const int Length, 
                         void* data)
{
  hsize_t dimsf = Length;

  if (!IsContained(GroupName))
    CreateGroup(GroupName.c_str());

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());

  //dset_id = H5Dopen(file_id_, Name.c_str());
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), type, 
                            filespace_id, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, type, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
// NEW ONE
void Galeri::HDF5::Write(const string& GroupName, const string& DataSetName,
                         int what)
{
  hid_t filespace_id = H5Screate(H5S_SCALAR);
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_INT, 
                      filespace_id, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                           H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& GroupName, const string& DataSetName,
                         double what)
{
  hid_t filespace_id = H5Screate(H5S_SCALAR);
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_DOUBLE, 
                            filespace_id, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, 
                           filespace_id, H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, int what)
{
  bool flag = IsContained(Name);

  hid_t dset_id;
  hid_t filespace_id = H5Screate(H5S_SCALAR);
  if (flag)
  {
    dset_id = H5Dopen(file_id_, Name.c_str());
  }
  else
    dset_id = H5Dcreate(file_id_, Name.c_str(), H5T_NATIVE_INT, filespace_id,
                        H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                           H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, double data)
{
  bool flag = IsContained(Name);

  hid_t dset_id;
  hid_t filespace_id = H5Screate(H5S_SCALAR);
  if (flag)
  {
    dset_id = H5Dopen(file_id_, Name.c_str());
  }
  else
    dset_id = H5Dcreate(file_id_, Name.c_str(), H5T_NATIVE_DOUBLE, 
                        filespace_id, H5P_DEFAULT);

  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, const Epetra_BlockMap& BlockMap)
{
  int NumMyPoints       = BlockMap.NumMyPoints();
  int NumMyElements     = BlockMap.NumMyElements();
  int NumGlobalElements = BlockMap.NumGlobalElements();
  int NumGlobalPoints   = BlockMap.NumGlobalPoints();
  int* MyGlobalElements = BlockMap.MyGlobalElements();
  int* ElementSizeList  = BlockMap.ElementSizeList();

  vector<int> NumMyElements_v(Comm_.NumProc());
  Comm_.GatherAll(&NumMyElements, &NumMyElements_v[0], 1);

  vector<int> NumMyPoints_v(Comm_.NumProc());
  Comm_.GatherAll(&NumMyPoints, &NumMyPoints_v[0], 1);

  Write(Name, "MyGlobalElements", NumMyElements, NumGlobalElements, 
        H5T_NATIVE_INT, MyGlobalElements);
  Write(Name, "ElementSizeList", NumMyElements, NumGlobalElements, 
        H5T_NATIVE_INT, ElementSizeList);
  Write(Name, "NumMyPoints", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyPoints_v[0]);

  Write(Name, "NumGlobalPoints", 1, Comm_.NumProc(), H5T_NATIVE_INT, &NumGlobalPoints);
  Write(Name, "NumGlobalElements", 1, Comm_.NumProc(), H5T_NATIVE_INT, &NumGlobalElements);
  Write(Name, "NumProc", Comm_.NumProc());
  Write(Name, "IndexBase", BlockMap.IndexBase());
  Write(Name, "Epetra_BlockMap");
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, Epetra_BlockMap*& BlockMap)
{
  int NumGlobalElements, NumGlobalPoints, IndexBase, NumProc;

  ReadBlockMapProperties(GroupName, NumGlobalElements, NumGlobalPoints,
                        IndexBase, NumProc);

  vector<int> NumMyPoints_v(Comm_.NumProc());
  vector<int> NumMyElements_v(Comm_.NumProc());

  Read(GroupName, "NumMyElements", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyElements_v[0]);
  Read(GroupName, "NumMyPoints", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyPoints_v[0]);
  int NumMyElements = NumMyElements_v[Comm_.MyPID()];
  int NumMyPoints   = NumMyPoints_v[Comm_.MyPID()];

  if (NumProc != Comm_.NumProc())
    throw(Exception(__FILE__, __LINE__,
                    "requested map not compatible with current number of",
                    "processors, " + toString(Comm_.NumProc()) + 
                    " vs. " + toString(NumProc)));

  vector<int> MyGlobalElements(NumMyElements);
  vector<int> ElementSizeList(NumMyElements);

  Read(GroupName, "MyGlobalElements", NumMyElements, NumGlobalElements, 
       H5T_NATIVE_INT, &MyGlobalElements[0]);

  Read(GroupName, "ElementSizeList", NumMyElements, NumGlobalElements, 
       H5T_NATIVE_INT, &ElementSizeList[0]);

  BlockMap = new Epetra_BlockMap(NumGlobalElements, NumMyElements, 
                                 &MyGlobalElements[0], &ElementSizeList[0],
                                 IndexBase, Comm_);
}

// ==========================================================================
void Galeri::HDF5::ReadBlockMapProperties(const string& GroupName, 
                                          int& NumGlobalElements,
                                          int& NumGlobalPoints,
                                          int& IndexBase,
                                          int& NumProc)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_BlockMap")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_BlockMap"));

  Read(GroupName, "NumGlobalElements", NumGlobalElements);
  Read(GroupName, "NumGlobalPoints", NumGlobalPoints);
  Read(GroupName, "IndexBase", IndexBase);
  Read(GroupName, "NumProc", NumProc);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, const Epetra_Map& Map)
{
  int MySize = Map.NumMyElements();
  int GlobalSize = Map.NumGlobalElements();
  int* MyGlobalElements = Map.MyGlobalElements();

  vector<int> NumMyElements(Comm_.NumProc());
  Comm_.GatherAll(&MySize, &NumMyElements[0], 1);

  Write(Name, "MyGlobalElements", MySize, GlobalSize, 
        H5T_NATIVE_INT, MyGlobalElements);
  Write(Name, "NumMyElements", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyElements[0]);
  Write(Name, "NumGlobalElements", 1, Comm_.NumProc(), H5T_NATIVE_INT, &GlobalSize);
  Write(Name, "NumProc", Comm_.NumProc());
  Write(Name, "IndexBase", Map.IndexBase());
  Write(Name, "Epetra_Map");
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, Epetra_Map*& Map)
{
  int NumGlobalElements, IndexBase, NumProc;

  ReadMapProperties(GroupName, NumGlobalElements, IndexBase, NumProc);

  vector<int> NumMyElements_v(Comm_.NumProc());

  Read(GroupName, "NumMyElements", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyElements_v[0]);
  int NumMyElements = NumMyElements_v[Comm_.MyPID()];

  if (NumProc != Comm_.NumProc())
    throw(Exception(__FILE__, __LINE__,
                    "requested map not compatible with current number of",
                    "processors, " + toString(Comm_.NumProc()) + 
                    " vs. " + toString(NumProc)));

  vector<int> MyGlobalElements(NumMyElements);

  Read(GroupName, "MyGlobalElements", NumMyElements, NumGlobalElements, 
       H5T_NATIVE_INT, &MyGlobalElements[0]);

  Map = new Epetra_Map(-1, NumMyElements, &MyGlobalElements[0], IndexBase, Comm_);
}

// ==========================================================================
void Galeri::HDF5::ReadMapProperties(const string& GroupName, 
                                     int& NumGlobalElements,
                                     int& IndexBase,
                                     int& NumProc)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_Map")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_Map"));

  Read(GroupName, "NumGlobalElements", NumGlobalElements);
  Read(GroupName, "IndexBase", IndexBase);
  Read(GroupName, "NumProc", NumProc);
}

// ================ //
// Epetra_IntVector //
// ================ //

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, const Epetra_IntVector& x)
{
  const Epetra_BlockMap& OriginalMap = x.Map();
  Epetra_Map LinearMap(OriginalMap.NumGlobalElements(), OriginalMap.IndexBase(), Comm_);
  Epetra_Import Importer(LinearMap, OriginalMap);

  // FIXME: To do only if the map is not already linear
  Epetra_IntVector LinearVector(LinearMap);
  LinearVector.Import(x, Importer, Insert);

  Write(Name, "GlobalLength", x.GlobalLength());
  Write(Name, "Values", LinearMap.NumMyElements(), LinearMap.NumGlobalElements(),
        H5T_NATIVE_INT, x.Values());
  Write(Name, "Epetra_IntVector");
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, Epetra_IntVector*& X)
{
  int GlobalLength;
   
  Epetra_Map LinearMap(GlobalLength, 0, Comm_);
  X = new Epetra_IntVector(LinearMap);

  Read(GroupName, "Values", LinearMap.NumMyElements(), 
       LinearMap.NumGlobalElements(), H5T_NATIVE_INT, X->Values());
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const Epetra_Map& Map,
                        Epetra_IntVector*& X)
{
  int GlobalLength;
  ReadIntVectorProperties(GroupName, GlobalLength);

  Epetra_Map LinearMap(GlobalLength, 0, Comm_);
  Epetra_IntVector LinearX(LinearMap);

  Read(GroupName, "Values", LinearMap.NumMyElements(), 
       LinearMap.NumGlobalElements(),
       H5T_NATIVE_INT, LinearX.Values());

  Epetra_Import Importer(Map, LinearMap);
  X = new Epetra_IntVector(Map);
  X->Import(LinearX, Importer, Insert);
}

// ==========================================================================
void Galeri::HDF5::ReadIntVectorProperties(const string& GroupName, 
                                           int& GlobalLength)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_IntVector")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_IntVector"));

  Read(GroupName, "GlobalLength", GlobalLength);
}

// =============== //
// Epetra_CrsGraph //
// =============== //

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, const Epetra_CrsGraph& Graph)
{
  // like RowMatrix, only without values
  int MySize = Graph.NumMyNonzeros();
  int GlobalSize = Graph.NumGlobalNonzeros();
  vector<int> ROW(MySize);
  vector<int> COL(MySize);

  int count = 0;
  int* RowIndices;
  int NumEntries;

  for (int i = 0; i < Graph.NumMyRows(); ++i)
  {
    Graph.ExtractMyRowView(i, NumEntries, RowIndices);
    for (int j = 0; j < NumEntries; ++j)
    {
      ROW[count] = Graph.GRID(i);
      COL[count] = Graph.GCID(RowIndices[j]);
      ++count;
    }
  }

  Write(Name, "ROW", MySize, GlobalSize, H5T_NATIVE_INT, &ROW[0]);
  Write(Name, "COL", MySize, GlobalSize, H5T_NATIVE_INT, &COL[0]);
  Write(Name, "MaxNumIndices", Graph.MaxNumIndices());
  Write(Name, "NumGlobalRows", Graph.NumGlobalRows());
  Write(Name, "NumGlobalCols", Graph.NumGlobalCols());
  Write(Name, "NumGlobalNonzeros", Graph.NumGlobalNonzeros());
  Write(Name, "NumGlobalDiagonals", Graph.NumGlobalDiagonals());
  Write(Name, "Epetra_CrsGraph");
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, Epetra_CrsGraph*& Graph)
{
  int NumGlobalRows, NumGlobalCols;
  int NumGlobalNonzeros, NumGlobalDiagonals, MaxNumIndices;

  ReadCrsGraphProperties(GroupName, NumGlobalRows,
                         NumGlobalCols, NumGlobalNonzeros,
                         NumGlobalDiagonals, MaxNumIndices);

  Epetra_Map RangeMap(NumGlobalRows, 0, Comm_);
  Epetra_Map DomainMap(NumGlobalCols, 0, Comm_);

  Read(GroupName, DomainMap, RangeMap, Graph);
}

// ==========================================================================
void Galeri::HDF5::ReadCrsGraphProperties(const string& GroupName, 
                                          int& NumGlobalRows,
                                          int& NumGlobalCols,
                                          int& NumGlobalNonzeros,
                                          int& NumGlobalDiagonals,
                                          int& MaxNumIndices)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_CrsGraph")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_CrsGraph"));

  Read(GroupName, "NumGlobalRows", NumGlobalRows);
  Read(GroupName, "NumGlobalCols", NumGlobalCols);
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const Epetra_Map& DomainMap, 
                        const Epetra_Map& RangeMap, Epetra_CrsGraph*& Graph)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found in database"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_CrsGraph")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_CrsGraph"));

  int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
  Read(GroupName, "NumGlobalNonzeros", NumGlobalNonzeros);
  Read(GroupName, "NumGlobalRows", NumGlobalRows);
  Read(GroupName, "NumGlobalCols", NumGlobalCols);

  // linear distribution for nonzeros
  int NumMyNonzeros = NumGlobalNonzeros / Comm_.NumProc();
  if (Comm_.MyPID() == 0)
    NumMyNonzeros += NumGlobalNonzeros % Comm_.NumProc();

  vector<int> ROW(NumMyNonzeros);
  Read(GroupName, "ROW", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &ROW[0]);

  vector<int> COL(NumMyNonzeros);
  Read(GroupName, "COL", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &COL[0]);

  Epetra_FECrsGraph* Graph2 = new Epetra_FECrsGraph(Copy, DomainMap, 0);

  for (int i = 0; i < NumMyNonzeros; )
  {
    int count = 1;
    while (ROW[i + count] == ROW[i])
      ++count;

    Graph2->InsertGlobalIndices(1, &ROW[i], count, &COL[i]);
    i += count;
  }

  Graph2->FillComplete(DomainMap, RangeMap);

  Graph = Graph2;
}

// ================ //
// Epetra_RowMatrix //
// ================ //

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, const Epetra_RowMatrix& Matrix)
{
  int MySize = Matrix.NumMyNonzeros();
  int GlobalSize = Matrix.NumGlobalNonzeros();
  vector<int> ROW(MySize);
  vector<int> COL(MySize);
  vector<double> VAL(MySize);

  int count = 0;
  int Length = Matrix.MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);
  int NumEntries;

  for (int i = 0; i < Matrix.NumMyRows(); ++i)
  {
    Matrix.ExtractMyRowCopy(i, Length, NumEntries, &Values[0], &Indices[0]);
    for (int j = 0; j < NumEntries; ++j)
    {
      ROW[count] = Matrix.RowMatrixRowMap().GID(i);
      COL[count] = Matrix.RowMatrixColMap().GID(Indices[j]);
      VAL[count] = Values[j];
      ++count;
    }
  }

  Write(Name, "ROW", MySize, GlobalSize, H5T_NATIVE_INT, &ROW[0]);
  Write(Name, "COL", MySize, GlobalSize, H5T_NATIVE_INT, &COL[0]);
  Write(Name, "VAL", MySize, GlobalSize, H5T_NATIVE_DOUBLE, &VAL[0]);
  Write(Name, "NumGlobalRows", Matrix.NumGlobalRows());
  Write(Name, "NumGlobalCols", Matrix.NumGlobalCols());
  Write(Name, "NumGlobalNonzeros", Matrix.NumGlobalNonzeros());
  Write(Name, "NumGlobalDiagonals", Matrix.NumGlobalDiagonals());
  Write(Name, "NormOne", Matrix.NormOne());
  Write(Name, "NormInf", Matrix.NormInf());
  Write(Name, "Epetra_RowMatrix");
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, Epetra_CrsMatrix*& A)
{
  int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
  int NumGlobalDiagonals, MaxNumEntries;
  double NormOne, NormInf;

  ReadCrsMatrixProperties(GroupName, NumGlobalRows, NumGlobalCols,
                          NumGlobalNonzeros, NumGlobalDiagonals, MaxNumEntries,
                          NormOne, NormInf);

  Epetra_Map RangeMap(NumGlobalRows, 0, Comm_);
  Epetra_Map DomainMap(NumGlobalCols, 0, Comm_);

  Read(GroupName, DomainMap, RangeMap, A);
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const Epetra_Map& DomainMap, 
                        const Epetra_Map& RangeMap, Epetra_CrsMatrix*& A)
{
  int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
  int NumGlobalDiagonals, MaxNumEntries;
  double NormOne, NormInf;

  ReadCrsMatrixProperties(GroupName, NumGlobalRows, NumGlobalCols,
                          NumGlobalNonzeros, NumGlobalDiagonals, MaxNumEntries,
                          NormOne, NormInf);

  int NumMyNonzeros = NumGlobalNonzeros / Comm_.NumProc();
  if (Comm_.MyPID() == 0)
    NumMyNonzeros += NumGlobalNonzeros % Comm_.NumProc();

  vector<int> ROW(NumMyNonzeros);
  Read(GroupName, "ROW", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &ROW[0]);

  vector<int> COL(NumMyNonzeros);
  Read(GroupName, "COL", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &COL[0]);

  vector<double> VAL(NumMyNonzeros);
  Read(GroupName, "VAL", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_DOUBLE, &VAL[0]);

  Epetra_FECrsMatrix* A2 = new Epetra_FECrsMatrix(Copy, DomainMap, 0);

  for (int i = 0; i < NumMyNonzeros; )
  {
    int count = 1;
    while (ROW[i + count] == ROW[i])
      ++count;

    A2->InsertGlobalValues(1, &ROW[i], count, &COL[i], &VAL[i]);
    i += count;
  }

  A2->FillComplete(DomainMap, RangeMap);

  A = A2;
}

// ==========================================================================
void Galeri::HDF5::ReadCrsMatrixProperties(const string& GroupName, 
                                           int& NumGlobalRows,
                                           int& NumGlobalCols,
                                           int& NumGlobalNonzeros,
                                           int& NumGlobalDiagonals,
                                           int& MaxNumEntries,
                                           double& NormOne,
                                           double& NormInf)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_RowMatrix")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_RowMatrix"));

  Read(GroupName, "NumGlobalRows", NumGlobalRows);
  Read(GroupName, "NumGlobalCols", NumGlobalCols);
  Read(GroupName, "NumGlobalNonzeros", NumGlobalNonzeros);
  Read(GroupName, "NumGlobalDiagonals", NumGlobalDiagonals);
  Read(GroupName, "MaxNumEntries", MaxNumEntries);
  Read(GroupName, "NormOne", NormOne);
  Read(GroupName, "NormInf", NormInf);
}


// ==========================================================================
void Galeri::HDF5::Write(const string& GroupName, const string& DataSetName,
                         int MySize, int GlobalSize, int type, const void* data)
{
  int Offset;
  Comm_.ScanSum(&MySize, &Offset, 1);
  Offset -= MySize;

  hsize_t MySize_t = MySize;
  hsize_t GlobalSize_t = GlobalSize;
  hsize_t Offset_t = Offset;

  hid_t filespace_id = H5Screate_simple(1, &GlobalSize_t, NULL); 

  // Create the dataset with default properties and close filespace.
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), type, filespace_id, 
                            H5P_DEFAULT);

  H5Sclose(filespace_id);

  // Each process defines dataset in memory and writes it to the hyperslab
  // in the file.

  hid_t memspace_id = H5Screate_simple(1, &MySize_t, NULL);

  // Select hyperslab in the file.
  filespace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &Offset_t, NULL, &MySize_t, NULL);

  status = H5Dwrite(dset_id, type, memspace_id, filespace_id,
                    plist_id_, data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
}

// ==========================================================================
void Galeri::HDF5::Write(const string& Name, const Teuchos::ParameterList& params)
{
  hid_t group_id = H5Gcreate(file_id_, Name.c_str(), 0);

  // Iterate through parameter list 
  WriteParameterListRecursive(params, group_id);

  // Finalize hdf5 file
  status = H5Gclose(group_id);
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const string& DataSetName,
                        const int type, const int Length, 
                        void* data)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  hsize_t dimsf = Length;

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());

  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str());

  herr_t status = H5Dread(dset_id, type, H5S_ALL, filespace_id,
                          H5P_DEFAULT, data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const string& DataSetName, int& data)
{
  // FIXME: GET COMMENTS
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str());

  herr_t status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const string& DataSetName, double& data)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str());

  herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

// ==========================================================================
void Galeri::HDF5::Read(const string& Name, Teuchos::ParameterList& params) 
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  // Open hdf5 file
  hid_t       group_id;  /* identifiers */
  herr_t      status;

  // Create group in the root group using absolute name.
  group_id = H5Gopen(file_id_, Name.c_str());

  // Iterate through parameter list 
  string FullName = "/" + Name;
  status = H5Giterate(group_id, FullName.c_str() , NULL, f_operator, &params);

  // Finalize hdf5 file
  status = H5Gclose(group_id);
}

// ================== //
// Epetra_MultiVector //
// ================== //

// ==========================================================================
void Galeri::HDF5::Write(const string& GroupName, const Epetra_MultiVector& X)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  hid_t       group_id, dset_id;
  hid_t       filespace_id, memspace_id;
  herr_t      status;

  // need a linear distribution to use hyperslabs
  Epetra_Map LinearMap(X.GlobalLength(), 0, Comm_);
  Epetra_MultiVector LinearX(LinearMap, X.NumVectors());
  Epetra_Import Importer(LinearMap, X.Map());
  LinearX.Import(X, Importer, Insert);

  int NumVectors = X.NumVectors();
  int GlobalLength = X.GlobalLength();

  // Create the dataspace for the dataset.
  hsize_t q_dimsf[] = {NumVectors, GlobalLength};
  filespace_id = H5Screate_simple(2, q_dimsf, NULL);

  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  group_id = H5Gopen(file_id_, GroupName.c_str());

  // Create the dataset with default properties and close filespace_id.
  dset_id = H5Dcreate(group_id, "Values", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT);

#if 0
  // Create property list for collective dataset write.
  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif

  // Select hyperslab in the file.
  hsize_t offset[] = {0, LinearX.Map().GID(0)};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {NumVectors, 1};
  hsize_t block[] = {1, LinearX.MyLength()};
  filespace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, stride, 
                      count, block);

  // Each process defines dataset in memory and writes it to the hyperslab in the file.
  hsize_t dimsm[] = {NumVectors * LinearX.MyLength()};
  memspace_id = H5Screate_simple(1, dimsm, NULL);

  // Write hyperslab
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, 
                    plist_id_, LinearX.Values());
  H5Gclose(group_id);
  H5Sclose(memspace_id);
  H5Sclose(filespace_id);
  H5Dclose(dset_id);

#if 0
  H5Pclose(plist_id_);
#endif

  Write(GroupName, "GlobalLength", GlobalLength);
  Write(GroupName, "NumVectors", NumVectors);
  Write(GroupName, "Epetra_MultiVector");
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const Epetra_Map& Map,
                        Epetra_MultiVector*& X)
{
  // first read it with linear distribution
  Epetra_MultiVector* LinearX;
  Read(GroupName, LinearX);
  
  // now build the importer to the actual one
  Epetra_Import Importer(Map, X->Map());
  X = new Epetra_MultiVector(Map, LinearX->NumVectors());
  X->Import(*LinearX, Importer, Insert);

  // delete memory
  delete LinearX;
}

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, Epetra_MultiVector*& LinearX)
{
  int GlobalLength, NumVectors;

  ReadMultiVectorProperties(GroupName, GlobalLength, NumVectors);

  cout << "GlobalLength = " << GlobalLength << endl;
  cout << "NumVectors = " << NumVectors << endl;

  hid_t group_id, dataset_id, space_id;
  hsize_t dims[2];
  herr_t status;

  group_id = H5Gopen(file_id_, GroupName.c_str());

  dataset_id = H5Dopen(group_id, "Values");
  space_id = H5Dget_space(dataset_id);
  if (H5Sget_simple_extent_ndims(space_id) != 2)
    throw(Exception(__FILE__, __LINE__,
                    "internal error, dimensions do not match"));

  // FIXME: this works, but I don't like it...
  status = H5Sget_simple_extent_dims(space_id, dims, NULL);
  cout << dims[0] << "   " << dims[1] << endl;
  if (dims[0] != NumVectors)
    throw(Exception(__FILE__, __LINE__,
                    "internal error, NumVectors does not match"));
  Epetra_Map LinearMap(-1, dims[1], 0, Comm());
  LinearX = new Epetra_MultiVector(LinearMap, NumVectors);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                   H5P_DEFAULT, LinearX->Values());
  H5Sclose(space_id);  
  H5Dclose(dataset_id);  
  H5Gclose(group_id);

  cout << "======================" << endl;
}

// ==========================================================================
void Galeri::HDF5::ReadMultiVectorProperties(const string& GroupName, 
                                             int& GlobalLength,
                                             int& NumVectors)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  string Label;
  Read(GroupName, Label);

  if (Label != "Epetra_MultiVector")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_MultiVector"));

  Read(GroupName, "GlobalLength", GlobalLength);
  Read(GroupName, "NumVectors", NumVectors);
}

// ============ //
// general data //
// ============ //

// ==========================================================================
void Galeri::HDF5::Read(const string& GroupName, const string& DataSetName,
                        int MySize, int GlobalSize,
                        const int type, void* data)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  // global size of the data to be read
  hsize_t MySize_t = MySize;
  hsize_t GlobalSize_t = GlobalSize;

  // offset
  int itmp;
  Comm_.ScanSum(&MySize, &itmp, 1);
  hsize_t Offset_t = itmp - MySize;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str());
  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str());
  //hid_t space_id = H5Screate_simple(1, &Offset_t, 0);

  // Select hyperslab in the file.
  hid_t filespace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &Offset_t, NULL, 
                      &MySize_t, NULL);

  hid_t mem_dataspace = H5Screate_simple (1, &MySize_t, NULL);

  herr_t status = H5Dread(dataset_id, type, mem_dataspace, filespace_id, 
                          H5P_DEFAULT, data);

  H5Sclose(mem_dataspace);
  H5Gclose(group_id);  
  //H5Sclose(space_id);  
  H5Dclose(dataset_id);  
//  H5Dclose(filespace_id);  
}

#endif
