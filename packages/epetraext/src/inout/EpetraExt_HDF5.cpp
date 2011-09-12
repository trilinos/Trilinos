/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "EpetraExt_ConfigDefs.h"


#ifdef HAVE_EPETRAEXT_HDF5


#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
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
#include "EpetraExt_Exception.h"
#include "EpetraExt_Utils.h"
#include "EpetraExt_HDF5.h"
#include "EpetraExt_DistArray.h"

#define CHECK_HID(hid_t) \
  { if (hid_t < 0) \
    throw(EpetraExt::Exception(__FILE__, __LINE__, \
                    "hid_t is negative")); }

#define CHECK_STATUS(status) \
  { if (status < 0) \
    throw(EpetraExt::Exception(__FILE__, __LINE__, \
                    "function H5Giterater returned a negative value")); }

// ==========================================================================
// data container and iterators to find a dataset with a given name
struct FindDataset_t
{
  std::string name;
  bool found;
};

static herr_t FindDataset(hid_t loc_id, const char *name, void *opdata)
{
  std::string& token = ((FindDataset_t*)opdata)->name;
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
    std::string key(params.name(it));
    if (params.isSublist(key)) 
    {
      // Sublist

      // Create subgroup for sublist.
      hid_t child_group_id = H5Gcreate(group_id, key.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
      if (params.isType<std::string>(key)) 
      {
        std::string value = params.get<std::string>(key);
        hsize_t len = value.size() + 1;
        dataspace_id = H5Screate_simple(1, &len, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, value.c_str());
        CHECK_STATUS(status);
        status = H5Dclose(dataset_id);
        CHECK_STATUS(status);
        status = H5Sclose(dataspace_id);
        CHECK_STATUS(status);
        found = true;
      } 

      if (params.isType<bool>(key)) 
      {
        // Use H5T_NATIVE_USHORT to store a bool value
        unsigned short value = params.get<bool>(key) ? 1 : 0;
        dataspace_id = H5Screate_simple(1, &one, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, &value);
        CHECK_STATUS(status);
        status = H5Dclose(dataset_id);
        CHECK_STATUS(status);
        status = H5Sclose(dataspace_id);
        CHECK_STATUS(status);
        found = true;
      } 
      
      if (params.isType<int>(key)) 
      {
        int value = params.get<int>(key);
        dataspace_id = H5Screate_simple(1, &one, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, &value);
        CHECK_STATUS(status);
        status = H5Dclose(dataset_id);
        CHECK_STATUS(status);
        status = H5Sclose(dataspace_id);
        CHECK_STATUS(status);
        found = true;
      } 

      if (params.isType<double>(key)) 
      {
        double value = params.get<double>(key);
        dataspace_id = H5Screate_simple(1, &one, NULL);
        dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                          H5P_DEFAULT, &value);
        CHECK_STATUS(status);
        status = H5Dclose(dataset_id);
        CHECK_STATUS(status);
        status = H5Sclose(dataspace_id);
        CHECK_STATUS(status);
        found = true;
      } 

      if (!found)
      {
        throw(EpetraExt::Exception(__FILE__, __LINE__, 
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
  if (strcmp(name, "__type__") == 0)
    return(0); // skip list identifier

  switch (statbuf.type) {
  case H5G_GROUP:
    sublist = &params->sublist(name);
    H5Giterate(loc_id, name , NULL, f_operator, sublist);
    break;
  case H5G_DATASET:
    hsize_t len;
    dataset_id = H5Dopen(loc_id, name, H5P_DEFAULT);
    space_id = H5Dget_space(dataset_id);
    if (H5Sget_simple_extent_ndims(space_id) != 1)
      throw(EpetraExt::Exception(__FILE__, __LINE__, 
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
      throw(EpetraExt::Exception(__FILE__, __LINE__,
                              "unsupported datatype")); // FIXME
    }
    H5Tclose(type_id);
    H5Sclose(space_id);  
    H5Dclose(dataset_id);  
    break;
  default:
    throw(EpetraExt::Exception(__FILE__, __LINE__,
                            "unsupported datatype")); // FIXME
  }
  return 0;
}

// ==========================================================================
EpetraExt::HDF5::HDF5(const Epetra_Comm& Comm) :
  Comm_(Comm),
  IsOpen_(false)
{}

// ==========================================================================
void EpetraExt::HDF5::Create(const std::string FileName)
{
  if (IsOpen())
    throw(Exception(__FILE__, __LINE__,
                    "an HDF5 is already open, first close the current one",
                    "using method Close(), then open/create a new one"));

  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
  // Create property list for collective dataset write.
  H5Pset_fapl_mpio(plist_id_, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif

#if 0
  unsigned int boh = H5Z_FILTER_MAX;
  H5Pset_filter(plist_id_, H5Z_FILTER_DEFLATE, H5Z_FILTER_MAX, 0, &boh);
#endif

  // create the file collectively and release property list identifier.
  file_id_ = H5Fcreate(FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, 
                      plist_id_);
  H5Pclose(plist_id_);

  IsOpen_ = true;
}

// ==========================================================================
void EpetraExt::HDF5::Open(const std::string FileName, int AccessType)
{
  if (IsOpen())
    throw(Exception(__FILE__, __LINE__,
                    "an HDF5 is already open, first close the current one",
                    "using method Close(), then open/create a new one"));

  FileName_ = FileName;

  // create the file collectively and release property list identifier.
  file_id_ = H5Fopen(FileName.c_str(), AccessType, H5P_DEFAULT);

#ifdef HAVE_MPI
// FIXME: DO I NEED THE MPIO_COLLECTIVE??
//  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
//  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif

  IsOpen_ = true;
}

// ==========================================================================
bool EpetraExt::HDF5::IsContained(std::string Name)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  FindDataset_t data;
  data.name = Name;
  data.found = false;

  //int idx_f =
  H5Giterate(file_id_, "/", NULL, FindDataset, (void*)&data);

  return(data.found);
}

// ============================ //
// Epetra_BlockMap / Epetra_Map //
// ============================ //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& Name, const Epetra_BlockMap& BlockMap)
{
  int NumMyPoints       = BlockMap.NumMyPoints();
  int NumMyElements     = BlockMap.NumMyElements();
  int NumGlobalElements = BlockMap.NumGlobalElements();
  int NumGlobalPoints   = BlockMap.NumGlobalPoints();
  int* MyGlobalElements = BlockMap.MyGlobalElements();
  int* ElementSizeList  = BlockMap.ElementSizeList();

  std::vector<int> NumMyElements_v(Comm_.NumProc());
  Comm_.GatherAll(&NumMyElements, &NumMyElements_v[0], 1);

  std::vector<int> NumMyPoints_v(Comm_.NumProc());
  Comm_.GatherAll(&NumMyPoints, &NumMyPoints_v[0], 1);

  Write(Name, "MyGlobalElements", NumMyElements, NumGlobalElements, 
        H5T_NATIVE_INT, MyGlobalElements);
  Write(Name, "ElementSizeList", NumMyElements, NumGlobalElements, 
        H5T_NATIVE_INT, ElementSizeList);
  Write(Name, "NumMyPoints", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyPoints_v[0]);

  // need to know how many processors currently host this map
  Write(Name, "NumProc", Comm_.NumProc());
  // write few more data about this map
  Write(Name, "NumGlobalPoints", 1, Comm_.NumProc(), H5T_NATIVE_INT, &NumGlobalPoints);
  Write(Name, "NumGlobalElements", 1, Comm_.NumProc(), H5T_NATIVE_INT, &NumGlobalElements);
  Write(Name, "IndexBase", BlockMap.IndexBase());
  Write(Name, "__type__", "Epetra_BlockMap");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Epetra_BlockMap*& BlockMap)
{
  int NumGlobalElements, NumGlobalPoints, IndexBase, NumProc;

  ReadBlockMapProperties(GroupName, NumGlobalElements, NumGlobalPoints,
                        IndexBase, NumProc);

  std::vector<int> NumMyPoints_v(Comm_.NumProc());
  std::vector<int> NumMyElements_v(Comm_.NumProc());

  Read(GroupName, "NumMyElements", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyElements_v[0]);
  Read(GroupName, "NumMyPoints", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyPoints_v[0]);
  int NumMyElements = NumMyElements_v[Comm_.MyPID()];
//  int NumMyPoints   = NumMyPoints_v[Comm_.MyPID()];

  if (NumProc != Comm_.NumProc())
    throw(Exception(__FILE__, __LINE__,
                    "requested map not compatible with current number of",
                    "processors, " + toString(Comm_.NumProc()) + 
                    " vs. " + toString(NumProc)));

  std::vector<int> MyGlobalElements(NumMyElements);
  std::vector<int> ElementSizeList(NumMyElements);

  Read(GroupName, "MyGlobalElements", NumMyElements, NumGlobalElements, 
       H5T_NATIVE_INT, &MyGlobalElements[0]);

  Read(GroupName, "ElementSizeList", NumMyElements, NumGlobalElements, 
       H5T_NATIVE_INT, &ElementSizeList[0]);

  BlockMap = new Epetra_BlockMap(NumGlobalElements, NumMyElements, 
                                 &MyGlobalElements[0], &ElementSizeList[0],
                                 IndexBase, Comm_);
}

// ==========================================================================
void EpetraExt::HDF5::ReadBlockMapProperties(const std::string& GroupName, 
                                          int& NumGlobalElements,
                                          int& NumGlobalPoints,
                                          int& IndexBase,
                                          int& NumProc)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_BlockMap")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_BlockMap",
                    "__type__ = " + Label));

  Read(GroupName, "NumGlobalElements", NumGlobalElements);
  Read(GroupName, "NumGlobalPoints", NumGlobalPoints);
  Read(GroupName, "IndexBase", IndexBase);
  Read(GroupName, "NumProc", NumProc);
}

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& Name, const Epetra_Map& Map)
{
  int MySize = Map.NumMyElements();
  int GlobalSize = Map.NumGlobalElements();
  int* MyGlobalElements = Map.MyGlobalElements();

  std::vector<int> NumMyElements(Comm_.NumProc());
  Comm_.GatherAll(&MySize, &NumMyElements[0], 1);

  Write(Name, "MyGlobalElements", MySize, GlobalSize, 
        H5T_NATIVE_INT, MyGlobalElements);
  Write(Name, "NumMyElements", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyElements[0]);
  Write(Name, "NumGlobalElements", 1, Comm_.NumProc(), H5T_NATIVE_INT, &GlobalSize);
  Write(Name, "NumProc", Comm_.NumProc());
  Write(Name, "IndexBase", Map.IndexBase());
  Write(Name, "__type__", "Epetra_Map");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Epetra_Map*& Map)
{
  int NumGlobalElements, IndexBase, NumProc;

  ReadMapProperties(GroupName, NumGlobalElements, IndexBase, NumProc);

  std::vector<int> NumMyElements_v(Comm_.NumProc());

  Read(GroupName, "NumMyElements", H5T_NATIVE_INT, Comm_.NumProc(), &NumMyElements_v[0]);
  int NumMyElements = NumMyElements_v[Comm_.MyPID()];

  if (NumProc != Comm_.NumProc())
    throw(Exception(__FILE__, __LINE__,
                    "requested map not compatible with current number of",
                    "processors, " + toString(Comm_.NumProc()) + 
                    " vs. " + toString(NumProc)));

  std::vector<int> MyGlobalElements(NumMyElements);

  Read(GroupName, "MyGlobalElements", NumMyElements, NumGlobalElements, 
       H5T_NATIVE_INT, &MyGlobalElements[0]);

  Map = new Epetra_Map(-1, NumMyElements, &MyGlobalElements[0], IndexBase, Comm_);
}

// ==========================================================================
void EpetraExt::HDF5::ReadMapProperties(const std::string& GroupName, 
                                     int& NumGlobalElements,
                                     int& IndexBase,
                                     int& NumProc)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_Map")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_Map",
                    "__type__ = " + Label));

  Read(GroupName, "NumGlobalElements", NumGlobalElements);
  Read(GroupName, "IndexBase", IndexBase);
  Read(GroupName, "NumProc", NumProc);
}

// ================ //
// Epetra_IntVector //
// ================ //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& Name, const Epetra_IntVector& x)
{
  if (x.Map().LinearMap())
  {
    Write(Name, "GlobalLength", x.GlobalLength());
    Write(Name, "Values", x.Map().NumMyElements(), x.Map().NumGlobalElements(),
          H5T_NATIVE_INT, x.Values());
  }
  else
  {
    // need to build a linear map first, the import data, then
    // finally write them 
    const Epetra_BlockMap& OriginalMap = x.Map();
    Epetra_Map LinearMap(OriginalMap.NumGlobalElements(), OriginalMap.IndexBase(), Comm_);
    Epetra_Import Importer(LinearMap, OriginalMap);

    Epetra_IntVector LinearX(LinearMap);
    LinearX.Import(x, Importer, Insert);

    Write(Name, "GlobalLength", x.GlobalLength());
    Write(Name, "Values", LinearMap.NumMyElements(), LinearMap.NumGlobalElements(),
          H5T_NATIVE_INT, LinearX.Values());
  }
  Write(Name, "__type__", "Epetra_IntVector");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Epetra_IntVector*& X)
{
  // gets the length of the std::vector
  int GlobalLength;
  ReadIntVectorProperties(GroupName, GlobalLength);

  // creates a new linear map and use it to read data
  Epetra_Map LinearMap(GlobalLength, 0, Comm_);
  X = new Epetra_IntVector(LinearMap);

  Read(GroupName, "Values", LinearMap.NumMyElements(), 
       LinearMap.NumGlobalElements(), H5T_NATIVE_INT, X->Values());
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const Epetra_Map& Map,
                        Epetra_IntVector*& X)
{
  // gets the length of the std::vector
  int GlobalLength;
  ReadIntVectorProperties(GroupName, GlobalLength);

  if (Map.LinearMap())
  {
    X = new Epetra_IntVector(Map);
    // simply read stuff and go home
    Read(GroupName, "Values", Map.NumMyElements(), Map.NumGlobalElements(),
         H5T_NATIVE_INT, X->Values());

  }
  else
  {
    // we need to first create a linear map, read the std::vector,
    // then import it to the actual nonlinear map
    Epetra_Map LinearMap(GlobalLength, Map.IndexBase(), Comm_);
    Epetra_IntVector LinearX(LinearMap);

    Read(GroupName, "Values", LinearMap.NumMyElements(), 
         LinearMap.NumGlobalElements(),
         H5T_NATIVE_INT, LinearX.Values());

    Epetra_Import Importer(Map, LinearMap);
    X = new Epetra_IntVector(Map);
    X->Import(LinearX, Importer, Insert);
  }
}

// ==========================================================================
void EpetraExt::HDF5::ReadIntVectorProperties(const std::string& GroupName, 
                                           int& GlobalLength)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_IntVector")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_IntVector",
                    "__type__ = " + Label));

  Read(GroupName, "GlobalLength", GlobalLength);
}

// =============== //
// Epetra_CrsGraph //
// =============== //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& Name, const Epetra_CrsGraph& Graph)
{
  if (!Graph.Filled())
    throw(Exception(__FILE__, __LINE__,
                    "input Epetra_CrsGraph is not FillComplete()'d"));

  // like RowMatrix, only without values
  int MySize = Graph.NumMyNonzeros();
  int GlobalSize = Graph.NumGlobalNonzeros();
  std::vector<int> ROW(MySize);
  std::vector<int> COL(MySize);

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
  Write(Name, "__type__", "Epetra_CrsGraph");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Epetra_CrsGraph*& Graph)
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
void EpetraExt::HDF5::ReadCrsGraphProperties(const std::string& GroupName, 
                                          int& NumGlobalRows,
                                          int& NumGlobalCols,
                                          int& NumGlobalNonzeros,
                                          int& NumGlobalDiagonals,
                                          int& MaxNumIndices)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_CrsGraph")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_CrsGraph",
                    "__type__ = " + Label));

  Read(GroupName, "NumGlobalRows",      NumGlobalRows);
  Read(GroupName, "NumGlobalCols",      NumGlobalCols);
  Read(GroupName, "NumGlobalNonzeros",  NumGlobalNonzeros);
  Read(GroupName, "NumGlobalDiagonals", NumGlobalDiagonals);
  Read(GroupName, "MaxNumIndices",      MaxNumIndices);
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const Epetra_Map& DomainMap, 
                        const Epetra_Map& RangeMap, Epetra_CrsGraph*& Graph)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found in database"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_CrsGraph")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_CrsGraph",
                    "__type__ = " + Label));

  int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
  Read(GroupName, "NumGlobalNonzeros", NumGlobalNonzeros);
  Read(GroupName, "NumGlobalRows", NumGlobalRows);
  Read(GroupName, "NumGlobalCols", NumGlobalCols);

  // linear distribution for nonzeros
  int NumMyNonzeros = NumGlobalNonzeros / Comm_.NumProc();
  if (Comm_.MyPID() == 0)
    NumMyNonzeros += NumGlobalNonzeros % Comm_.NumProc();

  std::vector<int> ROW(NumMyNonzeros);
  Read(GroupName, "ROW", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &ROW[0]);

  std::vector<int> COL(NumMyNonzeros);
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
void EpetraExt::HDF5::Write(const std::string& Name, const Epetra_RowMatrix& Matrix)
{
  if (!Matrix.Filled())
    throw(Exception(__FILE__, __LINE__,
                    "input Epetra_RowMatrix is not FillComplete()'d"));

  int MySize = Matrix.NumMyNonzeros();
  int GlobalSize = Matrix.NumGlobalNonzeros();
  std::vector<int> ROW(MySize);
  std::vector<int> COL(MySize);
  std::vector<double> VAL(MySize);

  int count = 0;
  int Length = Matrix.MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);
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
  Write(Name, "MaxNumEntries", Matrix.MaxNumEntries());
  Write(Name, "NormOne", Matrix.NormOne());
  Write(Name, "NormInf", Matrix.NormInf());
  Write(Name, "__type__", "Epetra_RowMatrix");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Epetra_CrsMatrix*& A)
{
  int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
  int NumGlobalDiagonals, MaxNumEntries;
  double NormOne, NormInf;

  ReadCrsMatrixProperties(GroupName, NumGlobalRows, NumGlobalCols,
                          NumGlobalNonzeros, NumGlobalDiagonals, MaxNumEntries,
                          NormOne, NormInf);

  // build simple linear maps for domain and range space
  Epetra_Map RangeMap(NumGlobalRows, 0, Comm_);
  Epetra_Map DomainMap(NumGlobalCols, 0, Comm_);

  Read(GroupName, DomainMap, RangeMap, A);
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const Epetra_Map& DomainMap, 
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

  std::vector<int> ROW(NumMyNonzeros);
  Read(GroupName, "ROW", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &ROW[0]);

  std::vector<int> COL(NumMyNonzeros);
  Read(GroupName, "COL", NumMyNonzeros, NumGlobalNonzeros, H5T_NATIVE_INT, &COL[0]);

  std::vector<double> VAL(NumMyNonzeros);
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
void EpetraExt::HDF5::ReadCrsMatrixProperties(const std::string& GroupName, 
                                           int& NumGlobalRows,
                                           int& NumGlobalCols,
                                           int& NumGlobalNonzeros,
                                           int& NumGlobalDiagonals,
                                           int& MaxNumEntries,
                                           double& NormOne,
                                           double& NormInf)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_RowMatrix")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_RowMatrix",
                    "__type__ = " + Label));

  Read(GroupName, "NumGlobalRows", NumGlobalRows);
  Read(GroupName, "NumGlobalCols", NumGlobalCols);
  Read(GroupName, "NumGlobalNonzeros", NumGlobalNonzeros);
  Read(GroupName, "NumGlobalDiagonals", NumGlobalDiagonals);
  Read(GroupName, "MaxNumEntries", MaxNumEntries);
  Read(GroupName, "NormOne", NormOne);
  Read(GroupName, "NormInf", NormInf);
}

// ============= //
// ParameterList //
// ============= //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const Teuchos::ParameterList& params)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Iterate through parameter list 
  WriteParameterListRecursive(params, group_id);

  // Finalize hdf5 file
  status = H5Gclose(group_id);
  CHECK_STATUS(status);

  Write(GroupName, "__type__", "Teuchos::ParameterList");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Teuchos::ParameterList& params) 
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Teuchos::ParameterList")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not a Teuchos::ParameterList",
                    "__type__ = " + Label));

  // Open hdf5 file
  hid_t       group_id;  /* identifiers */
  herr_t      status;

  // open group in the root group using absolute name.
  group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  CHECK_HID(group_id);

  // Iterate through parameter list 
  std::string xxx = "/" + GroupName;
  status = H5Giterate(group_id, xxx.c_str() , NULL, f_operator, &params);
  CHECK_STATUS(status);

  // Finalize hdf5 file
  status = H5Gclose(group_id);
  CHECK_STATUS(status);
}

// ================== //
// Epetra_MultiVector //
// ================== //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const Epetra_MultiVector& X, bool writeTranspose)
{

  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  hid_t       group_id, dset_id;
  hid_t       filespace_id, memspace_id;
  herr_t      status;

  // need a linear distribution to use hyperslabs
  Teuchos::RefCountPtr<Epetra_MultiVector> LinearX;

  if (X.Map().LinearMap())
    LinearX = Teuchos::rcp(const_cast<Epetra_MultiVector*>(&X), false);
  else
    {
      Epetra_Map LinearMap(X.GlobalLength(), X.Map().IndexBase(), Comm_);
      LinearX = Teuchos::rcp(new Epetra_MultiVector(LinearMap, X.NumVectors()));
      Epetra_Import Importer(LinearMap, X.Map());
      LinearX->Import(X, Importer, Insert);
    }

  int NumVectors = X.NumVectors();
  int GlobalLength = X.GlobalLength();

  // Whether or not we do writeTranspose or not is
  // handled by one of the components of q_dimsf, offset and count. 
  // They are determined by indexT
  int indexT(0);
  if (writeTranspose) indexT = 1;

  hsize_t q_dimsf[] = {GlobalLength, GlobalLength};
  q_dimsf[indexT] = NumVectors;

  filespace_id = H5Screate_simple(2, q_dimsf, NULL);

  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  // Create the dataset with default properties and close filespace_id.
  dset_id = H5Dcreate(group_id, "Values", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create property list for collective dataset write.
  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif


  // Select hyperslab in the file.
  hsize_t offset[] = {LinearX->Map().GID(0)-X.Map().IndexBase(),
		      LinearX->Map().GID(0)-X.Map().IndexBase()};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {LinearX->MyLength(),
		     LinearX->MyLength()};
  hsize_t block[] = {1, 1};
    
  // write vectors one by one
  for (int n(0); n < NumVectors; ++n)
    {
      // Select hyperslab in the file.
      offset[indexT] = n;
      count [indexT] = 1;

      filespace_id = H5Dget_space(dset_id);
      H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, stride,
			  count, block);

      // Each process defines dataset in memory and writes it to the hyperslab in the file.
      hsize_t dimsm[] = {LinearX->MyLength()};
      memspace_id = H5Screate_simple(1, dimsm, NULL);

      // Write hyperslab
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id,
			plist_id_, LinearX->operator[](n));
      CHECK_STATUS(status);
    }
  H5Gclose(group_id);
  H5Sclose(memspace_id);
  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id_);

  Write(GroupName, "GlobalLength", GlobalLength);
  Write(GroupName, "NumVectors", NumVectors);
  Write(GroupName, "__type__", "Epetra_MultiVector");
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const Epetra_Map& Map,
                        Epetra_MultiVector*& X, bool writeTranspose)
{
  // first read it with linear distribution
  Epetra_MultiVector* LinearX;
  Read(GroupName, LinearX, writeTranspose, Map.IndexBase());
  
  // now build the importer to the actual one
  Epetra_Import Importer(Map, LinearX->Map());
  X = new Epetra_MultiVector(Map, LinearX->NumVectors());
  X->Import(*LinearX, Importer, Insert);

  // delete memory
  delete LinearX;
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Epetra_MultiVector*& LinearX,
                           bool readTranspose, const int& indexBase)
{
  int GlobalLength, NumVectors;

  ReadMultiVectorProperties(GroupName, GlobalLength, NumVectors);

  hid_t group_id;
  hid_t memspace_id;

  // Whether or not we do readTranspose or not is
  // handled by one of the components of q_dimsf, offset and count.
  // They are determined by indexT
  int indexT(0);
  if (readTranspose) indexT = 1;

  hsize_t q_dimsf[] = {GlobalLength, GlobalLength};
  q_dimsf[indexT] = NumVectors;

  hid_t filespace_id = H5Screate_simple(2, q_dimsf, NULL);

  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  // Create the dataset with default properties and close filespace_id.
  hid_t dset_id = H5Dopen(group_id, "Values", H5P_DEFAULT);

  // Create property list for collective dataset write.
  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif
  H5Pclose(plist_id_);

  Epetra_Map LinearMap(GlobalLength, indexBase, Comm());
  LinearX = new Epetra_MultiVector(LinearMap, NumVectors);

  // Select hyperslab in the file.
  hsize_t offset[] = {LinearMap.GID(0) - indexBase, LinearMap.GID(0) - indexBase};
  hsize_t stride[] = {1, 1};

  // If readTranspose is false, we can read the data in one shot.
  // It would actually be possible to skip this first part and
  if (!readTranspose)
  {
  // Select hyperslab in the file.
  hsize_t count[] = {1, 1};
  hsize_t block[] = {LinearX->MyLength(), LinearX->MyLength()};

  offset[indexT]  = 0;
  count [indexT]  = NumVectors;
  block [indexT]  = 1;

  filespace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, stride, 
                      count, block);

  // Each process defines dataset in memory and writes it to the hyperslab in the file.
  hsize_t dimsm[] = {NumVectors * LinearX->MyLength()};
  memspace_id = H5Screate_simple(1, dimsm, NULL);

  // Write hyperslab
  CHECK_STATUS(H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, 
                       H5P_DEFAULT, LinearX->Values()));

  } else {
    // doing exactly the same as in write

    // Select hyperslab in the file.
    hsize_t count[] = {LinearX->MyLength(),
		       LinearX->MyLength()};
    hsize_t block[] = {1, 1};

    // write vectors one by one
    for (int n(0); n < NumVectors; ++n)
    {
      // Select hyperslab in the file.
      offset[indexT] = n;
      count [indexT] = 1;

      filespace_id = H5Dget_space(dset_id);
      H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, stride,
			  count, block);

      // Each process defines dataset in memory and writes it to the hyperslab in the file.
      hsize_t dimsm[] = {LinearX->MyLength()};
      memspace_id = H5Screate_simple(1, dimsm, NULL);

      // Read hyperslab
      CHECK_STATUS(H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id,
			   H5P_DEFAULT, LinearX->operator[](n)));

    }
  }

  CHECK_STATUS(H5Gclose(group_id));
  CHECK_STATUS(H5Sclose(memspace_id));
  CHECK_STATUS(H5Sclose(filespace_id));
  CHECK_STATUS(H5Dclose(dset_id));
}

// ==========================================================================
void EpetraExt::HDF5::ReadMultiVectorProperties(const std::string& GroupName, 
                                             int& GlobalLength,
                                             int& NumVectors)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "Epetra_MultiVector")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an Epetra_MultiVector",
                    "__type__ = " + Label));

  Read(GroupName, "GlobalLength", GlobalLength);
  Read(GroupName, "NumVectors", NumVectors);
}

// ========================= //
// EpetraExt::DistArray<int> //
// ========================= //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const DistArray<int>& x)
{
  if (x.Map().LinearMap())
  {
    Write(GroupName, "Values", x.Map().NumMyElements() * x.RowSize(), 
          x.Map().NumGlobalElements() * x.RowSize(),
          H5T_NATIVE_INT, x.Values());
  }
  else
  {
    // need to build a linear map first, the import data, then
    // finally write them 
    const Epetra_BlockMap& OriginalMap = x.Map();
    Epetra_Map LinearMap(OriginalMap.NumGlobalElements(), OriginalMap.IndexBase(), Comm_);
    Epetra_Import Importer(LinearMap, OriginalMap);

    EpetraExt::DistArray<int> LinearX(LinearMap, x.RowSize());
    LinearX.Import(x, Importer, Insert);

    Write(GroupName, "Values", LinearMap.NumMyElements() * x.RowSize(), 
          LinearMap.NumGlobalElements() * x.RowSize(),
          H5T_NATIVE_INT, LinearX.Values());
  }

  Write(GroupName, "__type__", "EpetraExt::DistArray<int>");
  Write(GroupName, "GlobalLength", x.GlobalLength());
  Write(GroupName, "RowSize", x.RowSize());
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const Epetra_Map& Map,
                           DistArray<int>*& X)
{
  // gets the length of the std::vector
  int GlobalLength, RowSize;
  ReadIntDistArrayProperties(GroupName, GlobalLength, RowSize);

  if (Map.LinearMap())
  {
    X = new EpetraExt::DistArray<int>(Map, RowSize);
    // simply read stuff and go home
    Read(GroupName, "Values", Map.NumMyElements() * RowSize, 
         Map.NumGlobalElements() * RowSize, H5T_NATIVE_INT, X->Values());
  }
  else
  {
    // we need to first create a linear map, read the std::vector,
    // then import it to the actual nonlinear map
    Epetra_Map LinearMap(GlobalLength, Map.IndexBase(), Comm_);
    EpetraExt::DistArray<int> LinearX(LinearMap, RowSize);

    Read(GroupName, "Values", LinearMap.NumMyElements() * RowSize, 
         LinearMap.NumGlobalElements() * RowSize,
         H5T_NATIVE_INT, LinearX.Values());

    Epetra_Import Importer(Map, LinearMap);
    X = new EpetraExt::DistArray<int>(Map, RowSize);
    X->Import(LinearX, Importer, Insert);
  }
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, DistArray<int>*& X)
{
  // gets the length of the std::vector
  int GlobalLength, RowSize;
  ReadIntDistArrayProperties(GroupName, GlobalLength, RowSize);

  // creates a new linear map and use it to read data
  Epetra_Map LinearMap(GlobalLength, 0, Comm_);
  X = new EpetraExt::DistArray<int>(LinearMap, RowSize);

  Read(GroupName, "Values", LinearMap.NumMyElements() * RowSize, 
       LinearMap.NumGlobalElements() * RowSize, H5T_NATIVE_INT, X->Values());
}

// ==========================================================================
void EpetraExt::HDF5::ReadIntDistArrayProperties(const std::string& GroupName, 
                                                 int& GlobalLength,
                                                 int& RowSize)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "EpetraExt::DistArray<int>")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an EpetraExt::DistArray<int>",
                    "__type__ = " + Label));

  Read(GroupName, "GlobalLength", GlobalLength);
  Read(GroupName, "RowSize", RowSize);
}

// ============================ //
// EpetraExt::DistArray<double> //
// ============================ //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const DistArray<double>& x)
{
  if (x.Map().LinearMap())
  {
    Write(GroupName, "Values", x.Map().NumMyElements() * x.RowSize(), 
          x.Map().NumGlobalElements() * x.RowSize(),
          H5T_NATIVE_DOUBLE, x.Values());
  }
  else
  {
    // need to build a linear map first, the import data, then
    // finally write them 
    const Epetra_BlockMap& OriginalMap = x.Map();
    Epetra_Map LinearMap(OriginalMap.NumGlobalElements(), OriginalMap.IndexBase(), Comm_);
    Epetra_Import Importer(LinearMap, OriginalMap);

    EpetraExt::DistArray<double> LinearX(LinearMap, x.RowSize());
    LinearX.Import(x, Importer, Insert);

    Write(GroupName, "Values", LinearMap.NumMyElements() * x.RowSize(), 
          LinearMap.NumGlobalElements() * x.RowSize(),
          H5T_NATIVE_DOUBLE, LinearX.Values());
  }

  Write(GroupName, "__type__", "EpetraExt::DistArray<double>");
  Write(GroupName, "GlobalLength", x.GlobalLength());
  Write(GroupName, "RowSize", x.RowSize());
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const Epetra_Map& Map,
                           DistArray<double>*& X)
{
  // gets the length of the std::vector
  int GlobalLength, RowSize;
  ReadDoubleDistArrayProperties(GroupName, GlobalLength, RowSize);

  if (Map.LinearMap())
  {
    X = new EpetraExt::DistArray<double>(Map, RowSize);
    // simply read stuff and go home
    Read(GroupName, "Values", Map.NumMyElements() * RowSize, 
         Map.NumGlobalElements() * RowSize, H5T_NATIVE_DOUBLE, X->Values());
  }
  else
  {
    // we need to first create a linear map, read the std::vector,
    // then import it to the actual nonlinear map
    Epetra_Map LinearMap(GlobalLength, Map.IndexBase(), Comm_);
    EpetraExt::DistArray<double> LinearX(LinearMap, RowSize);

    Read(GroupName, "Values", LinearMap.NumMyElements() * RowSize, 
         LinearMap.NumGlobalElements() * RowSize,
         H5T_NATIVE_DOUBLE, LinearX.Values());

    Epetra_Import Importer(Map, LinearMap);
    X = new EpetraExt::DistArray<double>(Map, RowSize);
    X->Import(LinearX, Importer, Insert);
  }
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, DistArray<double>*& X)
{
  // gets the length of the std::vector
  int GlobalLength, RowSize;
  ReadDoubleDistArrayProperties(GroupName, GlobalLength, RowSize);

  // creates a new linear map and use it to read data
  Epetra_Map LinearMap(GlobalLength, 0, Comm_);
  X = new EpetraExt::DistArray<double>(LinearMap, RowSize);

  Read(GroupName, "Values", LinearMap.NumMyElements() * RowSize, 
       LinearMap.NumGlobalElements() * RowSize, H5T_NATIVE_DOUBLE, X->Values());
}
//
// ==========================================================================
void EpetraExt::HDF5::ReadDoubleDistArrayProperties(const std::string& GroupName, 
                                                    int& GlobalLength,
                                                    int& RowSize)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  std::string Label;
  Read(GroupName, "__type__", Label);

  if (Label != "EpetraExt::DistArray<double>")
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not an EpetraExt::DistArray<double>",
                    "__type__ = " + Label));

  Read(GroupName, "GlobalLength", GlobalLength);
  Read(GroupName, "RowSize", RowSize);
}

// ======================= //
// basic serial data types //
// ======================= //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const std::string& DataSetName,
                            int what)
{
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_INT, 
                      filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                           H5P_DEFAULT, &what);
  CHECK_STATUS(status);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const std::string& DataSetName,
                            double what)
{
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_DOUBLE, 
                            filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, 
                           filespace_id, H5P_DEFAULT, &what);
  CHECK_STATUS(status);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const std::string& DataSetName, int& data)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);
  CHECK_STATUS(status);

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const std::string& DataSetName, double& data)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);
  CHECK_STATUS(status);

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, 
                            const std::string& DataSetName, 
                            const std::string& data)
{
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hsize_t len = 1;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataspace_id = H5Screate_simple(1, &len, NULL);

  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, data.size() + 1);

  hid_t dataset_id = H5Dcreate(group_id, DataSetName.c_str(), atype, 
                               dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_STATUS(H5Dwrite(dataset_id, atype, H5S_ALL, H5S_ALL, 
                        H5P_DEFAULT, data.c_str()));

  CHECK_STATUS(H5Tclose(atype));

  CHECK_STATUS(H5Dclose(dataset_id));

  CHECK_STATUS(H5Sclose(dataspace_id));

  CHECK_STATUS(H5Gclose(group_id));
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, 
                           const std::string& DataSetName, 
                           std::string& data)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  hid_t datatype_id = H5Dget_type(dataset_id);
//  size_t typesize_id = H5Tget_size(datatype_id);
  H5T_class_t typeclass_id = H5Tget_class(datatype_id);

  if(typeclass_id != H5T_STRING) 
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " is not a std::string"));
  char data2[160];
  CHECK_STATUS(H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, data2));
  data = data2;

  H5Dclose(dataset_id);  
  H5Gclose(group_id);  
}

// ============= //
// serial arrays //
// ============= //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const std::string& DataSetName,
                         const int type, const int Length, 
                         void* data)
{
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hsize_t dimsf = Length;

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), type, 
                            filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, type, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, data);
  CHECK_STATUS(status);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const std::string& DataSetName,
                        const int type, const int Length, 
                        void* data)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  hsize_t dimsf = Length;

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, type, H5S_ALL, filespace_id,
                          H5P_DEFAULT, data);
  CHECK_STATUS(status);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ================== //
// distributed arrays //
// ================== //

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, const std::string& DataSetName,
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

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), type, filespace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Sclose(filespace_id);

  // Each process defines dataset in memory and writes it to the hyperslab
  // in the file.

  hid_t memspace_id = H5Screate_simple(1, &MySize_t, NULL);

  // Select hyperslab in the file.
  filespace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &Offset_t, NULL, &MySize_t, NULL);

  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif

  status = H5Dwrite(dset_id, type, memspace_id, filespace_id,
                    plist_id_, data);
  CHECK_STATUS(status);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
  H5Pclose(plist_id_);
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, const std::string& DataSetName,
                        int MySize, int GlobalSize,
                        const int type, void* data)
{
  if (!IsOpen())
    throw(Exception(__FILE__, __LINE__, "no file open yet"));

  hsize_t MySize_t = MySize;

  // offset
  int itmp;
  Comm_.ScanSum(&MySize, &itmp, 1);
  hsize_t Offset_t = itmp - MySize;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);
  //hid_t space_id = H5Screate_simple(1, &Offset_t, 0);

  // Select hyperslab in the file.
  hid_t filespace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &Offset_t, NULL, 
                      &MySize_t, NULL);

  hid_t mem_dataspace = H5Screate_simple (1, &MySize_t, NULL);

  herr_t status = H5Dread(dataset_id, type, mem_dataspace, filespace_id, 
                          H5P_DEFAULT, data);
  CHECK_STATUS(status);

  H5Sclose(mem_dataspace);
  H5Gclose(group_id);  
  //H5Sclose(space_id);  
  H5Dclose(dataset_id);  
//  H5Dclose(filespace_id);  
}


#endif

