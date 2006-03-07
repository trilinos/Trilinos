#include "Galeri_ConfigDefs.h"
#ifdef HAVE_GALERI_HDF5
#include "Galeri_HDF5.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_FECrsMatrix.h"


herr_t find_dataset(hid_t loc_id, const char *name, void *opdata)
{
  string& token = ((find_dataset_struct*)opdata)->name;
  if (token == name)
    ((find_dataset_struct*)opdata)->found = true;

  return(0);
}

void Galeri::HDF5::Create(const string FileName)
{
  FileName_ = FileName;

#ifdef HAVE_MPI
  // Set up file access property list with parallel I/O access
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  cerr << "Not yet implemented" << endl;
  exit(EXIT_FAILURE);
#endif
  // create the file collectively and release property list identifier.
  file_id = H5Fcreate(FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, 
                      plist_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
}

void Galeri::HDF5::Open(const string FileName, int AccessType)
{
  FileName_ = FileName;

  // create the file collectively and release property list identifier.
  file_id = H5Fopen(FileName.c_str(), AccessType, H5P_DEFAULT);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
}

bool Galeri::HDF5::IsDataSet(string Name)
{
  find_dataset_struct data;
  data.name = Name;
  data.found = false;

  int idx_f = H5Giterate(file_id, "/", NULL, find_dataset, (void*)&data);

  return(data.found);
}

void Galeri::HDF5::Write(const string& Name, int what)
{
  bool flag = IsDataSet(Name);

  filespace = H5Screate(H5S_SCALAR);
  if (flag)
  {
    dset_id = H5Dopen(file_id, Name.c_str());
  }
  else
    dset_id = H5Dcreate(file_id, Name.c_str(), H5T_NATIVE_INT, filespace,
                        H5P_DEFAULT);

  status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace,
                    H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace);
}

void Galeri::HDF5::Write(const string& Name, double what)
{
  bool flag = IsDataSet(Name);

  filespace = H5Screate(H5S_SCALAR);
  if (flag)
  {
    dset_id = H5Dopen(file_id, Name.c_str());
  }
  else
    dset_id = H5Dcreate(file_id, Name.c_str(), H5T_NATIVE_DOUBLE, filespace,
                        H5P_DEFAULT);

  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace,
                    H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace);
}

void Galeri::HDF5::Write(const string& Name, const Epetra_Map& Map)
{
  int MySize = Map.NumMyElements();
  int GlobalSize = Map.NumGlobalElements();
  int* MyGlobalElements = Map.MyGlobalElements();

  string FullName = Name + "-" + toString(Comm_.NumProc());
  Write(FullName, MySize, GlobalSize, H5T_NATIVE_INT, MyGlobalElements);
  FullName = Name + "-" + toString(Comm_.NumProc()) + "-dist";
  Write(FullName, 1, Comm_.NumProc(), H5T_NATIVE_INT, &MySize);
}

void Galeri::HDF5::Write(const string& Name, const Epetra_Vector& x)
{
  const Epetra_BlockMap& OriginalMap = x.Map();
  Epetra_Map LinearMap(OriginalMap.NumGlobalElements(), OriginalMap.IndexBase(), Comm_);
  Epetra_Import Importer(LinearMap, OriginalMap);

  Epetra_Vector LinearVector(LinearMap);
  LinearVector.Import(x, Importer, Insert);

  Write(Name, LinearMap.NumMyElements(), LinearMap.NumGlobalElements(),
        H5T_NATIVE_DOUBLE, x.Values());
}

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

  string FullName = Name + "-ROW";
  Write(FullName, MySize, GlobalSize, H5T_NATIVE_INT, &ROW[0]);
  FullName = Name + "-COL";
  Write(FullName, MySize, GlobalSize, H5T_NATIVE_INT, &COL[0]);
  FullName = Name + "-VAL";
  Write(FullName, MySize, GlobalSize, H5T_NATIVE_DOUBLE, &VAL[0]);
}

void Galeri::HDF5::Write(const string& Name, int MySize, int GlobalSize, int type, const void* data)
{
  int Offset;
  Comm_.ScanSum(&MySize, &Offset, 1);
  Offset -= MySize;

  hsize_t MySize_t = MySize;
  hsize_t GlobalSize_t = GlobalSize;
  hsize_t Offset_t = Offset;

  filespace = H5Screate_simple(1, &GlobalSize_t, NULL); 

  // Create the dataset with default properties and close filespace.
  if (!IsDataSet(Name))
    dset_id = H5Dcreate(file_id, Name.c_str(), type, filespace,
                        H5P_DEFAULT);
  else
    dset_id = H5Dopen(file_id, Name.c_str());

  H5Sclose(filespace);

  // Each process defines dataset in memory and writes it to the hyperslab
  // in the file.

  memspace = H5Screate_simple(1, &MySize_t, NULL);

  // Select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &Offset_t, NULL, &MySize_t, NULL);

  // Create property list for collective dataset write.

  status = H5Dwrite(dset_id, type, memspace, filespace,
                    plist_id, data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
}

// following lines taken from from Roman Geus' FEMAXX code
static void h5_write_param_list_recursive(const Teuchos::ParameterList& params, hid_t group_id) {
    Teuchos::ParameterList::ConstIterator it = params.begin();
    for (; it != params.end(); ++ it) {
        std::string key(params.name(it));
        if (params.isSublist(key)) {
            //
            // Sublist
            //
            
            // Create subgroup for sublist.
            hid_t child_group_id = H5Gcreate(group_id, key.c_str(), 0);
            h5_write_param_list_recursive(params.sublist(key), child_group_id);
            H5Gclose(child_group_id);
        } else {
            //
            // Regular parameter
            //
            
            // Create dataspace/dataset.
            herr_t status;
            hsize_t one = 1;
            hid_t dataspace_id, dataset_id;

            // Write the dataset.
            if (params.isType<int>(key)) {
                int value = params.get<int>(key);
                dataspace_id = H5Screate_simple(1, &one, NULL);
                dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
                status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                                  H5P_DEFAULT, &value);
                status = H5Dclose(dataset_id);
                status = H5Sclose(dataspace_id);
            } else if (params.isType<double>(key)) {
                double value = params.get<double>(key);
                dataspace_id = H5Screate_simple(1, &one, NULL);
                dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
                status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                                  H5P_DEFAULT, &value);
                status = H5Dclose(dataset_id);
                status = H5Sclose(dataspace_id);
            } else if (params.isType<std::string>(key)) {
                std::string value = params.get<std::string>(key);
                hsize_t len = value.size() + 1;
                dataspace_id = H5Screate_simple(1, &len, NULL);
                dataset_id = H5Dcreate(group_id, key.c_str(), H5T_C_S1, dataspace_id, H5P_DEFAULT);
                status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, 
                                  H5P_DEFAULT, value.c_str());
                status = H5Dclose(dataset_id);
                status = H5Sclose(dataspace_id);
            } else if (params.isType<bool>(key)) {
                // Use H5T_NATIVE_USHORT to store a bool value
                unsigned short value = params.get<bool>(key) ? 1 : 0;
                dataspace_id = H5Screate_simple(1, &one, NULL);
                dataset_id = H5Dcreate(group_id, key.c_str(), H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT);
                status = H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, 
                                  H5P_DEFAULT, &value);
                status = H5Dclose(dataset_id);
                status = H5Sclose(dataspace_id);
            } else {
                throw("Unable to export parameter \"%s\" to hdf5 file: Unsupported type.", key.c_str());        
            }
        }
    }
}

void Galeri::HDF5::Write(const string& Name, const Teuchos::ParameterList& params)
{
  group_id = H5Gcreate(file_id, Name.c_str(), 0);

  // Iterate through parameter list 
  h5_write_param_list_recursive(params, group_id);

  // Finalize hdf5 file
  status = H5Gclose(group_id);
}

/**
 * Recursive Operator function called by H5Giterate for each entity in group.
 */
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
            throw("Dimensionality of parameters must be 1.");
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
            throw("Unsupported datatype"); // FIXME
        }
        H5Tclose(type_id);
        H5Sclose(space_id);  
        H5Dclose(dataset_id);  
        break;
    default:
        throw("Unsupported type"); // FIXME
    }
    return 0;
}

void Galeri::HDF5::Read(const string& Name, Teuchos::ParameterList& params) 
{
  // Open hdf5 file
  hid_t       group_id;  /* identifiers */
  herr_t      status;

  // Create group in the root group using absolute name.
  group_id = H5Gopen(file_id, Name.c_str());

  // Iterate through parameter list 
  string FullName = "/" + Name;
  status = H5Giterate(group_id, FullName.c_str() , NULL, f_operator, &params);

  // Finalize hdf5 file
  status = H5Gclose(group_id);
}


void Galeri::HDF5::Write(const string& Name, const Epetra_MultiVector& X)
{
  hid_t       group_id, dset_id;       /* file and dataset identifiers */
  hid_t       filespace, memspace;              /* file and memory dataspace identifiers */
  hid_t	    plist_id;                         /* property list identifier */
  herr_t      status;

  // Redistribute eigenvectors to force linear map
  Epetra_Map target_map(X.GlobalLength(), 0, Comm_);
  Epetra_MultiVector X0(target_map, X.NumVectors());
  X0.Export(X, Epetra_Export(X.Map(), target_map), Add);

  // Create the dataspace for the dataset.
  hsize_t q_dimsf[] = {X0.NumVectors(), X0.GlobalLength()};
  filespace = H5Screate_simple(2, q_dimsf, NULL);

  // Create the dataset with default properties and close filespace.
  dset_id = H5Dcreate(file_id, Name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT);

  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // Select hyperslab in the file.
  hsize_t offset[] = {0, X0.Map().GID(0)};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {X0.NumVectors(), 1};
  hsize_t block[] = {1, X0.MyLength()};
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

  // Each process defines dataset in memory and writes it to the hyperslab in the file.
  hsize_t dimsm[] = {X0.NumVectors() * X0.MyLength()};
  memspace = H5Screate_simple(1, dimsm, NULL);

  // Write hyperslab
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, X0.Values());
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);

  H5Pclose(plist_id);
}

void Galeri::HDF5::Read(const string& Name, Epetra_Map*& Map)
{
    hid_t dataset_id, space_id;
    hsize_t num_eigenpairs, dims[2];
    herr_t status;
    
    string FullName = "/" + Name + "-" + toString(Comm_.NumProc()) + "-dist";
    dataset_id = H5Dopen(file_id, FullName.c_str()); 
    space_id = H5Dget_space(dataset_id);
    status = H5Sget_simple_extent_dims(space_id, dims, NULL);
    vector<int> dist(Comm_.NumProc());
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                     &dist[0]);
    H5Sclose(space_id);  
    H5Dclose(dataset_id);  

    int MySize = dims[0];
    int GlobalSize;
    Comm_.ScanSum(&MySize, &GlobalSize, 1);
    GlobalSize -= MySize;
    hsize_t offset = GlobalSize;


    FullName = "/" + Name + "-" + toString(Comm_.NumProc());

    dataset_id = H5Dopen(file_id, FullName.c_str()); 
    space_id = H5Dget_space(dataset_id);
    if (H5Sget_simple_extent_ndims(space_id) != 1)
        throw("Dimensionality of map must be 1.");

    status = H5Sget_simple_extent_dims(space_id, dims, NULL);
    vector<int> MyGlobalElements(dims[0]);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                     &MyGlobalElements[0]);

    // FIXME: ADD IndexBase
    Map = new Epetra_Map(-1, dims[0], &MyGlobalElements[0], 0, Comm_);
    H5Sclose(space_id);  
    H5Dclose(dataset_id);  
}

void Galeri::HDF5::Read(const string& Name, const Epetra_Map& Map, 
                        Epetra_MultiVector*& X)
{
  hid_t dataset_id, space_id;
  hsize_t num_eigenpairs, dims[2];
  herr_t status;

  string FullName = "/" + Name;
  dataset_id = H5Dopen(file_id, FullName.c_str());
  space_id = H5Dget_space(dataset_id);
  if (H5Sget_simple_extent_ndims(space_id) != 2)
    throw("Dimensionality of modes/eigenvectors must be 2.");
  status = H5Sget_simple_extent_dims(space_id, dims, NULL);

  X = new Epetra_MultiVector(Map, dims[0]);

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   X->Values());
  H5Sclose(space_id);  
  H5Dclose(dataset_id);  
}

void Galeri::HDF5::Read(const string& Name, const Epetra_Map& Map, 
                        Epetra_CrsMatrix*& A)
{
  cout << "as;dlkfjas;lfjaslkfjaklsdfjalks;jdfalks;jdfkla;sjfd;" << endl;
  hid_t dataset_id, space_id;
  hsize_t num_eigenpairs, dims[2];
  herr_t status;

  string FullName = "/" + Name + "-ROW";
  dataset_id = H5Dopen(file_id, FullName.c_str());
  space_id = H5Dget_space(dataset_id);
  if (H5Sget_simple_extent_ndims(space_id) != 1)
    throw("Dimensionality of modes/eigenvectors must be 2.");
  status = H5Sget_simple_extent_dims(space_id, dims, NULL);

  cout << "ROW = " << dims[0] << endl;

  vector<int> ROW(dims[0]);

  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   &ROW[0]);
  H5Sclose(space_id);  
  H5Dclose(dataset_id);  

  // ===
  FullName = "/" + Name + "-COL";
  dataset_id = H5Dopen(file_id, FullName.c_str());
  space_id = H5Dget_space(dataset_id);
  if (H5Sget_simple_extent_ndims(space_id) != 1)
    throw("Dimensionality of modes/eigenvectors must be 2.");
  status = H5Sget_simple_extent_dims(space_id, dims, NULL);

  cout << "COL = " << dims[0] << endl;

  vector<int> COL(dims[0]);

  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   &COL[0]);
  H5Sclose(space_id);  
  H5Dclose(dataset_id);  

  // --- 
  FullName = "/" + Name + "-VAL";
  dataset_id = H5Dopen(file_id, FullName.c_str());
  space_id = H5Dget_space(dataset_id);
  if (H5Sget_simple_extent_ndims(space_id) != 1)
    throw("Dimensionality of modes/eigenvectors must be 2.");
  status = H5Sget_simple_extent_dims(space_id, dims, NULL);

  cout << "VAL = " << dims[0] << endl;

  vector<double> VAL(dims[0]);

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   &VAL[0]);
  H5Sclose(space_id);  
  H5Dclose(dataset_id);  

  Epetra_FECrsMatrix* A2 = new Epetra_FECrsMatrix(Copy, Map, 0); // FIXME: RECT MATRICES

  for (int i = 0; i < dims[0]; ++i)
    A2->InsertGlobalValues(1, &ROW[i], 1, &COL[i], &VAL[i]);

  A2->FillComplete();

  A = A2;
}

#endif
