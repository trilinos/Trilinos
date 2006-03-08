#include "Galeri_ConfigDefs.h"
#ifdef HAVE_GALERI_HDF5

#include "hdf5.h"
#include <iostream>
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <vector>
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Galeri_Utils.h"
#include "Teuchos_ParameterList.hpp"

namespace Galeri 
{
class HDF5 
{
  public: 
    HDF5(const Epetra_Comm& Comm) :
      Comm_(Comm)
    {}

    ~HDF5() {}

    void Create(const string FileName);

    void Open(const string FileName, int AccessType = H5F_ACC_RDWR);

    bool IsDataSet(const string Name);

    // OLD
    void Write(const string& Name, int what);

    void Write(const string& GroupName, const string& DataSetName, int what);

    void Write(const string& GroupName, const string& DataSetName, double what);

    void Write(const string& Name, double what);

    void Write(const string& Name, string Comment)
    {
      H5Gset_comment(file_id, Name.c_str(), Comment.c_str());
    }

    void Write(const string& Name, const int type, const int Length, void* data);

    void Write(const string& Name, const Epetra_Map& Map);

    void Write(const string& Name, const Epetra_Vector& x);

    void Write(const string& Name, const Epetra_RowMatrix& Matrix);

    void Write(const string& Name, int MySize, int GlobalSize, int type, const void* data);

    void Write(const string& GroupName, const string& DataSetName, int MySize, int GlobalSize, int type, const void* data);
    void Write(const string& Name, const Teuchos::ParameterList& List);

    void Read(const string& GroupName, const string& DataSetName, int& data);

    void Read(const string& GroupName, const string& DataSetName, double& data);

    void Read(const string& Name, Teuchos::ParameterList& List);

    void Read(const string& Name, Epetra_Map*& Map);

    void Read(const string& Name, Epetra_Vector*& X);

    void Read(const string& Name, const Epetra_Map& Map, Epetra_Vector*& X);

    void Read(const string& Name, Epetra_CrsMatrix*& A);

    void Read(const string& Name, const Epetra_Map& DomainMap, 
              const Epetra_Map& RangeMap, Epetra_CrsMatrix*& A);

    void Read(const string& GroupName, const string& DataSetName,
              int MySize, int GlobalSize,
              const int type, void* data);

    void Close()
    {
      H5Pclose(plist_id);
      H5Fclose(file_id);
    }

    void CreateGroup(const string& Name)
    {
      group_id = H5Gcreate(file_id, Name.c_str(), 0);
      H5Gclose(group_id);
    }

    void Read(const string& Name, string& Comment)
    {
      char comment[128];
      H5Gget_comment(file_id, Name.c_str(), 128, comment);
      Comment = comment;
    }

  private:
    const Epetra_Comm& Comm_;
    string FileName_;

    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t	count;	          /* hyperslab selection parameters */
    hsize_t	offset;
    hid_t	plist_id;                 /* property list identifier */
    int         i;
    herr_t	status;
    hid_t       group_id;
};
};

#endif
