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
    // @{ \name Constructor and destructor.
    //! ctor
    HDF5(const Epetra_Comm& Comm) :
      Comm_(Comm),
      IsOpen_(false)
    {}

    //! dtor
    ~HDF5() 
    {
      if (IsOpen())
        Close();
    }

    // @}
    // @{ \name Basic operations
    
    //! Creates a new file.
    void Create(const string FileName);

    //! Opens specified file with given access type.
    void Open(const string FileName, int AccessType = H5F_ACC_RDWR);

    //! Closes the file.
    void Close()
    {
      H5Pclose(plist_id_);
      H5Fclose(file_id_);
      IsOpen_ = false;
    }

    //! Returns \c true is a file has already been open using Open()/Create()
    bool IsOpen() const
    {
      return(IsOpen_);
    }

    //! Creates group \c GroupName.
    void CreateGroup(const string& GroupName)
    {
      hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), 0);
      H5Gclose(group_id);
    }

    //! Returns \c true is \c Name is contained in the database.
    bool IsContained(const string Name);

    // @}
    // @{ \name basic non-distributed data types
    
    //! Writes an integer. FIXME
    void Write(const string& Name, int data);

    //! Writes an double. FIXME
    void Write(const string& Name, double data);

    //! Writes an integer in group \c GroupName using intentified \c DataSetName.
    void Write(const string& GroupName, const string& DataSetName, int data);

    //! Writes a double in group \c GroupName using intentified \c DataSetName.
    void Write(const string& GroupName, const string& DataSetName, double data);

    //! Associates string \c Comment with group \c GroupName.
    void Write(const string& GroupName, string Comment)
    {
      H5Gset_comment(file_id_, GroupName.c_str(), Comment.c_str());
    }

    //! Reads the string associated with group \c GroupName.
    void Read(const string& GroupName, string& Comment)
    {
      char comment[128];
      H5Gget_comment(file_id_, GroupName.c_str(), 128, comment);
      Comment = comment;
    }

    // Writes serial array \c data, of length \c Length and type \c type, to \c Name.
    void Write(const string& Name, const int type, const int Length, void* data)
    {
      cout << "TO BE DONE " << endl;
    
    }

    // Reads serial array \c data, of length \c Length and type \c type, from \c Name.
    void Read(const string& Name, const int type, const int Length, void* data);

    //! Reads an integer from group \c /GroupName/DataSetName
    void Read(const string& GroupName, const string& DataSetName, int& data);

    //! Reads a double from group \c /GroupName/DataSetName
    void Read(const string& GroupName, const string& DataSetName, double& data);

    //! Reads serial array \c data, of type \c type, from group \c GroupNameusing dataset name \c DataSetName
    void Read(const string& GroupName, const string& DataSetName,
              const int type, const int Length, void* data);

    //! Writes serial array \c data, of type \c type, to group \c GroupNameusing dataset name \c DataSetName
    void Write(const string& GroupName, const string& DataSetName,
                         const int type, const int Length, 
                         void* data);

    // @}
    // @{ \name Epetra_Map
    
    //! Reads distributed array \c data, of type \c type, from group \c GroupNameusing dataset name \c DataSetName
    void Read(const string& GroupName, const string& DataSetName,
              int MySize, int GlobalSize,
              const int type, void* data);

    //! Writes distributed array \c data, of type \c type, to group \c GroupNameusing dataset name \c DataSetName
    void Write(const string& GroupName, const string& DataSetName, int MySize, int GlobalSize, int type, const void* data);

    // @}
    // @{ \name Epetra_Map

    //! Reads a map from \c GroupName.
    void Read(const string& GroupName, Epetra_Map*& Map);

    //! Writes a distributed Map to group \c GroupName.
    void Write(const string& GroupName, const Epetra_Map& Map);

    // @}
    // @{ \name Epetra_Vector

    //! Reads a vector from group \c GroupName, assumes linear distribution.
    void Read(const string& GroupName, Epetra_Vector*& X);

    //! Reads a vector from group \c GroupName using given map.
    void Read(const string& GroupName, const Epetra_Map& Map, Epetra_Vector*& X);

    //! Writes a distributed vector to group \c GroupName.
    void Write(const string& GroupName, const Epetra_Vector& x);

    // @}
    // @{ \name Epetra_RowMatrix/Epetra_CrsMatrix

    //! Writes a distributed RowMatrix to group \c GroupName.
    void Write(const string& GroupName, const Epetra_RowMatrix& Matrix);

    //! Reads a square matrix from group \c GroupName, assumes linear distribution.
    void Read(const string& GroupName, Epetra_CrsMatrix*& A);

    //! Reads a matrix from group \c GroupName with given range and domain maps.
    void Read(const string& GroupName, 
              const Epetra_Map& DomainMap, 
              const Epetra_Map& RangeMap, 
              Epetra_CrsMatrix*& A);

    // @}
    // @{ \name Teuchos::ParameterList

    //! Writes a parameter list to group \c GroupName.
    void Write(const string& GroupName, const Teuchos::ParameterList& List);

    //! Reads a parameter list from group \c GroupName.
    void Read(const string& GroupName, Teuchos::ParameterList& List);

  private:
    // @{
    // @} \name Private Data

    const Epetra_Comm& Comm() const
    {
      return(Comm_);
    }

    const Epetra_Comm& Comm_;
    string FileName_;

    hid_t       file_id_;
    hid_t	plist_id_;
    herr_t	status;
    bool IsOpen_;
    // @}
};
};

#endif
