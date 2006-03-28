#ifndef EPETRAEXT_HDF5_HANDLE_H
#define EPETRAEXT_HDF5_HANDLE_H
#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_EPETRAEXT_HDF5

namespace EpetraExt {

class Handle 
{
  public:
    virtual ~Handle() {}

    //! Returns the local number of elements.
    virtual int NumMyElements() const = 0;

    //! Returns the global number of elements.
    virtual int NumGlobalElements() const = 0;

    //! Returns the identifier of the distributed object.
    virtual string Type() const = 0;

    virtual bool HasInt() const = 0;

    virtual bool HasDouble() const = 0;

    //! Returns the size of integer data for local element \c EID.
    virtual int IntSize(const int EID) const = 0;

    //! Returns the size of double data for local element \c EID.
    virtual int DoubleSize(const int EID) const = 0;

    //! Packs all global information.
    virtual int GetLabels(vector<string>& IntLabels, 
                  vector<string>& DoubleLabels) const = 0;

    //! Packs all global information.
    virtual int GetLabels(vector<string>& IntLabels, vector<int>& IntLabelsData,
                  vector<string>& DoubleLabels, vector<double>& DoubleLabelsData) const = 0;

    //! Sets global information
    virtual int SetLabels(const vector<int>& IntLabelsData,
                  const vector<double>& DoubleLabelsData) = 0;

    //! Packs all data for local element \c EID in the specified arrays.
    virtual int Pack(const int EID, int* IntData, double* DoubleData) const = 0;

    //! Unpacks all data for local element \c EID in the specified arrays.
    virtual int UnPack(const int EID, int IntSize, int* IntData, 
               int DoubleSize, double* DoubleData) = 0;

    //! Performs any initialization procedure \e before unpacking.
    virtual int Initialize() = 0;

    //! Performs any finalization procedure \e after unpacking.
    virtual int Finalize() = 0;
};

#include "Epetra_Vector.h"

class Epetra_Vector_Handle : public Handle
{
  public:
    Epetra_Vector_Handle(Epetra_Vector& obj) :
      obj_(&obj)
    {}

    //! Returns the local number of elements.
    int NumMyElements() const
    {
      return(obj_->Map().NumMyElements());
    }

    //! Returns the global number of elements.
    int NumGlobalElements() const
    {
      return(obj_->Map().NumGlobalElements());
    }

    //! Returns the identifier of the distributed object.
    string Type() const
    {
      return("Handle<Epetra_Vector>");
    }

    bool HasInt() const
    {
      return(false);
    }

    bool HasDouble() const
    {
      return(true);
    }

    //! Returns the size of integer data for local element \c EID.
    int IntSize(const int EID) const
    {
      return(0);
    }

    //! Returns the size of double data for local element \c EID.
    int DoubleSize(const int EID) const
    {
      return(1);
    }

    //! Packs all global information.
    int GetLabels(vector<string>& IntLabels, 
                  vector<string>& DoubleLabels) const
    {
      IntLabels.resize(0);
      DoubleLabels.resize(0);

      return(0);
    }

    //! Packs all global information.
    int GetLabels(vector<string>& IntLabels, vector<int>& IntLabelsData,
                  vector<string>& DoubleLabels, vector<double>& DoubleLabelsData) const
    {
      return(0);
    }

    //! Sets global information
    int SetLabels(const vector<int>& IntLabelsData,
                  const vector<double>& DoubleLabelsData)
    {
      return(0);
    }

    //! Packs all data for local element \c EID in the specified arrays.
    int Pack(const int EID, int* IntData, double* DoubleData) const
    {
      DoubleData[0] = (*obj_)[EID];

      return(0);
    }

    //! Unpacks all data for local element \c EID in the specified arrays.
    int UnPack(const int EID, int IntSize, int* IntData, 
               int DoubleSize, double* DoubleData)
    {
      (*obj_)[EID] = DoubleData[0];

      return(0);
    }

    //! Performs any initialization procedure \e before unpacking.
    int Initialize()
    {
      // do nothing here
      return(0);
    }

    //! Performs any finalization procedure \e after unpacking.
    int Finalize()
    {
      // do nothing here
      return(0);
    }

  private:
    Epetra_Vector* obj_;
};

} // namespace EpetraExt
#endif
#endif
