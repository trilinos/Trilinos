//@HEADER
/*
************************************************************************

              EpetraExt: Linear Algebra Services Package
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov)

************************************************************************
*/
//@HEADER

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
    virtual std::string Type() const = 0;

    virtual bool HasInt() const = 0;

    virtual bool HasDouble() const = 0;

    //! Returns the size of integer data for local element \c EID.
    virtual int IntSize(const int EID) const = 0;

    //! Returns the size of double data for local element \c EID.
    virtual int DoubleSize(const int EID) const = 0;

    //! Packs all global information.
    virtual int GetLabels(std::vector<std::string>& IntLabels, 
                  std::vector<std::string>& DoubleLabels) const = 0;

    //! Packs all global information.
    virtual int GetLabels(std::vector<std::string>& IntLabels, std::vector<int>& IntLabelsData,
                  std::vector<std::string>& DoubleLabels, std::vector<double>& DoubleLabelsData) const = 0;

    //! Sets global information
    virtual int SetLabels(const std::vector<int>& IntLabelsData,
                  const std::vector<double>& DoubleLabelsData) = 0;

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
    std::string Type() const
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
    int GetLabels(std::vector<std::string>& IntLabels, 
                  std::vector<std::string>& DoubleLabels) const
    {
      IntLabels.resize(0);
      DoubleLabels.resize(0);

      return(0);
    }

    //! Packs all global information.
    int GetLabels(std::vector<std::string>& IntLabels, std::vector<int>& IntLabelsData,
                  std::vector<std::string>& DoubleLabels, std::vector<double>& DoubleLabelsData) const
    {
      return(0);
    }

    //! Sets global information
    int SetLabels(const std::vector<int>& IntLabelsData,
                  const std::vector<double>& DoubleLabelsData)
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
#endif /* EPETRAEXT_HDF5_HANDLE_H */
