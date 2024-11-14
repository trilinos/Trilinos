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

#ifndef EPETRAEXT_HDF5_HANDLE_H
#define EPETRAEXT_HDF5_HANDLE_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif
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
