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
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Epetra_BlockMap.h"
#include "Epetra_DistObject.h"
#include "EpetraExt_Exception.h"
#include "EpetraExt_Utils.h"
#include "EpetraExt_HDF5.h"
#include "EpetraExt_HDF5_Handle.h"

#define CHECK_HID(hid_t) \
  { if (hid_t < 0) \
    throw(EpetraExt::Exception(__FILE__, __LINE__, \
                    "hid_t is negative")); }

#define CHECK_STATUS(status) \
  { if (status < 0) \
    throw(EpetraExt::Exception(__FILE__, __LINE__, \
                    "function H5Giterater returned a negative value")); }

// ==========================================================================
void EpetraExt::HDF5::Write(const std::string& GroupName, 
                            const Handle& obj)
{
  int NumMyElements = obj.NumMyElements();
  int NumGlobalElements = obj.NumGlobalElements();

  // ===================== //
  // first get global info //
  // ===================== //
  
  std::vector<std::string> IntLabels, DoubleLabels;
  std::vector<int> IntLabelsData;
  std::vector<double> DoubleLabelsData;

  obj.GetLabels(IntLabels, IntLabelsData, DoubleLabels, DoubleLabelsData);

  CreateGroup(GroupName);

  for (unsigned int i = 0; i < IntLabels.size(); ++i)
    Write(GroupName, IntLabels[i], IntLabelsData[i]);

  for (unsigned int i = 0; i < DoubleLabels.size(); ++i)
    Write(GroupName, DoubleLabels[i], DoubleLabelsData[i]);

  // ====================================================== //
  // not compute the storage required by each local element //
  // ====================================================== //
  
  std::vector<int> IntSize(NumMyElements);
  std::vector<int> DoubleSize(NumMyElements);

  int TotalIntSize = 0, TotalDoubleSize = 0;

  std::vector<int> IntData;
  std::vector<double> DoubleData;

  for (int i = 0; i < NumMyElements; ++i)
  {
    IntSize[i] = obj.IntSize(i);
    DoubleSize[i] = obj.DoubleSize(i);
    TotalIntSize += IntSize[i];
    TotalDoubleSize += DoubleSize[i];
  }

  IntData.resize(TotalIntSize);
  DoubleData.resize(TotalDoubleSize);

  int IntCount = 0;
  int DoubleCount = 0;

  // ================================== //
  // pack all data and write it to file //
  // ================================== //

  for (int i = 0; i < NumMyElements; ++i)
  {
    obj.Pack(i, &IntData[IntCount], &DoubleData[DoubleCount]);
    IntCount += IntSize[i];
    DoubleCount += DoubleSize[i];
  }

  if (!IsContained(GroupName))
    CreateGroup(GroupName);
  
  Write(GroupName, "__type__", obj.Type());
  Write(GroupName, "NumGlobalElements", NumGlobalElements);
  Write(GroupName, "has int", obj.HasInt());  
  Write(GroupName, "has double", obj.HasDouble());  

  if (obj.HasInt())
  {
    Write(GroupName, "int ptr", NumMyElements, 
          NumGlobalElements, H5T_NATIVE_INT, &IntSize[0]);
    Write(GroupName, "int data", NumMyElements, 
          NumGlobalElements, H5T_NATIVE_INT, &IntData[0]);
  }

  if (obj.HasDouble())
  {
    Write(GroupName, "double ptr", NumMyElements, 
          NumGlobalElements, H5T_NATIVE_INT, &DoubleSize[0]);
    Write(GroupName, "double data", NumMyElements, 
          NumGlobalElements, H5T_NATIVE_DOUBLE, &DoubleData[0]);
  }
}

// ==========================================================================
void EpetraExt::HDF5::Read(const std::string& GroupName, Handle& obj)
{
  int NumMyElements = obj.NumMyElements();
  int NumGlobalElements = obj.NumGlobalElements();

  std::vector<std::string> IntLabels, DoubleLabels;

  obj.GetLabels(IntLabels, DoubleLabels);
  std::vector<int> IntLabelsData(IntLabels.size());
  std::vector<double> DoubleLabelsData(DoubleLabels.size());

  for (unsigned int i = 0; i < IntLabels.size(); ++i)
    Read(GroupName, IntLabels[i], IntLabelsData[i]);

  for (unsigned int i = 0; i < DoubleLabels.size(); ++i)
    Read(GroupName, DoubleLabels[i], DoubleLabelsData[i]);

  std::vector<int> IntSize(NumMyElements);
  std::vector<int> DoubleSize(NumMyElements);

  int TotalIntSize = 0, TotalDoubleSize = 0;
  int GrandTotalIntSize = 0, GrandTotalDoubleSize = 0;

  // read linear distribution
  if (obj.HasInt())
  {
    Read(GroupName, "int ptr", NumMyElements, 
         NumGlobalElements, H5T_NATIVE_INT, &IntSize[0]);
    for (int i = 0; i < NumMyElements; ++i)
      TotalIntSize += IntSize[i];
    Comm().SumAll(&TotalIntSize, &GrandTotalIntSize, 1);
  }

  if (obj.HasDouble())
  {
    Read(GroupName, "double ptr", NumMyElements, 
         NumGlobalElements, H5T_NATIVE_INT, &DoubleSize[0]);
    for (int i = 0; i < NumMyElements; ++i)
    {
      TotalDoubleSize += DoubleSize[i];
    }
    Comm().SumAll(&TotalDoubleSize, &GrandTotalDoubleSize, 1);
  }

  std::vector<int> IntData(TotalIntSize + 1);
  std::vector<double> DoubleData(TotalDoubleSize + 1);

  // read actual data
  if (obj.HasInt())
    Read(GroupName, "int data", TotalIntSize, 
         GrandTotalIntSize, H5T_NATIVE_INT, &IntData[0]);
  if (obj.HasDouble())
    Read(GroupName, "double data", TotalDoubleSize, 
         GrandTotalDoubleSize, H5T_NATIVE_DOUBLE, &DoubleData[0]);
  
  // now unpack data
  obj.Initialize();

  obj.SetLabels(IntLabelsData, DoubleLabelsData);

  int IntCount = 0, DoubleCount = 0;
  for (int i = 0; i < NumMyElements; ++i)
  {
    obj.UnPack(i, IntSize[i], &(IntData[IntCount]),
               DoubleSize[i], &(DoubleData[DoubleCount]));
    IntCount += IntSize[i];
    DoubleCount += DoubleSize[i];
  }

  obj.Finalize();
}

// ==========================================================================
void EpetraExt::HDF5::ReadHandleProperties(const std::string& GroupName, 
                                           std::string& Type,
                                           int& NumGlobalElements)
{
  if (!IsContained(GroupName))
    throw(Exception(__FILE__, __LINE__,
                    "requested group " + GroupName + " not found"));

  Read(GroupName, "__type__", Type);

  Read(GroupName, "NumGlobalElements", NumGlobalElements);
}
#endif
