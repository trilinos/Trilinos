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
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "EpetraExt_XMLWriter.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

// ============================================================================
EpetraExt::XMLWriter::
XMLWriter(const Epetra_Comm& comm, const std::string& FileName) :
  Comm_(comm),
  FileName_(FileName),
  IsOpen_(false)
{}

// ============================================================================
void EpetraExt::XMLWriter::
Create(const std::string& Label)
{
  if (Comm_.MyPID() == 0) 
  {
    std::ofstream of(FileName_.c_str());
    of << "<ObjectCollection Label=\"" << Label << "\">" << std::endl;
    of.close();
  }

  IsOpen_ = true;
}
  
// ============================================================================
void EpetraExt::XMLWriter:: Close()
{
  if (Comm_.MyPID() == 0) 
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</ObjectCollection>" << std::endl;
    of.close();
  }

  IsOpen_ = false;
}
  
// ============================================================================
void EpetraExt::XMLWriter::
Write(const std::string& Label, const std::vector<std::string>& Content)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  if (Comm_.MyPID()) return;

  std::ofstream of(FileName_.c_str(), std::ios::app);

  of << "<Text Label=\"" << Label << "\">" << std::endl;
  int Csize = (int) Content.size();
  for (int i = 0; i < Csize; ++i)
    of << Content[i] << std::endl;

  of << "</Text>" << std::endl;

  of.close();
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const std::string& Label, const Epetra_RowMatrix& Matrix)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  int Rows = Matrix.NumGlobalRows();
  int Cols = Matrix.NumGlobalRows();
  int Nonzeros = Matrix.NumGlobalNonzeros();

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "<PointMatrix Label=\"" << Label << '"'
      << " Rows=\"" << Rows << '"'
      << " Columns=\"" << Cols<< '"'
      << " Nonzeros=\"" << Nonzeros << '"'
      << " Type=\"double\" StartingIndex=\"0\">" << std::endl;
  }

  int Length = Matrix.MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);

  for (int iproc = 0; iproc < Comm_.NumProc(); iproc++)
  {
    if (iproc == Comm_.MyPID())
    {
      std::ofstream of(FileName_.c_str(), std::ios::app);
      of.precision(15);

      for (int i = 0; i < Matrix.NumMyRows(); ++i)
      {
        int NumMyEntries;
        Matrix.ExtractMyRowCopy(i, Length, NumMyEntries, &Values[0], &Indices[0]);

        int GRID = Matrix.RowMatrixRowMap().GID(i);

        for (int j = 0; j < NumMyEntries; ++j)
          of << GRID << " " << Matrix.RowMatrixColMap().GID(Indices[j])
             << " " << std::setiosflags(std::ios::scientific) << Values[j] << std::endl;
      }
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</PointMatrix>" << std::endl;
    of.close();
  }
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const std::string& Label, const Epetra_MultiVector& MultiVector)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  int Length = MultiVector.GlobalLength();
  int NumVectors = MultiVector.NumVectors();

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);

    of << "<MultiVector Label=\"" << Label 
      << "\" Length=\"" << Length << '"'
      << " NumVectors=\"" << NumVectors << '"'
      << " Type=\"double\">" << std::endl;
  }


  for (int iproc = 0; iproc < Comm_.NumProc(); iproc++)
  {
    if (iproc == Comm_.MyPID())
    {
      std::ofstream of(FileName_.c_str(), std::ios::app);

      of.precision(15);
      for (int i = 0; i < MultiVector.MyLength(); ++i)
      {
        for (int j = 0; j < NumVectors; ++j)
          of << std::setiosflags(std::ios::scientific) << MultiVector[j][i] << " ";
        of << std::endl;
      }
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</MultiVector>" << std::endl;
    of.close();
  }
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const std::string& Label, const Epetra_Map& Map)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  int NumGlobalElements = Map.NumGlobalElements();
  int* MyGlobalElements = Map.MyGlobalElements();

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);

    of << "<Map Label=\"" << Label 
      << "\" NumElements=\"" << NumGlobalElements << '"'
      << " IndexBase=\"" << Map.IndexBase() << '"'
      << " NumProc=\"" << Comm_.NumProc() << '"';

    of.close();
  }

  for (int iproc = 0; iproc < Comm_.NumProc(); ++iproc)
  {
    if (iproc == Comm_.MyPID())
    {
      std::ofstream of(FileName_.c_str(), std::ios::app);

      of << " ElementsOnProc" << iproc << "=\"" << Map.NumMyElements() << '"';
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << '>' << std::endl;
    of.close();
  }

  for (int iproc = 0; iproc < Comm_.NumProc(); iproc++)
  {
    if (iproc == Comm_.MyPID())
    {
      std::ofstream of(FileName_.c_str(), std::ios::app);

      of << "<Proc ID=\"" << Comm_.MyPID() << "\">" << std::endl;

      for (int i = 0; i < Map.NumMyElements(); ++i)
      {
        of << MyGlobalElements[i] << std::endl;
      }

      of << "</Proc>" << std::endl;
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</Map>" << std::endl;
    of.close();
  }
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const std::string& Label, Teuchos::ParameterList& List)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  if (Comm_.MyPID()) return;

  std::ofstream of(FileName_.c_str(), std::ios::app);

  of << "<List Label=\"" << Label << "\">" << std::endl;

  XMLParameterListWriter Writer;
  XMLObject Obj = Writer.toXML(List);

  of << Obj.toString();

  of << "</List>" << std::endl;

  of.close();
}
