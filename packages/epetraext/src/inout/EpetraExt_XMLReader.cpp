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
#include "EpetraExt_XMLReader.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

// ============================================================================
static void Tokenize(const std::string& str, std::vector<std::string>& tokens,
              const std::string& delimiters = " ")
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    // Found a token, add it to the std::vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}
using namespace Teuchos;

// ============================================================================
EpetraExt::XMLReader::XMLReader(const Epetra_Comm& comm, const std::string& FileName) :
  Comm_(comm)
{
#ifdef HAVE_TEUCHOS_EXPAT
  FileInputSource fileSrc(FileName);
  fileXML_ = rcp(new XMLObject(fileSrc.getObject()));
  IsOpen_ = true;
#else
  std::cerr << "Teuchos was not configured with support for expat." << std::endl;
  std::cerr << "Please reconfigure teuchos with --enable-teuchos-expat." << std::endl;
  exit(EXIT_FAILURE);
#endif
}

// ============================================================================
void EpetraExt::XMLReader::
Read(const std::string& Label, Epetra_CrsGraph*& Graph)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  Graph = 0;

  for (int i = 0; i < fileXML_->numChildren(); ++i)
  {
    const XMLObject& child = fileXML_->getChild(i);
    std::string tag = child.getTag();

    if (tag == "Graph")
    {
      if (child.hasAttribute("Label") && child.getRequired("Label") == Label)
      {
	bool debug = false;
        int NumGlobalRows = child.getRequiredInt("Rows");
        int NumGlobalCols = child.getRequiredInt("Columns");
        int NumGlobalEntries = child.getRequiredInt("Entries");
        int Offset = child.getRequiredInt("StartingIndex");
	if (debug) std::cout << NumGlobalCols << NumGlobalEntries << Offset << std::endl;

        Epetra_Map map(NumGlobalRows, 0, Comm_);
        Graph = new Epetra_CrsGraph(Copy, map, 0);

        for (int j = 0; j < child.numContentLines(); ++j)
        {
          std::vector<std::string> tokens;
          const std::string& line = child.getContentLine(j);
          Tokenize(line, tokens, " \n\r\t");
          if (tokens.size() < 2) continue;

          int row, col;
          row = atoi((char*)tokens[0].c_str());
          col = atoi((char*)tokens[1].c_str());

          if (map.LID(row) != -1)
            Graph->InsertGlobalIndices(row, 1, &col);
        }
        Graph->FillComplete();
      }
    }
  }
}

// ============================================================================
void EpetraExt::XMLReader::
Read(const std::string& Label, Epetra_CrsMatrix*& matrix)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  matrix = 0;

  for (int i = 0; i < fileXML_->numChildren(); ++i)
  {
    const XMLObject& child = fileXML_->getChild(i);
    std::string tag = child.getTag();

    if (tag == "PointMatrix")
    {
      if (child.hasAttribute("Label") && child.getRequired("Label") == Label)
      {
	bool debug = false;
        int NumGlobalRows = child.getRequiredInt("Rows");
        int NumGlobalCols = child.getRequiredInt("Columns");
        int NumGlobalNonzeros = child.getRequiredInt("Nonzeros");
        int Offset = child.getRequiredInt("StartingIndex");
	if (debug) std::cout << NumGlobalCols << NumGlobalNonzeros << Offset << std::endl;

        Epetra_Map map(NumGlobalRows, 0, Comm_);
        matrix = new Epetra_CrsMatrix(Copy, map, 0);

        for (int j = 0; j < child.numContentLines(); ++j)
        {
          std::vector<std::string> tokens;
          const std::string& line = child.getContentLine(j);
          Tokenize(line, tokens, " \n\r\t");
          if (tokens.size() < 3) continue;

          int row, col;
          double val;
          row = atoi((char*)tokens[0].c_str());
          col = atoi((char*)tokens[1].c_str());
          sscanf((char*)tokens[2].c_str(), "%lg", &val);
          //val = atof((char*)tokens[2].c_str());

          if (map.LID(row) != -1)
            matrix->InsertGlobalValues(row, 1, &val, &col);
        }
        matrix->FillComplete();
      }
    }
  }
}

// ============================================================================
void EpetraExt::XMLReader::
Read(const std::string& Label, Epetra_MultiVector*& MultiVector)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  MultiVector = 0;

  // read all file and create all objects in memory.
  for (int i = 0; i < fileXML_->numChildren(); ++i)
  {
    const XMLObject& child = fileXML_->getChild(i);
    std::string tag = child.getTag();

    if (tag == "MultiVector")
    {
      if (child.hasAttribute("Label") && child.getRequired("Label") == Label)
      {
        int GlobalLength = child.getRequiredInt("Length");
        int NumVectors = child.getRequiredInt("NumVectors");

        Epetra_Map Map(GlobalLength, 0, Comm_);
        MultiVector = new Epetra_MultiVector(Map, NumVectors);

        int count = 0;
        double val;
        for (int j = 0; j < child.numContentLines(); ++j)
        {
          std::vector<std::string> tokens;

          const std::string& line = child.getContentLine(j);

          Tokenize(line, tokens, " \n\r\t");

          if (tokens.size() == 0) continue;

          TEST_FOR_EXCEPTION(tokens.size() != (unsigned) NumVectors, std::logic_error,
                             "wrong number of tokens in line; "
                             << "tokens.size() = " << tokens.size() 
                             << ", NumVectors = " << NumVectors);
	  int tsize = (int) tokens.size();
          for (int k = 0; k < tsize; ++k)
          {
            if (Map.LID(count) != -1)
            {
              sscanf((char*)(tokens[k].c_str()), "%lf", &val);

              (*MultiVector)[k][Map.LID(count)] = val;
            }
          }
          ++count;
        }
      }
    }
  }
}

// ============================================================================
void EpetraExt::XMLReader::
Read(const std::string& Label, Epetra_Map*& Map)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  Map = 0;

  // read all file and create all objects in memory.
  for (int i = 0; i < fileXML_->numChildren(); ++i)
  {
    const XMLObject& child = fileXML_->getChild(i);
    std::string tag = child.getTag();

    if (tag == "Map")
    {
      if (child.hasAttribute("Label") && child.getRequired("Label") == Label)
      {
        int NumGlobalElements = child.getRequiredInt("NumElements");
        int IndexBase = child.getRequiredInt("IndexBase");
        int NumProc = child.getRequiredInt("NumProc");

        TEST_FOR_EXCEPTION(NumProc != Comm_.NumProc(), std::logic_error,
                           "Requested map defined with different number of processors, "
                           << "NumProc = " << NumProc << " while "
                           << "Comm.NumProc() = " << Comm_.NumProc());

        char str[80];
        sprintf(str, "ElementsOnProc%d", Comm_.MyPID());
        int NumMyElements = child.getRequiredInt(str);

        sprintf(str, "ElementsOnProc%d", Comm_.MyPID());

        std::vector<int> MyGlobalElements(NumMyElements);

        for (int iproc = 0; iproc < child.numChildren(); ++iproc)
        {
          const XMLObject& newChild = child.getChild(iproc);
          int count = 0;

          if (newChild.hasAttribute("ID") && 
              newChild.getRequiredInt("ID") == Comm_.MyPID())
          {
            for (int j = 0; j < newChild.numContentLines(); ++j)
            {
              std::vector<std::string> tokens;

              const std::string& line = newChild.getContentLine(j);

              Tokenize(line, tokens, " \n\r\t");
	      int tsize = (int) tokens.size();
              for (int k = 0; k < tsize; ++k)
              {
                MyGlobalElements[count++] = atoi((char*)tokens[k].c_str());
              }
            }
          }
        }

        Map = new Epetra_Map(NumGlobalElements, NumMyElements,
                             &MyGlobalElements[0], IndexBase, Comm_);
      }
    }
  }
}

// ============================================================================
void EpetraExt::XMLReader::
Read(const std::string& Label, std::vector<std::string>& Content)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  for (int i = 0; i < fileXML_->numChildren(); ++i)
  {
    const XMLObject& child = fileXML_->getChild(i);
    std::string tag = child.getTag();

    if (tag == "Text")
    {
      if (child.hasAttribute("Label") && child.getRequired("Label") == Label)
      {
        for (int j = 0; j < child.numContentLines(); ++j)
        {
          const std::string& line = child.getContentLine(j);
          if (line == "\n") continue;
          Content.push_back(line);
        }
      }
    }
  }
}

// ============================================================================
void EpetraExt::XMLReader::
Read(const std::string& Label, Teuchos::ParameterList& List)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  for (int i = 0; i < fileXML_->numChildren(); ++i)
  {
    const XMLObject& child = fileXML_->getChild(i);
    std::string tag = child.getTag();

    if (tag == "List")
    {
      if (child.hasAttribute("Label") && child.getRequired("Label") == Label)
      {
        Teuchos::XMLParameterListReader ListReader;
        List = ListReader.toParameterList(child.getChild(0));
      }
    }
  }
}
