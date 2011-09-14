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
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_XMLReader.h"
#include "EpetraExt_XMLWriter.h"
#include "Teuchos_ParameterList.hpp"

// Showing the usage of XML I/O.
// This example can be run with any number of processors.
//
// \author Marzio Sala, D-INFK/ETHZ.
//
// \date Last modified on 10-May-06.

int main (int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // define some Epetra objects

  int n = Comm.NumProc() * 4;
  Epetra_Map Map(n, 0, Comm);
  Epetra_MultiVector x(Map, 2); x.Random();
  Epetra_MultiVector b(Map, 2); x.Random();
  Epetra_CrsMatrix Matrix(Copy, Map, 0);
  // diagonal matrix
  for (int i = 0; i < Map.NumMyElements(); ++i)
  {
    int ii = Map.GID(i);
    double one = 1.0;
    Matrix.InsertGlobalValues(ii, 1, &one, &ii);
  }
  Matrix.FillComplete();

  Teuchos::ParameterList List;
  List.set("int parameter", 10);
  List.set("double parameter", 10.0);
  List.set("std::string parameter", "std::string");

  // ========================= //
  // Part I: generate XML file //
  // ========================= //
  
  EpetraExt::XMLWriter XMLWriter(Comm, "data.xml");

  std::vector<std::string> Content;
  Content.push_back("This is an example of description");
  Content.push_back("The description is as long as desired,");
  Content.push_back("just put it in a std::vector of strings.");

  XMLWriter.Create("MyProblem");
  XMLWriter.Write("Author", "myself and others");
  XMLWriter.Write("Date", "May 2006");
  XMLWriter.Write("MyMap", Map);
  XMLWriter.Write("MyMatrix", Matrix);
  XMLWriter.Write("MyLHS", x);
  XMLWriter.Write("MyRHS", b);
  XMLWriter.Write("MyContent", Content);
  XMLWriter.Write("MyParameters", List);
  XMLWriter.Close();

  // ================== //
  // Part II: read data //
  // ================== //
  
  EpetraExt::XMLReader XMLReader(Comm, "data.xml");

  Epetra_Map* MyMap;
  Epetra_CrsMatrix* MyMatrix;
  Epetra_MultiVector* MyLHS;
  Epetra_MultiVector* MyRHS;
  Teuchos::ParameterList MyParameters;
  std::vector<std::string> Author;
  std::vector<std::string> Date;
  std::vector<std::string> MyContent;

  XMLReader.Read("Author", Author);
  XMLReader.Read("Date", Date);
  XMLReader.Read("MyMap", MyMap);
  XMLReader.Read("MyMatrix", MyMatrix);
  XMLReader.Read("MyLHS", MyLHS);
  XMLReader.Read("MyRHS", MyRHS);
  XMLReader.Read("MyContent", MyContent);
  XMLReader.Read("MyParameters", MyParameters);

  std::cout << *MyMap;
  std::cout << *MyMatrix;
  std::cout << *MyLHS;
  std::cout << *MyRHS;
  if (Comm.MyPID() == 0)
  {
    int Msize = (int) MyContent.size();
    for (int i = 0; i < Msize; ++i)
      std::cout << MyContent[i] << std::endl;

    std::cout << MyParameters;
    std::cout << "Author = " << Author[0] << std::endl;
    std::cout << "Date   = " << Date[0] << std::endl;
  }

  delete MyMap;
  delete MyMatrix;
  delete MyLHS;
  delete MyRHS;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
