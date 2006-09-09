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
  List.set("string parameter", "string");

  // ========================= //
  // Part I: generate XML file //
  // ========================= //
  
  EpetraExt::XMLWriter XMLWriter(Comm, "data.xml");

  vector<string> Content;
  Content.push_back("This is an example of description");
  Content.push_back("The description is as long as desired,");
  Content.push_back("just put it in a vector of strings.");

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
  vector<string> Author;
  vector<string> Date;
  vector<string> MyContent;

  XMLReader.Read("Author", Author);
  XMLReader.Read("Date", Date);
  XMLReader.Read("MyMap", MyMap);
  XMLReader.Read("MyMatrix", MyMatrix);
  XMLReader.Read("MyLHS", MyLHS);
  XMLReader.Read("MyRHS", MyRHS);
  XMLReader.Read("MyContent", MyContent);
  XMLReader.Read("MyParameters", MyParameters);

  cout << *MyMap;
  cout << *MyMatrix;
  cout << *MyLHS;
  cout << *MyRHS;
  if (Comm.MyPID() == 0)
  {
    int Msize = (int) MyContent.size();
    for (int i = 0; i < Msize; ++i)
      cout << MyContent[i] << endl;

    cout << MyParameters;
    cout << "Author = " << Author[0] << endl;
    cout << "Date   = " << Date[0] << endl;
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
