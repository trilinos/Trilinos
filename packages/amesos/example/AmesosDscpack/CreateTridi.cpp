//
//  This code populates an empty EpetraCrsMatrix with a tridiagonal with 
//  -1 on the off-diagonals and 2 on the diagonal.  
//
//  This code was plaguerized from epetra/example/petra_power_method/cxx_main.cpp
//  presumably written by Mike Heroux.
//
//  Adapted by Ken Stanley - Aug 2003 
//
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

int CreateTridi(Epetra_CrsMatrix& A)
{

  Epetra_Map Map = A.RowMap();
  int NumMyElements = Map.NumMyElements();
  int NumGlobalElements = Map.NumGlobalElements();

  int * MyGlobalElements = new int[NumMyElements];
    Map.MyGlobalElements(MyGlobalElements);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[3];
  int *Indices = new int[3];
  double two = 2.0;
  int NumEntries;
  
  for (int i=0; i<NumMyElements; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 0;
	Indices[1] = 1;
	Values[0] = 2.0;
	Values[1] = -1.0;
	NumEntries = 2;
      }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
      {
	Indices[0] = NumGlobalElements-1;
	Indices[1] = NumGlobalElements-2;
	Values[0] = 2.0;
	Values[1] = -1.0;
	NumEntries = 2;
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i];
	Indices[2] = MyGlobalElements[i]+1;
	Values[0] = -1.0; 
	Values[1] = 2.0;
	Values[2] = -1.0;
	NumEntries = 3;
      }
    
    cout << " NumEntries = " << endl ; 
    cout << " NumEntries = " << NumEntries << " Indices[] = " << endl ; 
    cout << " NumEntries = " << NumEntries << " Indices[] = " 
	 << Indices[0] << endl ; 
    cout << " NumEntries = " << NumEntries << " Indices[] = " 
	 << Indices[0] << " " << Indices[1] << endl ; 
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
     // Put in the diagonal entry
    cout << " i = " <<  i << "MyGlobalElements[" << i << "]=" << MyGlobalElements[i] << endl ; 
     //     assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, &MyGlobalElements[i])==0);
    }
  
  // Finish up
  assert(A.TransformToLocal()==0);


  delete MyGlobalElements;
  delete Values;
  delete Indices;
  return 0;
}

