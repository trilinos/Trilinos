
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_BlockMap& Map)
  : Epetra_MultiVector(Map,1) // Vector is just special case of MultiVector
{
}
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_Vector& Source)
  : Epetra_MultiVector(Source) // Vector is just special case of MultiVector
{
}
//=============================================================================
Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *V)
  : Epetra_MultiVector(CV, Map, V, Map.NumMyEquations(), 1) // Vector is just special case of MultiVector
{
}
//=============================================================================
Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_Vector& Source, int Index)
  : Epetra_MultiVector(CV, Source, Index, 1) // Vector is just special case of MultiVector
{
}
//=========================================================================
Epetra_Vector::~Epetra_Vector(){}

//=============================================================================
int Epetra_Vector::ExtractCopy(double *V)
{
  return(Epetra_MultiVector::ExtractCopy(V, 1));
}

//=============================================================================
int Epetra_Vector::ExtractView(double **V)
{
  int junk;
  return(Epetra_MultiVector::ExtractView(V, &junk));
}

//=========================================================================
double& Epetra_Vector::operator [] (int Index)  {

   return(Values_[Index - IndexBase_]);
}

//=========================================================================
const double& Epetra_Vector::operator [] (int Index) const  {

   return(Values_[Index - IndexBase_]);
}
// Non-member functions

ostream& operator << (ostream& os, const Epetra_Vector& V)
{
  int MyPID = V.Map().Comm().MyPID();
  int NumProc = V.Map().Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int MyLength = V.MyLength();
      double * V_Values = V.Values();
      int * MyGlobalElements = V.Map().MyGlobalElements();
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);

      if (MyPID==0) {
	os.width(14);
	os <<  "     MyPID"; os << "    ";
	os.width(14);
	os <<  "      Global Index "; os << " ";
	os.width(20);
	os <<  "Value  ";
	os << endl;
      }
      
      for (int i=0; i < MyLength; i++)
	{
	  os.width(14);
	  os <<  MyPID; os << "    ";
	  os.width(14);
	  os <<  MyGlobalElements[i]; os << "    ";
	  os.width(20);
	  os <<  V_Values[i];
	  os << endl;
	}

      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }

    V.Map().Comm().Barrier(); // Do a few barriers to allow I/O to complete
    V.Map().Comm().Barrier();
    V.Map().Comm().Barrier();
  }

  return os;
}
