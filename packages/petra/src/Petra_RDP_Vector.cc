
#include "Petra_RDP_Vector.h"

//=============================================================================
Petra_RDP_Vector::Petra_RDP_Vector(const Petra_BlockMap& Map)
  : Petra_RDP_MultiVector(Map,1) // Vector is just special case of MultiVector
{
}
//=============================================================================
Petra_RDP_Vector::Petra_RDP_Vector(const Petra_RDP_Vector& Source)
  : Petra_RDP_MultiVector(Source) // Vector is just special case of MultiVector
{
}
//=============================================================================
Petra_RDP_Vector::Petra_RDP_Vector(Petra_DataAccess CV, const Petra_BlockMap& Map, double *V)
  : Petra_RDP_MultiVector(CV, Map, V, Map.NumMyEquations(), 1) // Vector is just special case of MultiVector
{
}
//=============================================================================
Petra_RDP_Vector::Petra_RDP_Vector(Petra_DataAccess CV, const Petra_RDP_Vector& Source, int Index)
  : Petra_RDP_MultiVector(CV, Source, Index, 1) // Vector is just special case of MultiVector
{
}
//=========================================================================
Petra_RDP_Vector::~Petra_RDP_Vector(){}

//=============================================================================
int Petra_RDP_Vector::ExtractCopy(double *V)
{
  return(Petra_RDP_MultiVector::ExtractCopy(V, 1));
}

//=============================================================================
int Petra_RDP_Vector::ExtractView(double **V)
{
  int junk;
  return(Petra_RDP_MultiVector::ExtractView(V, &junk));
}

//=========================================================================
double& Petra_RDP_Vector::operator [] (int Index)  {

   return(Values_[Index - IndexBase_]);
}

//=========================================================================
const double& Petra_RDP_Vector::operator [] (int Index) const  {

   return(Values_[Index - IndexBase_]);
}
// Non-member functions

ostream& operator << (ostream& os, const Petra_RDP_Vector& V)
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
