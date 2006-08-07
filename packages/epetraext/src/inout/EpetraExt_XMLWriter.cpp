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
XMLWriter(const Epetra_Comm& comm, const string& FileName) :
  Comm_(comm),
  FileName_(FileName),
  IsOpen_(false)
{}

// ============================================================================
void EpetraExt::XMLWriter::
Create(const string& Label)
{
  if (Comm_.MyPID() == 0) 
  {
    std::ofstream of(FileName_.c_str());
    of << "<ObjectCollection Label=\"" << Label << "\">" << endl;
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
    of << "</ObjectCollection>" << endl;
    of.close();
  }

  IsOpen_ = false;
}
  
// ============================================================================
void EpetraExt::XMLWriter::
Write(const string& Label, const vector<string>& Content)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  if (Comm_.MyPID()) return;

  std::ofstream of(FileName_.c_str(), std::ios::app);

  of << "<Text Label=\"" << Label << "\">" << endl;

  for (int i = 0; i < Content.size(); ++i)
    of << Content[i] << endl;

  of << "</Text>" << endl;

  of.close();
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const string& Label, const Epetra_RowMatrix& Matrix)
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
      << " Type=\"double\" StartingIndex=\"0\">" << endl;
  }

  int Length = Matrix.MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);

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
             << " " << setiosflags(ios::scientific) << Values[j] << endl;
      }
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</PointMatrix>" << endl;
    of.close();
  }
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const string& Label, const Epetra_MultiVector& MultiVector)
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
      << " Type=\"double\">" << endl;
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
          of << setiosflags(ios::scientific) << MultiVector[j][i] << " ";
        of << endl;
      }
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</MultiVector>" << endl;
    of.close();
  }
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const string& Label, const Epetra_Map& Map)
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
    of << '>' << endl;
    of.close();
  }

  for (int iproc = 0; iproc < Comm_.NumProc(); iproc++)
  {
    if (iproc == Comm_.MyPID())
    {
      std::ofstream of(FileName_.c_str(), std::ios::app);

      of << "<Proc ID=\"" << Comm_.MyPID() << "\">" << endl;

      for (int i = 0; i < Map.NumMyElements(); ++i)
      {
        of << MyGlobalElements[i] << endl;
      }

      of << "</Proc>" << endl;
      of.close();
    }
    Comm_.Barrier();
  }

  if (Comm_.MyPID() == 0)
  {
    std::ofstream of(FileName_.c_str(), std::ios::app);
    of << "</Map>" << endl;
    of.close();
  }
}

// ============================================================================
void EpetraExt::XMLWriter::
Write(const string& Label, Teuchos::ParameterList& List)
{
  TEST_FOR_EXCEPTION(IsOpen_ == false, std::logic_error,
                     "No file has been opened");

  if (Comm_.MyPID()) return;

  std::ofstream of(FileName_.c_str(), std::ios::app);

  of << "<List Label=\"" << Label << "\">" << endl;

  XMLParameterListWriter Writer;
  XMLObject Obj = Writer.toXML(List);

  of << Obj.toString();

  of << "</List>" << endl;

  of.close();
}
