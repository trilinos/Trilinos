/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include <iostream>
#include "ml_include.h"
#include "Teuchos_ParameterList.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Operator.h"
#include "MLAPI_DistributedMatrix.h"
#if defined(HAVE_ML_GALERI)
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#endif

using namespace std;

namespace MLAPI {

// ====================================================================== 
Operator Gallery(const string ProblemType,
                 const Space& MySpace)
{
#if defined(HAVE_ML_GALERI)
  int NumGlobalElements = MySpace.GetNumGlobalElements();

  Teuchos::ParameterList GaleriList;
  GaleriList.set("n", NumGlobalElements);
  Epetra_Map* Map = Galeri::CreateMap("Linear", GetEpetra_Comm(), GaleriList);
  Epetra_CrsMatrix* EpetraA = Galeri::CreateCrsMatrix(ProblemType, Map,
                                                      GaleriList);
  
  Operator A(MySpace,MySpace,EpetraA);

  delete Map;

  return (A);
#else
  ML_THROW("Configure with --enable-galeri", -1);
#endif

}

// ====================================================================== 
Operator GetShiftedLaplacian1D(const int NumGlobalElements,
                               const double Factor)
{

  double Pi = atan(1.0) * 4;
  double LambdaMin = 4.0 * pow((sin(1.0 / (NumGlobalElements + 1) * Pi / 2)),
                               2.0) * Factor;

  Space FineSpace(NumGlobalElements);

  DistributedMatrix A(FineSpace, FineSpace);

  if (GetMyPID() == 0) {
    for (int i = 0 ; i < NumGlobalElements ; ++i) {
      A(i, i) = 2.0 - LambdaMin;
      if (i)
        A(i, i - 1) = - 1.0;
      if (i != NumGlobalElements - 1)
        A(i, i + 1) = - 1.0;
    }
  }

  A.FillComplete();

  return(A);
}

// ====================================================================== 
Operator GetShiftedLaplacian2D(const int NX, const int NY, 
                               const double Factor,
                               const bool RandomScale)
{

  int NumGlobalElements = NX * NY;

  double Pi = atan(1.0) * 4;
  double LambdaMin = (4.0 * pow((sin(1.0 / (NX + 1) * Pi / 2)), 2.0) + 4.0 * pow((sin(1.0 / (NY + 1) * Pi / 2)), 2.0)) * Factor;

  Space FineSpace(NumGlobalElements);
  int n = 0;
  if (GetMyPID() == 0)
    n = NumGlobalElements;
  Space LocalizedSpace(-1, n);

  MultiVector Scale(LocalizedSpace);

  DistributedMatrix A(FineSpace, FineSpace);

  // assemble the matrix on processor 0 only
  if (GetMyPID() == 0) {

    if (RandomScale) {

      Scale.Random();
      for (int i = 0 ; i < Scale.GetMyLength() ; ++i)
        Scale(i) = 1.0001 + Scale(i);

      double sr, sc;
      for (int i = 0 ; i < NX ; ++i) {
        for (int j = 0 ; j < NY ; ++j) {
          int row = i + j * NX;
          sr = Scale(row);
          A(row, row) = sr * sr * (4.0 - LambdaMin);
          if (i > 0) {
            sc = Scale(row - 1);
            A(row, row - 1) =  sr * sc * (-1.0);
          }
          if (i < NX - 1) {
            sc = Scale(row + 1);
            A(row, row + 1) =  sr * sc * (-1.0);
          }
          if (j > 0) {
            sc = Scale(row - NX);
            A(row, row - NX) =  sr * sc * (-1.0);
          }
          if (j < NY - 1) {
            sc = Scale(row + NX);
            A(row, row + NX) =  sr * sc * (-1.0);
          }
        }
      }
    }
    else {

      for (int i = 0 ; i < NX ; ++i) {
        for (int j = 0 ; j < NY ; ++j) {
          int row = i + j * NX;
          A(row, row) = 4.0 - LambdaMin;
          if (i > 0) 
            A(row, row - 1) = -1.0;
          if (i < NX - 1)
            A(row, row + 1) = -1.0;
          if (j > 0) 
            A(row, row - NX) = -1.0;
          if (j < NY - 1)
            A(row, row + NX) = -1.0;
        }
      }
    }

  } // mypid == 0

  A.FillComplete();

  return (A);
}

// ====================================================================== 
Operator ReadMatrix(const char* FileName)
{

  int NumGlobalElements;
  int NumGlobalNonzeros;
  std::ifstream File;
  File.open(FileName);

  if (!File.good() )
    ML_THROW("Error opening file", -1);

  File >> NumGlobalElements;
  File >> NumGlobalNonzeros;

  if (NumGlobalElements <= 0)
    ML_THROW("Invalid number of global elements (" + GetString(NumGlobalElements) + ")", -1);

  if (NumGlobalNonzeros <= 0)
    ML_THROW("Invalid number of global nonzeros" + GetString(NumGlobalNonzeros) + ")", -1);

  if (GetMyPID()) File.close();

  Space FineSpace(NumGlobalElements);

  DistributedMatrix A(FineSpace, FineSpace);

  if (GetMyPID() == 0) {
    int row, col;
    double val;
    for (int i = 0 ; i < NumGlobalNonzeros ; ++i) {
      File >> row;
      File >> col;
      File >> val;

      if (row < 0 || row >= NumGlobalElements)
        ML_THROW("Invalid row number (" + GetString(row) + ")", -1);
      if (col < 0 || col >= NumGlobalElements)
        ML_THROW("Invalid col number (" + GetString(col) + ")", -1);

      A(row, col) = val;
    }
  }

  A.FillComplete();

  return(A);
}

// ====================================================================== 
Operator GetRecirc2D(const int NX, const int NY, const double conv,
                     const double diff)
{
  double LX = 1.0;    // length of the X-axis
  double LY = 1.0;    // lenght of the Y-axis

  Space FineSpace(NX * NY);

  DistributedMatrix A(FineSpace, FineSpace);

  if (GetMyPID() == 0) {

    double HX = LX / (NX + 1);  
    double HY = LY / (NY + 1);

    for (int ix = 0 ; ix < NX ; ++ix) {
      for( int iy = 0 ; iy < NY ; ++iy) {

        double X = HX * (ix + 1);
        double Y = HY * (iy + 1);

        double ConvX =  conv * 4 * X * (X - 1.) *(1. - 2 * Y) / HX;
        double ConvY = -conv * 4 * Y * (Y - 1.) *(1. - 2 * X) / HY;

        double lower = 0.0, upper = 0.0;
        double left = 0.0,  right = 0.0, center = 0.0;

        if( ConvX<0 ) {
          right  += ConvX;
          center -= ConvX;
        } else {
          left   -= ConvX;
          center += ConvX;
        }

        if( ConvY<0 ) {
          upper  += ConvY;
          center -= ConvY;
        } else {
          lower  -= ConvY;
          center += ConvY;
        }

        center += diff * 2. / (HX * HX) + diff * 2. / (HY * HY);
        left   -= diff / (HX * HX);
        right  -= diff / (HX * HX);
        lower  -= diff / (HY * HY);
        upper  -= diff / (HY * HY);

        int row = iy * NX + ix;

        if (ix != 0)
          A(row, row - 1) = left;

        if (ix != NX - 1)
          A(row, row + 1) = right;

        if (iy != 0)
          A(row, row - NX) = lower;

        if (iy != NY - 1)
          A(row, row + NX) = upper;

        A(row, row) = center;
      }
    }
  }

  A.FillComplete();

  return(A);
}

// ====================================================================== 
Teuchos::ParameterList ReadParameterList(const char* FileName)
{
  std::ifstream fp;

  fp.open(FileName);
  if (!fp.good()) 
    ML_THROW("Error opening file " + string(FileName), -1);

  Teuchos::ParameterList List;

  string line;
  while (getline(fp, line, '\n'))
  {
    char type = line[0];
    if (type == '#')
      continue;

    int i = line.find(" = ");
    string name = line.substr(2, i - 2);
    string value = line.substr(i + 3);

    if (type == 'i') 
      List.set(name, atoi(value.c_str()));
    else if (type == 'f')
      List.set(name, (double)atof(value.c_str()));
    else if (type == 's')
      List.set(name, value);
    else if (type == 'b') {
      if (value == "true")
        List.set(name, true);
      else
        List.set(name, false);
    }
    else
      ML_THROW("Type not valid", -1);
  }

  fp.close();

  return(List);
}

// ====================================================================== 
} // namespace MLAPI

#endif // HAVE_ML_MLAPI
