#ifndef GALERI_UTILS_H
#define GALERI_UTILS_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"

class Epetra_BlockMap;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_LinearProblem;
namespace Teuchos {
  class ParameterList;
}

namespace Galeri {

Epetra_MultiVector* 
CreateCartesianCoordinates(const string CoordType,
                           const Epetra_BlockMap* BlockMap,
                           Teuchos::ParameterList& List);

void  Solve(const Epetra_LinearProblem Problem);

void Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS,
             const Epetra_MultiVector* RHS);

double ComputeNorm(const Epetra_MultiVector* LHS,
                   const Epetra_MultiVector* RHS);

double ComputeNorm(const Epetra_RowMatrix* A,
                   const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);

string toString(const int& x);

string toString(const unsigned int& x);

string toString(const double& x);

void 
GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                         int & left, int & right, int & lower, int & upper);

void 
GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                         int& left, int& right, int& lower, int& upper,
                         int& left2, int& right2, int& lower2, int& upper2);

void 
GetNeighboursCartesian3d(const int i, const int nx, const int ny, const int nz,
                         int& left, int& right, int& lower, int& upper,
                         int& below, int& above);
void
PrintStencil2D(const Epetra_CrsMatrix* Matrix,
               const int nx, const int ny, int GID = -1);
}; // namespace Galeri

#endif
