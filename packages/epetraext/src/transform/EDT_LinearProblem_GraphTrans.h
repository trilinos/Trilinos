
#ifndef EDT_LINEARPROBLEM_GRAPHTRANS_H
#define EDT_LINEARPROBLEM_GRAPHTRANS_H

#include <Epetra_Transform.h>

class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_Export;
class Epetra_Import;

namespace EpetraExt {

class LinearProblem_GraphTrans : public SameTypeTransform<Epetra_LinearProblem>
{
  StructuralSameTypeTransform<Epetra_CrsGraph> & graphTrans_;

  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

  Epetra_LinearProblem * OldProblem_;
  Epetra_CrsGraph * OldGraph_;
  Epetra_CrsMatrix * OldMatrix_;
  Epetra_MultiVector * OldLHS_;
  Epetra_MultiVector * OldRHS_;
  Epetra_Map * OldRowMap_;

  Epetra_LinearProblem * NewProblem_;
  Epetra_CrsMatrix * NewMatrix_;
  Epetra_MultiVector * NewLHS_;
  Epetra_MultiVector * NewRHS_;

 public:

  ~LinearProblem_GraphTrans();

  //Note: User is responsible for destroying input Graph Transform.
  LinearProblem_GraphTrans( StructuralSameTypeTransform<Epetra_CrsGraph> & graph_trans )
  : graphTrans_(graph_trans),
    OldProblem_(0),
    OldGraph_(0),
    OldMatrix_(0),
    OldLHS_(0),
    OldRHS_(0),
    OldRowMap_(0),
    NewProblem_(0),
    NewMatrix_(0),
    NewLHS_(0),
    NewRHS_(0),
    Importer_(0),
    Exporter_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();

};

} //namespace EpetraExt

#endif //EDT_LINEARPROBLEM_GRAPHTRANS_H

