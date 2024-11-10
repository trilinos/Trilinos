/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_ELEMENT_BY_ELEMENT_SINGLE_ELEMENT_H
#define ML_ELEMENT_BY_ELEMENT_SINGLE_ELEMENT_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"
#ifdef HAVE_ML_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Operator.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_FECrsGraph.h"
#include <vector>
#include "ml_epetra.h"

namespace ML_Epetra {

class ElementByElement_SingleElement: public Epetra_Operator
{
  public:
    ElementByElement_SingleElement(Epetra_Comm& Comm,
                                   const int NumMyFEs,
                                   const int NumVerticesPerFE,
                                   int* MyFEs,
                                   const int NumPDEEqns,
                                   Epetra_SerialDenseMatrix* FEMatrix,
                                   const int NumMyBoundaryRows,
                                   int* MyBoundaryRows,
                                   const double* MyBoundaryValues,
                                   const Epetra_Map& GraphMap,
                                   const int MaxEntriesPerGraphRow = 0) :
      NumMyFEs_(NumMyFEs),
      NumVerticesPerFE_(NumVerticesPerFE),
      MyFEs_(MyFEs),
      NumPDEEqns_(NumPDEEqns),
      FEMatrix_(FEMatrix),
      NumMyBoundaryRows_(NumMyBoundaryRows),
      MyBoundaryRows_(MyBoundaryRows),
      MyBoundaryValues_(MyBoundaryValues),
      Comm_(Comm),
      Graph_(0)
    {
      // build the graph using GraphMap

      Graph_ = new Epetra_FECrsGraph(Copy, GraphMap, MaxEntriesPerGraphRow);

      for (int ie = 0; ie < NumMyFEs; ++ie)
      {
        const int* ptr = &(MyFEs_[ie * NumVerticesPerFE]);
        Graph_->InsertGlobalIndices(NumVerticesPerFE, ptr,
                                    NumVerticesPerFE, ptr);
      }

      Graph_->GlobalAssemble();

      // convert MyFEs into local column map ordering

      for (int ie = 0; ie < NumMyFEs; ++ie)
      {
        int* ptr = &(MyFEs_[ie * NumVerticesPerFE]);
        for (int i = 0; i < NumVerticesPerFE; ++i)
        {
          ptr[i] = Graph_->ColMap().LID(ptr[i]);
          assert (ptr[i] != -1);
        }
      }

      // build the map for the operator, which is the "extended"
      // version of GraphMap

      std::vector<int> MyGlobalElements2(Graph_->ColMap().NumMyElements() * NumPDEEqns);
      int* MyGlobalElements = Graph_->RowMap().MyGlobalElements();
      int NumMyElements = Graph_->RowMap().NumMyElements();

      for (int i = 0; i < NumMyElements; ++i)
        for (int j = 0; j < NumPDEEqns; ++j)
          MyGlobalElements2[i * NumPDEEqns + j] = MyGlobalElements[i] * NumPDEEqns + j;
      OperatorMap_ = new Epetra_Map(-1, NumMyElements * NumPDEEqns,
                                    &MyGlobalElements2[0], 0, Comm_);

      // expand the column map as well

      MyGlobalElements = Graph_->ColMap().MyGlobalElements();
      NumMyElements = Graph_->ColMap().NumMyElements();

      for (int i = 0; i < NumMyElements; ++i)
        for (int j = 0; j < NumPDEEqns; ++j)
          MyGlobalElements2[i * NumPDEEqns + j] = MyGlobalElements[i] * NumPDEEqns + j;
      OperatorColMap_ = new Epetra_Map(-1, NumMyElements * NumPDEEqns,
                                       &MyGlobalElements2[0], 0, Comm_);

      ColImporter_ = new Epetra_Import(*OperatorColMap_, *OperatorMap_);

      // memory allocation for element-by-element multiplication

      DenseX.Reshape(NumPDEEqns * NumVerticesPerFE, NumMyFEs_);
      DenseY.Reshape(NumPDEEqns * NumVerticesPerFE, NumMyFEs_);

      // now arrange the boundary conditions. I make the assumption that
      // each BC row has been specified on a different processor.

      MyGlobalElements2.resize(NumMyBoundaryRows_);
      for (int i = 0; i < NumMyBoundaryRows_; ++i)
        MyGlobalElements2[i] = MyBoundaryRows_[i];

      Epetra_Map BoundaryMap(-1, NumMyBoundaryRows_, &MyGlobalElements2[0], 0, Comm_);

      Epetra_Vector BoundaryVector(BoundaryMap);

      for (int i = 0; i < NumMyBoundaryRows_; ++i)
        BoundaryVector[i] = MyBoundaryValues_[i];

      // now build a map containing the boundaries for ghost nodes only
      int count = 0;
      for (int i = 0; i < NumMyBoundaryRows_; ++i)
      {
        if (GraphMap.LID(MyBoundaryRows_[i]) == -1)
          MyGlobalElements2[++count] = MyBoundaryRows_[i];
      }

      ColBoundaryMap_ = new Epetra_Map(-1, count, &MyGlobalElements2[0], 0, Comm_);

      ColBoundaryVector_ = new Epetra_Vector(*ColBoundaryMap_);

      Epetra_Import Importer(*ColBoundaryMap_, BoundaryMap);
      ColBoundaryVector_->Import(BoundaryVector, Importer, Insert);

      ColBoundaryImporter_ = new Epetra_Import(*OperatorColMap_, *ColBoundaryMap_);

      for (int i = 0; i < NumMyBoundaryRows_; ++i)
        MyBoundaryRows_[i] = OperatorColMap_->LID(MyBoundaryRows_[i]);
    }

    ~ElementByElement_SingleElement()
    {
      delete Graph_;
      delete OperatorMap_;
      delete OperatorColMap_;
      delete ColImporter_;
      delete ColBoundaryMap_;
      delete ColBoundaryVector_;
      delete ColBoundaryImporter_;
    }

    const Epetra_CrsGraph& Graph() const
    {
      return(*Graph_);
    }

    int SetUseTranspose(bool UseTranspose)
    {
      if (UseTranspose)
        ML_CHK_ERR(-1);
      return(0);
    }

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      assert (X.NumVectors() == 1); // FIXME

      // import ghost nodes
      Epetra_MultiVector ColX(*OperatorColMap_, X.NumVectors());
      ML_CHK_ERR(ColX.Import(X, *ColImporter_, Insert));
      // fix boundary conditions
      ML_CHK_ERR(SetMyBoundaryRows(ColX));
      ML_CHK_ERR(ColX.Import(*ColBoundaryVector_, *ColBoundaryImporter_, Insert));

      // now redistribute the vector into Dense
      for (int ie = 0; ie < NumMyFEs_; ++ie)
      {
        const int* ptr = &(MyFEs_[ie * NumVerticesPerFE_]);
        for (int i = 0; i < NumVerticesPerFE_; ++i)
        {
          for (int j = 0; j < NumPDEEqns_; ++j)
            DenseX(i * NumPDEEqns_ + j, ie) = ColX[0][ptr[i] * NumPDEEqns_ + j];
        }
      }

      DenseY.Multiply('N', 'N', 1.0, *FEMatrix_, DenseX, 0.0);

      ColX.PutScalar(0.0);

      // put values back into Y
      for (int ie = 0; ie < NumMyFEs_; ++ie)
      {
        const int* ptr = &(MyFEs_[ie * NumVerticesPerFE_]);
        for (int i = 0; i < NumVerticesPerFE_; ++i)
        {
          for (int j = 0; j < NumPDEEqns_; ++j)
            ColX[0][ptr[i] * NumPDEEqns_ + j] += DenseY(i * NumPDEEqns_ + j, ie);
        }
      }

      Y.PutScalar(0.0);
      ML_CHK_ERR(Y.Export(ColX, *ColImporter_, Add));
      ML_CHK_ERR(ResetMyBoundaryRows(Y));
      return(0);
    }

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      ML_CHK_ERR(-1);
    }

    double NormInf() const
    {
      return(-1.0);
    }

    const char* Label() const
    {
      return("ML_Epetra::ElementByElementMatrix");
    }

    bool UseTranspose() const
    {
      return(false);
    }

    bool HasNormInf() const
    {
      return(false);
    }

    const Epetra_Comm& Comm() const
    {
      return(Comm_);
    }

    const Epetra_Map& OperatorDomainMap() const
    {
      return(*OperatorMap_);
    }

    const Epetra_Map& OperatorRangeMap() const
    {
      return(*OperatorMap_);
    }

    const Epetra_BlockMap& Map() const
    {
      return(*OperatorMap_);
    }

    int SetMyBoundaryRows(Epetra_MultiVector& Y) const
    {
      assert (Y.NumVectors() == 1);

      for (int i = 0; i < NumMyBoundaryRows_; ++i)
        Y[0][MyBoundaryRows_[i]] = MyBoundaryValues_[i];
      return(0);
    }

    int ResetMyBoundaryRows(Epetra_MultiVector& Y) const
    {
      assert (Y.NumVectors() == 1);

      for (int i = 0; i < NumMyBoundaryRows_; ++i)
        Y[0][MyBoundaryRows_[i]] = 0.0;
      return(0);
    }

  private:
    const int NumMyFEs_;
    const int NumVerticesPerFE_;
    int* MyFEs_;
    const int NumPDEEqns_;
    const Epetra_SerialDenseMatrix* FEMatrix_;
    const int NumMyBoundaryRows_;
    int* MyBoundaryRows_;
    const double* MyBoundaryValues_;
    const Epetra_Comm& Comm_;
    Epetra_FECrsGraph* Graph_;
    Epetra_Map* OperatorMap_;
    Epetra_Map* OperatorColMap_;
    Epetra_Import* ColImporter_;
    mutable Epetra_SerialDenseMatrix DenseX, DenseY;
    Epetra_Map* ColBoundaryMap_;
    Epetra_Vector* ColBoundaryVector_;
    Epetra_Import* ColBoundaryImporter_;

}; // class ElementByElementMatrix

} // namespace ML_Epetra

#endif // HAVE_ML_EPETRA
#endif // ML_ELEMENT_BY_ELEMENT_SINGLE_ELEMENT_H
