#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_EPETRAEXT)
#include "ml_MatrixFreePreconditioner.h"
#include "ml_aggregate.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_MapColoring.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "ml_RowMatrix.h"

using namespace EpetraExt;

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
MatrixFreePreconditioner(const Epetra_Operator& Operator,
                             const Epetra_CrsGraph& Graph,
                             Teuchos::ParameterList& List,
                             const Epetra_MultiVector& NullSpace) :
  Label_("ML matrix-free preconditioner"),
  Comm_(Operator.Comm()),
  Operator_(Operator),
  Graph_(Graph),
  R_(0),
  C_(0)
{
  // ML communicator, here based on MPI_COMM_WORLD
  ML_Comm_Create(&Comm_ML_);

  ML_CHK_ERRV(Compute(NullSpace));
}

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
~MatrixFreePreconditioner()
{
  ML_Comm_Destroy(&Comm_ML_);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
SetUseTranspose(bool UseTranspose)
{
  ML_RETURN(-1);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  ML_RETURN(-1);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Coarsen(ML_Operator*A, ML_Aggregate** MLAggr, ML_Operator** P, 
        ML_Operator** R, ML_Operator** C)
{
  // Aggregate object, with settings
  ML_Aggregate_Create(MLAggr);
  
  string CoarsenType = List_.get("aggregation: type", "Uncoupled");
  double Threshold   = List_.get("aggregation: threshold", 0.0);

  ML_Aggregate_Set_MaxLevels(*MLAggr, 2);
  ML_Aggregate_Set_StartLevel(*MLAggr, 0);
  ML_Aggregate_Set_Threshold(*MLAggr, Threshold);
  (*MLAggr)->cur_level = 0;
  ML_Aggregate_Set_Reuse(*MLAggr);
  (*MLAggr)->keep_agg_information = 1;
  
  *P = ML_Operator_Create(Comm_ML());

  if (CoarsenType == "Uncoupled") 
    (*MLAggr)->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "METIS")
    (*MLAggr)->coarsen_scheme = ML_AGGR_METIS;
  else 
  {
    ML_CHK_ERR(-1);
  }

  int NumAggregates = ML_Aggregate_Coarsen(*MLAggr, A, P, Comm_ML());
  if (NumAggregates == 0)
  {
    cerr << "Found 0 aggregates, perhaps the problem is too small." << endl;
    ML_CHK_ERR(-2);
  }

  *R = ML_Operator_Create(Comm_ML());

  ML_Operator_Transpose_byrow(*P, *R);

  *C = ML_Operator_Create(Comm_ML());

  // FIXME: try to build an Epetra_CrsMatrix directly to save memory
  // Note: I must create an Epetra_CrsMatrix object because I need the graph
  // to use EpetraExt coloring routines.
  ML_rap(*R, A, *P, *C, ML_MSR_MATRIX);

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Compute(const Epetra_MultiVector& NullSpace)
{
  int    NumPDEEqns  = List_.get("PDE equations", 1);

  // ML wrapper for Graph_
  ML_Operator* Graph_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraCrsGraph(const_cast<Epetra_CrsGraph*>(&Graph_), Graph_ML);

  // ========================================================= //
  // building P and R for block graph. This is done by working //
  // on the Graph_ object. Support is provided for local       //
  // aggregation schemes only so that all is basically local.  //
  // Then, build the block graph coarse problem.               //
  // ========================================================= //
  
  ML_Aggregate* MLAggr = 0;
  ML_Operator* BlockPtent_ML = 0, *BlockRtent_ML = 0,* CoarseGraph_ML = 0;
  ML_CHK_ERR(Coarsen(Graph_ML, &MLAggr, &BlockPtent_ML, &BlockRtent_ML, 
                     &CoarseGraph_ML));

  Epetra_CrsMatrix* BlockPtent,* GraphCoarse;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(BlockPtent_ML, BlockPtent));
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(CoarseGraph_ML, GraphCoarse));

  // ================================================== //
  // coloring of block graph:                           //
  // - color of block row `i' is given by `ColorMap[i]' //
  // - number of colors is ColorMap.NumColors().        //
  // ================================================== //
  
  CrsGraph_MapColoring MapColoringTransform(CrsGraph_MapColoring::JONES_PLASSMAN,
                                            0, false, true);

  Epetra_MapColoring& ColorMap = MapColoringTransform(const_cast<Epetra_CrsGraph&>(GraphCoarse->Graph()));

  delete GraphCoarse;
  ML_Operator_Destroy(&CoarseGraph_ML);

  // collect aggregate information, and mark all nodes that are
  // connected with each aggregate. These nodes will have a possible
  // nonzero entry after the matrix-matrix product between the Operator_
  // and the tentative prolongator.

  int NumAggregates = BlockPtent_ML->invec_leng;
  vector<vector<int> > aggregates(NumAggregates);
  vector<int> NodeList;

  for (int i = 0; i < Graph_.NumMyRows(); ++i)
  {
    int AID = MLAggr->aggr_info[0][i];

    int NumEntries;
    int* Indices;

    Graph_.ExtractMyRowView(i, NumEntries, Indices);

    vector<int>::iterator iter;
    for (int k = 0; k < NumEntries; ++k)
    {
      const int& GCID = Graph_.ColMap().GID(Indices[k]);
      iter = find(aggregates[AID].begin(), aggregates[AID].end(), GCID);
      if (iter == aggregates[AID].end())
        aggregates[AID].push_back(GCID);
      iter = find(NodeList.begin(), NodeList.end(), GCID);
      if (iter == NodeList.end())
        NodeList.push_back(GCID);
    }
  }
  
  ML_Aggregate_Destroy(&MLAggr);

  int NumColors = ColorMap.MaxNumColors();

  cout << "# colors on processor " << Comm().MyPID() << " = " 
       << ColorMap.NumColors() << endl;
  if (Comm().MyPID() == 0)
    cout << "Maximum # of colors = " << NumColors << endl;

  int NullSpaceDim = NullSpace.NumVectors();

  // grab map form the operator
  const Epetra_Map& FineMap = Operator_.OperatorDomainMap();
  Epetra_Map CoarseMap(-1, NumAggregates * NullSpaceDim, 0, Comm());
  Epetra_Map NodeListMap(-1, NodeList.size(), &NodeList[0], 0, Comm());

  Epetra_MultiVector* ColoredP = new Epetra_MultiVector(FineMap, NumColors * NullSpaceDim);
  ColoredP->PutScalar(0.0);
  for (int i = 0; i < BlockPtent->NumMyRows(); ++i)
  {
    int NumEntries;
    int* Indices;
    double* Values;
    BlockPtent->ExtractMyRowView(i, NumEntries, Values, Indices);
    assert (NumEntries == 1); // this is the block P
    const int& Color = ColorMap[Indices[0]] - 1;
    (*ColoredP)[Color][i] = Values[0];
#if 0
    for (int k = 0; k < NumPDEEqns; ++k)
      for (int j = 0; j < NullSpaceDim; ++j)
        ColoredP[(Color * NullSpaceDim + j)][i * NumPDEEqns + k] = 
          NullSpace[j][i * NumPDEEqns + k];
#endif
  }

  Epetra_MultiVector* ColoredAP = new Epetra_MultiVector(Graph_.RangeMap(), 
                                                         NumColors * NullSpaceDim);
  Operator_.Apply(*ColoredP, *ColoredAP);
  delete ColoredP;

  Epetra_MultiVector ExtColoredAP(NodeListMap, ColoredAP->NumVectors());
  Epetra_Import Importer(NodeListMap, Graph_.RangeMap());
  ExtColoredAP.Import(*ColoredAP, Importer, Insert);

  // populate the actual AP operator.

  Epetra_FECrsMatrix* AP = new Epetra_FECrsMatrix(Copy, FineMap, 0);

  for (int i = 0; i < NumAggregates; ++i)
  {
    for (int j = 0; j < aggregates[i].size(); ++j)
    {
      int GRID = aggregates[i][j];
      int LRID = NodeListMap.LID(GRID);
      assert (LRID != -1);
      int GCID = CoarseMap.GID(i);
      assert (GCID != -1);
      int color = ColorMap[i] - 1;
      double val = ExtColoredAP[color][LRID];
      if (val != 0.0)
        AP->InsertGlobalValues(1, &GRID, 1, &GCID, &val);
    }
  }

  aggregates.resize(0);
  delete ColoredAP;

  AP->GlobalAssemble(false);
  AP->FillComplete(CoarseMap, FineMap);

  ML_Operator* AP_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(AP, AP_ML);

  ML_Operator* C_ML;
  ML_matmat_mult(BlockRtent_ML, AP_ML, &C_ML);

  ML_Operator_Destroy(&AP_ML);
  delete AP;

  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(C_ML, C_));
  ML_Operator_Destroy(&C_ML);

  // At this point I have:
  // - C_ is the coarse matrix
  // - P_ is the prolongator operator
  // - R_ is the restriction operator

  // free memory
  ML_Operator_Destroy(&BlockPtent_ML);
  delete BlockPtent;
  ML_Operator_Destroy(&BlockRtent_ML);
  ML_Operator_Destroy(&Graph_ML);

  return(0);
}
#endif
