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
#include "Epetra_VbrMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_MapColoring.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "ml_RowMatrix.h"
#include "ml_lapack.h"

using namespace EpetraExt;

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
MatrixFreePreconditioner(const Epetra_Operator& Operator,
                             const Epetra_CrsGraph& Graph,
                             Teuchos::ParameterList& List,
                             Epetra_MultiVector& NullSpace) :
  IsComputed_(false),
  Label_("ML matrix-free preconditioner"),
  Comm_(Operator.Comm()),
  Operator_(Operator),
  Graph_(Graph),
  R_(0),
  C_(0),
  Time_(0),
  NewNullSpace_(0)
{
  List_ = List;

  ML_Set_PrintLevel(10); // FIXME

  // ML communicator, here based on MPI_COMM_WORLD
  ML_Comm_Create(&Comm_ML_);

  Time_ = new Epetra_Time(Comm());

  ML_CHK_ERRV(Compute(NullSpace));
}

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
~MatrixFreePreconditioner()
{
  if (Time_) delete Time_;

  if (NewNullSpace_) delete NewNullSpace_;

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
        ML_Operator** R, ML_Operator** C, int NumPDEEqns, int NullSpaceDim,
        double* NullSpace)
{
  Time().ResetStartTime();

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

  ML_Aggregate_Set_NullSpace(*MLAggr, NumPDEEqns, NullSpaceDim, NullSpace,
                             A->invec_leng);

  int NumAggregates = ML_Aggregate_Coarsen(*MLAggr, A, P, Comm_ML());
  if (NumAggregates == 0)
  {
    cerr << "Found 0 aggregates, perhaps the problem is too small." << endl;
    ML_CHK_ERR(-2);
  }

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5)
    cout << "Coarsening time = " << Time().ElapsedTime() << " (s)" << endl;

  Time().ResetStartTime();

  *R = ML_Operator_Create(Comm_ML());

  ML_Operator_Transpose_byrow(*P, *R);

  *C = ML_Operator_Create(Comm_ML());

  // FIXME: try to build an Epetra_CrsMatrix directly to save memory
  // Note: I must create an Epetra_CrsMatrix object because I need the graph
  // to use EpetraExt coloring routines.
  ML_rap(*R, A, *P, *C, ML_MSR_MATRIX);

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5)
    cout << "(block) R, P, C construction time = " << Time().ElapsedTime() << " (s)" << endl;

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Compute(Epetra_MultiVector& NullSpace)
{
  int NullSpaceDim = NullSpace.NumVectors();
  // get parameters from the list
  int    NumPDEEqns  = List_.get("PDE equations", 1);

  // =============================== //
  // basic checkings and some output //
  // =============================== //
  
  int OperatorDomainPoints =  Operator_.OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int GraphBlockRows = Graph_.NumGlobalBlockRows();
  int GraphNnz = Graph_.NumGlobalNonzeros();

  if (OperatorDomainPoints != OperatorRangePoints)
    ML_CHK_ERR(-1); // only square matrices

  if (GraphBlockRows != OperatorRangePoints / NumPDEEqns)
    ML_CHK_ERR(-2); // num PDEs seems not constant

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) 
  {
    ML_print_line("-",78);
    cout << "*** " << endl;
    cout << "*** ML_Epetra::MatrixFreePreconditioner" << endl;
    cout << "***" << endl;
    cout << "The operator domain map has " << OperatorDomainPoints;
    cout << " points, the operator range map "  << OperatorRangePoints << endl;
    cout << "The graph has " << GraphBlockRows << " rows and " << GraphNnz << endl;
    cout << "Processors used in computation = " << Comm().NumProc() << endl;
    cout << "Number of PDE equations = " << NumPDEEqns << endl;
    cout << "Null space dimension = " << NullSpaceDim << endl;
  }

  // ML wrapper for Graph_
  ML_Operator* Graph_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraCrsGraph(const_cast<Epetra_CrsGraph*>(&Graph_), Graph_ML);

  // ========================================================= //
  // building P and R for block graph. This is done by working //
  // on the Graph_ object. Support is provided for local       //
  // aggregation schemes only so that all is basically local.  //
  // Then, build the block graph coarse problem.               //
  // ========================================================= //
  
  ML_Aggregate* BlockAggr_ML = 0;
  ML_Operator* BlockPtent_ML = 0, *BlockRtent_ML = 0,* CoarseGraph_ML = 0;

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) cout << endl;

  ML_CHK_ERR(Coarsen(Graph_ML, &BlockAggr_ML, &BlockPtent_ML, &BlockRtent_ML, 
                     &CoarseGraph_ML));

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) cout << endl;

  Epetra_CrsMatrix* BlockPtent,* GraphCoarse,* BlockRtent;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(BlockPtent_ML, BlockPtent));
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(CoarseGraph_ML, GraphCoarse));
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(BlockRtent_ML, BlockRtent));

  // ================================================== //
  // coloring of block graph:                           //
  // - color of block row `i' is given by `ColorMap[i]' //
  // - number of colors is ColorMap.NumColors().        //
  // ================================================== //
  
  if (MyPID() == 0)
    cout << "Coloring with JONES_PLASSMAN..." << endl;

  Time().ResetStartTime();

  CrsGraph_MapColoring MapColoringTransform(CrsGraph_MapColoring::JONES_PLASSMAN,
                                            0, false, true);

  Epetra_MapColoring& ColorMap = MapColoringTransform(const_cast<Epetra_CrsGraph&>(GraphCoarse->Graph()));

  delete GraphCoarse;
  ML_Operator_Destroy(&CoarseGraph_ML);

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5)
    cout << "Coloring time = " << Time().ElapsedTime() << " (s)" << endl;

  // collect aggregate information, and mark all nodes that are
  // connected with each aggregate. These nodes will have a possible
  // nonzero entry after the matrix-matrix product between the Operator_
  // and the tentative prolongator.

  Time().ResetStartTime();

  int NumAggregates = BlockPtent_ML->invec_leng;
  vector<vector<int> > aggregates(NumAggregates);
  vector<int> BlockNodeList;

  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    int AID = BlockAggr_ML->aggr_info[0][i];

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
      iter = find(BlockNodeList.begin(), BlockNodeList.end(), GCID);
      if (iter == BlockNodeList.end())
        BlockNodeList.push_back(GCID);
    }
  }
  
  // get some other information about the aggregates, to be used
  // in the QR factorization of the null space. NodesOfAggregate
  // contains the local ID of block rows contained in each aggregate.

  vector< vector<int> > NodesOfAggregate(NumAggregates);

  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    int AID = BlockAggr_ML->aggr_info[0][i];
    NodesOfAggregate[AID].push_back(i);
  }

  // finally get rid of the ML_Aggregate structure.
  ML_Aggregate_Destroy(&BlockAggr_ML);

  const Epetra_Map& FineMap = Operator_.OperatorDomainMap();
  Epetra_Map CoarseMap(-1, NumAggregates * NullSpaceDim, 0, Comm());
  Epetra_BlockMap BlockNodeListMap(-1, BlockNodeList.size(), &BlockNodeList[0], 
                                   NumPDEEqns, 0, Comm());

  // ====================== //
  // process the null space //
  // ====================== //

  NewNullSpace_ = new Epetra_MultiVector(NullSpace.Map(), NullSpaceDim);

  int MaxAggrSize = 0;
  for (int i = 0; i < NumAggregates; ++i)
  {
    const int& MySize = NodesOfAggregate[i].size();
    if (MySize > MaxAggrSize) MaxAggrSize = MySize;
  }

  if (NullSpaceDim == 1)
  {
    double* ns_ptr = NullSpace.Values();

    for (int AID = 0; AID < NumAggregates; ++AID)
    {
      double dtemp = 0.0;
      for (int j = 0; j < NodesOfAggregate[AID].size(); j++)
        dtemp += (ns_ptr[NodesOfAggregate[AID][j]] * ns_ptr[NodesOfAggregate[AID][j]]); 
      dtemp = sqrt(dtemp);
      dtemp = 1.0 / dtemp;

      for (int j = 0; j < NodesOfAggregate[AID].size(); j++)
        ns_ptr[NodesOfAggregate[AID][j]] *= dtemp;
    }
  }
  else
  {
    // FIXME
    vector<double> q_ptr(100 * MaxAggrSize * NumPDEEqns * MaxAggrSize * NumPDEEqns);
    vector<double> r_ptr(100 * MaxAggrSize * NumPDEEqns * NullSpaceDim);

    vector<double> work(NullSpaceDim);
    int info;

    double* new_ptr = NewNullSpace_->Values();
    for (int AID = 0; AID < NumAggregates; ++AID)
    {
      int MySize = NodesOfAggregate[AID].size();
      int MyFullSize = NodesOfAggregate[AID].size() * NumPDEEqns;
      int lwork = NullSpaceDim;

      for (int k = 0; k < NullSpaceDim; ++k)
        for (int j = 0; j < MySize; ++j)
          for (int m = 0; m < NumPDEEqns; ++m)
            q_ptr[k * MyFullSize + j * NumPDEEqns + m] = 
              NullSpace[k][NodesOfAggregate[AID][j] * NumPDEEqns + m];

      DGEQRF_F77(&MyFullSize, &NullSpaceDim, &q_ptr[0], 
                 &MyFullSize, &r_ptr[0], &work[0], &lwork, &info);

      ML_CHK_ERR(info);

      if (work[0] > lwork) work.resize((int) work[0]);

      // FIXME: next null space here

      DORGQR_F77(&MyFullSize, &NullSpaceDim, &NullSpaceDim, 
                 &q_ptr[0], &MyFullSize, &r_ptr[0], &work[0], &lwork, &info);
      ML_CHK_ERR(info); // dgeqtr returned a non-zero

#if 0
         for (int k = 0; k < MyFullSize; ++k)
         {
           for (int kk = 0; kk < NullSpaceDim; ++kk)
             cout << "af3> " << q_ptr[kk * MyFullSize + k] << "  ";
           cout << endl;
         }
#endif

      if (work[0] > lwork) work.resize((int) work[0]);

      // insert the Q block into the null space

      for (int k = 0; k < NullSpaceDim; ++k)
        for (int j = 0; j < MySize; ++j)
          for (int m = 0; m < NumPDEEqns; ++m)
          {
            int LRID = NodesOfAggregate[AID][j] * NumPDEEqns + m;
            double& val = q_ptr[k * MyFullSize + j * NumPDEEqns + m];
            NullSpace[k][LRID] = val;
          }
    }
  }

  int NumColors = ColorMap.MaxNumColors();

  if (ML_Get_PrintLevel() > 5)
    cout << "# colors on processor " << Comm().MyPID() << " = "
        << ColorMap.NumColors() << endl;
  if (MyPID() == 0)
    cout << "Maximum # of colors = " << NumColors << endl;

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
#if 0
    (*ColoredP)[Color][i] = Values[0];
#else
    const int& MyLength = ColoredP->MyLength();
    for (int k = 0; k < NumPDEEqns; ++k)
      for (int j = 0; j < NullSpaceDim; ++j)
        (*ColoredP)[(Color * NullSpaceDim + j)][i * NumPDEEqns + k] = 
          NullSpace[j][i * NumPDEEqns + k];
#endif
  }

  Epetra_MultiVector* ColoredAP = new Epetra_MultiVector(Operator_.OperatorRangeMap(), 
                                                         NumColors * NullSpaceDim);
  Operator_.Apply(*ColoredP, *ColoredAP);
  delete ColoredP;

  Epetra_MultiVector ExtColoredAP(BlockNodeListMap, ColoredAP->NumVectors());
  Epetra_Import Importer(BlockNodeListMap, Operator_.OperatorRangeMap());
  ExtColoredAP.Import(*ColoredAP, Importer, Insert);

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5)
    cout << "colored A times P time = " << Time().ElapsedTime() << " (s)" << endl;

  Time().ResetStartTime();

  // populate the actual AP operator.

  Epetra_FECrsMatrix* AP = new Epetra_FECrsMatrix(Copy, FineMap, 0);

  for (int i = 0; i < NumAggregates; ++i)
  {
    for (int j = 0; j < aggregates[i].size(); ++j)
    {
      int GRID = aggregates[i][j];
      int LRID = BlockNodeListMap.LID(GRID); // this is the block ID
      assert (LRID != -1); // FIXME
      int GCID = CoarseMap.GID(i) / NullSpaceDim;
      assert (GCID != -1); // FIXME
      int color = ColorMap[i] - 1;
#if 0
      double val = ExtColoredAP[color][LRID];
      if (val != 0.0)
        AP->InsertGlobalValues(1, &GRID, 1, &GCID, &val);
#else
      for (int k = 0; k < NumPDEEqns; ++k)
        for (int j = 0; j < NullSpaceDim; ++j)
        {
          double val = ExtColoredAP[color * NullSpaceDim + j][LRID * NumPDEEqns + k];
          if (val != 0.0)
          {
            int GRID2 = GRID * NumPDEEqns + k;
            int GCID2 = GCID * NullSpaceDim + j;
            AP->InsertGlobalValues(1, &GRID2, 1, &GCID2, &val);
          }
        }
#endif
    }
  }

  aggregates.resize(0);
  delete ColoredAP;

  AP->GlobalAssemble(false);
  AP->FillComplete(CoarseMap, FineMap);

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5)
    cout << "Construction of final AP time = " << Time().ElapsedTime() << " (s)" << endl;

  ML_Operator* AP_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(AP, AP_ML);

  // ======== //
  // create R //
  // ======== //
  
  R_ = new Epetra_FECrsMatrix(Copy, CoarseMap, 0);

  for (int AID = 0; AID < NumAggregates; ++AID)
  {
    const int& MySize = NodesOfAggregate[AID].size();

    // FIXME: make it faster
    for (int j = 0; j < MySize; ++j)
      for (int m = 0; m < NumPDEEqns; ++m)
        for (int k = 0; k < NullSpaceDim; ++k)
        {
          int LCID = NodesOfAggregate[AID][j] * NumPDEEqns + m;
          int GCID = FineMap.GID(LCID);
          assert (GCID != -1);

          double& val = NullSpace[k][LCID];

          int GRID = CoarseMap.GID(AID * NullSpaceDim + k);
          assert (GRID != -1); // FIXME
          R_->InsertGlobalValues(1, &GRID, 1, &GCID, &val);
        }
  }

  R_->GlobalAssemble(false);
  R_->FillComplete(FineMap, CoarseMap);
  R_->OptimizeStorage(); // FIXME: TO BE DONE FOR ALL OBJECTS

  ML_Operator* R_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(R_, R_ML);

  // ======== //
  // Create C //
  // ======== //

  ML_Operator* C_ML;
  ML_matmat_mult(R_ML, AP_ML, &C_ML);

  ML_Operator_Destroy(&AP_ML);
  ML_Operator_Destroy(&R_ML);
  delete AP;

  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(C_ML, C_));
  ML_Operator_Destroy(&C_ML);

  // FIXME: NOT TRUE
  // At this point I have:
  // - C_ is the coarse matrix
  // - P_ is the prolongator operator
  // - R_ is the restriction operator

#if 0
  // free memory
  ML_Operator_Destroy(&BlockPtent_ML);
  delete BlockPtent;
  ML_Operator_Destroy(&BlockRtent_ML);
  ML_Operator_Destroy(&Graph_ML);
#endif

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) ML_print_line("-",78);

  IsComputed_ = true;

  return(0);
}
#endif
