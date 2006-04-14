#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT)
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
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

#include "Epetra_MapColoring.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "ml_RowMatrix.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_lapack.h"

using namespace EpetraExt;

#define ML_MFP_ADDITIVE 0
#define ML_MFP_HYBRID   1

// ============================================================================ 
ML_Epetra::MatrixFreePreconditioner::
MatrixFreePreconditioner(const Epetra_Operator& Operator,
                             const Epetra_CrsGraph& Graph,
                             Teuchos::ParameterList& List,
                             Epetra_MultiVector& NullSpace) :
  Comm_ML_(0),
  Comm_(Operator.Comm()),
  Label_("ML matrix-free preconditioner"),
  IsComputed_(false),
  PrecType_(ML_MFP_HYBRID),
  omega_(1.00),
  Operator_(Operator),
  Graph_(Graph),
  R_(0),
  C_(0),
  C_ML_(0),
  MLP_(0),
  Time_(0)
{
  List_ = List;

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

  if (R_) delete R_;

  if (MLP_) delete MLP_;

  if (C_) delete C_;

  if (C_ML_) ML_Operator_Destroy(&C_ML_);

  ML_Comm_Destroy(&Comm_ML_);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyJacobi(Epetra_MultiVector& X, const double omega) const
{
  ML_CHK_ERR(ApplyInvBlockDiag(omega, X, 0.0, X));

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyJacobi(Epetra_MultiVector& X, const Epetra_MultiVector& B,
            const double omega, Epetra_MultiVector& tmp) const
{
  Operator_.Apply(X, tmp);
  tmp.Update(1.0, B, -1.0);
  ML_CHK_ERR(ApplyInvBlockDiag(omega, X, 1.0, tmp));
  ///ML_CHK_ERR(X.Multiply('T', 'N', omega, *InvBlockDiag_, tmp, 1.0));

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyInverse(const Epetra_MultiVector& Xinput, Epetra_MultiVector& Y) const
{
  ResetStartTime();

  // recall:: Y is a view of Xinput for AztecOO
  Epetra_MultiVector X(Xinput);

  if (!Y.Map().SameAs(R_->OperatorDomainMap())) ML_CHK_ERR(-1);
  if (X.NumVectors() != Y.NumVectors()) ML_CHK_ERR(-1);

  Epetra_MultiVector X_c(R_->OperatorRangeMap(), X.NumVectors());
  Epetra_MultiVector Y_c(R_->OperatorRangeMap(), X.NumVectors());

  // ================================= //
  // ADDITIVE TWO-LEVEL PRECONDITIONER //
  // ================================= //

  if (PrecType_ == ML_MFP_ADDITIVE)
  {
    ML_CHK_ERR(R_->Multiply(false, X, X_c));

    ML_CHK_ERR(MLP_->ApplyInverse(X_c, Y_c));

    ML_CHK_ERR(R_->Multiply(true, Y_c, Y));

    ML_CHK_ERR(ApplyInvBlockDiag(1.0, X, 1.0, X));
    ////ML_CHK_ERR(Y.Multiply('T', 'N', 1.0, *InvBlockDiag_, X, 1.0));
  }
  else if (PrecType_ == ML_MFP_HYBRID)
  {
    Epetra_MultiVector Ytmp(X.Map(), X.NumVectors());

    // apply pre-smoother
    ML_CHK_ERR(ApplyJacobi(X, omega_));
    // new residual
    ML_CHK_ERR(Operator_.Apply(X, Ytmp));
    ML_CHK_ERR(Ytmp.Update(1.0, Y, -1.0));
    // restrict to coarse
    ML_CHK_ERR(R_->Multiply(false, Ytmp, X_c));
    // solve coarse problem
    ML_CHK_ERR(MLP_->ApplyInverse(X_c, Y_c));
    // prolongate back
    ML_CHK_ERR(R_->Multiply(true, Y_c, Ytmp));
    // add to solution, X now has the correction
    ML_CHK_ERR(X.Update(1.0, Ytmp, 1.0));
    // apply post-smoother
    ML_CHK_ERR(ApplyJacobi(X, Y, omega_, Ytmp));
    Y = X;
  }
  else
    ML_CHK_ERR(-3); // type not recognized

  AddAndResetStartTime("ApplyInverse()");

  return(0);
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

  *R = ML_Operator_Create(Comm_ML());

  ML_Operator_Transpose_byrow(*P, *R);

  *C = ML_Operator_Create(Comm_ML());

  // FIXME: try to build an Epetra_CrsMatrix directly to save memory
  // Note: I must create an Epetra_CrsMatrix object because I need the graph
  // to use EpetraExt coloring routines.
  ////ML_rap(*R, A, *P, *C, ML_EpetraCRS_MATRIX);
  ML_rap(*R, A, *P, *C, ML_MSR_MATRIX);

  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
Compute(Epetra_MultiVector& NullSpace)
{
  Epetra_Time TotalTime(Comm());

  const int NullSpaceDim = NullSpace.NumVectors();
  // get parameters from the list
  string PrecType = List_.get("prec: type", "hybrid");
  string ColoringType = List_.get("coloring: type", "JONES_PLASSMAN");
  string DiagonalColoringType = List_.get("diagonal coloring: type", "JONES_PLASSMAN");
  int OutputLevel = List_.get("output", 10);
  omega_ = List_.get("smoother: damping", omega_);
  ML_Set_PrintLevel(OutputLevel);

  // ================ //
  // check parameters //
  // ================ //

  if (PrecType == "hybrid")
    PrecType_ = ML_MFP_HYBRID;
  else if (PrecType == "additive")
    PrecType_ = ML_MFP_ADDITIVE;
  else
    ML_CHK_ERR(-3); // not recognized

  // =============================== //
  // basic checkings and some output //
  // =============================== //
  
  int OperatorDomainPoints =  Operator_.OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int GraphBlockRows = Graph_.NumGlobalBlockRows();
  int GraphNnz = Graph_.NumGlobalNonzeros();
  int NumPDEEqns = OperatorRangePoints / GraphBlockRows;

  if (OperatorDomainPoints != OperatorRangePoints)
    ML_CHK_ERR(-1); // only square matrices

  if (OperatorRangePoints % NumPDEEqns != 0)
    ML_CHK_ERR(-2); // num PDEs seems not constant

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) 
  {
    ML_print_line("=",78);
    cout << "*** " << endl;
    cout << "*** ML_Epetra::MatrixFreePreconditioner" << endl;
    cout << "***" << endl;
    cout << "The operator domain map has " << OperatorDomainPoints;
    cout << " points, the operator range map "  << OperatorRangePoints << endl;
    cout << "The graph has " << GraphBlockRows << " rows and " << GraphNnz << " nonzeros" << endl;
    cout << "Processors used in computation = " << Comm().NumProc() << endl;
    cout << "Number of PDE equations        = " << NumPDEEqns << endl;
    cout << "Null space dimension           = " << NullSpaceDim << endl;
    cout << "Preconditioner type            = " << PrecType << endl;
    cout << "Coloring type                  = " << ColoringType << endl;
    cout << "Diagonal coloring type         = " << DiagonalColoringType << endl;
  }

  ResetStartTime();

  // probes for the block diagonal of the matrix.
  ML_CHK_ERR(GetBlockDiagonal(DiagonalColoringType));

  AddAndResetStartTime("block diagonal construction", true);

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

  if (ML_Get_PrintLevel() > 5 && MyPID() == 0) cout << endl;

  ML_CHK_ERR(Coarsen(Graph_ML, &BlockAggr_ML, &BlockPtent_ML, &BlockRtent_ML, 
                     &CoarseGraph_ML));

  if (ML_Get_PrintLevel() > 5 && MyPID() == 0) cout << endl;

  Epetra_CrsMatrix* GraphCoarse;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(CoarseGraph_ML, GraphCoarse));

  int NumAggregates = BlockPtent_ML->invec_leng;
  ML_Operator_Destroy(&BlockRtent_ML);
  ML_Operator_Destroy(&CoarseGraph_ML);

  AddAndResetStartTime("construction of block C, R, and P", true);
  if (ML_Get_PrintLevel() > 5) cout << endl;

  // ================================================== //
  // coloring of block graph:                           //
  // - color of block row `i' is given by `ColorMap[i]' //
  // - number of colors is ColorMap.NumColors().        //
  // ================================================== //
  
  ResetStartTime();

  CrsGraph_MapColoring* MapColoringTransform;
  
  if (ColoringType == "JONES_PLASSMAN")
    MapColoringTransform = new CrsGraph_MapColoring (CrsGraph_MapColoring::JONES_PLASSMAN,
                                                     0, false, true);
  else if (ColoringType == "PSEUDO_PARALLEL")
    MapColoringTransform = new CrsGraph_MapColoring (CrsGraph_MapColoring::PSEUDO_PARALLEL,
                                                     0, false, true);
  else 
    ML_CHK_ERR(-1);

  Epetra_MapColoring& ColorMap = (*MapColoringTransform)(const_cast<Epetra_CrsGraph&>(GraphCoarse->Graph()));

  const int NumColors = ColorMap.MaxNumColors();

  delete GraphCoarse;

  AddAndResetStartTime("coarse graph coloring", true);
  if (ML_Get_PrintLevel() > 5) cout << endl;

  // get some other information about the aggregates, to be used
  // in the QR factorization of the null space. NodesOfAggregate
  // contains the local ID of block rows contained in each aggregate.

  // FIXME: make it faster
  vector< vector<int> > NodesOfAggregate(NumAggregates);

  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    int AID = BlockAggr_ML->aggr_info[0][i];
    NodesOfAggregate[AID].push_back(i);
  }

  int MaxAggrSize = 0;
  for (int i = 0; i < NumAggregates; ++i)
  {
    const int& MySize = NodesOfAggregate[i].size();
    if (MySize > MaxAggrSize) MaxAggrSize = MySize;
  }

  // collect aggregate information, and mark all nodes that are
  // connected with each aggregate. These nodes will have a possible
  // nonzero entry after the matrix-matrix product between the Operator_
  // and the tentative prolongator.

  vector<vector<int> > aggregates(NumAggregates);
  vector<int>::iterator iter;

  for (int i = 0; i < NumAggregates; ++i)
    aggregates[i].reserve(MaxAggrSize);

  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    int AID = BlockAggr_ML->aggr_info[0][i];

    int NumEntries;
    int* Indices;

    Graph_.ExtractMyRowView(i, NumEntries, Indices);

    for (int k = 0; k < NumEntries; ++k)
    {
      // FIXME: use hash??
      const int& GCID = Graph_.ColMap().GID(Indices[k]);

      iter = find(aggregates[AID].begin(), aggregates[AID].end(), GCID);
      if (iter == aggregates[AID].end())
        aggregates[AID].push_back(GCID);
    }
  }
  
  int* BlockNodeList = Graph_.ColMap().MyGlobalElements();

  // finally get rid of the ML_Aggregate structure.
  ML_Aggregate_Destroy(&BlockAggr_ML);

  const Epetra_Map& FineMap = Operator_.OperatorDomainMap();
  Epetra_Map CoarseMap(-1, NumAggregates * NullSpaceDim, 0, Comm());
  Epetra_Map* BlockNodeListMap = new Epetra_Map(-1, Graph_.ColMap().NumMyElements(),
                                                BlockNodeList, 0, Comm());

  vector<int> NodeList(Graph_.ColMap().NumMyElements() * NumPDEEqns);
  for (int i = 0; i < Graph_.ColMap().NumMyElements(); ++i)
    for (int m = 0; m < NumPDEEqns; ++m)
      NodeList[i * NumPDEEqns + m] = BlockNodeList[i] * NumPDEEqns + m;
  Epetra_Map NodeListMap(-1, NodeList.size(), &NodeList[0], 0, Comm());

  AddAndResetStartTime("data structures", true);

  // ====================== //
  // process the null space //
  // ====================== //

  Epetra_MultiVector NewNullSpace(NullSpace.Map(), NullSpaceDim);
  NewNullSpace.PutScalar(0.0);

  if (NullSpaceDim == 1)
  {
    double* ns_ptr = NullSpace.Values();

    for (int AID = 0; AID < NumAggregates; ++AID)
    {
      double dtemp = 0.0;
      for (int j = 0; j < NodesOfAggregate[AID].size(); j++)
        for (int m = 0; m < NumPDEEqns; ++m)
        {
          const int& pos = NodesOfAggregate[AID][j] * NumPDEEqns + m;
          dtemp += (ns_ptr[pos] * ns_ptr[pos]);
        }
      dtemp = sqrt(dtemp);

      NewNullSpace[0][AID] = dtemp;

      dtemp = 1.0 / dtemp;

      for (int j = 0; j < NodesOfAggregate[AID].size(); j++)
        for (int m = 0; m < NumPDEEqns; ++m)
          ns_ptr[NodesOfAggregate[AID][j] * NumPDEEqns + m] *= dtemp;
    }
  }
  else
  {
    // FIXME
    vector<double> qr_ptr(MaxAggrSize * NumPDEEqns * MaxAggrSize * NumPDEEqns);
    vector<double> tmp_ptr(MaxAggrSize * NumPDEEqns * NullSpaceDim);

    vector<double> work(NullSpaceDim);
    int info;

    for (int AID = 0; AID < NumAggregates; ++AID)
    {
      int MySize = NodesOfAggregate[AID].size();
      int MyFullSize = NodesOfAggregate[AID].size() * NumPDEEqns;
      int lwork = NullSpaceDim;

      for (int k = 0; k < NullSpaceDim; ++k)
        for (int j = 0; j < MySize; ++j)
          for (int m = 0; m < NumPDEEqns; ++m)
            qr_ptr[k * MyFullSize + j * NumPDEEqns + m] = 
              NullSpace[k][NodesOfAggregate[AID][j] * NumPDEEqns + m];

      DGEQRF_F77(&MyFullSize, (int*)&NullSpaceDim, &qr_ptr[0], 
                 &MyFullSize, &tmp_ptr[0], &work[0], &lwork, &info);

      ML_CHK_ERR(info);

      if (work[0] > lwork) work.resize((int) work[0]);

      // the upper triangle of qr_tmp is now R, so copy that into the 
      //  new nullspace

      for (int j = 0; j < NullSpaceDim; j++)
        for (int k = j; k < NullSpaceDim; k++)
          NewNullSpace[k][AID * NullSpaceDim + j] = qr_ptr[j + MyFullSize * k];
		 
      // to get this block of P, need to run qr_tmp through another LAPACK 
      // function:

      DORGQR_F77(&MyFullSize, (int*)&NullSpaceDim, (int*)&NullSpaceDim, 
                 &qr_ptr[0], &MyFullSize, &tmp_ptr[0], &work[0], &lwork, &info);
      ML_CHK_ERR(info); // dgeqtr returned a non-zero

      if (work[0] > lwork) work.resize((int) work[0]);

      // insert the Q block into the null space

      for (int k = 0; k < NullSpaceDim; ++k)
        for (int j = 0; j < MySize; ++j)
          for (int m = 0; m < NumPDEEqns; ++m)
          {
            int LRID = NodesOfAggregate[AID][j] * NumPDEEqns + m;
            double& val = qr_ptr[k * MyFullSize + j * NumPDEEqns + m];
            NullSpace[k][LRID] = val;
          }
    }
  }

  AddAndResetStartTime("null space setup", true);

  if (ML_Get_PrintLevel() > 5)
    cout << "# colors on processor " << Comm().MyPID() << " = "
        << ColorMap.NumColors() << endl;
  if (ML_Get_PrintLevel() > 5 && MyPID() == 0)
    cout << "Maximum # of colors = " << NumColors << endl;

  Epetra_MultiVector* ColoredP = new Epetra_MultiVector(FineMap, NumColors * NullSpaceDim);
  ColoredP->PutScalar(0.0);

  for (int i = 0; i < BlockPtent_ML->outvec_leng; ++i)
  {
    int allocated = 1;
    int NumEntries;
    int Indices;
    double Values;
    int ierr = ML_Operator_Getrow(BlockPtent_ML, 1 ,&i, allocated,
                                    &Indices,&Values,&NumEntries);
    if (ierr < 0)
      ML_CHK_ERR(-1);
      
    assert (NumEntries == 1); // this is the block P
    const int& Color = ColorMap[Indices] - 1;
    const int& MyLength = ColoredP->MyLength();
    for (int k = 0; k < NumPDEEqns; ++k)
      for (int j = 0; j < NullSpaceDim; ++j)
        (*ColoredP)[(Color * NullSpaceDim + j)][i * NumPDEEqns + k] = 
          NullSpace[j][i * NumPDEEqns + k];
  }

  ML_Operator_Destroy(&BlockPtent_ML);

  Epetra_MultiVector* ColoredAP = new Epetra_MultiVector(Operator_.OperatorRangeMap(), 
                                                         NumColors * NullSpaceDim);
  Operator_.Apply(*ColoredP, *ColoredAP);
  delete ColoredP;

  // FIXME: only if NumProc > 1
  // FIXME: delete after usage
  Epetra_MultiVector* ExtColoredAP = 0;
  ExtColoredAP = new Epetra_MultiVector(NodeListMap, ColoredAP->NumVectors());
  Epetra_Import* Importer;
  Importer = new Epetra_Import(NodeListMap, Operator_.OperatorRangeMap());
  ExtColoredAP->Import(*ColoredAP, *Importer, Insert);
  delete Importer;
  delete ColoredAP;

  AddAndResetStartTime("computation of AP", true); 

  // populate the actual AP operator, skip some controls to make it faster

  Epetra_FECrsMatrix* AP = new Epetra_FECrsMatrix(Copy, FineMap, MaxAggrSize * NumPDEEqns);

  for (int i = 0; i < NumAggregates; ++i)
  {
    for (int j = 0; j < aggregates[i].size(); ++j)
    {
      int GRID = aggregates[i][j];
      int LRID = BlockNodeListMap->LID(GRID); // this is the block ID
      //assert (LRID != -1);
      int GCID = CoarseMap.GID(i * NullSpaceDim);
      //assert (GCID != -1); 
      int color = ColorMap[i] - 1;
      for (int k = 0; k < NumPDEEqns; ++k)
        for (int j = 0; j < NullSpaceDim; ++j)
        {
          double val = (*ExtColoredAP)[color * NullSpaceDim + j][LRID * NumPDEEqns + k];
          if (val != 0.0)
          {
            int GRID2 = GRID * NumPDEEqns + k;
            int GCID2 = GCID + j;
            AP->InsertGlobalValues(1, &GRID2, 1, &GCID2, &val);
            //if (ierr < 0) ML_CHK_ERR(ierr);
          }
        }
    }
  }

  aggregates.resize(0);
  delete ExtColoredAP;
  delete BlockNodeListMap;
  delete MapColoringTransform;

  AP->GlobalAssemble(false);
  AP->FillComplete(CoarseMap, FineMap);
  AP->OptimizeStorage();

  AddAndResetStartTime("computation of the final AP", true); 

  ML_Operator* AP_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(AP, AP_ML);

  // ======== //
  // create R //
  // ======== //
  
  vector<int> REntries(NumAggregates * NullSpaceDim);
  for (int AID = 0; AID < NumAggregates; ++AID)
  {
    for (int m = 0; m < NullSpaceDim; ++m)
      REntries[AID * NullSpaceDim + m] = NodesOfAggregate[AID].size() * NumPDEEqns;
  }

  R_ = new Epetra_CrsMatrix(Copy, CoarseMap, &REntries[0], true);
  REntries.resize(0);

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
          int ierr = R_->InsertGlobalValues(GRID, 1, &val, &GCID);
          if (ierr < 0)
            ML_CHK_ERR(-1);
        }
  }

  NodesOfAggregate.resize(0);

  R_->FillComplete(FineMap, CoarseMap);
  R_->OptimizeStorage(); // FIXME: TO BE DONE FOR ALL OBJECTS

  ML_Operator* R_ML = ML_Operator_Create(Comm_ML());
  ML_Operator_WrapEpetraMatrix(R_, R_ML);

  AddAndResetStartTime("computation of the R", true); 

  // ======== //
  // Create C //
  // ======== //

  C_ML_ = ML_Operator_Create(Comm_ML());
  ML_2matmult(R_ML, AP_ML, C_ML_, ML_MSR_MATRIX);

  ML_Operator_Destroy(&AP_ML);
  ML_Operator_Destroy(&R_ML);
  delete AP;

  C_ = new ML_Epetra::RowMatrix(C_ML_, &Comm(), false);
  assert (R_->OperatorRangeMap().SameAs(C_->OperatorDomainMap()));

  double SetupTime = TotalTime.ElapsedTime();
  TotalTime.ResetStartTime();

  AddAndResetStartTime("computation of the C", true); 

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) 
  {
    cout << "Matrix-free preconditioner built. Now building solver for C..." << endl; 
  }

  Teuchos::ParameterList& sublist = List_.sublist("ML list");
  sublist.set("PDE equations", NullSpaceDim);
  sublist.set("null space: type", "pre-computed");
  sublist.set("null space: dimension", NewNullSpace.NumVectors());
  sublist.set("null space: vectors", NewNullSpace.Values());

  MLP_ = new MultiLevelPreconditioner(*C_, sublist, true);

  IsComputed_ = true;

  assert (MLP_ != 0);

  AddAndResetStartTime("computation of the preconditioner for C", true); 

  if (MyPID() == 0 && ML_Get_PrintLevel() > 5) 
  {
    cout << endl;
    cout << "Total CPU time for construction (all included) = ";
    cout << TotalCPUTime() << endl;
    ML_print_line("=",78);
  }

  return(0);
}

// ============================================================================ 
double ML_Epetra::MatrixFreePreconditioner::
TotalCPUTime() const
{
  double TotalCPUTime = 0.0;
  map<string, double>::iterator iter2;

  for (iter2 = TimeTable.begin(); iter2 != TimeTable.end(); ++iter2)
  {
    TotalCPUTime += iter2->second;
  }

  return(TotalCPUTime);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
GetBlockDiagonal(string DiagonalColoringType)
{
  int OperatorDomainPoints = Operator_.OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int GraphBlockRows = Graph_.NumGlobalBlockRows();
  int GraphNnz = Graph_.NumGlobalNonzeros();
  int NumPDEEqns = OperatorRangePoints / GraphBlockRows;

  CrsGraph_MapColoring MapColoringTransform(CrsGraph_MapColoring::JONES_PLASSMAN,
                                            0, true, true);

  Epetra_MapColoring& ColorMap = MapColoringTransform(const_cast<Epetra_CrsGraph&>(Graph_));

  const int NumColors = ColorMap.MaxNumColors();

  Epetra_MultiVector X(Operator_.OperatorDomainMap(), NumPDEEqns * NumColors);
  X.PutScalar(0.0);

  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    int color = ColorMap[i] - 1;
    for (int j = 0; j < NumPDEEqns; ++j)
    {
      X[color * NumPDEEqns + j][i * NumPDEEqns + j] = 1.0;
    }
  }

  Epetra_MultiVector AX(Operator_.OperatorRangeMap(), NumPDEEqns * NumColors);

  Operator_.Apply(X, AX);

  InvBlockDiag_.resize(Operator_.OperatorRangeMap().NumMyElements() * NumPDEEqns);
  
  char job = 'A';
  vector<double> S(NumPDEEqns), U(NumPDEEqns * NumPDEEqns), VT(NumPDEEqns * NumPDEEqns);
  vector<double> WORK(5 * NumPDEEqns);
  int LWORK = 5 * NumPDEEqns, INFO;

  // extract the diagonals

  Epetra_SerialDenseMatrix V(NumPDEEqns, NumPDEEqns);
  Epetra_SerialDenseSVD SVD;
  SVD.SetMatrix(V);

  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    int color = ColorMap[i] - 1;
    int offset = i * NumPDEEqns * NumPDEEqns;

    // extract the block
    for (int j = 0; j < NumPDEEqns; ++j)
    {
      for (int k = 0; k < NumPDEEqns; ++k)
      {
        V(j, k) = AX[color * NumPDEEqns + j][i * NumPDEEqns + k];
      }
    }

    // invert the block
    SVD.Invert();
    
    // set the inverted block
    for (int j = 0; j < NumPDEEqns; ++j)
    {
      for (int k = 0; k < NumPDEEqns; ++k)
      {
        InvBlockDiag_[offset + j * NumPDEEqns + k] = (*SVD.InvertedMatrix())(j, k);
      }
    }
  }

  /* some possible output for debugging
  Epetra_MultiVector XXX(Copy, Operator_.OperatorRangeMap(), &InvBlockDiag_[0],
                         Operator_.OperatorRangeMap().NumMyElements(), NumPDEEqns);
  */
  return(0);
}

// ============================================================================ 
int ML_Epetra::MatrixFreePreconditioner::
ApplyInvBlockDiag(const double alpha, Epetra_MultiVector& X,
                  const double beta, const Epetra_MultiVector& B) const
{
  assert (X.NumVectors() == 1); // to be fixed
  // FIXME
  int OperatorRangePoints =  Operator_.OperatorRangeMap().NumGlobalPoints();
  int GraphBlockRows = Graph_.NumGlobalBlockRows();
  int NumPDEEqns = OperatorRangePoints / GraphBlockRows;
  int NumPDEEqns2 = NumPDEEqns * NumPDEEqns;

  char trans = 'N';
  int NumVectorsX = X.NumVectors();
  vector<double> tmp(NumPDEEqns);

  size_t len = sizeof(double) * NumPDEEqns;
  for (int i = 0; i < Graph_.NumMyBlockRows(); ++i)
  {
    memcpy(&tmp[0], &(B[0][i * NumPDEEqns]), len);

    int offset = i * NumPDEEqns2;
#if 0
    cout << InvBlockDiag_[offset] << " " << InvBlockDiag_[offset + 1] << endl;
    cout << InvBlockDiag_[offset + 2] << " " << InvBlockDiag_[offset + 3] << endl;
    cout << endl;
#endif

    DGEMM_F77(&trans, &trans, &NumPDEEqns, &NumVectorsX, &NumPDEEqns,
              (double*)&alpha, (double*)&InvBlockDiag_[offset], &NumPDEEqns, 
              &tmp[0], &NumPDEEqns, (double*)&beta, 
              (double*)&X[0][i * NumPDEEqns], &NumPDEEqns);
  }

  return(0);
}

#endif
