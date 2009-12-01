#include "TPICrsMatrix.h"

TPICrsMatrix::TPICrsMatrix(const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow)
: Epetra_BasicRowMatrix(RowMap.Comm())
, NumMyRows_(RowMap.NumMyPoints())
, IndexOffset_(NumMyRows_+1)
, Indices_(0)
, Values_(0)
{
  SetMaps(RowMap,ColMap,ColMap,RowMap);
  NumMyEntries_ = 0;
  for (int i=0; i<NumMyRows_; ++i) 
  {
    IndexOffset_[i] = NumMyEntries_;
    NumMyEntries_ += NumEntriesPerRow[i];
  }
  IndexOffset_[NumMyRows_] = NumMyEntries_;
  Indices_.Size(NumMyEntries_);  // this sizes and sets to zero
  Values_.Size(NumMyEntries_);
}

TPICrsMatrix::~TPICrsMatrix()
{}

int TPICrsMatrix::RNNZ(int MyRow) const
{ 
  // assumes that MyRow is in [0,NumMyRows_)
  return IndexOffset_[MyRow+1] - IndexOffset_[MyRow];
}

int TPICrsMatrix::InsertMyValues(int MyRow, int NumEntries, const double *values, const int *indices)
{
  if (MyRow < 0 || MyRow >= NumMyRows_) EPETRA_CHK_ERR(-1);
  if (NumEntries < 0 || NumEntries > RNNZ(MyRow)) EPETRA_CHK_ERR(-2);
  // copy into arrays, filter using column map
  const Epetra_Map &cmap = RowMatrixColMap();
  int entered = 0;
  bool filtered = false;
  double * rvals = &Values_[IndexOffset_[MyRow]];
  int    * rinds = &Indices_[IndexOffset_[MyRow]];
  for (int e=0; e<NumEntries; ++e) {
    if (cmap.MyLID(indices[e])) {
      rinds[entered] = indices[e];
      rvals[entered] = values[e];
      ++entered;  
    }
    else {
      filtered = true;
    }
  }
  // set the tail to zero,
  const int numinrow = RNNZ(MyRow);
  bool zeros = (entered < numinrow);
  while (entered < numinrow) {
    rinds[entered] = 0;
    rvals[entered] = 0.0;
    ++entered;
  }
  if (filtered) return 1;
  if (zeros)    return 2;
  return 0;
}

int TPICrsMatrix::FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap)
{
  SetMaps(RowMatrixRowMap(),RowMatrixColMap(),DomainMap,RangeMap);
  return 0;
}

int TPICrsMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  if (MyRow < 0 || MyRow >= NumMyRows_) EPETRA_CHK_ERR(-1);
  NumEntries = RNNZ(MyRow);
  return 0;
}

int TPICrsMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  if (MyRow < 0 || MyRow >= NumMyRows_) EPETRA_CHK_ERR(-1);
  const int numinrow = RNNZ(MyRow);
  if (Length < numinrow)                EPETRA_CHK_ERR(-2);
  NumEntries = numinrow;
  const double * rvals = &Values_[IndexOffset_[MyRow]];
  const int    * rinds = &Indices_[IndexOffset_[MyRow]];
  for (int i=0; i<numinrow; ++i) {
    Indices[i] = rinds[i];
    Values[i] = rvals[i];
  }
  return 0;
}

int TPICrsMatrix::ExtractMyEntryView(int CurEntry, double * & Value, int & RowIndex, int & ColIndex)
{
  if (CurEntry < 0 || CurEntry >= NumMyEntries_) EPETRA_CHK_ERR(-1);
  Value = &Values_[CurEntry];
  ColIndex = Indices_[CurEntry];
  for (int j=0; j<NumMyRows_; ++j) {
    if (CurEntry < IndexOffset_[j+1]) {
      RowIndex = j;
      break;
    }
  }
  return 0;
}

int TPICrsMatrix::ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const
{
  if (CurEntry < 0 || CurEntry >= NumMyEntries_) EPETRA_CHK_ERR(-1);
  Value = &Values_[CurEntry];
  ColIndex = Indices_[CurEntry];
  for (int j=0; j<NumMyRows_; ++j) {
    if (CurEntry < IndexOffset_[j+1]) {
      RowIndex = j;
      break;
    }
  }
  return 0;
}

int TPICrsMatrix::SetUseTranspose(bool use_transpose) 
{ return -1; }


int TPICrsMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{ 
  int ierr = Multiply(TransA,1.0,X,0.0,Y);
  EPETRA_CHK_ERR(ierr); 
  return ierr;
}


int TPICrsMatrix::Multiply(bool TransA, double alpha, const Epetra_MultiVector& X, double beta, Epetra_MultiVector& Y) const
{
  // TODO: currently does not support X==Y, must add copy of X
  if (TransA == true) EPETRA_CHK_ERR(-1);
  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-1); // Need same number of vectors in each MV
  }

  UpdateImportVector(NumVectors); // Make sure Import and Export Vectors are compatible
  UpdateExportVector(NumVectors);

  const double ** Xp = (const double**) X.Pointers();
  double ** Yp = (double**) Y.Pointers();

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if (Importer()!=0) {
    EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
    Xp = (const double**)ImportVector_->Pointers();
  }

  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  // Thus, computation is done into ExportVector_, instead of into Y
  double beta_lcl = ( Comm().MyPID() == 0 ? beta : 0.0 );
  if (Exporter()!=0) {
    Yp = (double**)ExportVector_->Pointers();
    beta_lcl = 0.0;
  }

  // Do actual computation
  TPICrsMatrix::tpi_crs_matrix_apply( NumMyRows_, IndexOffset_.Values(), Indices_.Values(), (const double *)Values_.Values(), 
                                      NumVectors, alpha, Xp, beta_lcl, Yp );

  if (Exporter()!=0) {
    Y.Scale(beta);
    Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
  }
  // Handle case of rangemap being a local replicated map
  if (!OperatorRangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

void TPICrsMatrix::tpi_work_crs_matrix_apply( TPI_Work * work )
{
  const struct tpi_crs_matrix * const h =
    (struct tpi_crs_matrix *) work->info ;

  const int    * const A_pc = h->A_pc ;
  const int    * const A_ia = h->A_ia ;
  const double * const A_a  = h->A_a ;
  const int     NumVectors = h->numVectors ;
  const double * const * const X = h->x ;
        double * const * const Y = h->y ;
  const double alpha = h->alpha;
  const double  beta = h->beta;

  const int nRow  = h->nRow ;
  const int chunk = ( nRow + work->count - 1 ) / work->count ;

  int rowEnd = chunk * ( work->rank + 1 );
  int row    = chunk * work->rank ;

  if ( nRow < rowEnd ) { rowEnd = nRow ; }

  for ( ; row < rowEnd ; ++row ) {
    const int jEnd = A_pc[ row + 1 ];
    for (int k = 0; k < NumVectors; ++k) {
      const double * const x = X[k];
      double * const y = Y[k];
      int j = A_pc[ row ];
      double tmp = 0 ;
      for ( ; j < jEnd ; ++j ) { tmp += A_a[j] * x[ A_ia[j] ]; }
      tmp = alpha*tmp + beta*y[row];
      y[ row ] = tmp;
    }
  }
}

void TPICrsMatrix::tpi_crs_matrix_apply(
  const int      nRow ,
  const int    * A_pc ,
  const int    * A_ia ,
  const double * A_a ,
  const int      numVectors ,
  double alpha ,
  const double ** x ,
  double beta ,
        double ** y )
{
  struct tpi_crs_matrix h = { 0 , NULL , NULL , NULL , 0, NULL , NULL };
  h.nRow = nRow ;
  h.A_pc = A_pc ;
  h.A_ia = A_ia ;
  h.A_a  = A_a ;
  h.numVectors = numVectors;
  h.alpha = alpha;
  h.x    = x ;
  h.beta = beta;
  h.y    = y ;
  TPI_Run_threads( tpi_work_crs_matrix_apply , & h , 0 );
}
