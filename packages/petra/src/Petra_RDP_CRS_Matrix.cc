#include "Petra_RDP_CRS_Matrix.h"

//==============================================================================
Petra_RDP_CRS_Matrix::Petra_RDP_CRS_Matrix(Petra_DataAccess CV, const Petra_Map& RowMap, int *NumEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyEquations()),
    CV_(CV)
{
  Graph_ = new Petra_CRS_Graph(CV, RowMap, NumEntriesPerRow);
  InitializeDefaults();
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_CRS_Matrix::Petra_RDP_CRS_Matrix(Petra_DataAccess CV, const Petra_Map& RowMap, int NumEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyEquations()),
    CV_(CV)
{
  Graph_ = new Petra_CRS_Graph(CV, RowMap, NumEntriesPerRow);
  InitializeDefaults();
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_CRS_Matrix::Petra_RDP_CRS_Matrix(Petra_DataAccess CV, const Petra_Map& RowMap, 
					   const Petra_Map& ColMap, int *NumEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyEquations()),
    CV_(CV)
{
  
  cerr << "Not Implemented." << endl; abort();
  Graph_ = new Petra_CRS_Graph(CV, RowMap, ColMap, NumEntriesPerRow);
  InitializeDefaults();
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_CRS_Matrix::Petra_RDP_CRS_Matrix(Petra_DataAccess CV, const Petra_Map& RowMap, 
					   const Petra_Map& ColMap, int NumEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
   NumMyRows_(RowMap.NumMyEquations()),
    CV_(CV)
{
  cerr << "Not Implemented." << endl; abort();
  Graph_ = new Petra_CRS_Graph(CV, RowMap, ColMap, NumEntriesPerRow);
  InitializeDefaults();
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_CRS_Matrix::Petra_RDP_CRS_Matrix(Petra_DataAccess CV, const Petra_CRS_Graph & Graph) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_((Petra_CRS_Graph*) &Graph),
    Allocated_(false),
    StaticGraph_(true),
    NumMyRows_(Graph.NumMyRows_),
    CV_(CV)
{
  InitializeDefaults();
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_CRS_Matrix::Petra_RDP_CRS_Matrix(const Petra_RDP_CRS_Matrix & Matrix) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(Matrix.Allocated_),
    StaticGraph_(false),
    Values_(0),
    All_Values_(0),
    NormInf_(-1.0),
    NormOne_(-1.0),
    NumMyRows_(Matrix.NumMyRows_),
    ImportVector_(0),
    LenImports_(0),
    LenExports_(0),
    Imports_(0),
    Exports_(0),
    IntImports_(0),
    IntExports_(0),
    CV_(Copy)
{
  Graph_ = new Petra_CRS_Graph(Matrix.Graph());
  int ierr = Allocate();
  for (int i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    for (int j=0; j< NumEntries; j++) Values_[i][j] = Matrix.Values_[i][j];
  }
}

//==============================================================================
void Petra_RDP_CRS_Matrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  Values_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  ImportVector_ = 0;
  LenImports_ = 0;
  LenExports_ = 0;
  Imports_ = 0;
  Exports_ = 0;
  IntImports_ = 0;
  IntExports_ = 0;

  NumEntriesPerRow_  = 0;
  NumAllocatedEntriesPerRow_ = 0;
  Indices_ = 0;

  return;
}

//==============================================================================
int Petra_RDP_CRS_Matrix::Allocate() {

  int i, j;
  
  // Set direct access pointers to graph info (needed for speed)
  NumEntriesPerRow_ = Graph_->NumIndicesPerRow();
  NumAllocatedEntriesPerRow_ = Graph_->NumAllocatedIndicesPerRow();
  Indices_ = Graph_->Indices();

  // Allocate Values array
  Values_ = new double*[NumMyRows_];

  // Allocate and initialize entries if we are copying data
  if (CV_==Copy) {
    for (i=0; i<NumMyRows_; i++) {
      int NumAllocatedEntries = NumAllocatedEntriesPerRow_[i];

      if (NumAllocatedEntries > 0) Values_[i] = new double[NumAllocatedEntries];
      else Values_[i] = 0;

      for (j=0; j< NumAllocatedEntries; j++) Values_[i][j] = 0.0; // Fill values with zero
    }
  }	 
  else {
    for (i=0; i<NumMyRows_; i++) {
      Values_[i] = 0;
    }
  }
    SetAllocated(true);
    return(0);
}
//==============================================================================
Petra_RDP_CRS_Matrix::~Petra_RDP_CRS_Matrix(){

  int i;

  if (CV_==Copy) {
    if (All_Values_!=0) delete [] All_Values_;
    else for (i=0; i<NumMyRows_; i++) if (NumAllocatedEntriesPerRow_[i] >0) delete Values_[i];
  }

  if (ImportVector_!=0) delete ImportVector_;
  ImportVector_=0;
    
  if (Imports_!=0) delete Imports_;
  Imports_=0;
  if (Exports_!=0) delete Exports_;
  Exports_=0;
  if (IntImports_!=0) delete IntImports_;
  IntImports_=0;
  if (IntExports_!=0) delete IntExports_;
  IntExports_=0;
    
  delete [] Values_;
  if (!StaticGraph()) delete Graph_; // We created the graph, so must delete it.

  NumMyRows_ = 0;
  
  Allocated_ = false;
}

//==============================================================================
int Petra_RDP_CRS_Matrix::PutScalar(double Scalar) 
{
  for (int i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    for (int j=0; j< NumEntries; j++) Values_[i][j] = Scalar;
  }
  return(0);
}
//==========================================================================
int Petra_RDP_CRS_Matrix::InsertGlobalValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (IndicesAreLocal()) return(-2); // Cannot insert global values into local graph
  if (IndicesAreContiguous()) return(-3); // Indices cannot be individually deleted and newed
  Graph_->SetIndicesAreGlobal(true);
  Row = Graph_->LRID(Row); // Find local row number for this global row index

  return(InsertValues(Row, NumEntries, Values, Indices));

}

//==========================================================================
int Petra_RDP_CRS_Matrix::InsertMyValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (IndicesAreGlobal()) return(-2); // Cannot insert global values into filled graph
  if (IndicesAreContiguous()) return(-3); // Indices cannot be individually deleted and new
  Graph_->SetIndicesAreLocal(true);

  return(InsertValues(Row, NumEntries, Values, Indices));

}

//==========================================================================
int Petra_RDP_CRS_Matrix::InsertValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (StaticGraph()) return(-2); // If the matrix graph is fully constructed, we cannot insert new values

  int j;
  double * tmp_Values;
  int ierr = 0;

  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range
    
  if (CV_==View) {
    if (Values_[Row]!=0) ierr = 2; // This row has be defined already.  Issue warning.
    Values_[Row] = Values;
  }
  else {
    
    int start = NumEntriesPerRow_[Row];
    int stop = start + NumEntries;
    int NumAllocatedEntries = NumAllocatedEntriesPerRow_[Row];
    if (stop > NumAllocatedEntries){
      if (NumAllocatedEntries==0) Values_[Row] = new double[NumEntries]; // Row was never allocated, so do it
      else {
	ierr = 1; // Out of room.  Must delete and allocate more space...
	tmp_Values = new double[stop];
	for (j=0; j< start; j++) tmp_Values[j] = Values_[Row][j]; // Copy existing entries
	delete [] Values_[Row]; // Delete old storage
	Values_[Row] = tmp_Values; // Set pointer to new storage
      }
    }
        
    for (j=start; j<stop; j++) Values_[Row][j] = Values[j-start];
  }

  Graph_->InsertIndices(Row, NumEntries, Indices); // Update graph

  return(ierr);

}

//==========================================================================
int Petra_RDP_CRS_Matrix::ReplaceGlobalValues(int Row, int NumEntries, double * Values, int *Indices) {

  int j;
  int ierr = 0;
  int Loc;

  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.

   Row = Graph_->LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindGlobalIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] = Values[j];
    else return(-2); // Value not found
  }

  return(ierr);
}

//==========================================================================
int Petra_RDP_CRS_Matrix::ReplaceMyValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (!IndicesAreLocal()) return(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc;

  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.

    
  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindMyIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] = Values[j];
    else return(-2); // Value not found
  }

  return(ierr);
}

//==========================================================================
int Petra_RDP_CRS_Matrix::SumIntoGlobalValues(int Row, int NumEntries, double * Values, int *Indices) {

  int j;
  int ierr = 0;
  int Loc;

  Row = Graph_->LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindGlobalIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] += Values[j];
    else return(-2); // Value not found
  }

  return(ierr);
}

//==========================================================================
int Petra_RDP_CRS_Matrix::SumIntoMyValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (!IndicesAreLocal()) return(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc;

  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindMyIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] += Values[j];
    else return(-2); // Value not found
  }

  return(ierr);
}

//==========================================================================
int Petra_RDP_CRS_Matrix::TransformToLocal() {
  return(TransformToLocal((Petra_BlockMap *) (&RowMap()), (Petra_BlockMap *) (&ColMap())));
}

//==========================================================================
int Petra_RDP_CRS_Matrix::TransformToLocal(Petra_BlockMap *DomainMap, Petra_BlockMap *RangeMap) {
  
  if (!StaticGraph()) Graph_->MakeIndicesLocal(*DomainMap, *RangeMap);
  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values
  if (!StaticGraph()) Graph_->TransformToLocal(DomainMap, RangeMap);


  return(0);
}

//==========================================================================
int Petra_RDP_CRS_Matrix::SortEntries() {

  if (!IndicesAreLocal()) return(-1);
  if (Sorted()) return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.

  
  for (int i=0; i<NumMyRows_; i++){

    double * Values = Values_[i];
    int NumEntries = NumEntriesPerRow_[i];
    int * Indices = Indices_[i];

    int n = NumEntries;
    int m = n/2;
    
    while (m > 0) {
      int max = n - m;
      for (int j=0; j<max; j++)
        {
	  for (int k=j; k>=0; k-=m)
            {
	      if (Indices[k+m] >= Indices[k])
		break;
	      double dtemp = Values[k+m];
	      Values[k+m] = Values[k];
	      Values[k] = dtemp;
	      int itemp = Indices[k+m];
	      Indices[k+m] = Indices[k];
	      Indices[k] = itemp;
            }
        }
      m = m/2;
    }
  }
  Graph_->SetSorted(true); // This also sorted the graph
  return(0);
}

//==========================================================================
int Petra_RDP_CRS_Matrix::MergeRedundantEntries() {

  int i, j, k;

  if (NoRedundancies()) return(0);
  if (!Sorted()) return(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal
  // Note:  This function assumes that SortEntries was already called.

  for (i=0; i<NumMyRows_; i++){
    int NumEntries = NumEntriesPerRow_[i];
    if (NumEntries>0) {
      double * const Values = Values_[i];
      int * const Indices = Indices_[i];
      int j0 = 0;
      int jj0 = Indices[j0];
      for (j=1; j<NumEntries; j++) {
	int jj = Indices[j];
	if (jj==jj0) {// Check if index is repeated
	  Values[j0] += Values[j];
	  for (k=j; k<NumEntries-1; k++) Values[k] = Values[k+1]; // Sum up values
	  NumEntries--;
	}
	else {
	  j0=j; // Redefine comparison index value
	  jj0=Indices[j0];
	}
      }
    }
  }
  
  return(Graph_->RemoveRedundantIndices()); // Remove redundant indices and then return


}

//==========================================================================
int Petra_RDP_CRS_Matrix::OptimizeStorage() {

  int i, j;

  if (StorageOptimized()) return(0); // Have we been here before?

  bool Contiguous = true; // Assume contiguous is true
  for (i=1; i<NumMyRows_; i++){
    int NumEntries = NumEntriesPerRow_[i];
    int NumAllocatedEntries = NumAllocatedEntriesPerRow_[i];
      
    // Check if NumEntries is same as NumAllocatedEntries and 
    // check if end of beginning of current row starts immediately after end of previous row.
    if ((NumEntries!=NumAllocatedEntries) || (Values_[i]!=Values_[i-1]+NumEntries)) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumEntries and NumAllocatedEntries
  //        for the last row could be different, but I don't think it matters.


  if ((CV_==View) && !Contiguous) return(-1);  // This is user data, it's not contiguous and we can't make it so.

  int ierr = Graph_->OptimizeStorage(); // Make sure graph has optimized storage
  if (ierr) return(ierr);

  if (Contiguous) return(0); // Everything is done.  Return

 // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
  int NumMyNonzeros = Graph_->NumMyNonzeros();

  // Allocate one big integer array for all index values
  All_Values_ = new double[NumMyNonzeros];
  
  // Set Entries_ to point into All_Entries_
  
  double * tmp = All_Values_;
  for (i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    for (j=0; j<NumEntries; j++) tmp[j] = Values_[i][j];
    if (Values_[i] !=0) delete [] Values_[i];
    Values_[i] = tmp;
    tmp += NumEntries;
  }
  
  
    return(0);
}
//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * Values,
					 int * Indices) const 
{

  int ierr = Graph_->ExtractGlobalRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) return(ierr);

  return(ExtractGlobalRowCopy(Row, Length, NumEntries, Values));
}
//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * Values,
					 int * Indices) const 
{

  int ierr = Graph_->ExtractMyRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) return(ierr);

  return(ExtractMyRowCopy(Row, Length, NumEntries, Values));
}
//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * Values) const 
{

  int Row0 = Graph_->RowMap().LID(Row); // Normalize row range

  return(ExtractMyRowCopy(Row0, Length, NumEntries, Values));
}


//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * Values) const 
{
  int j;

  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range

  NumEntries = NumEntriesPerRow_[Row];
  if (Length < NumEntries) return(-2); // Not enough space for copy. Needed size is passed back in NumEntries


  for(j=0; j<NumEntries; j++)Values[j] = Values_[Row][j];
  
  return(0);
}


//==============================================================================
int Petra_RDP_CRS_Matrix::ExtractDiagonalCopy(Petra_RDP_Vector & Diagonal) const {
	
  if (!Filled()) return(-1); // Can't get diagonal unless matrix is filled
  if (!RowMap().SameAs(Diagonal.Map())) return(-2); // Maps must be the same

  int Base = IndexBase();
  for(int i=0; i<NumMyRows_; i++){
    int Row = i + Base;
    int NumEntries = NumEntriesPerRow_[i];
    int * Indices = Indices_[i];
    Diagonal[i] = 0.0;
    for (int j=0; j<NumEntries; j++) {
      int Col = Indices[j];
      if (Row==Col) {
	Diagonal[j] = Values_[i][j];
	break;
      }
    }
  }
  return(0);
}
//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& Values, int *& Indices) const 
{

  int ierr = Graph_->ExtractGlobalRowView(Row, NumEntries, Indices);
  if (ierr) return(ierr);

  return(ExtractGlobalRowView(Row, NumEntries, Values));
}
//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractMyRowView(int Row, int & NumEntries, double *& Values, int *& Indices) const 
{

  int ierr = Graph_->ExtractMyRowView(Row, NumEntries, Indices);
  if (ierr) return(ierr);

  return(ExtractMyRowView(Row, NumEntries, Values));
}
//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& Values) const 
{

  int Row0 = Graph_->RowMap().LID(Row); // Normalize row range

  return(ExtractMyRowView(Row0, NumEntries, Values));
}

//==========================================================================
int Petra_RDP_CRS_Matrix::ExtractMyRowView(int Row, int & NumEntries, double *& Values) const 
{

  if (Row < 0 || Row >= NumMyRows_) return(-1); // Not in Row range

  NumEntries = NumEntriesPerRow_[Row];

  Values = Values_[Row];
  
  return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::Multiply(bool TransA, const Petra_RDP_Vector& x, Petra_RDP_Vector& y) const {
//
// This function forms the product y = A * x or y = A' * x
//

  if (!Filled()) return (-1); // Matrix must be filled.

  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double * xp = (double*)x.Values();
  double *yp = (double*)y.Values();
  int NumMyCols_ = NumMyCols();


  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),1); // Create import vector if needed
      ImportVector_->Import(x, *Importer(), Insert);
      xp = (double*)ImportVector_->Values();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*)ExportVector_->Values();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      double sum = 0.0;
      for (j=0; j < NumEntries; j++) sum += RowValues[j] * xp[RowIndices[j]];

      yp[i] = sum;

    }
    if (Exporter()!=0) y.Export(*ExportVector_, *Exporter(), Add); // Fill y with Values from export vector
  }

  else { // Transpose operation


    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),1); // Create Export vector if needed
      ExportVector_->Import(x, *Exporter(), Insert);
      xp = (double*)ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),1); // Create import vector if needed
      yp = (double*)ImportVector_->Values();
    }

    // Do actual computation

    for (i=0; i < NumMyCols_; i++) yp[i] = 0.0; // Initialize y for transpose multiply
        
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (j=0; j < NumEntries; j++) yp[RowIndices[j]] += RowValues[j] * xp[i];
    }
    if (Importer()!=0) y.Export(*ImportVector_, *Importer(), Add); // Fill y with Values from export vector
  }

    UpdateFlops(2*NumGlobalNonzeros());
    return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::Multiply(bool TransA, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const {
//
// This function forms the product Y = A * Y or Y = A' * X
//

  if (!Filled()) return (-1); // Matrix must be filled.

  int i, j, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();
  int NumMyCols_ = NumMyCols();


  // Need to better manage the Import and Export vectors:
  // - Need accessor functions
  // - Need to make the NumVector match (use a View to do this)
  // - Need to look at RightScale and ColSum routines too.

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()<NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),NumVectors); // Create import vector if needed
      ImportVector_->Import(X, *Importer(), Insert);
      Xp = (double**)ImportVector_->Pointers();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()<NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      Yp = (double**)ExportVector_->Pointers();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (k=0; k<NumVectors; k++) {
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
	Yp[k][i] = sum;
      }
    }
    if (Exporter()!=0) Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
  }
  else { // Transpose operation


    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()<NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      ExportVector_->Import(X, *Exporter(), Insert);
      Xp = (double**)ExportVector_->Pointers();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()<NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),NumVectors); // Create import vector if needed
      Yp = (double**)ImportVector_->Pointers();
    }

    // Do actual computation



        for (k=0; k<NumVectors; k++) 
	  for (i=0; i < NumMyCols_; i++) Yp[k][i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (k=0; k<NumVectors; k++) {
	for (j=0; j < NumEntries; j++) Yp[k][RowIndices[j]] += RowValues[j] * Xp[k][i];
      }
    }
    if (Importer()!=0) Y.Export(*ImportVector_, *Importer(), Add); // Fill Y with Values from export vector
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Petra_RDP_Vector& x, Petra_RDP_Vector& y) const {
//
// This function find y such that Ly = x or Uy = x or the transpose cases.
//

  if (!Filled()) return (-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) return (-2);
  if ((!Upper) && (!LowerTriangular())) return (-3);
  if ((!UnitDiagonal) && (NoDiagonal())) return (-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) return(-5); // Need each row to have a diagonal
      

  int i, j, j0;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  int NumMyCols_ = NumMyCols();

  // If upper, point to last row
  if (Upper) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Values += NumMyRows_-1;
  }
    
    double *xp = (double*)x.Values();
    double *yp = (double*)y.Values();

  if (!Trans) {

    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	double sum = 0.0;
	for (j=j0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[i] = xp[i] - sum;
	else yp[i] = (xp[i] - sum)/RowValues[0];

      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[i] = xp[i] - sum;
	else yp[i] = (xp[i] - sum)/RowValues[NumEntries];

      }
    }
  }

  // ***********  Transpose case *******************************

  else {

    if (xp!=yp) for (i=0; i < NumMyCols_; i++) yp[i] = xp[i]; // Initialize y for transpose solve
    
    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
    
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) yp[i] = yp[i]/RowValues[0];
	for (j=j0; j < NumEntries; j++) yp[RowIndices[j]] -= RowValues[j] * yp[i];
      }
    }
    else {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
    
      for (i=NumMyRows_-1; i >= 0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal)  yp[i] = yp[i]/RowValues[NumEntries];
	for (j=0; j < NumEntries; j++) yp[RowIndices[j]] -= RowValues[j] * yp[i];
      }
    }

  }
    UpdateFlops(2*NumGlobalNonzeros());
    return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const {
//
// This function find Y such that LY = X or UY = X or the transpose cases.
//

  if (!Filled()) return (-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) return (-2);
  if ((!Upper) && (!LowerTriangular())) return (-3);
  if ((!UnitDiagonal) && (NoDiagonal())) return (-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) return(-5); // Need each row to have a diagonal

  int i, j, j0, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double diag;
  int NumMyCols_ = NumMyCols();

  // If upper, point to last row
  if (Upper) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Values += NumMyRows_-1;
  }

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();

  if (!Trans) {
    

    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal) diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=j0; j < NumEntries; j++) sum += RowValues[j] * Yp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][i] = Xp[k][i] - sum;
	  else Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=0; j < NumEntries; j++) sum += RowValues[j] * Yp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][i] = Xp[k][i] - sum;
	  else Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
  }
  // ***********  Transpose case *******************************

  else {

    for (k=0; k<NumVectors; k++) 
      if (Yp[k]!=Xp[k]) for (i=0; i < NumMyRows_; i++) Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) diag = 1.0/RowValues[j0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  if (!UnitDiagonal) Yp[k][i] = Yp[k][i]*diag;
	  for (j=j0; j < NumEntries; j++) Yp[k][RowIndices[j]] -= RowValues[j] * Yp[k][i];
	}
      }
    }
    else {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      
      for (i=NumMyRows_-1; i>=0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal)  Yp[k][i] = Yp[k][i]/Xp[k][i];
	for (j=0; j < NumEntries; j++) Yp[k][RowIndices[j]] -= RowValues[j] * Yp[k][i];
      }
    }
  }
  
  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::InvRowSums(Petra_RDP_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//

  if (!Filled()) return (-1); // Matrix must be filled.
  if (!Graph().RangeMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the range of A
  int ierr = 0;
  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double * xp = (double*)x.Values();


  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    double * RowValues  = *Values++;
    double scale = 0.0;
    for (j=0; j < NumEntries; j++) scale += fabs(RowValues[j]);
    if (scale<Petra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Petra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  return(ierr);
}

//=============================================================================
//=============================================================================
int Petra_RDP_CRS_Matrix::InvColSums(Petra_RDP_Vector& x) const {
//
// Put inverse of the sum of absolute values of the jth column of A in x[j].
//

  if (!Filled()) return (-1); // Matrix must be filled.
  if (!Graph().DomainMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the domain of A
  
  double * xp = (double*)x.Values();
  Petra_RDP_Vector * x_tmp = 0;
  int NumMyCols_ = NumMyCols();
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (Importer()!=0) {
    x_tmp = new Petra_RDP_Vector(ImportMap()); // Create import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int ierr = 0;
  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  for (i=0; i < NumMyCols_; i++) xp[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    int *    ColIndices = *Indices++;
    double * RowValues  = *Values++;
    for (j=0; j < NumEntries; j++) xp[ColIndices[j]] += fabs(RowValues[j]);
  }

  if (Importer()!=0){
    x.Export(*x_tmp, *Importer(), Add); // Fill x with Values from import vector
    delete x_tmp;
    xp = (double*) x.Values();
  }
  // Invert values, don't allow them to get too large
  for (i=0; i < NumMyRows_; i++) {
    double scale = xp[i];
    if (scale<Petra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Petra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  return(ierr);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::LeftScale(const Petra_RDP_Vector& x) {
//
// This function scales the ith row of A by x[i].
//

  if (!Filled()) return (-1); // Matrix must be filled.
  if (!Graph().RangeMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the range of A

  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double * xp = (double*)x.Values();


  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    double * RowValues  = *Values++;
    double scale = xp[i];
    for (j=0; j < NumEntries; j++)  RowValues[j] *= scale;
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::RightScale(const Petra_RDP_Vector& x) {
//
// This function scales the jth row of A by x[j].
//

  if (!Filled()) return (-1); // Matrix must be filled.
  if (!Graph().DomainMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the domain of A

  double *xp = (double*)x.Values();
  Petra_RDP_MultiVector * x_tmp = 0;

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if (Importer()!=0) {
    x_tmp = new Petra_RDP_Vector(ImportMap()); // Create import vector if needed
    x_tmp->Import(x,*Importer(), Insert); // x_tmp will have all the values we need
    xp = (double*)x_tmp->Values();
  }

  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    int *    ColIndices = *Indices++;
    double * RowValues  = *Values++;
    double scale = xp[i];
    for (j=0; j < NumEntries; j++)  RowValues[j] *=  xp[ColIndices[j]];
  }
  if (x_tmp!=0) delete x_tmp;
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
double Petra_RDP_CRS_Matrix::NormInf() const {

  if (NormInf_>-1.0) return(NormInf_);

  int * NumEntriesPerRow = NumEntriesPerRow_;
  double ** Values = Values_;
  double Local_NormInf = 0.0;
  for (int i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++ ;
    double * RowValues  = *Values++;
    double sum = 0.0;
    for (int j=0; j < NumEntries; j++) sum += fabs(RowValues[j]);
    
    Local_NormInf = maxfn(Local_NormInf, sum);
  }
  Comm().MaxAll(&Local_NormInf, &NormInf_, 1);
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf_);
}
//=============================================================================
double Petra_RDP_CRS_Matrix::NormOne() const {

  if (NormOne_>-1.0) return(NormOne_);

  if (!Filled()) return (-1); // Matrix must be filled.

  Petra_RDP_Vector * x = new Petra_RDP_Vector(RowMap()); // Need temp vector for column sums
  
  double * xp = (double*)x->Values();
  Petra_RDP_MultiVector * x_tmp = 0;
  int NumMyCols_ = NumMyCols();
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (Importer()!=0) {
    x_tmp = new Petra_RDP_Vector(ImportMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int ierr = 0;
  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  for (i=0; i < NumMyCols_; i++) xp[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    int *    ColIndices = *Indices++;
    double * RowValues  = *Values++;
    for (j=0; j < NumEntries; j++) xp[ColIndices[j]] += fabs(RowValues[j]);
  }
  if (Importer()!=0) x->Export(*x_tmp, *Importer(), Add); // Fill x with Values from temp vector
  x->MaxValue(&NormOne_); // Find max
  if (x_tmp!=0) delete x_tmp;
  delete x;
  UpdateFlops(NumGlobalNonzeros());
  return(NormOne_);
}
//=========================================================================
int Petra_RDP_CRS_Matrix::Import(const Petra_RDP_CRS_Matrix& SourceMatrix, 
				 const Petra_Import & Importer, Petra_CombineMode CombineMode) {


  if (!RowMap().SameAs(Importer.TargetMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Importer.SourceMap())) return(-3);

  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumRemoteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int *ExportLIDs = Importer.ExportLIDs();
  int *RemoteLIDs = Importer.RemoteLIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumEntries() + 2;
  int Nsend = SizeOfPacket * Importer.NumSend();
  int Nrecv = SizeOfPacket * Importer.NumRecv();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenExports_, Exports_, IntExports_, LenImports_, Imports_, IntImports_,
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_RDP_CRS_Matrix::Export(const Petra_RDP_CRS_Matrix& SourceMatrix, 
				 const Petra_Export & Exporter, Petra_CombineMode CombineMode) {


  if (!RowMap().SameAs(Exporter.TargetMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Exporter.SourceMap())) return(-3);

  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumRemoteIDs();
  int NumExportIDs = Exporter.NumExportIDs();
  int *ExportLIDs = Exporter.ExportLIDs();
  int *RemoteLIDs = Exporter.RemoteLIDs();
  int *PermuteToLIDs = Exporter.PermuteToLIDs();
  int *PermuteFromLIDs = Exporter.PermuteFromLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumEntries() + 2;
  int Nsend = SizeOfPacket * Exporter.NumSend();
  int Nrecv = SizeOfPacket * Exporter.NumRecv();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenExports_, Exports_, IntExports_, LenImports_, Imports_, IntImports_,
#ifdef PETRA_MPI
		    Exporter.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_RDP_CRS_Matrix::Import(const Petra_RDP_CRS_Matrix& SourceMatrix, 
				 const Petra_Export & Exporter, Petra_CombineMode CombineMode) {


  if (!RowMap().SameAs(Exporter.SourceMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Exporter.TargetMap())) return(-3);

  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumExportIDs();
  int NumExportIDs = Exporter.NumRemoteIDs();
  int *ExportLIDs = Exporter.RemoteLIDs();
  int *RemoteLIDs = Exporter.ExportLIDs();
  int *PermuteToLIDs = Exporter.PermuteFromLIDs();
  int *PermuteFromLIDs = Exporter.PermuteToLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumEntries() + 2;
  int Nsend = SizeOfPacket * Exporter.NumRecv();
  int Nrecv = SizeOfPacket * Exporter.NumSend();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenImports_, Imports_, IntImports_, LenExports_, Exports_, IntExports_,
#ifdef PETRA_MPI
		    Exporter.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_RDP_CRS_Matrix::Export(const Petra_RDP_CRS_Matrix& SourceMatrix, 
				 const Petra_Import & Importer, Petra_CombineMode CombineMode) {


  if (!RowMap().SameAs(Importer.SourceMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Importer.TargetMap())) return(-3);

  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumExportIDs();
  int NumExportIDs = Importer.NumRemoteIDs();
  int *ExportLIDs = Importer.RemoteLIDs();
  int *RemoteLIDs = Importer.ExportLIDs();
  int *PermuteToLIDs = Importer.PermuteFromLIDs();
  int *PermuteFromLIDs = Importer.PermuteToLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumEntries() + 2;
  int Nsend = SizeOfPacket * Importer.NumRecv();
  int Nrecv = SizeOfPacket * Importer.NumSend();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenImports_, Imports_, IntImports_, LenExports_, Exports_, IntExports_,
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_RDP_CRS_Matrix::DoTransfer(const Petra_RDP_CRS_Matrix& SourceMatrix, 
				     Petra_CombineMode CombineMode,
				     int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, 
				     int NumExportIDs, 
				     int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, 
				     int * ExportLIDs,
				     int Nsend, int Nrecv, int SizeOfPacket,
				     int & LenExports, double * & Exports, int * & IntExports,
				     int & LenImports, double * & Imports, int * & IntImports,
#ifdef PETRA_MPI
				     GSComm_Plan & Plan, 
#endif
				     bool DoReverse){
  int ierr = 0;
  ierr = 
    CopyAndPermute(*this, SourceMatrix, 
		       NumSameIDs, 
		       NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs);

  if (ierr!=0) return(ierr);
  if (CombineMode==Zero) return(0); // All done if CombineMode only involves copying and permuting


  // Need to send global row ID, NumIndices for the row, then the indices.
  if (Nsend>LenExports) {
    if (LenExports>0) {
      delete [] Exports;
      delete [] IntExports;
    }
    Exports = new double[Nsend];
    IntExports = new int[Nsend];
    LenExports = Nsend;
  }

  ierr =  Pack(SourceMatrix, NumExportIDs, ExportLIDs, Exports, IntExports);

  if (ierr!=0) return(ierr);
  
  if (Nrecv>LenImports) {
    if (LenImports>0) {
      delete Imports;
      delete IntImports;
    }
    Imports = new double[Nrecv];
    IntImports = new int[Nrecv];
    LenImports = Nrecv;
  }
  
#ifdef PETRA_MPI
  if (RowMap().DistributedGlobal()) {
  int msgtag = 32765;
  bool GSComm_OK;
  GSComm_Comm GSComm;

  if (DoReverse)
  // Do the exchange of remote values data
    GSComm_OK = GSComm.DoReverse( Plan, msgtag, 
				  reinterpret_cast<char *> (Exports), 
				  SizeOfPacket * sizeof( double ),
				  reinterpret_cast<char *> (Imports) );
  else
    GSComm_OK = GSComm.Do( Plan, msgtag, 
			   reinterpret_cast<char *> (Exports), 
			   SizeOfPacket * sizeof( double ),
			   reinterpret_cast<char *> (Imports) );

  if (!GSComm_OK) return(-1);

  if (DoReverse)
    // Do the exchange of remote index data
  GSComm_OK = GSComm.DoReverse( Plan, msgtag, 
				reinterpret_cast<char *> (IntExports), 
				SizeOfPacket * sizeof( int ),
				reinterpret_cast<char *> (IntImports) );
  else
    GSComm_OK = GSComm.Do( Plan, msgtag, 
			   reinterpret_cast<char *> (IntExports), 
			   SizeOfPacket * sizeof( int ),
			   reinterpret_cast<char *> (IntImports) );
  
  if (!GSComm_OK) return(-1);

  ierr = 
    UnpackAndCombine( *this, SizeOfPacket, NumRemoteIDs, RemoteLIDs, Imports, IntImports, CombineMode);
  }
#endif

  return(ierr);
  
}
 
//=========================================================================
int Petra_RDP_CRS_Matrix::CopyAndPermute(Petra_RDP_CRS_Matrix & Target, 
					 const Petra_RDP_CRS_Matrix & Source,
					 int NumSameIDs, 
					 int NumPermuteIDs, int * PermuteToLIDs,
					 int *PermuteFromLIDs) {
  
  int i;
  
  int Row, NumEntries;
  int * Indices;
  double * Values;
  int FromRow, ToRow;
  
  // Do copy first
  if (NumSameIDs>0) {
    if (Source.IndicesAreLocal()) {
      int MaxNumEntries = Source.MaxNumEntries();
      Indices = new int[MaxNumEntries];  // Need some temporary space
      Values = new double[MaxNumEntries];  // Need some temporary space
      
      for (i=0; i<NumSameIDs; i++) {
	Row = Target.GRID(i);
	assert(Source.ExtractGlobalRowCopy(Row, MaxNumEntries, NumEntries, Values, Indices)==0); // Set pointers
	// Place into target matrix.  Depends on Petra_DataAccess copy/view and static/dynamic graph.
	if (Target.StaticGraph())
	  assert(Target.ReplaceGlobalValues(Row, NumEntries, Values, Indices)==0);
	else
	  assert(Target.InsertGlobalValues(Row, NumEntries, Values, Indices)==0); 
      }
      delete [] Values;
      delete [] Indices;
    }
    else { // Source.IndiceAreGlobal()
      for (i=0; i<NumSameIDs; i++) {
	Row = Target.GRID(i);
	assert(Source.ExtractGlobalRowView(Row, NumEntries, Values, Indices)==0); // Set pointers
	// Place into target matrix.  Depends on Petra_DataAccess copy/view and static/dynamic graph.
	if (Target.StaticGraph())
	  assert(Target.ReplaceGlobalValues(Row, NumEntries, Values, Indices)==0); 
	else
	  assert(Target.InsertGlobalValues(Row, NumEntries, Values, Indices)==0); 
      }
    }	
  }

  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (Source.IndicesAreLocal()) {
      int MaxNumEntries = Source.MaxNumEntries();
      Indices = new int[MaxNumEntries];  // Need some temporary space
      Values = new double[MaxNumEntries];  // Need some temporary space
      
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = Source.GRID(PermuteFromLIDs[i]);
	ToRow = Target.GRID(PermuteToLIDs[i]);
	assert(Source.ExtractGlobalRowCopy(FromRow, MaxNumEntries, NumEntries, Values, Indices)==0); // Set pointers
	// Place into target matrix.  Depends on Petra_DataAccess copy/view and static/dynamic graph.
	if (Target.StaticGraph())
	  assert(Target.ReplaceGlobalValues(ToRow, NumEntries, Values, Indices)==0);
	else
	  assert(Target.InsertGlobalValues(ToRow, NumEntries, Values, Indices)==0); 
      }
      delete [] Values;
      delete [] Indices;
    }
    else { // Source.IndiceAreGlobal()
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = Source.GRID(PermuteFromLIDs[i]);
	ToRow = Target.GRID(PermuteToLIDs[i]);
	assert(Source.ExtractGlobalRowView(FromRow, NumEntries, Values, Indices)==0); // Set pointers
	// Place into target matrix.  Depends on Petra_DataAccess copy/view and static/dynamic graph.
	if (Target.StaticGraph())
	  assert(Target.ReplaceGlobalValues(ToRow, NumEntries, Values, Indices)==0); 
	else
	  assert(Target.InsertGlobalValues(ToRow, NumEntries, Values, Indices)==0); 
      }
    }
  }	
    
  return(0);
}

//=========================================================================
int Petra_RDP_CRS_Matrix::Pack(const Petra_RDP_CRS_Matrix & Source,
			       int NumSendIDs, int * SendLIDs, double * Sends, 
			       int * IntSends) {
  
  int i;
  
  int NumEntries;
  int * Indices;
  double * Values;
  int FromRow;
  double * valptr;
  int * intptr;
  
  if (NumSendIDs>0) {

    // Each segment of Sends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source matrix
    // next int:  NumEntries, Number of indices in row.
    // next NumEntries: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumEntries - NumEntries ints) will be wasted but we need fixed
    //   sized segments for current communication routines.

    valptr = Sends;
    intptr = IntSends;
    int GlobalMaxNumEntries = Source.GlobalMaxNumEntries();
    int Inc = GlobalMaxNumEntries + 2;
    for (i=0; i<NumSendIDs; i++) {
      FromRow = Source.GRID(SendLIDs[i]);
      *intptr = FromRow;
      Values = valptr + 2; 
      Indices = intptr + 2;
      assert(Source.ExtractGlobalRowCopy(FromRow, GlobalMaxNumEntries, NumEntries, Values, Indices)==0); // Set pointers
      intptr[1] = NumEntries; // Load second slot of segment
      valptr += Inc; // Point to next segment
      intptr += Inc; // Point to next segment
    }
    // assert(ptr-Sends==Nsend); // Sanity check on send count

  }
    
  return(0);
}

//=========================================================================
int Petra_RDP_CRS_Matrix::UnpackAndCombine(Petra_RDP_CRS_Matrix & Target, int SizeOfPacket,
					   int NumRecvIDs, int * RecvLIDs, 
					   double * Recvs, int * IntRecvs, Petra_CombineMode CombineMode) {
  int NumEntries;
  int * Indices;
  double * Values;
  int ToRow;
  int i;
  
  double * valptr;
  int * intptr;
  // Unpack it...

  if (NumRecvIDs>0) {

    // Each segment of Sends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source matrix
    // next int:  NumEntries, Number of indices in row.
    // next NumEntries: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumEntries - NumEntries ints) will be 
    //  wasted but we need fixed sized segments for current communication routines.

    valptr = Recvs;
    intptr = IntRecvs;
    
    int Inc = SizeOfPacket;
    for (i=0; i<NumRecvIDs; i++) {
      ToRow = Target.GRID(RecvLIDs[i]);
      assert((intptr[0])==ToRow); // Sanity check
      NumEntries = intptr[1];
      Values = valptr + 2; 
      Indices = intptr + 2; 
      if (CombineMode==Add) {
	if (Target.StaticGraph())
	  // Replace any current values
	  assert(Target.SumIntoGlobalValues(ToRow, NumEntries, Values, Indices)==0);
	else
	  // Insert values
	  assert(Target.InsertGlobalValues(ToRow, NumEntries, Values, Indices)>=0);
      }
      else if (CombineMode==Insert) {
	if (Target.StaticGraph())
	  // Replace any current values
	  assert(Target.ReplaceGlobalValues(ToRow, NumEntries, Values, Indices)==0);
	else
	  // Insert values
	  assert(Target.InsertGlobalValues(ToRow, NumEntries, Values, Indices)>=0);
      }
      valptr += Inc; // Point to next segment
      intptr += Inc; // Point to next segment
    }
    // assert(ptr-Recvs==Nrecv); // Sanity check on receive count
  }
  
  return(0);
}

#ifdef PETRA_LEVELSCHEDULING

//=============================================================================
int Petra_RDP_CRS_Matrix::LevelSolve(bool Upper, bool Trans, bool UnitDiagonal, 
					   const Petra_RDP_Vector& x, Petra_RDP_Vector& y) {
//
// This function find y such that Ly = x or Uy = x or the transpose cases.
//

  if (!Filled()) return (-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) return (-2);
  if ((!Upper) && (!LowerTriangular())) return (-3);
  if ((!UnitDiagonal) && (NoDiagonal())) return (-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows())) return(-5); // Need each row to have a diagonal
     
  int i, j, j0, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  double *xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  int MyThreadID = Comm().MyThreadID(); // Get my thread ID

  int NumThreads = Graph_->NumThreads_;
  int ** ThreadStartRows = Graph_->ThreadStartRows_;
  int NumLevels = Graph_->NumLevels_;
  int * LevelOrder = Graph_->LevelOrder_;

  
  if (Trans) return(-10); // Level Scheduling can only work with Notranspose case

  for (k=0; k<NumLevels; k++) {
#pragma omp parallel for default(shared) private(kk,MyThreadStart,MyThreadStop,j0,i,ii,NumEntries,RowIndices,RowValues,sum,j)
    for (int kk = 0; kk<NumThreads; kk++) { // kludge for now
      MyThreadID = kk;
    int MyThreadStart = ThreadStartRows[k][MyThreadID];
    int MyThreadStop = ThreadStartRows[k][MyThreadID+1] - 1;
    
    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=MyThreadStart; i <= MyThreadStop; i++) {
	int ii = LevelOrder[i];
	int      NumEntries = NumEntriesPerRow[ii];
	int *    RowIndices = Indices[ii];
	double * RowValues  = Values[ii];
	double sum = 0.0;
	for (j=j0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[ii] = xp[ii] - sum;
	else yp[ii] = (xp[ii] - sum)/RowValues[0];
	
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=MyThreadStart; i<= MyThreadStop; i++) {
	int ii = LevelOrder[i];
	int      NumEntries = NumEntriesPerRow[ii] - j0;
	int *    RowIndices = Indices[ii];
	double * RowValues  = Values[ii];
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[ii] = xp[ii] - sum;
	else yp[ii] = (xp[ii] - sum)/RowValues[NumEntries];
	
      }
    }
    
    // Comm().NodeBarrier();
    }
  }
  UpdateFlops(2*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Petra_RDP_CRS_Matrix::LevelSolve(bool Upper, bool Trans, bool UnitDiagonal, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) {
//
// This function find Y such that LY = X or UY = X or the transpose cases.
//

  if (!Filled()) return (-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) return (-2);
  if ((!Upper) && (!LowerTriangular())) return (-3);
  if ((!UnitDiagonal) && (NoDiagonal())) return (-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) return(-5); // Need each row to have a diagonal

  int i, j, j0, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double diag;

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();

  int MyThreadID = Comm().MyThreadID(); // Get my thread ID

  int ** ThreadStartRows = Graph_->ThreadStartRows_;
  int NumLevels = Graph_->NumLevels_;
  int * LevelOrder = Graph_->LevelOrder_;

  if (Trans) return(-10); // Level Scheduling can only work with Notranspose case

  for (k=0; k<NumLevels; k++) {
    int MyThreadStart = ThreadStartRows[k][MyThreadID];
    int MyThreadStop = ThreadStartRows[k][MyThreadID+1] - 1;
    
    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=MyThreadStart; i <= MyThreadStop ; i++) {
	int ii = LevelOrder[i];
	int      NumEntries = NumEntriesPerRow[ii];
	int *    RowIndices = Indices[ii];
	double * RowValues  = Values[ii];
	if (!UnitDiagonal) diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=j0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][ii] = Xp[k][ii] - sum;
	  else Yp[k][ii] = (Xp[k][ii] - sum)*diag;
	}
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=MyThreadStart; i <= MyThreadStop ; i++) {
	int ii = LevelOrder[i];
	int      NumEntries = NumEntriesPerRow[ii] - j0;
	int *    RowIndices = Indices[ii];
	double * RowValues  = Values[ii];
	if (!UnitDiagonal) diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][ii] = Xp[k][ii] - sum;
	  else Yp[k][ii] = (Xp[k][ii] - sum)*diag;
	}
      }
    }
    
    Comm().NodeBarrier(); // At the end of each level all threads on a node should sync
  }
  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}
#endif

// Non-member functions

ostream& operator << (ostream& os, const Petra_RDP_CRS_Matrix& A) {
  int MyPID = A.RowMap().Comm().MyPID();
  int NumProc = A.RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);
      if (MyPID==0) {
	os <<  "\nNumber of Global Rows        = "; os << A.NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << A.NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << A.NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << A.NumGlobalNonzeros(); os << endl;
	os <<    "Global Maximum Num Entries   = "; os << A.GlobalMaxNumEntries(); os << endl;
	if (A.LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (A.UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (A.NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }

      os <<  "\nNumber of My Rows        = "; os << A.NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << A.NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << A.NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << A.NumMyNonzeros(); os << endl;
      os <<    "My Maximum Num Entries   = "; os << A.MaxNumEntries(); os << endl; os << endl;

      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    A.Comm().Barrier();
    A.Comm().Barrier();
    A.Comm().Barrier();
  }

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);
      int NumMyRows = A.NumMyRows();
      int MaxNumIndices = A.MaxNumEntries();
      int * Indices  = new int[MaxNumIndices];
      double * Values  = new double[MaxNumIndices];
      int NumIndices;
      int i, j;

      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Row Index ";
	os.width(10);
	os <<  "   Col Index ";
	os.width(20);
	os <<  "   Value     ";
	os << endl;
      }
      for (i=0; i<NumMyRows; i++) {
	int Row = A.GRID(i); // Get global row number
	A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Values, Indices);
	
	for (j = 0; j < NumIndices ; j++) {   
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  Row ; os << "    ";	
	  os.width(10);
	  os <<  Indices[j]; os << "    ";
	  os.width(20);
	  os <<  Values[j]; os << "    ";
	  os << endl;
	}
      }

      delete [] Indices;
      delete [] Values;
      
      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
  }

  return os;
}
