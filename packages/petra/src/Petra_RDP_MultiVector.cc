
#include "Petra_RDP_MultiVector.h"


//=============================================================================

// Petra_BlockMap Constructor

Petra_RDP_MultiVector::Petra_RDP_MultiVector(const Petra_BlockMap& Map, int NumVectors)
  : IndexBase_(Map.IndexBase()),
    Values_(0),
    Pointers_(0),
    Map_(&Map),
    Flops_(0.0),
    MyLength_(Map.NumMyEquations()),
    GlobalLength_(Map.NumGlobalEquations()),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(Map.NumMyEquations()),
    DistributedGlobal_(Map.DistributedGlobal()),
    Allocated_(false),
    Seed_(1),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Comm_(&Map.Comm())
{
    AllocateForCopy();
    
    for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Values_+i*Stride_;

  PutScalar(0.0); // Fill all vectors with zero.
}
//==========================================================================

// Copy Constructor

Petra_RDP_MultiVector::Petra_RDP_MultiVector(const Petra_RDP_MultiVector& Source)
  : IndexBase_(Source.IndexBase_),
    Values_(0),
    Pointers_(0),
    Map_(Source.Map_),
    Flops_(0.0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(Source.NumVectors_),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(0),
    DistributedGlobal_(Source.DistributedGlobal_),
    Allocated_(false),
    Seed_(Source.Seed_),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Comm_(Source.Comm_)
{
  AllocateForCopy();
  
  double ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[i];
  
  DoCopy();
  
}
//==========================================================================

// This constructor copies in or makes view of a standard Fortran array

Petra_RDP_MultiVector::Petra_RDP_MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
					     double *A, int MyLDA, int NumVectors)
  : IndexBase_(Map.IndexBase()),
    Values_(0),
    Pointers_(0),
    Map_(&Map),
    Flops_(0.0),
    MyLength_(Map.NumMyEquations()),
    GlobalLength_(Map.NumGlobalEquations()),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(Map.NumMyEquations()),
    DistributedGlobal_(Map.DistributedGlobal()),
    Allocated_(false),
    Seed_(1),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Comm_(&Map.Comm())
{
  
  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = A + i*MyLDA;
  
   if (CV==Copy) DoCopy();
   else DoView();
  
}

//==========================================================================

// This constructor copies in or makes view of a C/C++ array of pointer

Petra_RDP_MultiVector::Petra_RDP_MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
					      double **ArrayOfPointers, int NumVectors)
  : IndexBase_(Map.IndexBase()),
    Values_(0),
    Pointers_(0),
    Map_(&Map),
    Flops_(0.0),
    MyLength_(Map.NumMyEquations()),
    GlobalLength_(Map.NumGlobalEquations()),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(Map.NumMyEquations()),
    DistributedGlobal_(Map.DistributedGlobal()),
    Allocated_(false),
    Seed_(1),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Comm_(&Map.Comm())
{
  if (CV==Copy) AllocateForCopy();
  else AllocateForView();
  
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = ArrayOfPointers[i];
  
   if (CV==Copy) DoCopy();
   else DoView();
  
}

//==========================================================================

// This constructor copies or makes view of selected vectors, specified in Indices, 
// from an existing MultiVector

Petra_RDP_MultiVector::Petra_RDP_MultiVector(Petra_DataAccess CV, const Petra_RDP_MultiVector& Source, 
					     int *Indices, int NumVectors)
  : IndexBase_(Source.IndexBase_),
    Values_(0),
    Pointers_(0),
    Map_(Source.Map_),
    Flops_(0.0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(0),
    DistributedGlobal_(Source.DistributedGlobal_),
    Allocated_(false),
    Seed_(1),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Comm_(Source.Comm_)
{
  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  double ** Source_Pointers = Source.Pointers();  
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[Indices[i]];
  
   if (CV==Copy) DoCopy();
   else DoView();
  
}

//==========================================================================

// This interface copies or makes view of a range of vectors from an existing MultiVector

Petra_RDP_MultiVector::Petra_RDP_MultiVector(Petra_DataAccess CV, const Petra_RDP_MultiVector& Source, 
					     int StartIndex, int NumVectors)
  : IndexBase_(Source.IndexBase_),
    Values_(0),
    Pointers_(0),
    Map_(Source.Map_),
    Flops_(0.0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(0),
    DistributedGlobal_(Source.DistributedGlobal_),
    Allocated_(false),
    Seed_(1),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Comm_(Source.Comm_)
{
  
  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  double ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[StartIndex+i];

   if (CV==Copy) DoCopy();
   else DoView();
}

//=========================================================================
Petra_RDP_MultiVector::~Petra_RDP_MultiVector(){

  if (!Allocated_) return;
  if (!UserAllocated_) delete [] Values_;

  delete [] Pointers_;

  delete [] DoubleTemp_;
  delete [] IntTemp_;

   //cout << Map().Comm().MyPID() << " LenExports_ = " << LenExports_<< endl;
   //cout << Map().Comm().MyPID() << " LenImports_ = " << LenImports_<< endl;
   //cout << Map().Comm().MyPID() << " Exports_ = " << Exports_<< endl;
   //cout << Map().Comm().MyPID() << " Imports_ = " << Imports_<< endl;
  if (Exports_!=0) delete [] Exports_;
  if (Imports_!=0) delete [] Imports_;

}

//=========================================================================
int Petra_RDP_MultiVector::AllocateForCopy(void)
{
  
  if (Allocated_) return(0);
    
  if (NumVectors_<=0) 
    Petra_error("Number of vectors must be greater than zero", NumVectors_);

  Stride_ = Map_->NumMyEquations();
  Values_ = new double[Stride_ * NumVectors_];
  Pointers_ = new double *[NumVectors_];

  DoubleTemp_ = new double[NumVectors_];
  IntTemp_ = new int[NumVectors_];
  
  if (DistributedGlobal_)
    Seed_ = 2*Comm_->MyPID() + 1;
  else
    Seed_ = 1; // Replicated Local MultiVectors get the same seed on all PEs

  Allocated_ = true;
  UserAllocated_ = false;
  return(0);
}

//=========================================================================
int Petra_RDP_MultiVector::DoCopy(void)
{
  // On entry Pointers_ contains pointers to the incoming vectors.  These
  // pointers are the only unique piece of information for each of the 
  // constructors.

  // \interal { Optimization of this function can impact performance since it 
  //           involves a fair amount of memory traffic.}

  // On exit, Pointers_ is redefined to point to its own MultiVector vectors.

  for (int i = 0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = Values_+i*Stride_;
      Pointers_[i] = to;
      for (int j=0; j<MyLength_; j++) *to++ = *from++;
    }

  return(0);
}
//=========================================================================
int Petra_RDP_MultiVector::AllocateForView(void)
{
  
  if (NumVectors_<=0) 
    Petra_error("Number of vectors must be greater than zero", NumVectors_);
  
  Pointers_ = new double *[NumVectors_];
  
  DoubleTemp_ = new double[NumVectors_];
  IntTemp_ = new int[NumVectors_];

  if (DistributedGlobal_)
    Seed_ = 2*Comm_->MyPID() + 1;
  else
    Seed_ = 1; // Replicated Local MultiVectors get the same seed on all PEs

  Allocated_ = true;
  UserAllocated_ = true;
  
  return(0);
}

//=========================================================================
int Petra_RDP_MultiVector::DoView(void)
{
  // On entry Pointers_ contains pointers to the incoming vectors.  These
  // pointers are the only unique piece of information for each of the 
  // constructors.


  Values_ = Pointers_[0];

  if (NumVectors_==1) {
    Stride_ = Map_->NumMyEquations();
    ConstantStride_ = true;
    return(0);
  }

  // Remainder of code checks if MultiVector has regular stride

  Stride_ = Pointers_[1] - Pointers_[0];
  ConstantStride_ = false;

  for (int i = 1; i < NumVectors_-1; i++) if (Pointers_[i+1] - Pointers_[i] != Stride_) return(0);

  ConstantStride_ = true;

  return(0);
}
//=========================================================================
int Petra_RDP_MultiVector::Random(void)
{
  // Generate random numbers drawn from a uniform distribution on
  // the interval (-1,1) using a multiplicative congruential generator
  // with modulus 2^31 - 1.
  
  int i,j;
  const double a = 16807.0, BigInt=2147483647.0, DbleOne=1.0, DbleTwo=2.0;
  
  for(i=0; i < NumVectors_; i++)
    for (j=0; j<MyLength_; j++){
      Seed_ = fmod( a*Seed_, BigInt );
      Pointers_[i][j] = DbleTwo*(Seed_/BigInt)-DbleOne;
    }
  
  return(0);
  
}
 
//=========================================================================

// Extract a copy of a Petra_MultiVector.  Put in a user's Fortran-style array

int Petra_RDP_MultiVector::ExtractCopy(double *A, int MyLDA) const {
  if (NumVectors_>1 && Stride_ > MyLDA) return(-1); // LDA not big enough
  
  for (int i=0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = A + i*MyLDA;
      for (int j=0; j<MyLength_; j++) *to++ = *from++;
    }
  
  return(0);
}
      
      

//=========================================================================

// Extract a copy of a Petra_MultiVector.  Put in a user's array of pointers

int Petra_RDP_MultiVector::ExtractCopy(double **ArrayOfPointers) const {
  for (int i=0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = ArrayOfPointers[i];
      for (int j=0; j<MyLength_; j++) *to++ = *from++;
    }
  
  return(0);
}
      
      

//=========================================================================

// Extract a view of a Petra_MultiVector.  Set up a user's Fortran-style array

int Petra_RDP_MultiVector::ExtractView(double **A, int *MyLDA) const {
  if (!ConstantStride_) return(-1);  // Can't make a Fortran-style view if not constant stride
  *MyLDA = Stride_; // Set user's LDA
  *A = Values_; // Set user's value pointer
  return(0);
}
      
      

//=========================================================================

// Extract a view of a Petra_MultiVector.  Put in a user's array of pointers

int Petra_RDP_MultiVector::ExtractView(double ***ArrayOfPointers) const {
  *ArrayOfPointers = Pointers_;
  
  return(0);
}
      
      
//=========================================================================
int Petra_RDP_MultiVector::PutScalar(double Scalar) {

  // Fills RDP_MultiVector with the value Scalar **/

  for (int i = 0; i < NumVectors_; i++)
    for (int j=0; j<MyLength_; j++) Pointers_[i][j] = Scalar;
  return(0);
}

//=========================================================================
int Petra_RDP_MultiVector::Import(const Petra_RDP_MultiVector& A, const Petra_Import & Importer,
				  Petra_CombineMode CombineMode) {

  if (NumVectors_ != A.NumVectors()) return(-1);
  if (!Map_->SameAs(Importer.TargetMap())) return(-2);
  if ((A.Pointers()!=Pointers_) && (!A.Map().SameAs(Importer.SourceMap()))) return(-3);
  
  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumRemoteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int * ExportLIDs = Importer.ExportLIDs();
  int *RemoteLIDs = Importer.RemoteLIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int Nsend = NumVectors_ * Importer.NumSend();
  int Nrecv = NumVectors_ * Importer.NumRecv();

  return(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		  PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, 
		    LenExports_, Exports_, LenImports_, Imports_, 
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_RDP_MultiVector::Export(const Petra_RDP_MultiVector& A, const Petra_Export & Exporter,
				  Petra_CombineMode CombineMode) {

  if (NumVectors_ != A.NumVectors()) return(-1);
  if (!Map_->SameAs(Exporter.TargetMap())) return(-2);
  if ((A.Pointers()!=Pointers_) && (!A.Map().SameAs(Exporter.SourceMap()))) return(-3);
  
  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumRemoteIDs();
  int NumExportIDs = Exporter.NumExportIDs();
  int * ExportLIDs = Exporter.ExportLIDs();
  int *RemoteLIDs = Exporter.RemoteLIDs();
  int *PermuteToLIDs = Exporter.PermuteToLIDs();
  int *PermuteFromLIDs = Exporter.PermuteFromLIDs();
  int Nsend = NumVectors_ * Exporter.NumSend();
  int Nrecv = NumVectors_ * Exporter.NumRecv();

  return(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		  PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, LenExports_, Exports_, 
		  LenImports_, Imports_, 
#ifdef PETRA_MPI
            Exporter.GSPlan(), 
#endif
            false));
}
//=========================================================================
int Petra_RDP_MultiVector::Import(const Petra_RDP_MultiVector& A, const Petra_Export & Exporter,
				  Petra_CombineMode CombineMode) {

  if (NumVectors_ != A.NumVectors()) return(-1);
  if (!Map_->SameAs(Exporter.SourceMap())) return(-2);
  if ((A.Pointers()!=Pointers_) && (!A.Map().SameAs(Exporter.TargetMap()))) return(-3);
  
  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumExportIDs();
  int NumExportIDs = Exporter.NumRemoteIDs();
  int * ExportLIDs = Exporter.RemoteLIDs();
  int *RemoteLIDs = Exporter.ExportLIDs();
  int *PermuteToLIDs = Exporter.PermuteFromLIDs();
  int *PermuteFromLIDs = Exporter.PermuteToLIDs();
  int Nsend = NumVectors_ * Exporter.NumRecv();
  int Nrecv = NumVectors_ * Exporter.NumSend();

  return(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs, 
		  PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, LenImports_, Imports_, 
		  LenExports_, Exports_, 
#ifdef PETRA_MPI
            Exporter.GSPlan(), 
#endif
            true));
}
//=========================================================================
int Petra_RDP_MultiVector::Export(const Petra_RDP_MultiVector& A, const Petra_Import & Importer,
				  Petra_CombineMode CombineMode) {

  if (NumVectors_ != A.NumVectors()) return(-1);
  if (!Map_->SameAs(Importer.SourceMap())) return(-2);
  if ((A.Pointers()!=Pointers_) && (!A.Map().SameAs(Importer.TargetMap()))) return(-3);
  
  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumExportIDs();
  int NumExportIDs = Importer.NumRemoteIDs();
  int * ExportLIDs = Importer.RemoteLIDs();
  int *RemoteLIDs = Importer.ExportLIDs();
  int *PermuteToLIDs = Importer.PermuteFromLIDs();
  int *PermuteFromLIDs = Importer.PermuteToLIDs();
  int Nsend = NumVectors_ * Importer.NumRecv();
  int Nrecv = NumVectors_ * Importer.NumSend();

  return(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,  
		    PermuteToLIDs, PermuteFromLIDs,  RemoteLIDs, ExportLIDs, Nsend, Nrecv, 
		    LenImports_, Imports_, LenExports_, Exports_, 
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_RDP_MultiVector::DoTransfer(const Petra_RDP_MultiVector& A, Petra_CombineMode CombineMode, 
				      int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, int NumExportIDs, 
				      int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, 
				      int * ExportLIDs,
				      int Nsend, int Nrecv,
				      int & LenExports, double * & Exports, int & LenImports, 
				      double * & Imports, 
#ifdef PETRA_MPI
				      GSComm_Plan & Plan, 
#endif
				      bool DoReverse){
  int ierr = 0;
  double **A_Pointers = A.Pointers();

  int * FirstElementEntryList = 0;
  int * A_FirstElementEntryList = 0;
  int * ElementSizeList = 0;
  int * A_ElementSizeList = 0;
  if (!  Map().ConstantElementSize()) {
    FirstElementEntryList =   Map().FirstElementEntryList();
    ElementSizeList =   Map().ElementSizeList();
  }
  if (!A.Map().ConstantElementSize()) {
    A_FirstElementEntryList = A.Map().FirstElementEntryList();
    A_ElementSizeList = A.Map().ElementSizeList();
  }

  ierr = 
    CopyAndPermute(Pointers_, A_Pointers, NumVectors_, A.Map().MaxElementSize(), 
		   A.Map().ConstantElementSize(),
		   A_ElementSizeList, A_FirstElementEntryList,
		   FirstElementEntryList, 
		   NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs);

  if (ierr!=0) return(ierr);
  if (CombineMode==Zero) return(0); // All done if CombineMode only involves copying and permuting
  
  if (Nsend>LenExports) {
    if (LenExports>0) delete Exports;
    Exports = new double[Nsend];
    LenExports = Nsend;
  }
  
  ierr = 
    Pack(A_Pointers, NumVectors_, A.Map().MaxElementSize(), A.Map().ConstantElementSize(),
	 A_ElementSizeList, A_FirstElementEntryList, 
	 NumExportIDs, ExportLIDs, Exports);

  if (ierr!=0) return(ierr);
  
  if (Nrecv>LenImports) {
    if (LenImports>0) delete [] Imports;
    Imports = new double[Nrecv];
    LenImports = Nrecv;
  }

#ifdef PETRA_MPI
  if (DistributedGlobal_) {
    int msgtag = 32765;
    
    GSComm_Comm GSComm;
    
    int SizeOfPacket = NumVectors_ * Map().MaxElementSize(); 
    bool GSComm_OK;
    if (DoReverse)
      // Do the exchange of remote data
      GSComm_OK = GSComm.DoReverse( Plan, msgtag, 
				    reinterpret_cast<char *> (Exports), 
				    SizeOfPacket * sizeof( double ),
				    reinterpret_cast<char *> (Imports) );
    
    else
      // Do the exchange of remote data
      GSComm_OK = GSComm.Do( Plan, msgtag, 
			     reinterpret_cast<char *> (Exports), 
			     SizeOfPacket * sizeof( double ),
			     reinterpret_cast<char *> (Imports) );
    if (!GSComm_OK) return(-1);
    
    ierr = 
      UnpackAndCombine(Pointers_, NumVectors_, Map().MaxElementSize(), Map().ConstantElementSize(),
		       ElementSizeList, FirstElementEntryList, NumRemoteIDs, RemoteLIDs, 
		       Imports, CombineMode);
  }
#endif
  return(ierr);
}

//=========================================================================
int Petra_RDP_MultiVector::CopyAndPermute(double ** To, double ** From, int NumVectors, 
					      int MaxElementSize, bool ConstantElementSize, 
					      int * FromElementSizeList, 
					      int * FromFirstElementEntryList, int * ToFirstElementEntryList,
					      int NumSameIDs, 
					      int NumPermuteIDs, int * PermuteToLIDs, int *PermuteFromLIDs) {

  int i, j, jj, jjj, k;
  
  int NumSameEntries;

  bool Case1 = false;
  bool Case2 = false;
  // bool Case3 = false;

  if (MaxElementSize==1) {
    Case1 = true;
    NumSameEntries = NumSameIDs;
  }
  else if (ConstantElementSize) {
    Case2 = true;
    NumSameEntries = NumSameIDs * MaxElementSize;
  }
  else {
    // Case3 = true;
    NumSameEntries = FromFirstElementEntryList[NumSameIDs];
  }

  // Short circuit for the case where the source and target vector is the same.
  if (To==From) NumSameEntries = 0;
  
  // Do copy first
  if (NumSameIDs>0)
    if (To!=From) {
      for (i=0; i < NumVectors; i++)
	for (j=0; j<NumSameEntries; j++)
	  To[i][j] = From[i][j];
    }
  // Do local permutation next
  if (NumPermuteIDs>0) {
  
    // Point entry case
    if (Case1) {
      
      if (NumVectors==1)
	for (j=0; j<NumPermuteIDs; j++) 
	  To[0][PermuteToLIDs[j]] = From[0][PermuteFromLIDs[j]];
      
      else {
	for (j=0; j<NumPermuteIDs; j++) {
	  jj = PermuteToLIDs[j];
	  jjj = PermuteFromLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    To[i][jj] = From[i][jjj];
	}
      }
    }
    // constant block size case
    else if (Case2) {
      
      for (j=0; j<NumPermuteIDs; j++) {
	jj = MaxElementSize*PermuteToLIDs[j];
	jjj = MaxElementSize*PermuteFromLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    To[i][jj+k] = From[i][jjj+k];
      }
    }
    
    // variable block size case
    else {
      
      for (j=0; j<NumPermuteIDs; j++) {
	jj = ToFirstElementEntryList[PermuteToLIDs[j]];
	jjj = FromFirstElementEntryList[PermuteFromLIDs[j]];
	int ElementSize = FromElementSizeList[PermuteFromLIDs[j]];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    To[i][jj+k] = From[i][jjj+k];
      }
    }
  }
  return(0);
}

//=========================================================================
int Petra_RDP_MultiVector::Pack(double ** From, int NumVectors, 
				int MaxElementSize, bool ConstantElementSize, 
				int * FromElementSizeList, 
				int * FromFirstElementEntryList,
				int NumSendIDs, int * SendLIDs, double * Sends) {

  int i, j, jj, k;
  // int Nsend;
  
  double * ptr;

  if (NumSendIDs>0) {
    ptr = Sends;
    
    
    // Point entry case
    if (MaxElementSize==1) {
      
      if (NumVectors==1) for (j=0; j<NumSendIDs; j++) *ptr++ = From[0][SendLIDs[j]];

      else {
	for (j=0; j<NumSendIDs; j++) {
	  jj = SendLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    *ptr++ = From[i][jj];
	}
      }
    }

    // constant block size case
    else if (ConstantElementSize) {
      
      for (j=0; j<NumSendIDs; j++) {
	jj = MaxElementSize*SendLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    *ptr++ = From[i][jj+k];
      }
    }
    
    // variable block size case
    else {
      
      int SizeOfPacket = NumVectors*MaxElementSize;
      for (j=0; j<NumSendIDs; j++) {
	ptr = Sends + j*SizeOfPacket;
	jj = FromFirstElementEntryList[SendLIDs[j]];
	int ElementSize = FromElementSizeList[SendLIDs[j]];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    *ptr++ = From[i][jj+k];
      }
    }
    // assert(ptr-Sends==Nsend); // Sanity check on send count
  }

  return(0);
}

//=========================================================================
int Petra_RDP_MultiVector::UnpackAndCombine(double ** To, int NumVectors, 
					int MaxElementSize, bool ConstantElementSize, 
					int * ToElementSizeList, int * ToFirstElementEntryList,
					int NumRecvIDs, int * RecvLIDs, 
					double * Recvs, Petra_CombineMode CombineMode ) {
  int i, j, jj, k;
  
  double * ptr;
  // Unpack it...

  if (NumRecvIDs>0) {

    ptr = Recvs;
    
    // Point entry case
    if (MaxElementSize==1) {
      
      if (NumVectors==1) {
	if (CombineMode==Add)
	  for (j=0; j<NumRecvIDs; j++) To[0][RecvLIDs[j]] += *ptr++; // Add to existing value
	else if(CombineMode==Insert)
	  for (j=0; j<NumRecvIDs; j++) To[0][RecvLIDs[j]] = *ptr++;
	// Note:  The following form of averaging is not a true average if more that one value is combined.
	//        This might be an issue in the future, but we leave this way for now.
	else // if(CombineMode==Average)
	  for (j=0; j<NumRecvIDs; j++) {To[0][RecvLIDs[j]] += *ptr++; To[0][RecvLIDs[j]] *= 0.5;}
      }

      else {  // NumVectors>1

	if (CombineMode==Add) {
	  for (j=0; j<NumRecvIDs; j++) {
	    jj = RecvLIDs[j];
	    for (i=0; i<NumVectors; i++)
	      To[i][jj] += *ptr++; // Add to existing value
	  }
	}
	else if(CombineMode==Insert) {
	  for (j=0; j<NumRecvIDs; j++) {
	    jj = RecvLIDs[j];
	    for (i=0; i<NumVectors; i++)
	      To[i][jj] = *ptr++;
	  }
	}
	// Note:  The following form of averaging is not a true average if more that one value is combined.
	//        This might be an issue in the future, but we leave this way for now.
	else { // if (CombineMode==Average)
	  for (j=0; j<NumRecvIDs; j++) {
	    jj = RecvLIDs[j];
	    for (i=0; i<NumVectors; i++)
	      { To[i][jj] += *ptr++;  To[i][jj] *= 0.5;}
	  }
	}
      }
    }

    // constant block size case

    else if (ConstantElementSize) {
   
      if (CombineMode==Add) {
	for (j=0; j<NumRecvIDs; j++) {
	  jj = MaxElementSize*RecvLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    for (k=0; k<MaxElementSize; k++)
	      To[i][jj+k] += *ptr++; // Add to existing value
	}
      }
      else if(CombineMode==Insert) {
	for (j=0; j<NumRecvIDs; j++) {
	  jj = MaxElementSize*RecvLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    for (k=0; k<MaxElementSize; k++)
	      To[i][jj+k] = *ptr++;
	}
      }
	// Note:  The following form of averaging is not a true average if more that one value is combined.
	//        This might be an issue in the future, but we leave this way for now.
	else { // if (CombineMode==Average)
	for (j=0; j<NumRecvIDs; j++) {
	  jj = MaxElementSize*RecvLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    for (k=0; k<MaxElementSize; k++)
	      { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}
	}
      }
    }
    
    // variable block size case

    else {
      
      int SizeOfPacket = NumVectors*MaxElementSize;

      if (CombineMode==Add) {
	for (j=0; j<NumRecvIDs; j++) {
	  ptr = Recvs + j*SizeOfPacket;
	  jj = ToFirstElementEntryList[RecvLIDs[j]];
	  int ElementSize = ToElementSizeList[RecvLIDs[j]];
	  for (i=0; i<NumVectors; i++)
	    for (k=0; k<ElementSize; k++)
	      To[i][jj+k] += *ptr++; // Add to existing value
	}
      }
      else  if(CombineMode==Insert){
	for (j=0; j<NumRecvIDs; j++) {
	  ptr = Recvs + j*SizeOfPacket;
	  jj = ToFirstElementEntryList[RecvLIDs[j]];
	  int ElementSize = ToElementSizeList[RecvLIDs[j]];
	  for (i=0; i<NumVectors; i++)
	    for (k=0; k<ElementSize; k++)
	      To[i][jj+k] = *ptr++;
	}
      }
	// Note:  The following form of averaging is not a true average if more that one value is combined.
	//        This might be an issue in the future, but we leave this way for now.
	else { // if (CombineMode==Average)
	for (j=0; j<NumRecvIDs; j++) {
	  ptr = Recvs + j*SizeOfPacket;
	  jj = ToFirstElementEntryList[RecvLIDs[j]];
	  int ElementSize = ToElementSizeList[RecvLIDs[j]];
	  for (i=0; i<NumVectors; i++)
	    for (k=0; k<ElementSize; k++)
	      { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}
	}
      }
    }
  }
  
  return(0);
}

//=========================================================================
int Petra_RDP_MultiVector::Dot(const Petra_RDP_MultiVector& A, double *Result) const {

  // Dot product of two RDP_MultiVectors 

  int i;
  if (NumVectors_ != A.NumVectors()) return(-1);
  if (MyLength_ != A.MyLength()) return(-2);
  
  Petra_BLAS & Blas = *new Petra_BLAS();
  
  double **A_Pointers = A.Pointers();

  for (i=0; i < NumVectors_; i++) DoubleTemp_[i] = Blas.DOT(MyLength_, Pointers_[i], A_Pointers[i]);
  
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);

  delete &Blas;
  return(0);
}
//=========================================================================
int Petra_RDP_MultiVector::Abs(const Petra_RDP_MultiVector& A) {

  // this[i][j] = fabs(A[i][j])

  int i, j;
  if (NumVectors_ != A.NumVectors()) return(-1);
  if (MyLength_ != A.MyLength()) return(-2);

  double **A_Pointers = A.Pointers();

  for (i=0; i < NumVectors_; i++) 
    for (j=0; j < MyLength_; j++)
      Pointers_[i][j] = fabs(A_Pointers[i][j]);

  return(0);
}
//=========================================================================
int Petra_RDP_MultiVector::Reciprocal(const Petra_RDP_MultiVector& A) {

  // this[i][j] = 1.0/(A[i][j])

  int ierr = 0;
  int i, j;
  if (NumVectors_ != A.NumVectors()) return(-1);
  if (MyLength_ != A.MyLength()) return(-2);

  double **A_Pointers = A.Pointers();

  for (i=0; i < NumVectors_; i++) 
    for (j=0; j < MyLength_; j++) {
      double value = A_Pointers[i][j];
      // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)     
      if (fabs(value)<Petra_MinDouble) {
	if (value==0.0) ierr = 1;
	else if (ierr!=1) ierr = 2;
	Pointers_[i][j] = PETRA_SGN(value) * Petra_MaxDouble;
      }
      else
	Pointers_[i][j] = 1.0/value;
    }
  
  return(ierr);
}
  //=========================================================================
  int Petra_RDP_MultiVector::Scale (double Scalar) {

    // scales a RDP_MultiVector in place by a scalar

    Petra_BLAS & Blas = *new Petra_BLAS();
  
    for (int i = 0; i < NumVectors_; i++)
      Blas.SCAL(MyLength_, Scalar, Pointers_[i]);

    UpdateFlops(GlobalLength_*NumVectors_);
  delete &Blas;

    return(0);
  }

  //=========================================================================
  int Petra_RDP_MultiVector::Scale (double ScalarA, const Petra_RDP_MultiVector& A) {

    // scales a RDP_MultiVector by a scalar and put in the this:
    // this = ScalarA * A

  if (NumVectors_ != A.NumVectors()) return(-1);
  if (MyLength_ != A.MyLength()) return(-2);

    double **A_Pointers = (double**)A.Pointers();
    
    for (int i = 0; i < NumVectors_; i++)
      for (int j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j];

    UpdateFlops(GlobalLength_*NumVectors_);

    return(0);
  }

  //=========================================================================
  int Petra_RDP_MultiVector::Update(double ScalarA, const Petra_RDP_MultiVector& A, double Scalar) {

    int i, j;

    // linear combination of two RDP_MultiVectors: this = Scalar * this + ScalarA * A

  if (NumVectors_ != A.NumVectors()) return(-1);
  if (MyLength_ != A.MyLength()) return(-2);

    double **A_Pointers = (double**)A.Pointers();

    if (Scalar==0.0)
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j];
	UpdateFlops(GlobalLength_*NumVectors_);
      }
    else if (Scalar==1.0)
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = Pointers_[i][j] + ScalarA * A_Pointers[i][j];
	UpdateFlops(2*GlobalLength_*NumVectors_);
      }
    else if (ScalarA==1.0)
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = Scalar * Pointers_[i][j] + A_Pointers[i][j];
	UpdateFlops(2*GlobalLength_*NumVectors_);
      }
    else
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = Scalar * Pointers_[i][j] +
					                    ScalarA *  A_Pointers[i][j];
	UpdateFlops(3*GlobalLength_*NumVectors_);
      }

    return(0);
  }

//=========================================================================
int Petra_RDP_MultiVector::Update(double ScalarA, const Petra_RDP_MultiVector& A, 
				  double ScalarB, const Petra_RDP_MultiVector& B, double Scalar) {
  
  int i, j;
  
  // linear combination of three RDP_MultiVectors: 
  // this = Scalar * this + ScalarA * A + ScalarB * B
  
  if (ScalarA==0.0) return(Update(ScalarB, B, Scalar));
  if (ScalarB==0.0) return(Update(ScalarA, A, Scalar));
			   
  if (NumVectors_ != A.NumVectors() || NumVectors_ != B.NumVectors()) return(-1);
  if (MyLength_ != A.MyLength() || MyLength_ != B.MyLength()) return(-2);
  
    double **A_Pointers = (double**)A.Pointers();
    double **B_Pointers = (double**)B.Pointers();

    if (Scalar==0.0)
      {
	if (ScalarA==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =           A_Pointers[i][j] + 
						                ScalarB * B_Pointers[i][j];
	    UpdateFlops(2*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j] +
						                          B_Pointers[i][j];
	    UpdateFlops(2*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j] + 
						                ScalarB * B_Pointers[i][j];
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
      }
    else if (Scalar==1.0)
      {
	if (ScalarA==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] +=           A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarA * A_Pointers[i][j] +
						                           B_Pointers[i][j];
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarA * A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
      }
    else
      {
	if (ScalarA==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  Scalar *    Pointers_[i][j] +
						                           A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  Scalar *    Pointers_[i][j] +
						                 ScalarA * A_Pointers[i][j] +
						                           B_Pointers[i][j];
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  Scalar *    Pointers_[i][j] +
						                 ScalarA * A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(5*GlobalLength_*NumVectors_);
	  }
      }


    return(0);
  }

//=========================================================================
int  Petra_RDP_MultiVector::Norm1 (double* Result)  const {
  
  // 1-norm of each vector in RDP_MultiVector 
  Petra_BLAS & Blas = *new Petra_BLAS();
    
  int i;
  for (i=0; i < NumVectors_; i++) DoubleTemp_[i] = Blas.ASUM(MyLength_, Pointers_[i]);
  
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);

  delete &Blas;
  return(0);
}

//=========================================================================
int  Petra_RDP_MultiVector::Norm2 (double* Result)  const {
  
  // 2-norm of each vector in RDP_MultiVector 
  
  int i, j;
  for (i=0; i < NumVectors_; i++) 
    {
      double sum = 0.0;
      for (j=0; j < MyLength_; j++) sum += Pointers_[i][j] * Pointers_[i][j];
      DoubleTemp_[i] = sum;
    }
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  for (i=0; i < NumVectors_; i++) Result[i] = sqrt(Result[i]);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);
  
  return(0);
}

//=========================================================================
int  Petra_RDP_MultiVector::NormInf (double* Result)  const {
  
  // Inf-norm of each vector in RDP_MultiVector 
  
  Petra_BLAS & Blas = *new Petra_BLAS();
  
  int i, j;
  for (i=0; i < NumVectors_; i++) 
    {
      j = Blas.IAMAX(MyLength_, Pointers_[i]);
      DoubleTemp_[i] = fabs(Pointers_[i][j]);
    }
  Comm_->MaxAll(DoubleTemp_, Result, NumVectors_);
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  delete &Blas;
  
  return(0);
}

//=========================================================================
int  Petra_RDP_MultiVector::NormWeighted (const Petra_RDP_MultiVector& Weights, double* Result)  const {
  
  // Weighted 2-norm of each vector in RDP_MultiVector 

  // If only one vector in Weights, we assume it will be used as the weights for all vectors

  int i, j;
  bool OneW = false;
  if (Weights.NumVectors()==1) OneW = true;
  else if (NumVectors_ != Weights.NumVectors()) return(-1);

  if (MyLength_ != Weights.MyLength()) return(-2);

  double *W = Weights.Values();
  double **W_Pointers = Weights.Pointers();
  
  for (i=0; i < NumVectors_; i++) 
    {
      if (!OneW) W = W_Pointers[i]; // If Weights has the same number of vectors as this, use each weight vector
      double sum = 0.0;
      for (j=0; j < MyLength_; j++) {
        double tmp = Pointers_[i][j]/W[j];
        sum += tmp * tmp;
      }
      DoubleTemp_[i] = sum;
    }
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  double OneOverN = 1.0 / (double) GlobalLength_;
  for (i=0; i < NumVectors_; i++) Result[i] = sqrt(Result[i]*OneOverN);
  
  UpdateFlops(3*GlobalLength_*NumVectors_);
  
  return(0);
  }

//=========================================================================
int  Petra_RDP_MultiVector::MinValue (double* Result)  const {
  
  // Minimum value of each vector in RDP_MultiVector 
  
  int i, j;
  for (i=0; i < NumVectors_; i++) 
    {
      double MinVal = Pointers_[i][0];
      for (j=1; j< MyLength_; j++) MinVal = minfn(MinVal,Pointers_[i][j]); 
      DoubleTemp_[i] = MinVal;
    }
  Comm_->MinAll(DoubleTemp_, Result, NumVectors_);
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  
  return(0);
}

//=========================================================================
int  Petra_RDP_MultiVector::MaxValue (double* Result)  const {
  
  // Maximum value of each vector in RDP_MultiVector 
  
  int i, j;
  for (i=0; i < NumVectors_; i++) 
    {
      double MaxVal = Pointers_[i][0];
      for (j=1; j< MyLength_; j++) MaxVal = maxfn(MaxVal,Pointers_[i][j]); 
      DoubleTemp_[i] = MaxVal;
    }
  Comm_->MaxAll(DoubleTemp_, Result, NumVectors_);
  
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  
  return(0);
}

//=========================================================================
int  Petra_RDP_MultiVector::MeanValue (double* Result)  const {
  
  // Mean value of each vector in RDP_MultiVector 
  
  int i, j;
  double fGlobalLength = 1.0/maxfn((double) GlobalLength_, 1.0);
  
  for (i=0; i < NumVectors_; i++) 
    {
      double sum = 0.0;
      for (j=0; j < MyLength_; j++) sum += Pointers_[i][j];
      DoubleTemp_[i] = sum;
    }
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  for (i=0; i < NumVectors_; i++) Result[i] = Result[i]*fGlobalLength;
  
  UpdateFlops(GlobalLength_*NumVectors_);

  return(0);
}

  //=========================================================================
  int  Petra_RDP_MultiVector::Multiply (char TransA, char TransB, double ScalarAB, 
					const Petra_RDP_MultiVector& A, 
					const Petra_RDP_MultiVector& B,
					double Scalar ) {

    // This routine performs a variety of matrix-matrix multiply operations, interpreting
    // the Petra_MultiVector (this-aka C , A and B) as 2D matrices.  Variations are due to
    // the fact that A, B and C can be local replicated or global distributed
    // Petra_MultiVectors and that we may or may not operate with the transpose of 
    // A and B.  Possible cases are:

    //                                       Num
    //     OPERATIONS                        case  Notes
    // 1) C(local) = A^X(local) * B^X(local)  4   (X=Trans or Not, No comm needed) 
    // 2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)
    // 3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no comm needed)

    // Note that the following operations are not meaningful for 
    // 1D distributions:

    // 1) C(local) = A^T(distr) * B^T(distr)  1
    // 2) C(local) = A  (distr) * B^X(distr)  2
    // 3) C(distr) = A^X(local) * B^X(local)  4
    // 4) C(distr) = A^X(local) * B^X(distr)  4
    // 5) C(distr) = A^T(distr) * B^X(local)  2
    // 6) C(local) = A^X(distr) * B^X(local)  4
    // 7) C(distr) = A^X(distr) * B^X(local)  4
    // 8) C(local) = A^X(local) * B^X(distr)  4

    // Total of 32 case (2^5).

    
    //if (!ConstantStride_    ||
    //!A.ConstantStride() ||
    //!B.ConstantStride()    ) return(-1); // Return error

    // Check for compatible dimensions

    int A_nrows = (TransA=='T') ? A.NumVectors() : A.MyLength();
    int A_ncols = (TransA=='T') ? A.MyLength() : A.NumVectors();
    int B_nrows = (TransB=='T') ? B.NumVectors() : B.MyLength();
    int B_ncols = (TransB=='T') ? B.MyLength() : B.NumVectors();

    double Scalar_local = Scalar; // local copy of Scalar

    if (MyLength_        != A_nrows     ||
	A_ncols        != B_nrows     ||
	NumVectors_    != B_ncols  ) return(-2); // Return error

    bool A_is_local = (!A.DistributedGlobal());
    bool B_is_local = (!B.DistributedGlobal());
    bool C_is_local = (!DistributedGlobal_);
    bool Case1 = ( A_is_local &&  B_is_local &&  C_is_local);  // Case 1 above
    bool Case2 = (!A_is_local && !B_is_local &&  C_is_local && TransA=='T' );// Case 2
    bool Case3 = (!A_is_local &&  B_is_local && !C_is_local && TransA=='N');// Case 3
  
    // Test for meaningful cases

    if (Case1 || Case2 || Case3)
      {
	if (Scalar!=0.0 && Case2)
	  {
	    const int MyPID = Comm_->MyPID();
	    if (MyPID!=0) Scalar_local = 0.0;
	  }

        // Check if A, B, C have constant stride, if not then make temp copy (strided)

        Petra_RDP_MultiVector * A_tmp, * B_tmp, *C_tmp;
        if (!ConstantStride_) C_tmp = new Petra_RDP_MultiVector(*this);
        else C_tmp = this;
          
        if (!A.ConstantStride()) A_tmp = new Petra_RDP_MultiVector(A);
        else A_tmp = (Petra_RDP_MultiVector *) &A;
    
        if (!B.ConstantStride()) B_tmp = new Petra_RDP_MultiVector(B);
        else B_tmp = (Petra_RDP_MultiVector *) &B;
    	
    
	int m = MyLength_;
	int n = NumVectors_;
	int k = A_ncols;
	int lda = A_tmp->Stride();
	int ldb = B_tmp->Stride();
	int ldc = C_tmp->Stride();
	double *Ap = A_tmp->Values();
	double *Bp = B_tmp->Values();
	double *Cp = C_tmp->Values();
   
	Petra_BLAS & Blas = *new Petra_BLAS();
	Blas.GEMM(TransA, TransB,  m, n, k, ScalarAB,
		  Ap, lda, Bp, ldb, Scalar_local, Cp, ldc);
	delete &Blas;

      
	// FLOP Counts
	//                                       Num
	//     OPERATIONS                        case  Notes 
	// 1) C(local) = A^X(local) * B^X(local)  4   (X=Trans or Not, No comm needed)
	// 2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)      
	// 3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no comm needed)

	// For Case 1 we only count the local operations, since we are interested in serial
	// cost.  Computation on other processors is redundant.
	if (Case1)
	  {	    
	    UpdateFlops(2*m*n*k);
	    if (ScalarAB!=1.0) UpdateFlops(m*n);
	    if (Scalar==1.0) UpdateFlops(m*n); else if (Scalar!=0.0) UpdateFlops(2*m*n);
	  }
	else if (Case2)
	  {	    
	    UpdateFlops(2*m*n*A.GlobalLength());
	    if (ScalarAB!=1.0) UpdateFlops(m*n);
	    if (Scalar==1.0) UpdateFlops(m*n); else if (Scalar!=0.0) UpdateFlops(2*m*n);
	  }
	else
	  {
	    UpdateFlops(2*GlobalLength_*n*k);
	    if (ScalarAB!=1.0) UpdateFlops(GlobalLength_*n);
	    if (Scalar==1.0) UpdateFlops(GlobalLength_*n);
	    else if (Scalar!=0.0) UpdateFlops(2*GlobalLength_*n);
	  }

	// If A or B were not strided, dispose of extra copies.
	if (!A.ConstantStride()) delete A_tmp;
	if (!B.ConstantStride()) delete B_tmp;

	// If C was not strided, copy from strided version and delete
	if (!ConstantStride_) 
	  {
	    C_tmp->ExtractCopy(Pointers_);
	    delete C_tmp;
	  }

	// If Case 2 then sum up C and distribute it to all processors.

	if (Case2) return(Reduce());

	return(0);

      }
    else return(-3); // Return error: not a supported operation

  }


//=========================================================================
int Petra_RDP_MultiVector::Multiply(double ScalarAB, const Petra_RDP_MultiVector& A, const Petra_RDP_MultiVector& B,
		       double Scalar) {
  
  int i, j;
  
  // Hadamard product of two RDP_MultiVectors: 
  // this = Scalar * this + ScalarAB * A * B (element-wise)
  
  if (ScalarAB==0.0) return(Scale(Scalar));
			   
  if (A.NumVectors() != 1 && A.NumVectors() != B.NumVectors()) return(-1); // A must have one column or be the same as B.
  if (NumVectors_ != B.NumVectors()) return(-2);
  if (MyLength_ != A.MyLength() || MyLength_ != B.MyLength()) return(-3);
  
  double **A_Pointers = (double**)A.Pointers();
  double **B_Pointers = (double**)B.Pointers();

  int IncA = 1;
  if (A.NumVectors() == 1 ) IncA = 0;

    if (Scalar==0.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  A[j] * B_Pointers[i][j];
	  }
	  UpdateFlops(GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarAB * A[j] *
					      B_Pointers[i][j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
    }
    else if (Scalar==1.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] +=  A[j] * B_Pointers[i][j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarAB * A[j] *
					      B_Pointers[i][j];
	    }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
    }
    else { // if (Scalar!=1.0 && Scalar !=0 ) 
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  Scalar * Pointers_[i][j] + A[j] * B_Pointers[i][j];
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = Scalar * Pointers_[i][j] + ScalarAB * A[j] *
					      B_Pointers[i][j];
	    }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
    }
  return(0);
}
//=========================================================================
int Petra_RDP_MultiVector::ReciprocalMultiply(double ScalarAB, const Petra_RDP_MultiVector& A, const Petra_RDP_MultiVector& B,
		       double Scalar) {
  
  int i, j;
  
  // Hadamard product of two RDP_MultiVectors: 
  // this = Scalar * this + ScalarAB * B / A (element-wise)
  
  if (ScalarAB==0.0) return(Scale(Scalar));
			   
  if (A.NumVectors() != 1 && A.NumVectors() != B.NumVectors()) return(-1); // A must have one column or be the same as B.
  if (NumVectors_ != B.NumVectors()) return(-2);
  if (MyLength_ != A.MyLength() || MyLength_ != B.MyLength()) return(-3);
  
  double **A_Pointers = (double**)A.Pointers();
  double **B_Pointers = (double**)B.Pointers();

  int IncA = 1;
  if (A.NumVectors() == 1 ) IncA = 0;

    if (Scalar==0.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarAB * 
					      B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
    }
    else if (Scalar==1.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] +=  B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarAB *
					      B_Pointers[i][j] / A[j];
	    }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
    }
    else { // if (Scalar!=1.0 && Scalar !=0 ) 
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  Scalar * Pointers_[i][j] + B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = Scalar * Pointers_[i][j] + ScalarAB * 
					      B_Pointers[i][j] / A[j];
	    }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
    }
  return(0);
}
//=======================================================================
double*& Petra_RDP_MultiVector::operator [] (int index)  {
  
  //  Petra_RDP_MultiVector::operator [] --- return non-const reference 
  
  return(Pointers_[index]);
}

//========================================================================
Petra_RDP_MultiVector& Petra_RDP_MultiVector::operator = (const Petra_RDP_MultiVector& Source) {
  
  // Check for special case of this=Source
  if (this != &Source) Assign(Source);
  
  return(*this);
}

//=========================================================================
void Petra_RDP_MultiVector::Assign(const Petra_RDP_MultiVector& A) {
  
  if (NumVectors_ != A.NumVectors())
    Petra_error("Number of vectors incompatible in Assign", NumVectors_);
  if (MyLength_ != A.MyLength())
    Petra_error("Length of vectors incompatible in Assign", MyLength_);
  
  double ** A_Pointers = A.Pointers();
  for (int i = 0; i< NumVectors_; i++)
      for (int j=0; j<MyLength_; j++) Pointers_[i][j] = A_Pointers[i][j];
    return;    
  }

  //=========================================================================
  int  Petra_RDP_MultiVector::Reduce() {

    // Global reduction on each entry of a Replicated Local MultiVector

    int i, j;
    double * tmp = new double[MyLength_];
    for (i = 0; i < NumVectors_; i++)
      {
	for (j=0; j< MyLength_; j++) tmp[j] = Pointers_[i][j];
	Comm_->SumAll(tmp, Pointers_[i], MyLength_);
      }
    delete [] tmp;

    // UpdateFlops(0);  No serial Flops in this function
    return(0);
  }
// Non-member functions

ostream& operator << (ostream& os, const Petra_RDP_MultiVector& A)
{
  int MyPID = A.Map().Comm().MyPID();
  int NumProc = A.Map().Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumVectors = A.NumVectors();
      int MyLength = A.MyLength();
      int * MyGlobalElements = A.Map().MyGlobalElements();
      double ** A_Pointers = A.Pointers();
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);

      if (MyPID==0) {
	os.width(14);
	os <<  "     MyPID"; os << "    ";
	os.width(14);
	os <<  "      Global Index "; os << " ";
	for (int j = 0; j < NumVectors ; j++)
	  {   
	    os.width(20);
	    os <<  "Value  ";
	  }
	os << endl;
      }
      
      for (int i=0; i < MyLength; i++)
	{
	  os.width(14);
	  os <<  MyPID; os << "    ";
	  os.width(14);
	  os <<  MyGlobalElements[i]; os << "    ";
	  for (int j = 0; j < NumVectors ; j++)
	    {   
	      os.width(20);
	      os <<  A_Pointers[j][i];
	    }
	  os << endl;
	}
      os << flush;

      // Reset os flags

      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }

    // Do a few global ops to give I/O a chance to complete
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
  }
  return os;
}
