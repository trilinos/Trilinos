
#include "Petra_BlockMap.h"


//==============================================================================
// Petra_BlockMap constructor for a Petra-defined uniform linear distribution of constant block size elements.
Petra_BlockMap::Petra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Petra_Comm& Comm)
  : LID_(0),
    NumGlobalElements_(NumGlobalElements),
    MyGlobalElements_(0),
    FirstElementEntryList_(0),
    ElementSize_(ElementSize),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    MinMyElementSize_(ElementSize),
    MaxMyElementSize_(ElementSize),
    MinElementSize_(ElementSize),
    MaxElementSize_(ElementSize),
    ConstantElementSize_(true),
    LinearMap_(true),
    DistributedGlobal_(true)
{
  // Each processor gets roughly numGlobalEquations/p equations
  // This routine automatically defines a linear partitioning of a
  // map with numGlobalEquations across the processors
  // specified in the given Petra_Comm

    if (NumGlobalElements_ < 0) 
      {
	cerr << "Petra_BlockMap: ERROR, NumGlobalElements < 0. Aborting." << endl;
	abort();
      }
   if (ElementSize_ <= 0) 
    {
      cerr << "Petra_BlockMap: ERROR, ElementSize <= 0. Aborting." << endl;
      abort();
    }
  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  NumMyElements_ = NumGlobalElements_/NumProc;
  int remainder = NumGlobalElements_%NumProc;
  int start_index = MyPID*(NumMyElements_ + 1);

  if (MyPID<remainder) NumMyElements_++;
  else start_index -= (MyPID-remainder);

  NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
  NumMyEquations_ = NumMyElements_ * ElementSize_;

  MinMyElementSize_ = ElementSize_;
  MaxMyElementSize_ = ElementSize_;
  MinElementSize_ = ElementSize_;
  MaxElementSize_ = ElementSize_;

  MinAllGID_ = IndexBase_;
  MaxAllGID_ = MinAllGID_ + NumGlobalElements_ - 1;
  MinMyGID_ = start_index + IndexBase_;
  MaxMyGID_ = MinMyGID_ + NumMyElements_ - 1;
  MinLID_ = 0;
  MaxLID_ = MinLID_ + NumMyElements_ - 1;
  if (NumProc==1) DistributedGlobal_ = false;

  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}
//==============================================================================
// Petra_BlockMap constructor for a user-defined linear distribution of constant block size elements.
Petra_BlockMap::Petra_BlockMap(int NumGlobalElements, int NumMyElements, 
			       int ElementSize, int IndexBase, const Petra_Comm& Comm)
  : LID_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(NumMyElements),
    MyGlobalElements_(0),
    FirstElementEntryList_(0),
    ElementSize_(ElementSize),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    MinMyElementSize_(ElementSize),
    MaxMyElementSize_(ElementSize),
    MinElementSize_(ElementSize),
    MaxElementSize_(ElementSize),
    ConstantElementSize_(true),
    LinearMap_(true),
    DistributedGlobal_(true)
{
  // Each processor gets NumMyElements equations

  if (NumGlobalElements_ < -1) 
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements < -1. Aborting." << endl;
      abort();
    }
  if (NumMyElements_ < 0) 
    {
      cerr << "Petra_BlockMap: ERROR, NumMyElements < 0. Aborting." << endl;
      abort();
    }
  if (ElementSize_ <= 0) 
    {
      cerr << "Petra_BlockMap: ERROR, ElementSize <= 0. Aborting." << endl;
      abort();
    }

  // Get processor information

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();

   // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (NumGlobalElements==NumMyElements || NumProc==1)
    {
      NumGlobalElements_ = NumMyElements_;
      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumGlobalElements_ <<" or should be set to -1 for me to compute it. Aborting." << endl;
      abort();
    }
      if (NumProc==1) DistributedGlobal_ = false;
      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;
      
      MinAllGID_ = IndexBase_;
      MaxAllGID_ = MinAllGID_ + NumGlobalElements_ - 1;
      MinMyGID_ = IndexBase_;
      MaxMyGID_ = MinMyGID_ + NumMyElements_ - 1;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;
    }
  else if (NumProc > 1)
    {
      // Sum up all local element counts to get global count
     Comm_->SumAll(&NumMyElements_, &NumGlobalElements_, 1);

      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumGlobalElements_ <<" . Aborting." << endl;
      abort();
    }

      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;

      MinAllGID_ = IndexBase_;
      MaxAllGID_ = MinAllGID_ + NumGlobalElements_ - 1;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;

      // Use the ScanSum function to compute a prefix sum of the number of equations
      Comm_->ScanSum(&NumMyElements_, &MaxMyGID_, 1);

      int start_index = MaxMyGID_ - NumMyElements_;
      MinMyGID_ = start_index + IndexBase_;
      MaxMyGID_ = MinMyGID_ + NumMyElements_ - 1;
    }
  else
    {
      cerr << "Petra_BlockMap: INTERNAL ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumMyElements 
	   << "and NumProc = " << NumProc << ".  Should = 1. Aborting." << endl;
      abort();
    }
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}
//==============================================================================
// Petra_BlockMap constructor for a user-defined arbitrary distribution of constant block size elements.
Petra_BlockMap::Petra_BlockMap(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
			       int ElementSize, int IndexBase, const Petra_Comm& Comm)
  : LID_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(NumMyElements),
    MyGlobalElements_(MyGlobalElements),
    FirstElementEntryList_(0),
    ElementSize_(ElementSize),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    MinMyElementSize_(ElementSize),
    MaxMyElementSize_(ElementSize),
    MinElementSize_(ElementSize),
    MaxElementSize_(ElementSize),
    ConstantElementSize_(true),
    LinearMap_(false),
    DistributedGlobal_(true)
{
  int i;
  // Each processor gets NumMyElements equations

  if (NumGlobalElements_ < -1) 
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements < -1. Aborting." << endl;
      abort();
    }
  if (NumMyElements_ < 0) 
    {
      cerr << "Petra_BlockMap: ERROR, NumMyElements < 0. Aborting." << endl;
      abort();
    }
  if (ElementSize_ <= 0) 
    {
      cerr << "Petra_BlockMap: ERROR, ElementSize <= 0. Aborting." << endl;
      abort();
    }

  // Allocate storage for global index list information

  if (NumMyElements>0) MyGlobalElements_ = new int[NumMyElements];

  // Get processor information

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  if (NumMyElements>0) {
    // Compute min/max GID on this processor
    MinMyGID_ = MyGlobalElements[0];
    MaxMyGID_ = MyGlobalElements[0];
    for (i = 0; i < NumMyElements; i++)
      {
	MyGlobalElements_[i] = MyGlobalElements[i];
	MinMyGID_ = minfn(MinMyGID_,MyGlobalElements[i]);
	MaxMyGID_ = maxfn(MaxMyGID_,MyGlobalElements[i]);
      }
  }
  else {
    MinMyGID_ = IndexBase_;
    MaxMyGID_ = IndexBase_;
  }
    
    // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (NumGlobalElements==NumMyElements || NumProc==1)
    {
      NumGlobalElements_ = NumMyElements_;
      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumGlobalElements_ <<" or should be set to -1 for me to compute it. Aborting." << endl;
      abort();
    }
      if (NumProc==1) DistributedGlobal_ = false;
      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;
      
      MinAllGID_ = MinMyGID_;
      MaxAllGID_ = MaxMyGID_;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;
    }
  else if (NumProc > 1)
    {
      // Sum up all local element counts to get global count
      Comm_->SumAll(&NumMyElements_, &NumGlobalElements_, 1);
      // Check to see if user's value for NumGlobalElements is either -1
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	{
	  cerr << "Petra_BlockMap: ERROR, NumGlobalElements = " << NumGlobalElements
	       << ". Should = "<< NumGlobalElements_ <<" . Aborting." << endl;
	  abort();
	}
      
      NumGlobalEquations_ = NumGlobalElements_ * ElementSize_;
      NumMyEquations_ = NumMyElements_ * ElementSize_;
      
      MinLID_ = 0;
      MaxLID_ = maxfn(MinLID_ + NumMyElements_ - 1,MinLID_);
      
      // Use the Allreduce function to find min/max GID 
      int *tmp_send = new int[2];
      int *tmp_recv = new int[2];
      tmp_send[0] = - MinMyGID_; // Negative sign lets us do one reduction
      tmp_send[1] =   MaxMyGID_;
      Comm_->MaxAll(tmp_send, tmp_recv, 2);
      MinAllGID_ = - tmp_recv[0];
      MaxAllGID_ =   tmp_recv[1];
      delete [] tmp_send;
      delete [] tmp_recv;
      if (MinAllGID_ < IndexBase_)
	{
	  cerr << "Petra_BlockMap: ERROR, Minimum element index = " << MinAllGID_
	       << ". which is less than the Index Base = "<< IndexBase_ <<" . Aborting." << endl;
	  abort();
	}

    }
  else
    {
      cerr << "Petra_BlockMap: INTERNAL ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumMyElements 
	   << "and NumProc = " << NumProc << ".  Should = 1. Aborting." << endl;
      abort();
    }

  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}

//==============================================================================
// Petra_BlockMap constructor for a user-defined arbitrary distribution of variable block size elements.
Petra_BlockMap::Petra_BlockMap(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
			       int *ElementSizeList, int IndexBase, const Petra_Comm& Comm)
  : LID_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(NumMyElements),
    FirstElementEntryList_(0),
    ElementSize_(0),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(IndexBase),
    Comm_(&Comm),
    Directory_(0),
    ConstantElementSize_(false),
    LinearMap_(false),
    DistributedGlobal_(true)
{
  int i;
  // Each processor gets NumMyElements equations

  if (NumGlobalElements_ < -1) 
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements < -1. Aborting." << endl;
      abort();
    }
  if (NumMyElements_ < 0) 
    {
      cerr << "Petra_BlockMap: ERROR, NumMyElements < 0. Aborting." << endl;
      abort();
    }
  for (i=0; i<NumMyElements_; i++)
    if (ElementSizeList[i] <= 0) 
      {
	cerr << "Petra_BlockMap: ERROR, ElementSizeList[" <<i<<" <= 0. Aborting." << endl;
	abort();
      }
  
  // Allocate storage for global index list and element size information

  if (NumMyElements>0) {
    MyGlobalElements_ = new int[NumMyElements];
    ElementSizeList_ = new int[NumMyElements];
  }
  // Get processor information

  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  
  if (NumMyElements>0) {
    // Compute min/max GID and element size, number of equations on this processor
    MinMyGID_ = MyGlobalElements[0];
    MaxMyGID_ = MyGlobalElements[0];
    MinMyElementSize_ = ElementSizeList[0];
    MaxMyElementSize_ = ElementSizeList[0];
    NumMyEquations_ = 0;
    for (i = 0; i < NumMyElements; i++)
      {
	MyGlobalElements_[i] = MyGlobalElements[i];
	ElementSizeList_[i] = ElementSizeList[i];
	MinMyGID_ = minfn(MinMyGID_,MyGlobalElements[i]);
	MaxMyGID_ = maxfn(MaxMyGID_,MyGlobalElements[i]);
	MinMyElementSize_ = minfn(MinMyElementSize_,ElementSizeList[i]);
	MaxMyElementSize_ = maxfn(MaxMyElementSize_,ElementSizeList[i]);
	NumMyEquations_ += ElementSizeList[i];
      }
  }
  else {
    MinMyGID_ = IndexBase_;
    MaxMyGID_ = IndexBase_;
    MinMyElementSize_ = 1;
    MaxMyElementSize_ = 1;
    NumMyEquations_ = 0;
  }

  
  
  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (NumGlobalElements==NumMyElements || NumProc==1)
    {
      NumGlobalElements_ = NumMyElements_;
      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
    {
      cerr << "Petra_BlockMap: ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumGlobalElements_ <<" or should be set to -1 for me to compute it. Aborting." << endl;
      abort();
    }
      if (NumProc==1) DistributedGlobal_ = false;
      NumGlobalEquations_ = NumMyEquations_;
      
      MinAllGID_ = MinMyGID_;
      MaxAllGID_ = MaxMyGID_;
      MinLID_ = 0;
      MaxLID_ = MinLID_ + NumMyElements_ - 1;
      MinElementSize_ = MinMyElementSize_;
      MaxElementSize_ = MaxMyElementSize_;
    }
  else if (NumProc > 1)
    {
      // Sum up all local element and equation counts to get global counts
      int *tmp_send = new int[4];
      int *tmp_recv = new int[4];
      tmp_send[0] = NumMyElements_;
      tmp_send[1] = NumMyEquations_;
      Comm_->SumAll(tmp_send, tmp_recv, 2);
      NumGlobalElements_ =  tmp_recv[0];
      NumGlobalEquations_ = tmp_recv[1];

      // Check to see if user's value for NumGlobalElements is either -1 
      // (in which case we use our computed value) or matches ours.
      if (NumGlobalElements!=-1 && NumGlobalElements!=NumGlobalElements_)
	{
	  cerr << "Petra_BlockMap: ERROR, NumGlobalElements = " << NumGlobalElements
	       << ". Should = "<< NumGlobalElements_ <<" . Aborting." << endl;
	  abort();
	}
      
      MinLID_ = 0;
      MaxLID_ = maxfn(MinLID_ + NumMyElements_ - 1,MinLID_);
      
      // Use the MaxAll function to find min/max GID 
      tmp_send[0] = - MinMyGID_; // Negative signs lets us do one reduction
      tmp_send[1] =   MaxMyGID_;
      tmp_send[2] = - MinMyElementSize_;
      tmp_send[3] =   MaxMyElementSize_;

      Comm_->MaxAll(tmp_send, tmp_recv, 4);

      MinAllGID_ =      - tmp_recv[0];
      MaxAllGID_ =        tmp_recv[1];
      MinElementSize_ = - tmp_recv[2];
      MaxElementSize_ =   tmp_recv[3];

      delete [] tmp_send;
      delete [] tmp_recv;

      if (MinAllGID_ < IndexBase_)
	{
	  cerr << "Petra_BlockMap: ERROR, Minimum element index = " << MinAllGID_
	       << ". which is less than the Index Base = "<< IndexBase_ <<" . Aborting." << endl;
	  abort();
	}

    }
  else
    {
      cerr << "Petra_BlockMap: INTERNAL ERROR, NumGlobalElements = " << NumGlobalElements
	   << ". Should = "<< NumMyElements 
	   << "and NumProc = " << NumProc << ".  Should = 1. Aborting." << endl;
      abort();
    }

  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}
//==============================================================================
Petra_BlockMap::Petra_BlockMap(const Petra_BlockMap& map)
  : LID_(0),
    NumGlobalElements_(map.NumGlobalElements_),
    NumMyElements_(map.NumMyElements_),
    MyGlobalElements_(0),
    FirstElementEntryList_(0),
    ElementSize_(map.ElementSize_),
    ElementSizeList_(0),
    EquationToBlockList_(0),
    IndexBase_(map.IndexBase_),
    Comm_(map.Comm_),
    Directory_(0),
    NumGlobalEquations_(map.NumGlobalEquations_),
    NumMyEquations_(map.NumMyEquations_),
    MinAllGID_(map.MinAllGID_),
    MaxAllGID_(map.MaxAllGID_),
    MinMyGID_(map.MinMyGID_),
    MaxMyGID_(map.MaxMyGID_),
    MinLID_(map.MinLID_),
    MaxLID_(map.MaxLID_),
    MinMyElementSize_(map.MinMyElementSize_),
    MaxMyElementSize_(map.MaxMyElementSize_),
    MinElementSize_(map.MinElementSize_),
    MaxElementSize_(map.MaxElementSize_),
    ConstantElementSize_(map.ConstantElementSize_),
    LinearMap_(map.LinearMap_),
    DistributedGlobal_(map.DistributedGlobal_)
{
  int i;

  if (map.MyGlobalElements_!=0)
    {
      MyGlobalElements_ = new int[NumMyElements_];
      
      for(i=0; i<NumMyElements_; i++)
	MyGlobalElements_[i] = map.MyGlobalElements_[i];
    }
  if (map.FirstElementEntryList_!=0)
    {
      FirstElementEntryList_ = new int[NumMyElements_+1];
      
      for(i=0; i<NumMyElements_+1; i++)
	FirstElementEntryList_[i] = map.FirstElementEntryList_[i];
    }
  if (map.ElementSizeList_!=0)
    {
      ElementSizeList_ = new int[NumMyElements_];
      
      for(i=0; i<NumMyElements_; i++)
	ElementSizeList_[i] = map.ElementSizeList_[i];
    }
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.

}

//==============================================================================
Petra_BlockMap::~Petra_BlockMap(void)  {

  if (MyGlobalElements_ != 0 && NumMyElements_>0) delete [] MyGlobalElements_;
  MyGlobalElements_ = 0;
  
  if (FirstElementEntryList_ != 0 && NumMyElements_>0) delete [] FirstElementEntryList_;
  FirstElementEntryList_ = 0;
  
  if (ElementSizeList_ != 0 && NumMyElements_>0) delete [] ElementSizeList_;
  ElementSizeList_ = 0;

  if (EquationToBlockList_ != 0 && NumMyEquations_>0) delete [] EquationToBlockList_;
  EquationToBlockList_ = 0;

  if (Directory_ !=0) delete Directory_;
  Directory_ = 0;

  if (LID_ !=0 && NumMyElements_>0) delete [] LID_;
  LID_ = 0;
}


//==============================================================================
bool Petra_BlockMap::SameAs(const Petra_BlockMap & Map) const
{

  if (this == &Map) return(true);
    
  if (MinAllGID_ != Map.MinAllGID() ||
      MaxAllGID_ != Map.MaxAllGID() ||
      NumGlobalElements_!=Map.NumGlobalElements() ||
      IndexBase_!=Map.IndexBase() ) return(false);
  
  if (ConstantElementSize_) {
    if (ElementSize_!=Map.ElementSize()) return(false);
    else return(true);
  }
  else {

    // If we get this far, we need to check local properties and then check across
    // all processors to see if local properties are all true

    int MySameMap = 1; // Assume not needed
    if (NumMyElements_!=Map.NumMyElements()) MySameMap = 0;
    
    if (MySameMap==1) 
      for (int i=0; i<NumMyElements_; i++) 
	if (GID(i) != Map.GID(i)) MySameMap = 0;

    // Now get min of MySameMap across all processors

    int GlobalSameMap = 0;
    assert(Comm().MinAll(&MySameMap, &GlobalSameMap, 1)==0);
    
    return(GlobalSameMap==1);
  }
  return(false);
}


//==============================================================================
int Petra_BlockMap::MyGlobalElements(int * MyGlobalElements) const
{
  // If the global element list is not create, then do so.  This can only happen when
  // a linear distribution has been specified.  Thus we can easily construct the update
  // list in this case.

  int i;
  if (MyGlobalElements_==0)
      for (i = 0; i<NumMyElements_; i++)
	MyGlobalElements[i] = MinMyGID_ + i;
  else
    for (i = 0; i<NumMyElements_; i++)
      MyGlobalElements[i] = MyGlobalElements_[i];
  
  return(0);
}

//==============================================================================
int * Petra_BlockMap::MyGlobalElements() const {
  
  // If ElementSizeList not built, do so
  if (MyGlobalElements_==0 && NumMyElements_>0) {
    int * tmp = new int[NumMyElements_+1];
    MyGlobalElements(tmp);
    (int * &) MyGlobalElements_ = tmp;
  }
  return(MyGlobalElements_);
  
}
//==============================================================================
int Petra_BlockMap::FirstElementEntryList(int * FirstElementEntryList) const
{
  // If the first element entry list is not create, then do so.  

  // Note: This array is of length NumMyElement+1

  int i;

  if (FirstElementEntryList_==0) {
    FirstElementEntryList[0] = 0; // First element of first entry is always zero
    
    if (ConstantElementSize_)
      for (i = 0; i<NumMyElements_; i++)
	FirstElementEntryList[i+1] = FirstElementEntryList[i] + ElementSize_;
    else
      for (i = 0; i<NumMyElements_; i++)
	FirstElementEntryList[i+1] = FirstElementEntryList[i] + ElementSizeList_[i];
  }
  else 
    for (i = 0; i<=NumMyElements_; i++)
      FirstElementEntryList[i] = FirstElementEntryList_[i];
  return(0);
}

//==============================================================================
int * Petra_BlockMap::FirstElementEntryList() const {

  // If ElementSizeList not built, do so
  if (FirstElementEntryList_==0 && NumMyElements_>0) {
    int * tmp = new int[NumMyElements_+1];
    FirstElementEntryList(tmp);
    (int * &) FirstElementEntryList_ = tmp;
 }
  return(FirstElementEntryList_);
  
}
//==============================================================================
int Petra_BlockMap::ElementSizeList(int * ElementSizeList) const
{
  // If the element size list is not create, then do so.  This can only happen when
  // a constant element size has been specified.  Thus we can easily construct the block size
  // list in this case.

  int i;
  if (ElementSizeList_==0)
    for (i = 0; i<NumMyElements_; i++)
      ElementSizeList[i] = ElementSize_;
  else
    for (i = 0; i<NumMyElements_; i++)
      ElementSizeList[i] = ElementSizeList_[i];
  
  return(0);
}
//==============================================================================
int * Petra_BlockMap::ElementSizeList() const {

  // If ElementSizeList not built, do so
  if (ElementSizeList_==0 && NumMyElements_>0) {
    int * tmp = new int[NumMyElements_];
    ElementSizeList(tmp);
    (int * &) ElementSizeList_ = tmp;
 }
  return(ElementSizeList_);
  
}
//==============================================================================
int Petra_BlockMap::EquationToBlockList(int * EquationToBlockList) const
{
  // Build an array such that the local block ID is stored for each equation

  int i;
  if (EquationToBlockList_==0) {
    int * ptr = EquationToBlockList;
    for (i = 0; i<NumMyElements_; i++) {
      int Size = ElementSize(i);
      for (int j=0; j<Size; j++) *ptr++ = i;
    }
  }
  else
    for (i = 0; i<NumMyEquations_; i++)
      EquationToBlockList[i] = EquationToBlockList_[i];
  
  return(0);
}
//==============================================================================
int * Petra_BlockMap::EquationToBlockList() const {

  // If EquationToBlockList not built, do so
  if (EquationToBlockList_==0 && NumMyEquations_>0) {
    int * tmp = new int[NumMyEquations_];
    EquationToBlockList(tmp);
    (int * &) EquationToBlockList_ = tmp;
 }
  return(EquationToBlockList_);
  
}
//==============================================================================
int Petra_BlockMap::ElementSize(int LID) const {

  if (ConstantElementSize()) return(ElementSize_);
  else
    return(ElementSizeList_[LID]);
  
}
//==============================================================================
void Petra_BlockMap::GlobalToLocalSetup() {

  int i;

  if (NumGlobalElements_==0) return; // Nothing to do

  else if (LinearMap() || (!DistributedGlobal()) || NumMyElements_==0) {
    if (Directory_ ==0) Directory_ = new Petra_Directory(this); // Make directory
    return; // Nothing else to do
  }
  else {
    // Build LID_ vector to make look up of local index values fast
    
    int SpanGID = MaxMyGID_ - MinMyGID_ + 1;
    LID_ = new int[SpanGID];
    
    for (i=0; i<SpanGID; i++) LID_[i] = -1; // Fill all locations with -1
    
    for (i=0; i<NumMyElements_; i++) {
      int tmp = MyGlobalElements_[i]-MinMyGID_;
      assert(tmp>=0); assert(tmp <SpanGID);
      LID_[MyGlobalElements_[i]-MinMyGID_] = i; // Spread local indices
    }
    
    if (Directory_ ==0) Directory_ = new Petra_Directory(this); // Make directory
  }
}

//==============================================================================
int Petra_BlockMap::LID(int GID) const {

  if (GID<MinMyGID_ || GID > MaxMyGID_) return(-1); // Out of range
  else if (!DistributedGlobal()) return(GID-IndexBase_); // I own all indices
  else if (LinearMap()) return(GID-MinMyGID_); // Can compute with an offset
  else return(LID_[GID-MinMyGID_]); // Find it in LID array
}
//==============================================================================
int Petra_BlockMap::GID(int LID) const {

  if (LID<MinLID_ || LID>MaxLID_) return(IndexBase_-1); // Out of range
  else if (!DistributedGlobal()) return(LID+IndexBase_); // I own all indices
  else if (LinearMap()) return(LID+MinMyGID_); // Can compute with an offset
  else return(MyGlobalElements_[LID]); // Find it in MyGlobalElements array
}
//==============================================================================
int Petra_BlockMap::FindLocalBlockID(int EquationID, int & BlockID, int & BlockOffset) const {

  int ierr = 0;

  if (EquationID>=NumMyEquations_) return(-1); // Equation is out of range

  if (ConstantElementSize()) {
    BlockID = EquationID/MaxElementSize_;
    BlockOffset = EquationID%MaxElementSize_;
    return(0);
  }
  else {
    if (EquationToBlockList_==0) {
      int * tmp;
      if (EquationToBlockList_==0) tmp = new int[NumMyEquations_];
      ierr = EquationToBlockList(tmp);
      if (ierr!=0) return(ierr);
      (int * &) EquationToBlockList_ = tmp;  // Allows assignment in a const method
    }
    if (FirstElementEntryList_==0) {
      int * tmp;
      if (FirstElementEntryList_==0) tmp = new int[NumMyElements_];
      ierr = FirstElementEntryList(tmp);
      if (ierr!=0) return(ierr);
      (int * &) FirstElementEntryList_ = tmp;  // Allows assignment in a const method
    }

    BlockID = EquationToBlockList_[EquationID];
    BlockOffset = EquationID - FirstElementEntryList_[BlockID];
    return(0);
  }
}
//==============================================================================
int Petra_BlockMap::RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList, int * SizeList) const {

  return(Directory_->GetDirectoryEntries(NumIDs, GIDList, PIDList, LIDList, SizeList));
}
// Non-member functions

ostream& operator << (ostream& os, const Petra_BlockMap & Map)
{
  int * MyGlobalElements = Map.MyGlobalElements();
  int * FirstElementEntryList = 0;
  int * ElementSizeList = 0;
  if (!Map.ConstantElementSize()) {
    FirstElementEntryList = Map.FirstElementEntryList();
    ElementSizeList = Map.ElementSizeList();
  }
  int MyPID = Map.Comm().MyPID();
  int NumProc = Map.Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);
      if (MyPID==0) {
	os <<  "\nNumber of Global Elements  = "; os << Map.NumGlobalElements(); os << endl;
	os <<    "Number of Global Equations = "; os << Map.NumGlobalEquations(); os << endl;
	os <<    "Maximum of all GIDs        = "; os << Map.MaxAllGID(); os << endl;
	os <<    "Minimum of all GIDs        = "; os << Map.MinAllGID(); os << endl;
	os <<    "Index Base                 = "; os << Map.IndexBase(); os << endl;
	if (Map.ConstantElementSize())
	  os <<  "Constant Element Size      = "; os << Map.ElementSize(); os << endl;
      }
      os << endl;

	os <<    "Number of Local Elements   = "; os << Map.NumMyElements(); os << endl;
	os <<    "Number of Local Equations  = "; os << Map.NumMyEquations(); os << endl;
	os <<    "Maximum of my GIDs         = "; os << Map.MaxMyGID(); os << endl;
	os <<    "Minimum of my GIDs         = "; os << Map.MinMyGID(); os << endl;
      os << endl;

      os.width(14);
      os <<  "     MyPID"; os << "    ";
      os.width(14);
      os <<  "       Local Index "; os << " ";
      os.width(14);
      os <<  "      Global Index "; os << " ";
      if (!Map.ConstantElementSize()) {
	os.width(14);
	os <<" FirstElementEntry "; os << " ";
	os.width(14);
	os <<"   ElementSize "; os << " ";
      }
      os << endl;
    
      for (int i=0; i < Map.NumMyElements(); i++) {
	os.width(14);
	os <<  MyPID; os << "    ";
	os.width(14);
	os <<  i; os << "    ";
	os.width(14);
	os <<  MyGlobalElements[i]; os << "    ";
	if (!Map.ConstantElementSize()) {	  
	  os.width(14);
	  os << FirstElementEntryList[i]; os << "    ";
	  os.width(14);
	  os << ElementSizeList[i]; os << "    ";
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
    Map.Comm().Barrier();
    Map.Comm().Barrier();
    Map.Comm().Barrier();
  }
  return os;
}
