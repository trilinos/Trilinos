
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_HashTable.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblemRedistor.h"
#include "Epetra_Util.h"
//=============================================================================
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(const Epetra_LinearProblem& Problem, 
																													 const Epetra_Map & RedistMap)
  : Epetra_Object(Map, "Epetra::LinearProblemRedistor"),
    OrigProblem_(Problem),
    RedistProblem_(0),
    RedistMap_((Epetra_Map *) &RedistMap),
    Replicated_(false),
    ConstructTranspose_(false),
    MakeDataContiguous_(false),
    MapGenerated_(false),
		ptr_(0)
{
}
//=============================================================================
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(const Epetra_LinearProblem& Problem, 
																													 const bool Replicate)
  : Epetra_Object(Map, "Epetra::LinearProblemRedistor"),
    OrigProblem_(Problem),
    RedistProblem_(0),
    RedistMap_((Epetra_Map *) &RedistMap),
    Replicated_(Replicate),
    ConstructTranspose_(false),
    MakeDataContiguous_(false),
    MapGenerated_(false),
		ptr_(0)
{

	if (Problem.GetMatrix()!=0)
		throw ReportError("Matrix must be defined for the input problem in order to generate RedistMap", -1);

	int ierr = GenerateRedistMap();

	if (ierr!=0) 
		throw ReportError("Error generating redistribution map.  Input problem map is not one to one.", ierr);
}
//=============================================================================
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(const Epetra_LinearProblemRedistor& Source)
  : Epetra_Object(Source)
    OrigProblem_(Source.OrigProblem_),
    RedistProblem_(0),
    RedistMap_(Source.RedistMap_),
    Replicated_(Source.Replicated_),
    ConstructTranspose_(Source.ConstructTranspose_),
    MakeDataContiguous_(Source.MakeDataContiguous_),
		ptr_(0)
{

	RedistProblem_ = Epetra_LinearProblem(*Source.RedistProblem_);
}
//=========================================================================
Epetra_LinearProblemRedistor::~Epetra_LinearProblemRedistor(){

	if (ptr_!=0) {delete [] ptr_; ptr_=0;}
	if (MapGenerated_ && RedistMap_!=0) {delete RedistMap_; RedistMap_=0;}
	if (RedistProblem_!=0) {delete RedistProblem_; RedistProblem_=0;}

}

//=========================================================================
int Epetra_LinearProblemRedistor::GenerateRedistMap() {


  if (MapGenerated_) return(0);

	const Epetra_Map & SourceMap = OrigProblem_->GetMatrix().RowMatrixRowMap();
	const Epetra_Comm & SourceMap.Comm();
	int IndexBase = SourceMap.IndexBase();

	// Build a list of contiguous GIDs that will be used in either case below.
	int NumMyRedistElements = 0;
	if ((Comm.MyPID()==0) || Replicated_) NumMyRedistElements = SourceMap.NumGlobalElements();
	
	// Now build the GID list to broadcast SourceMapGIDs to all processors that need it (or just to PE 0).
	int * ContigIDs = 0;
	if (NumMyRedistElements>0) ContigIDs = new int[NumMyRedistElements];
	for (int i=0; i<NumMyRedistElements; i++) ContigIDs[i] = IndexBase + i;
	
	// Case 1: If the map for the input matrix is not a linear contiguous map, then we have to collect
	// the map indices in order to construct a compatible map containing all GIDs.
	if (!SourceMap.LinearMap()) {

		// First generate a linear map of the same distribution as RowMatrixRowMap
		Epetra_Map SourceLinearMap(-1, SourceMap.NumMyElements(), IndexBase, Comm);
		// Generate Int vector containing GIDs of SourceMap
		Epetra_IntVector SourceMapGIDs(View, SourceLinearMap, SourceMap.MyGlobalElements());
		// Now Build target map for SourceMapGIDs and Importer to get IDs
		Epetra_Map GIDsTargetMap(-1, NumMyRedistElements, ContigIDs, IndexBase, Comm);
		if (NumMyRedistElements>0) delete [] ContigIDs;

		Epetra_Import GIDsImporter(GIDsTargetMap, SourceMap);
		Epetra_IntVector TargetMapGIDs(GIDsTargetMap);

		// Now Send SourceMap GIDs to PE 0, and all processors if Replicated is true
		EPETRA_CHK_ERR(TargetMapGIDs.Import(SourceMapGIDs, GIDsImporter, Insert));

		// Finally, create RedistMap containing all GIDs of SourceMap on PE 0, or all PEs of Replicate is true
		RedistMap_ = new Epetra_Map(-1, NumMyRedistElements, TargetMapGIDs.Values, IndexBase, Comm);

	}
	// Case 2: If the map has contiguous IDs then we can simply build the map right away using the list
	// of contiguous GIDs directly
	else
		RedistMap_ = new Epetra_Map(-1, NumMyRedistElements, ContigIDs, IndexBase, Comm);

	MapGenerated_ = true;

  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::CreateRedistProblem(const bool ConstructTranspose, const bool MakeDataContiguous, 
																											Epetra_LinearProblem *& RedistProblem)
{
  
	if (Allocated_) return(0);
  
  int NumMyElements = Map().NumMyElements();
  ElementColors_ = new int[NumMyElements];
  for (int i=0; i< NumMyElements; i++) ElementColors_[i] = ElementColors[i*Increment];
  Allocated_ = true;
  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::GenerateLists() const {
  int NumMyElements = Map().NumMyElements();
  if (NumMyElements==0) return(0); // Nothing to do

  if (ListsAreValid_) return(0); // Already been here

  if (ListsAreGenerated_) DeleteLists();  // Delete any existing lists


  // Scan the ElementColors to determine how many colors we have
  NumColors_ = 1;
  FirstColor_ = new ListItem(ElementColors_[0]); // Initialize First color in list
  for (int i=1; i<NumMyElements; i++) if (!InItemList(ElementColors_[i])) NumColors_++;

  // Create hash table that maps color IDs to the integers 0,...NumColors_
  ColorIDs_ = new Epetra_HashTable(NumColors_);
  ListOfColors_ = new int[NumColors_];
  ListItem * CurItem = FirstColor_;
  {for (int i=0; i<NumColors_; i++) {
    ColorIDs_->Add(CurItem->ItemValue, i); // Create hash table entry
    ListOfColors_[i] = CurItem->ItemValue; // Put color value in a list of colors
    CurItem = CurItem->NextItem;
  }}
  Epetra_Util util;
  util.Sort(true, NumColors_, ListOfColors_, 0, 0, 0, 0); // Sort List of colors in ascending order
  // Count the number of IDs of each color
  ColorCount_ = new int[NumColors_];
  {for (int i=0; i<NumColors_; i++) ColorCount_[i] = 0;}
  {for (int i=0; i<NumMyElements; i++) ColorCount_[ColorIDs_->Get(ElementColors_[i])]++;}

  // Finally build list of IDs grouped by color
  ColorLists_ = new int *[NumColors_];
  {for (int i=0; i<NumColors_; i++) ColorLists_[i] = new int[ColorCount_[i]];}
  {for (int i=0; i<NumColors_; i++) ColorCount_[i] = 0;} // Reset so we can use for counting
  {for (int i=0; i<NumMyElements; i++) {
    int j = ColorIDs_->Get(ElementColors_[i]);
    ColorLists_[j][ColorCount_[j]++] = i;
  }}
  ListsAreValid_ = true;
  ListsAreGenerated_ = true;

  return(0);
}
//=========================================================================
bool Epetra_LinearProblemRedistor::InItemList(int ColorValue) const {
  bool ColorFound = false;
  ListItem * CurColor = 0;
  ListItem * NextColor = FirstColor_;
  while (!ColorFound && NextColor!=0) {
    CurColor = NextColor;
    NextColor = CurColor->NextItem;
    if (ColorValue==CurColor->ItemValue) ColorFound = true;
  }

  if (!ColorFound) CurColor->NextItem = new ListItem(ColorValue);

  return(ColorFound);
}
//=========================================================================
int Epetra_LinearProblemRedistor::NumElementsWithColor(int Color) const  {
  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = ColorIDs_->Get(Color);
  if (arrayIndex>-1) return(ColorCount_[arrayIndex]);
  else return(0);
}
//=========================================================================
int * Epetra_LinearProblemRedistor::ColorLIDList(int Color) const  {
  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = ColorIDs_->Get(Color);
  if (arrayIndex>-1) return(ColorLists_[arrayIndex]);
  else return(0);
}
//=========================================================================
Epetra_Map * Epetra_LinearProblemRedistor::GenerateMap(int Color) const {

  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = ColorIDs_->Get(Color);
  int NumElements = ColorCount_[arrayIndex];
  int * ColorElementLIDs = ColorLIDList(Color);
  int * ColorElementGIDs = new int[NumElements];
  for (int i=0; i<NumElements; i++) ColorElementGIDs[i] = Map().GID(ColorElementLIDs[i]);
  Epetra_Map * map = new Epetra_Map(-1, NumElements, ColorElementGIDs, 
				    Map().IndexBase(), Map().Comm());
  delete [] ColorElementGIDs;
  return(map);
}
//=========================================================================
Epetra_BlockMap * Epetra_LinearProblemRedistor::GenerateBlockMap(int Color) const {

  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = ColorIDs_->Get(Color);
  int NumElements = ColorCount_[arrayIndex];
  int * ColorElementLIDs = ColorLIDList(Color);
  int * ColorElementSizes = new int[NumElements];
  int * ColorElementGIDs = new int[NumElements];
  for (int i=0; i<NumElements; i++) ColorElementGIDs[i] = Map().GID(ColorElementLIDs[i]);
  int * MapElementSizes = Map().ElementSizeList();

  {for (int i=0; i<NumElements; i++) 
    ColorElementSizes[i] = MapElementSizes[ColorElementLIDs[i]];}

  Epetra_BlockMap * map = new Epetra_BlockMap(-1, NumElements, ColorElementGIDs, 
					      ColorElementSizes,
					      Map().IndexBase(), Map().Comm());

  delete [] ColorElementGIDs;
  delete [] ColorElementSizes;

  return(map);
}
//=========================================================================
void Epetra_LinearProblemRedistor::Print(ostream& os) const {
  int MyPID = Map().Comm().MyPID();
  int NumProc = Map().Comm().NumProc();
  
  if (MyPID==0) os 
    << endl 
    << " *****************************************" << endl
    << " Coloring information arranged map element" << endl 
    << " *****************************************" << endl
    << endl;
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyElements1 =Map(). NumMyElements();
      int * MyGlobalElements1 = Map().MyGlobalElements();

      if (MyPID==0) {
	os.width(8);
	os <<  "     MyPID"; os << "    ";
	os.width(12);
	os <<  "GID  ";
	os.width(20);
	os <<  "Color  ";
	os << endl;
      }
      for (int i=0; i < NumMyElements1; i++) {
	os.width(10);
	os <<  MyPID; os << "    ";
	os.width(10);
	os << MyGlobalElements1[i] << "    ";
	os.width(20);
	os <<  ElementColors_[i];
	os << endl;
      }
      os << flush; 
    }

    // Do a few global ops to give I/O a chance to complete
    Map().Comm().Barrier();
    Map().Comm().Barrier();
    Map().Comm().Barrier();
  }

  if (MyPID==0) os 
    << endl 
    << " **************************************" << endl
    << " Coloring information arranged by color" << endl 
    << " **************************************" << endl
    << endl;
  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      if (NumColors()==0) os << " No colored elements on processor " << MyPID << endl;
      else {
        os << "Number of colors in map = " << NumColors() << endl
	         << "Default color           = " << DefaultColor() << endl << endl;
        if (MyPID==0) {
	        os.width(8);
	        os <<  "     MyPID"; os << "    ";
	        os.width(12);
	        os <<  "LID  ";
	        os.width(20);
	        os <<  "Color  ";
	        os << endl;
        }
	      int * ColorValues = ListOfColors();
	      for (int ii=0; ii<NumColors(); ii++) {
	        int CV = ColorValues[ii];
	  int ColorCount = NumElementsWithColor(CV);
	  int * LIDList = ColorLIDList(CV);
	  
	  
	  for (int i=0; i < ColorCount; i++) {
	    os.width(10);
	    os <<  MyPID; os << "    ";
	    os.width(10);
	    os << LIDList[i] << "    ";
	    os.width(20);
	    os << CV;
	    os << endl;
	  }
	  os << flush; 
	}
      }
    }
    // Do a few global ops to give I/O a chance to complete
    Map().Comm().Barrier();
    Map().Comm().Barrier();
    Map().Comm().Barrier();
  }}
  return;
}
//=========================================================================
int Epetra_LinearProblemRedistor::CheckSizes(const Epetra_DistObject& Source) {
  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::CopyAndPermute(const Epetra_DistObject& Source, int NumSameIDs, 
				       int NumPermuteIDs, int * PermuteToLIDs, 
				       int *PermuteFromLIDs) {

  const Epetra_LinearProblemRedistor & A = dynamic_cast<const Epetra_LinearProblemRedistor &>(Source);

  int * From = A.ElementColors();
  int *To = ElementColors_;

  // Do copy first
  if (NumSameIDs>0)
    if (To!=From) {
      for (int j=0; j<NumSameIDs; j++)
	To[j] = From[j];
    }
  // Do local permutation next
  if (NumPermuteIDs>0)
    for (int j=0; j<NumPermuteIDs; j++) 
      To[PermuteToLIDs[j]] = From[PermuteFromLIDs[j]];
  
  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::PackAndPrepare(const Epetra_DistObject & Source, int NumExportIDs, int * ExportLIDs,
				      int Nsend, int Nrecv,
				      int & LenExports, char * & Exports, int & LenImports, 
				      char * & Imports, 
				      int & SizeOfPacket, Epetra_Distributor & Distor) {



  const Epetra_LinearProblemRedistor & A = dynamic_cast<const Epetra_LinearProblemRedistor &>(Source);

  int  * From = A.ElementColors();
  int * IntExports = 0;
  int * IntImports = 0;

  if (Nsend>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = Nsend;
    IntExports = new int[LenExports];
    Exports = (char *) IntExports;
  }

  if (Nrecv>LenImports) {
    if (LenImports>0) delete [] Imports;
    LenImports = Nrecv;
    IntImports = new int[LenImports];
    Imports = (char *) IntImports;
  }

  SizeOfPacket = sizeof(int); 

  int * ptr;

  if (NumExportIDs>0) {
    ptr = (int *) Exports;    
    for (int j=0; j<NumExportIDs; j++) *ptr++ = From[ExportLIDs[j]];
  }
  
  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::UnpackAndCombine(const Epetra_DistObject & Source,
					 int NumImportIDs, int * ImportLIDs, 
					char * Imports, int & SizeOfPacket, 
					 Epetra_Distributor & Distor, 
					 Epetra_CombineMode CombineMode ) {
  int j;
  
  if(    CombineMode != Add
      && CombineMode != Zero
      && CombineMode != Insert
      && CombineMode != AbsMax )
    EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

  if (NumImportIDs<=0) return(0);

  int * To = ElementColors_;

  int * ptr;
  // Unpack it...

  ptr = (int *) Imports;
    
  if (CombineMode==Add)
    for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] += *ptr++; // Add to existing value
  else if(CombineMode==Insert)
    for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = *ptr++;
  else if(CombineMode==AbsMax)
    for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = EPETRA_MAX( To[ImportLIDs[j]],abs(*ptr++));
  
  return(0);
}

