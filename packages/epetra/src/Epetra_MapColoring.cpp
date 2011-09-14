
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_HashTable.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Util.h"
//=============================================================================
Epetra_MapColoring::Epetra_MapColoring(const Epetra_BlockMap& map, int * elementColors, 
				       const int defaultColor)
  : Epetra_DistObject(map, "Epetra::MapColoring"),
    DefaultColor_(defaultColor),
    ColorIDs_(0),
    FirstColor_(0),
    NumColors_(0),
    ListOfColors_(0),
    ColorCount_(0),
    ElementColors_(0),
    ColorLists_(0),
    Allocated_(false),
    ListsAreGenerated_(false),
    ListsAreValid_(false)
{
  Allocate(elementColors, 1);
}
//=============================================================================
Epetra_MapColoring::Epetra_MapColoring(const Epetra_BlockMap& map,
				       const int defaultColor)
  : Epetra_DistObject(map, "Epetra::MapColoring"),
    DefaultColor_(defaultColor),
    ColorIDs_(0),
    FirstColor_(0),
    NumColors_(0),
    ListOfColors_(0),
    ColorCount_(0),
    ElementColors_(0),
    ColorLists_(0),
    Allocated_(false),
    ListsAreGenerated_(false),
    ListsAreValid_(false)
{
  Allocate(&DefaultColor_, 0);
}
//=============================================================================
Epetra_MapColoring::Epetra_MapColoring(const Epetra_MapColoring& Source)
  : Epetra_DistObject(Source),
    DefaultColor_(Source.DefaultColor_),
    ColorIDs_(0),
    FirstColor_(0),
    NumColors_(0),
    ListOfColors_(0),
    ColorCount_(0),
    ElementColors_(0),
    ColorLists_(0),
    Allocated_(false),
    ListsAreGenerated_(false),
    ListsAreValid_(false)
{
  Allocate(Source.ElementColors_, 1);
}
//=========================================================================
Epetra_MapColoring::~Epetra_MapColoring(){


  if (Allocated_ && Map().NumMyElements()>0) delete [] ElementColors_;
  if (ListsAreGenerated_) DeleteLists();
}

//=========================================================================
int Epetra_MapColoring::DeleteLists() const {


  if (ListsAreGenerated_) {
    for (int i=0; i<NumColors_; i++) if (ColorLists_[i]!=0) delete [] ColorLists_[i];
    delete [] ColorLists_;
    delete [] ColorCount_;
    delete [] ListOfColors_;
    delete ColorIDs_;
    ListItem * CurItem = FirstColor_;
    while (CurItem!=0) {
      ListItem * NextItem = CurItem->NextItem;
      delete CurItem;
      CurItem = NextItem;
    }
  }
  ListsAreValid_ = false;
  return(0);
}

//=========================================================================
int Epetra_MapColoring::Allocate(int * elementColors, int Increment)
{
  
  if (Allocated_) return(0);
  
  int NumMyElements = Map().NumMyElements();
  if (NumMyElements>0) ElementColors_ = new int[NumMyElements];
  for (int i=0; i< NumMyElements; i++) ElementColors_[i] = elementColors[i*Increment];
  Allocated_ = true;
  return(0);
}

//=========================================================================
int Epetra_MapColoring::GenerateLists() const {
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
bool Epetra_MapColoring::InItemList(int ColorValue) const {
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
int Epetra_MapColoring::NumElementsWithColor(int Color) const  {
  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = -1;
  if( ColorIDs_ )
    arrayIndex = ColorIDs_->Get(Color);
  if (arrayIndex>-1) return(ColorCount_[arrayIndex]);
  else return(0);
}
//=========================================================================
int * Epetra_MapColoring::ColorLIDList(int Color) const  {
  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = -1;
  if( ColorIDs_ )
    arrayIndex = ColorIDs_->Get(Color);
  if (arrayIndex>-1) return(ColorLists_[arrayIndex]);
  else return(0);
}
//=========================================================================
Epetra_Map * Epetra_MapColoring::GenerateMap(int Color) const {

  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = -1;
  if( ColorIDs_ )
    arrayIndex = ColorIDs_->Get(Color);
  int NumElements = 0;
  int * ColorElementLIDs = 0;
  int * ColorElementGIDs =0;
  if (arrayIndex>-1) NumElements = ColorCount_[arrayIndex];
  if (NumElements>0) {
    ColorElementLIDs = ColorLIDList(Color);
    ColorElementGIDs = new int[NumElements];
    for (int i=0; i<NumElements; i++) ColorElementGIDs[i] = Map().GID(ColorElementLIDs[i]);
  }
  Epetra_Map * map = new Epetra_Map(-1, NumElements, ColorElementGIDs, 
				    Map().IndexBase(), Map().Comm());
  if (ColorElementGIDs!=0) delete [] ColorElementGIDs;
  return(map);
}
//=========================================================================
Epetra_BlockMap * Epetra_MapColoring::GenerateBlockMap(int Color) const {

  if (!ListsAreValid_) GenerateLists(); 
  int arrayIndex = -1;
  if( ColorIDs_ )
    arrayIndex = ColorIDs_->Get(Color);
  int NumElements = 0;
  int * ColorElementLIDs = 0;
  int * ColorElementSizes = 0;
  int * ColorElementGIDs = 0;
  if (arrayIndex>-1) NumElements = ColorCount_[arrayIndex];
  if (NumElements>0) {
    ColorElementLIDs = ColorLIDList(Color);
    ColorElementSizes = new int[NumElements];
    ColorElementGIDs = new int[NumElements];
    for (int i=0; i<NumElements; i++) ColorElementGIDs[i] = Map().GID(ColorElementLIDs[i]);
  }
  int * MapElementSizes = Map().ElementSizeList();

  {for (int i=0; i<NumElements; i++) 
    ColorElementSizes[i] = MapElementSizes[ColorElementLIDs[i]];}

  Epetra_BlockMap * map = new Epetra_BlockMap(-1, NumElements, ColorElementGIDs, 
					      ColorElementSizes,
					      Map().IndexBase(), Map().Comm());

  if (ColorElementGIDs!=0) delete [] ColorElementGIDs;
  if (ColorElementSizes!=0) delete [] ColorElementSizes;

  return(map);
}
//=========================================================================
void Epetra_MapColoring::Print(ostream& os) const {
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
int Epetra_MapColoring::MaxNumColors() const {

  if (!ListsAreValid_) GenerateLists(); 
  int tmp1 = NumColors_, tmp2;
  Map().Comm().MaxAll(&tmp1, &tmp2, 1);
  return(tmp2);
}
//=========================================================================
int Epetra_MapColoring::CheckSizes(const Epetra_SrcDistObject& Source) {
  (void)Source;//prevents unused variable compiler warning
  return(0);
}

//=========================================================================
int Epetra_MapColoring::CopyAndPermute(const Epetra_SrcDistObject& Source,
                                       int NumSameIDs, 
				       int NumPermuteIDs,
                                       int * PermuteToLIDs, 
				       int *PermuteFromLIDs,
                                       const Epetra_OffsetIndex * Indexor)
{
  (void)Indexor;
  const Epetra_MapColoring & A = dynamic_cast<const Epetra_MapColoring &>(Source);

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
int Epetra_MapColoring::PackAndPrepare(const Epetra_SrcDistObject & Source,
                                       int NumExportIDs,
                                       int * ExportLIDs,
				       int & LenExports,
                                       char * & Exports,
				       int & SizeOfPacket,
				       int * Sizes,
				       bool & VarSizes,
                                       Epetra_Distributor & Distor)
{
  (void)Sizes;
  (void)VarSizes;
  (void)Distor;
  const Epetra_MapColoring & A = dynamic_cast<const Epetra_MapColoring &>(Source);

  int  * From = A.ElementColors();
  int * IntExports = 0;

  SizeOfPacket = (int)sizeof(int); 

  if (NumExportIDs*SizeOfPacket>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = NumExportIDs*SizeOfPacket;
    IntExports = new int[LenExports];
    Exports = (char *) IntExports;
  }

  int * ptr;

  if (NumExportIDs>0) {
    ptr = (int *) Exports;    
    for (int j=0; j<NumExportIDs; j++) ptr[j] = From[ExportLIDs[j]];
  }
  
  return(0);
}

//=========================================================================
int Epetra_MapColoring::UnpackAndCombine(const Epetra_SrcDistObject & Source,
					 int NumImportIDs,
                                         int * ImportLIDs, 
                                         int LenImports,
					 char * Imports,
                                         int & SizeOfPacket, 
					 Epetra_Distributor & Distor, 
					 Epetra_CombineMode CombineMode,
                                         const Epetra_OffsetIndex * Indexor )
{
  (void)Source;
  (void)LenImports;
  (void)Imports;
  (void)SizeOfPacket;
  (void)Distor;
  (void)Indexor;
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
    for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] += ptr[j]; // Add to existing value
  else if(CombineMode==Insert)
    for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = ptr[j];
  else if(CombineMode==AbsMax) {
    for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = 0;
    for (j=0; j<NumImportIDs; j++)  To[ImportLIDs[j]] = EPETRA_MAX( To[ImportLIDs[j]],std::abs(ptr[j]));
  }
  
  return(0);
}

