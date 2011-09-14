
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

#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"



//=============================================================================
// Constructor

Epetra_DistObject::Epetra_DistObject(const Epetra_BlockMap& map)
  : Epetra_Object("Epetra::DistObject"),
    Map_(map),
    Comm_(&Map_.Comm()),
    DistributedGlobal_(map.DistributedGlobal()),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Sizes_(0)
{}

//=============================================================================
// Constructor (label given)

Epetra_DistObject::Epetra_DistObject(const Epetra_BlockMap& map, const char* const label)
  : Epetra_Object(label),
    Map_(map),
    Comm_(&Map_.Comm()),
    DistributedGlobal_(map.DistributedGlobal()),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Sizes_(0)
{}

//=============================================================================
// Copy Constructor

Epetra_DistObject::Epetra_DistObject(const Epetra_DistObject& Source)
  : Epetra_Object(Source),
    Map_(Source.Map_),
    Comm_(&Map_.Comm()),
    DistributedGlobal_(Source.DistributedGlobal_),
    Exports_(0),
    Imports_(0),
    LenExports_(0),
    LenImports_(0),
    Sizes_(0)
{}

//=============================================================================
Epetra_DistObject::~Epetra_DistObject(){


  if (LenExports_!=0) {
    delete[] Exports_;
    Exports_ = 0;
    LenExports_ = 0;
  }
  if (LenImports_!=0) {
    delete[] Imports_;
    Imports_ = 0;
    LenImports_ = 0;
  }

  if (Sizes_!=0) delete [] Sizes_;
  Sizes_ = 0;
}

//=============================================================================
int Epetra_DistObject::Import(const Epetra_SrcDistObject& A, 
			      const Epetra_Import& Importer,
			      Epetra_CombineMode CombineMode,
                              const Epetra_OffsetIndex * Indexor) 
{

  if (!Map_.SameAs(Importer.TargetMap())) EPETRA_CHK_ERR(-2);
  if (!A.Map().SameAs(Importer.SourceMap())) EPETRA_CHK_ERR(-3);
  
  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumRemoteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int* ExportLIDs = Importer.ExportLIDs();
  int* RemoteLIDs = Importer.RemoteLIDs();
  int* PermuteToLIDs = Importer.PermuteToLIDs();
  int* PermuteFromLIDs = Importer.PermuteFromLIDs();

  EPETRA_CHK_ERR(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
			    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, 
			    LenExports_, Exports_, LenImports_, Imports_, Importer.Distributor(), 
			    false, Indexor));
  return(0);
}

//=============================================================================
int Epetra_DistObject::Export(const Epetra_SrcDistObject& A, 
			      const Epetra_Export& Exporter,
			      Epetra_CombineMode CombineMode,
                              const Epetra_OffsetIndex * Indexor) 
{

  if (!Map_.SameAs(Exporter.TargetMap())) EPETRA_CHK_ERR(-2);
  if (!A.Map().SameAs(Exporter.SourceMap())) EPETRA_CHK_ERR(-3);
  
  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumRemoteIDs();
  int NumExportIDs = Exporter.NumExportIDs();
  int* ExportLIDs = Exporter.ExportLIDs();
  int* RemoteLIDs = Exporter.RemoteLIDs();
  int* PermuteToLIDs = Exporter.PermuteToLIDs();
  int* PermuteFromLIDs = Exporter.PermuteFromLIDs();

  EPETRA_CHK_ERR(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
			    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs,
			    LenExports_, Exports_,LenImports_, Imports_, Exporter.Distributor(), 
			    false, Indexor));
  return(0);
}

//=============================================================================
int Epetra_DistObject::Import(const Epetra_SrcDistObject& A, 
			      const Epetra_Export& Exporter,
			      Epetra_CombineMode CombineMode,
                              const Epetra_OffsetIndex * Indexor) 
{

  if (!Map_.SameAs(Exporter.SourceMap())) EPETRA_CHK_ERR(-2);
  if (!A.Map().SameAs(Exporter.TargetMap())) EPETRA_CHK_ERR(-3);
  
  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumExportIDs();
  int NumExportIDs = Exporter.NumRemoteIDs();
  int* ExportLIDs = Exporter.RemoteLIDs();
  int* RemoteLIDs = Exporter.ExportLIDs();
  int* PermuteToLIDs = Exporter.PermuteFromLIDs();
  int* PermuteFromLIDs = Exporter.PermuteToLIDs();

  EPETRA_CHK_ERR(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs, 
			    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs,
			    LenImports_, Imports_, LenExports_, Exports_, Exporter.Distributor(), 
			    true, Indexor));
  return(0);
}

//=============================================================================
int Epetra_DistObject::Export(const Epetra_SrcDistObject& A, 
			      const Epetra_Import& Importer,
			      Epetra_CombineMode CombineMode,
                              const Epetra_OffsetIndex * Indexor) 
{

  if (!Map_.SameAs(Importer.SourceMap())) EPETRA_CHK_ERR(-2);
  if (!A.Map().SameAs(Importer.TargetMap())) EPETRA_CHK_ERR(-3);
  
  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumExportIDs();
  int NumExportIDs = Importer.NumRemoteIDs();
  int* ExportLIDs = Importer.RemoteLIDs();
  int* RemoteLIDs = Importer.ExportLIDs();
  int* PermuteToLIDs = Importer.PermuteFromLIDs();
  int* PermuteFromLIDs = Importer.PermuteToLIDs();

  EPETRA_CHK_ERR(DoTransfer(A, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,  
			    PermuteToLIDs, PermuteFromLIDs,  RemoteLIDs, ExportLIDs,
			    LenImports_, Imports_, LenExports_, Exports_, 
			    Importer.Distributor(), true, Indexor));
  return(0);
}

//=============================================================================
int Epetra_DistObject::DoTransfer(const Epetra_SrcDistObject& A, 
				  Epetra_CombineMode CombineMode, 
				  int NumSameIDs, 
				  int NumPermuteIDs, 
				  int NumRemoteIDs, 
				  int NumExportIDs, 
				  int* PermuteToLIDs, 
				  int* PermuteFromLIDs, 
				  int* RemoteLIDs, 
				  int* ExportLIDs,
				  int& LenExports, 
				  char*& Exports,
				  int& LenImports, 
				  char*& Imports, 
				  Epetra_Distributor& Distor, 
				  bool DoReverse,
                                  const Epetra_OffsetIndex * Indexor)
{

  EPETRA_CHK_ERR(CheckSizes(A));

  if (NumSameIDs + NumPermuteIDs > 0) {
    EPETRA_CHK_ERR(CopyAndPermute(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs,Indexor));
  }

  if (CombineMode==Zero) 
    return(0); // All done if CombineMode only involves copying and permuting
  
  int SizeOfPacket; 
  bool VarSizes = false;
  if( NumExportIDs > 0) {
    delete [] Sizes_;
    Sizes_ = new int[NumExportIDs];
  }
  EPETRA_CHK_ERR(PackAndPrepare(A, NumExportIDs, ExportLIDs,
                 LenExports, Exports, SizeOfPacket, Sizes_, VarSizes, Distor));

  if ((DistributedGlobal_ && DoReverse) || (A.Map().DistributedGlobal() && !DoReverse)) {
    if (DoReverse) {
      // Do the exchange of remote data
      if( VarSizes ) {
        EPETRA_CHK_ERR(Distor.DoReverse(Exports, SizeOfPacket, Sizes_, LenImports, Imports));
      }
      else {
        EPETRA_CHK_ERR(Distor.DoReverse(Exports, SizeOfPacket, LenImports, Imports));
      }
    }
    else {
      // Do the exchange of remote data
      if( VarSizes ) {
        EPETRA_CHK_ERR(Distor.Do(Exports, SizeOfPacket, Sizes_, LenImports, Imports));
      }
      else {
        EPETRA_CHK_ERR(Distor.Do(Exports, SizeOfPacket, LenImports, Imports));
      }
    }
    EPETRA_CHK_ERR(UnpackAndCombine(A, NumRemoteIDs, RemoteLIDs, LenImports, Imports, SizeOfPacket, Distor, CombineMode, Indexor));
  }

  return(0);
}

//=============================================================================
void Epetra_DistObject::Print(ostream& os) const {
  int MyPID = Comm().MyPID();
  int NumProc = Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      Comm().PrintInfo(os);
      os << "Length of Export buffer (in chars) = " << LenExports_ << endl;
      os << "Length of Import buffer (in chars) = " << LenImports_ << endl;
      os << flush;
    }
  }
  return;
}

//------------------------------------------------------------------------------
Epetra_DistObject& Epetra_DistObject::operator=(const Epetra_DistObject& src)
{
  (void)src;
  //not currently supported
  bool throw_error = true;
  if (throw_error) {
    throw ReportError("Epetra_DistObject::operator= is not supported.",-1);
  }

  return(*this);
}
