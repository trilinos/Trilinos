//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include <EpetraExt_ConfigDefs.h>
#include <EpetraExt_MatrixMatrix.h>

#include <EpetraExt_MMHelpers.h>

#include <EpetraExt_Transpose_RowMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Directory.h>
#include <Epetra_HashTable.h>
#include <Epetra_Distributor.h>
#include <Epetra_IntSerialDenseVector.h>

#ifdef HAVE_VECTOR
#include <vector>
#endif

#ifdef HAVE_MPI
#include <Epetra_MpiDistributor.h>
#endif


#include <Teuchos_TimeMonitor.hpp>

namespace EpetraExt {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static inline int C_estimate_nnz(const Epetra_CrsMatrix & A, const Epetra_CrsMatrix &B){
  // Follows the NZ estimate in ML's ml_matmatmult.c
  int Aest=(A.NumMyRows()>0)? A.NumMyNonzeros()/A.NumMyRows():100;
  int Best=(B.NumMyRows()>0)? B.NumMyNonzeros()/B.NumMyRows():100;

  int nnzperrow=(int)(sqrt((double)Aest) + sqrt((double)Best) - 1);
  nnzperrow*=nnzperrow;

  return (int)(A.NumMyRows()*nnzperrow*0.75 + 100);
}

// Commented out unused, file-local function.
#if 0
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static inline int auto_resize(std::vector<int> &x,int num_new){
  int newsize=x.size() + EPETRA_MAX((int)x.size(),num_new);
  x.resize(newsize);
  return newsize;
}
#endif // 0

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int aztecoo_and_ml_compatible_map_union(const Epetra_CrsMatrix &B, const LightweightCrsMatrix &Bimport, Epetra_Map*& unionmap, std::vector<int>& Cremotepids,
                                        std::vector<int> &Bcols2Ccols, std::vector<int> &Icols2Ccols)
{
#ifdef HAVE_MPI

#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM M5 CMap 1")));
#endif

 // So we need to merge the ColMap of B and the TargetMap of Bimport in an AztecOO/Ifpack/ML compliant order.
  Epetra_Util util;
  int i,j,MyPID = B.Comm().MyPID(), NumProc = B.Comm().NumProc();
  int Bstart=0, Istart=0, Cstart=0,Pstart=0;

  const Epetra_Map & BColMap       = B.ColMap();
  const Epetra_Map & DomainMap     = B.DomainMap();
  const LightweightMap & IColMap   = Bimport.ColMap_;

  int Nb         = BColMap.NumMyElements();
  int_type * Bgids = 0;
  BColMap.MyGlobalElementsPtr(Bgids);
  int Ni         = IColMap.NumMyElements();
  int_type * Igids    = 0;
  if(Ni>0)
    IColMap.MyGlobalElementsPtr(Igids);

  if((int)Bcols2Ccols.size() != Nb) Bcols2Ccols.resize(Nb);
  if((int)Icols2Ccols.size() != Ni) Icols2Ccols.resize(Ni);

  // Since we're getting called, we know we have to be using an MPI implementation of Epetra.
  // Which means we should have an MpiDistributor for both B and Bimport.
  // Unless all of B's columns are owned by the calling proc (e.g. MueLu for A*Ptent w/ uncoupled aggregation)
  Epetra_MpiDistributor *Distor=0;
  if(B.Importer()) {
    Distor=dynamic_cast<Epetra_MpiDistributor*>(&B.Importer()->Distributor());
    if(!Distor) EPETRA_CHK_ERR(-2);
  }

  // **********************
  // Stage 1: Get the owning PIDs
  // **********************
  // Note: if B doesn't have an importer, the calling proc owns all its colids...
  std::vector<int> Bpids(Nb), Ipids(Ni);
  if(B.Importer()) {EPETRA_CHK_ERR(util.GetPids(*B.Importer(),Bpids,true));}
  else Bpids.assign(Nb,-1);

  if(Ni != (int) Bimport.ColMapOwningPIDs_.size()) {
    EPETRA_CHK_ERR(-21);
  }
  for(i=0;i<Ni;i++){
    Ipids[i] = (Bimport.ColMapOwningPIDs_[i]==MyPID)?(-1):(Bimport.ColMapOwningPIDs_[i]);
  }

  // **********************
  // Stage 2: Allocate memory (make things too big)
  // **********************
  int Csize=Nb+Ni;
  int Psize=Nb+Ni;
  std::vector<int_type> Cgids(Csize);
  Cremotepids.resize(Psize);

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM M5 CMap 2")));
#endif

  // **********************
  // Stage 3: Local Unknowns
  // **********************
  if(!B.Importer() || (B.Importer()->NumSameIDs() == DomainMap.NumMyElements())) {
    // B's colmap has all of the domain elements.  We can just copy those into the start of the array.
    DomainMap.MyGlobalElements(Cgids.size() ? &Cgids[0] : 0);
    Cstart=DomainMap.NumMyElements();
    Bstart=DomainMap.NumMyElements();

    for(i=0; i<DomainMap.NumMyElements(); i++) Bcols2Ccols[i] = i;
  }
  else {
    // There are more entries in the DomainMap than B's ColMap.  So we stream through both B and Bimport for the copy.
    int NumDomainElements     = DomainMap.NumMyElements();
    for(i = 0; i < NumDomainElements; i++) {
      int_type GID = (int_type) DomainMap.GID64(i);
      int LID = BColMap.LID(GID);
      // B has this guy
      if(LID!=-1) {
        Bcols2Ccols[LID]=Cstart;
        Cgids[Cstart] = GID;
        Cstart++;
        Bstart++;
      }
      else {
        // B import has this guy
        LID = IColMap.LID(GID);
        if(LID!=-1) {
          Icols2Ccols[LID]=Cstart;
          Cgids[Cstart] = GID;
          Cstart++;
        }
      }
    }
  }

  // Now advance Ilast to the last owned unknown in Bimport
  for(i=0,j=0; i<Ni && Ipids[i]==-1; i++) {
    while(Cgids[j]!=Igids[i]) j++;
    Icols2Ccols[i]=j;
  }
  Istart=i;


#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM M5 CMap 3")));
#endif


  // **********************
  // Stage 4: Processor-by-processor set_union
  // **********************
  // NOTE: Intial sizes for Btemp/Itemp from B's distributor.  This should be exact for Btemp and a decent guess for Itemp
  int initial_temp_length = 0;
  const int * lengths_from=0;
  if(Distor) {
    lengths_from= Distor->LengthsFrom();
    for(i=0; i < Distor->NumReceives(); i++) initial_temp_length += lengths_from[i];
  }
  else initial_temp_length=100;

  std::vector<int_type> Btemp(initial_temp_length),  Itemp(initial_temp_length);
  std::vector<int> Btemp2(initial_temp_length), Itemp2(initial_temp_length);


  while (Bstart < Nb || Istart < Ni) {
    int Bproc=NumProc+1, Iproc=NumProc+1, Cproc;

    // Find the next active processor ID
    if(Bstart < Nb) Bproc=Bpids[Bstart];
    if(Istart < Ni) Iproc=Ipids[Istart];

    Cproc = (Bproc < Iproc)?Bproc:Iproc;

    if(Bproc == Cproc && Iproc != Cproc) {
      // Only B has this processor.  Copy the data.
      // B: Find the beginning of the next processor
      for(i=Bstart; i<Nb && Bpids[i]==Bproc; i++) {}
      int Bnext=i;

      // Copy data to C
      int tCsize = Bnext-Bstart;
      if(Btemp.size() < (size_t)tCsize) {Btemp2.resize(tCsize);}

      for(i=Bstart; i<Bnext; i++) {
        Cremotepids[i-Bstart+Pstart] = Cproc;
        Cgids[i-Bstart+Cstart]       = Bgids[i];
        Btemp2[i-Bstart]             = i;
      }

      // Sort & record reindexing
      int *Bptr2 = Btemp2.size() ? &Btemp2[0] : 0;
      util.Sort(true, tCsize, &Cgids[Cstart], 0, 0, 1, &Bptr2, 0, 0);

      for(i=0, j=Cstart; i<tCsize; i++){
        while(Cgids[j] != Bgids[Btemp2[i]]) j++;
        Bcols2Ccols[Btemp2[i]] =  j;
      }
      Cstart+=tCsize;
      Pstart+=tCsize;
      Bstart=Bnext;
    }
    else if(Bproc != Cproc && Iproc == Cproc) {
      // Only I has this processor.  Copy the data.
      // I: Find the beginning of the next processor
      for(i=Istart; i<Ni && Ipids[i]==Iproc; i++) {}
      int Inext=i;

      // Copy data to C
      int tCsize = Inext-Istart;
      if(Itemp.size() < (size_t)tCsize) {Itemp2.resize(tCsize);}

      for(i=Istart; i<Inext; i++) {
        Cremotepids[i-Istart+Pstart] = Cproc;
        Cgids[i-Istart+Cstart]       = Igids[i];
        Itemp2[i-Istart]             = i;
      }

      // Sort & record reindexing
      int *Iptr2 = Itemp2.size() ? &Itemp2[0] : 0;
      util.Sort(true, tCsize, &Cgids[Cstart], 0, 0, 1, &Iptr2, 0, 0);

      for(i=0, j=Cstart; i<tCsize; i++){
        while(Cgids[j] != Igids[Itemp2[i]]) j++;
        Icols2Ccols[Itemp2[i]] =  j;
      }
      Cstart+=tCsize;
      Pstart+=tCsize;
      Istart=Inext;
    }
    else {
      // Both B and I have this processor, so we need to do a set_union.  So we need to sort.
      int Bnext, Inext;
      // B: Find the beginning of the next processor
      for(i=Bstart; i<Nb && Bpids[i]==Bproc; i++) {}
      Bnext=i;

      // I: Find the beginning of the next processor
      for(i=Istart; i<Ni && Ipids[i]==Iproc; i++) {}
      Inext=i;

      // Copy data to temp
      int tBsize = Bnext-Bstart;
      int tIsize = Inext-Istart;
      //      int tCsize = tBsize+tIsize;
      if(Btemp.size() < (size_t)tBsize) {Btemp.resize(tBsize); Btemp2.resize(tBsize);}
      if(Itemp.size() < (size_t)tIsize) {Itemp.resize(tIsize); Itemp2.resize(tIsize);}

      for(i=Bstart; i<Bnext; i++) {Btemp[i-Bstart]=Bgids[i]; Btemp2[i-Bstart]=i;}
      for(i=Istart; i<Inext; i++) {Itemp[i-Istart]=Igids[i]; Itemp2[i-Istart]=i;}

      // Sort & set_union
      int *Bptr2 = Btemp2.size() ? &Btemp2[0] : 0; int *Iptr2 = Itemp2.size() ? &Itemp2[0] : 0;
      util.Sort(true, tBsize, Btemp.size() ? &Btemp[0] : 0, 0, 0, 1, &Bptr2, 0, 0);
      util.Sort(true, tIsize, Itemp.size() ? &Itemp[0] : 0, 0, 0, 1, &Iptr2, 0, 0);
      typename std::vector<int_type>::iterator mycstart = Cgids.begin()+Cstart;
      typename std::vector<int_type>::iterator last_el=std::set_union(Btemp.begin(),Btemp.begin()+tBsize,Itemp.begin(),Itemp.begin()+tIsize,mycstart);

      for(i=0, j=Cstart; i<tBsize; i++){
        while(Cgids[j] != Bgids[Btemp2[i]]) j++;
        Bcols2Ccols[Btemp2[i]] =  j;
      }

      for(i=0, j=Cstart; i<tIsize; i++){
        while(Cgids[j] != Igids[Itemp2[i]]) j++;
        Icols2Ccols[Itemp2[i]] =  j;
      }

      for(i=Pstart; i<(last_el - mycstart) + Pstart; i++) Cremotepids[i]=Cproc;
      Cstart = (last_el - mycstart) + Cstart;
      Pstart = (last_el - mycstart) + Pstart;
      Bstart=Bnext;
      Istart=Inext;
    }
  } // end while

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM M5 CMap 4")));
#endif

  // Resize the RemotePIDs down
  Cremotepids.resize(Pstart);

  // **********************
  // Stage 5: Call constructor
  // **********************
  // Make the map
  unionmap=new Epetra_Map((int_type) -1,Cstart,Cgids.size() ? &Cgids[0] : 0, (int_type) B.ColMap().IndexBase64(),
          B.Comm(),B.ColMap().DistributedGlobal(),(int_type) B.ColMap().MinAllGID64(),(int_type) B.ColMap().MaxAllGID64());
  return 0;
#else
  return -1;
#endif
  }

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
// Provide a "resize" operation for double*'s.
inline void resize_doubles(int nold,int nnew,double*& d){
  if(nnew > nold){
    double *tmp = new double[nnew];
    for(int i=0; i<nold; i++)
      tmp[i]=d[i];
    delete [] d;
    d=tmp;
  }
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int  mult_A_B_newmatrix(const Epetra_CrsMatrix & A,
                        const Epetra_CrsMatrix & B,
                        const CrsMatrixStruct& Bview,
                        std::vector<int> & Bcol2Ccol,
                        std::vector<int> & Bimportcol2Ccol,
                        std::vector<int>& Cremotepids,
                        Epetra_CrsMatrix& C,
			bool keep_all_hard_zeros
){
#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Newmatrix SerialCore")));
#endif

  // *****************************
  // Improved Parallel Gustavson in Local IDs
  // *****************************
  const Epetra_Map * colmap_C = &(C.ColMap());


  int NumMyDiagonals=0; // Counter to speed up ESFC

  int m=A.NumMyRows();
  int n=colmap_C->NumMyElements();
  int i,j,k;

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));

  // DataPointers for B, Bimport
  int *Browptr, *Bcolind;
  double *Bvals;
  EPETRA_CHK_ERR(B.ExtractCrsDataPointers(Browptr,Bcolind,Bvals));

  int *Irowptr=0, *Icolind=0;
  double *Ivals=0;
  if(Bview.importMatrix){
    Irowptr = &Bview.importMatrix->rowptr_[0];
    Icolind = (Bview.importMatrix->colind_.size()>0)?(&Bview.importMatrix->colind_[0]):0;
    Ivals   = (Bview.importMatrix->vals_.size()>0)?(&Bview.importMatrix->vals_[0]):0;
  }

  // MemorySetup: If somebody else is sharing this C's graphdata, make a new one.
  // This is needed because I'm about to walk all over the CrsGrapData...
  C.ExpertMakeUniqueCrsGraphData();

  // The status array will contain the index into colind where this entry was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with -1's indicating that there are no entries yet.
  std::vector<int> c_status(n, -1);

  // Classic csr assembly (low memory edition)
  int CSR_alloc=C_estimate_nnz(A,B);
  if(CSR_alloc < n) CSR_alloc = n;
  int CSR_ip=0,OLD_ip=0;
  Epetra_IntSerialDenseVector & CSR_rowptr = C.ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector & CSR_colind = C.ExpertExtractIndices();
  double *&                     CSR_vals   = C.ExpertExtractValues();

  CSR_rowptr.Resize(m+1);
  CSR_colind.Resize(CSR_alloc);
  resize_doubles(0,CSR_alloc,CSR_vals);

  // Static Profile stuff
  std::vector<int> NumEntriesPerRow(m);

  // For each row of A/C
  for(i=0; i<m; i++){
    bool found_diagonal=false;
    CSR_rowptr[i]=CSR_ip;

    for(k=Arowptr[i]; k<Arowptr[i+1]; k++){
      int Ak      = Acolind[k];
      double Aval = Avals[k];
      if(!keep_all_hard_zeros && Aval==0) continue;

      if(Bview.targetMapToOrigRow[Ak] != -1){
        // Local matrix
        int Bk = Bview.targetMapToOrigRow[Ak];
        for(j=Browptr[Bk]; j<Browptr[Bk+1]; ++j) {
          int Cj=Bcol2Ccol[Bcolind[j]];

          if(Cj==i && !found_diagonal) {found_diagonal=true; NumMyDiagonals++;}

          if(c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]=Aval*Bvals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]+=Aval*Bvals[j];
        }
      }
      else{
        // Remote matrix
        int Ik = Bview.targetMapToImportRow[Ak];
        for(j=Irowptr[Ik]; j<Irowptr[Ik+1]; ++j) {
          int Cj=Bimportcol2Ccol[Icolind[j]];

          if(Cj==i && !found_diagonal) {found_diagonal=true; NumMyDiagonals++;}

          if(c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]=Aval*Ivals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]+=Aval*Ivals[j];
        }
      }
    }
    NumEntriesPerRow[i]=CSR_ip-CSR_rowptr[i];

    // Resize for next pass if needed
    if(CSR_ip + n > CSR_alloc){
      resize_doubles(CSR_alloc,2*CSR_alloc,CSR_vals);
      CSR_alloc*=2;
      CSR_colind.Resize(CSR_alloc);
    }
    OLD_ip=CSR_ip;
  }

  CSR_rowptr[m]=CSR_ip;

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Newmatrix Final Sort")));
#endif

  // Sort the entries
  Epetra_Util::SortCrsEntries(m, &CSR_rowptr[0], &CSR_colind[0], &CSR_vals[0]);

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Newmatrix Fast IE")));
#endif

  // Do a fast build of C's importer
  Epetra_Import * Cimport=0;
  int *RemotePIDs = Cremotepids.size()?&Cremotepids[0]:0;
  int NumExports=0;
  int *ExportLIDs=0, *ExportPIDs=0;
  if(Bview.importMatrix) {
    NumExports = Bview.importMatrix->ExportLIDs_.size();
    ExportLIDs = Bview.importMatrix->ExportLIDs_.size()?&Bview.importMatrix->ExportLIDs_[0]:0;
    ExportPIDs = Bview.importMatrix->ExportPIDs_.size()?&Bview.importMatrix->ExportPIDs_[0]:0;
  }
  else if(B.Importer()) {
    // Grab the exports from B proper
    NumExports = B.Importer()->NumExportIDs();
    ExportLIDs = B.Importer()->ExportLIDs();
    ExportPIDs = B.Importer()->ExportPIDs();
  }


  if(B.Importer() && C.ColMap().SameAs(B.ColMap()))
    Cimport = new Epetra_Import(*B.Importer()); // Because the domain maps are the same
  else if(!C.ColMap().SameAs(B.DomainMap()))
    Cimport = new Epetra_Import(C.ColMap(),B.DomainMap(),Cremotepids.size(),RemotePIDs,NumExports,ExportLIDs,ExportPIDs);


#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Newmatrix ESFC")));
#endif

  // Update the CrsGraphData
  C.ExpertStaticFillComplete(B.DomainMap(),A.RangeMap(),Cimport,0,NumMyDiagonals);

  return 0;
}




/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int mult_A_B_reuse(const Epetra_CrsMatrix & A,
                   const Epetra_CrsMatrix & B,
                   CrsMatrixStruct& Bview,
                   std::vector<int> & Bcol2Ccol,
                   std::vector<int> & Bimportcol2Ccol,
                   Epetra_CrsMatrix& C,
		   bool keep_all_hard_zeros){

#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Reuse SerialCore")));
#endif

  // *****************************
  // Improved Parallel Gustavson in Local IDs
  // *****************************
  const Epetra_Map * colmap_C = &(C.ColMap());

  int m=A.NumMyRows();
  int n=colmap_C->NumMyElements();
  int i,j,k;

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));

  // DataPointers for B, Bimport
  int *Browptr, *Bcolind;
  double *Bvals;
  EPETRA_CHK_ERR(B.ExtractCrsDataPointers(Browptr,Bcolind,Bvals));

  int *Irowptr=0, *Icolind=0;
  double *Ivals=0;
  if(Bview.importMatrix){
    Irowptr = &Bview.importMatrix->rowptr_[0];
    Icolind = &Bview.importMatrix->colind_[0];
    Ivals   = &Bview.importMatrix->vals_[0];
  }

  // DataPointers for C
  int *CSR_rowptr, *CSR_colind;
  double *CSR_vals;
  EPETRA_CHK_ERR(C.ExtractCrsDataPointers(CSR_rowptr,CSR_colind,CSR_vals));


  // The status array will contain the index into colind where this dude was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with -1's indicating that there are no entries yet.
  std::vector<int> c_status(n, -1);

  // Classic csr assembly
  int CSR_alloc=CSR_rowptr[m] - CSR_rowptr[0];
  int CSR_ip=0,OLD_ip=0;

 // For each row of A/C
  for(i=0; i<m; i++){
    for(k=Arowptr[i]; k<Arowptr[i+1]; k++){
      int Ak=Acolind[k];
      double Aval = Avals[k];
      if(!keep_all_hard_zeros && Aval==0) continue;

      if(Bview.targetMapToOrigRow[Ak] != -1){
        // Local matrix
        int Bk = Bview.targetMapToOrigRow[Ak];
        for(j=Browptr[Bk]; j<Browptr[Bk+1]; ++j) {
          int Cj=Bcol2Ccol[Bcolind[j]];

          if(c_status[Cj]<OLD_ip){
            // New entry
            if(CSR_ip >= CSR_alloc) EPETRA_CHK_ERR(-13);
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]=Aval*Bvals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]+=Aval*Bvals[j];
        }
      }
      else{
        // Remote matrix
        int Ik = Bview.targetMapToImportRow[Ak];
        for(j=Irowptr[Ik]; j<Irowptr[Ik+1]; ++j) {
          int Cj=Bimportcol2Ccol[Icolind[j]];

          if(c_status[Cj]<OLD_ip){
            // New entry
            if(CSR_ip >= CSR_alloc) EPETRA_CHK_ERR(-14);
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]=Aval*Ivals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]+=Aval*Ivals[j];
        }
      }
    }
    OLD_ip=CSR_ip;
  }

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Reuse Final Sort")));
#endif

  // Sort the entries
  Epetra_Util::SortCrsEntries(m, &CSR_rowptr[0], &CSR_colind[0], &CSR_vals[0]);

  return 0;
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
//kernel method for computing the local portion of C = A*B
template<typename int_type>
  int mult_A_B_general(const Epetra_CrsMatrix & A,
                       CrsMatrixStruct & Aview,
                       const Epetra_CrsMatrix & B,
                       CrsMatrixStruct& Bview,
                       Epetra_CrsMatrix& C,
                       bool call_FillComplete_on_result,
		       bool keep_all_hard_zeros)
{
  int C_firstCol = Bview.colMap->MinLID();
  int C_lastCol = Bview.colMap->MaxLID();

  int C_firstCol_import = 0;
  int C_lastCol_import = -1;

  int_type* bcols = 0;
  Bview.colMap->MyGlobalElementsPtr(bcols);
  int_type* bcols_import = NULL;
  if (Bview.importMatrix != NULL) {
    C_firstCol_import = Bview.importMatrix->ColMap_.MinLID();
    C_lastCol_import = Bview.importMatrix->ColMap_.MaxLID();
    Bview.importMatrix->ColMap_.MyGlobalElementsPtr(bcols_import);
  }

  int C_numCols = C_lastCol - C_firstCol + 1;
  int C_numCols_import = C_lastCol_import - C_firstCol_import + 1;

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  // Allocate workspace memory
  double* dwork = new double[C_numCols];
  int_type* iwork = new int_type[C_numCols];
  int_type *c_cols=iwork;
  double *c_vals=dwork;
  int *c_index=new int[C_numCols];

  int  i, j, k;
  bool C_filled=C.Filled();

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));

  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.

  // Mark indices as empty w/ -1
  for(k=0;k<C_numCols;k++) c_index[k]=-1;

  //loop over the rows of A.
  for(i=0; i<A.NumMyRows(); ++i) {

    int_type global_row = (int_type) Aview.rowMap->GID64(i);

    //loop across the i-th row of A and for each corresponding row
    //in B, loop across columns and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i. In other words,
    //as we stride across B(k,:) we're calculating updates for row i of the
    //result matrix C.

    /* Outline of the revised, ML-inspired algorithm

    C_{i,j} = \sum_k A_{i,k} B_{k,j}

    This algorithm uses a "middle product" formulation, with the loop ordering of
    i, k, j.  This means we compute a row of C at a time, but compute partial sums of
    each entry in row i until we finish the k loop.

    This algorithm also has a few twists worth documenting.

    1) The first major twist involves the c_index, c_cols and c_vals arrays.  The arrays c_cols
    and c_vals store the *local* column index and values accumulator respectively.  These
    arrays are allocated to a size equal to the max number of local columns in C, namely C_numcols.
    The value c_current tells us how many non-zeros we currently have in this row.

    So how do we take a LCID and find the right accumulator?  This is where the c_index array
    comes in.  At the start (and stop) and the i loop, c_index is filled with -1's.  Now
    whenever we find a LCID in the k loop, we first loop at c_index[lcid].  If this value is
    -1 we haven't seen this entry yet.  In which case we add the appropriate stuff to c_cols
    and c_vals and then set c_index[lcid] to the location of the accumulator (c_current before
    we increment it).  If the value is NOT -1, this tells us the location in the c_vals/c_cols
    arrays (namely c_index[lcid]) where our accumulator lives.

    This trick works because we're working with local ids.  We can then loop from 0 to c_current
    and reset c_index to -1's when we're done, only touching the arrays that have changed.
    While we're at it, we can switch to the global ids so we can call [Insert|SumInto]GlobalValues.
    Thus, the effect of this trick is to avoid passes over the index array.

    2) The second major twist involves handling the remote and local components of B separately.
    (ML doesn't need to do this, because its local ordering scheme is consistent between the local
    and imported components of B.)  Since they have different column maps, they have inconsistent
    local column ids.  This means the "second twist" won't work as stated on both matrices at the
    same time.  While this could be handled any number of ways, I have chosen to do the two parts
    of B separately to make the code easier to read (and reduce the memory footprint of the MMM).
    */

    // Local matrix: Zero Current counts for matrix
    int c_current=0;

    // Local matrix: Do the "middle product"
    for(k=Arowptr[i]; k<Arowptr[i+1]; ++k) {
      int Ak      = Acolind[k];
      double Aval = Avals[k];
      // We're skipping remote entries on this pass.
      if(Bview.remote[Ak] || (!keep_all_hard_zeros && Aval==0)) continue;

      int* Bcol_inds = Bview.indices[Ak];
      double* Bvals_k = Bview.values[Ak];

      for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
        int col=Bcol_inds[j];
        if(c_index[col]<0){
          // We haven't seen this entry before; add it.  (In ML, on
          // the first pass, you haven't seen any of the entries
          // before, so they are added without the check.  Not sure
          // how much more efficient that would be; depends on branch
          // prediction.  We've favored code readability here.)
          c_cols[c_current]=col;
          c_vals[c_current]=Aval*Bvals_k[j];
          c_index[col]=c_current;
          c_current++;
        }
        else{
          // We've already seen this entry; accumulate it.
          c_vals[c_index[col]]+=Aval*Bvals_k[j];
        }
      }
    }
    // Local matrix: Reset c_index and switch c_cols to GIDs
    for(k=0; k<c_current; k++){
      c_index[c_cols[k]]=-1;
      c_cols[k]=bcols[c_cols[k]]; // Switch from local to global IDs.
    }
    // Local matrix: Insert.
    //
    // We should check global error results after the algorithm is
    // through.  It's probably safer just to let the algorithm run all
    // the way through before doing this, since otherwise we have to
    // remember to free all allocations carefully.
    //
    // FIXME (mfh 27 Mar 2015) This code collects error codes, but
    // doesn't do anything with them.  This results in build warnings
    // (set but unused variable).  Thus, I'm commenting out error code
    // collection for now.
#if 0
    int err = C_filled ?
      C.SumIntoGlobalValues(global_row,c_current,c_vals,c_cols)
      :
      C.InsertGlobalValues(global_row,c_current,c_vals,c_cols);
#else
    if (C_filled) {
      C.SumIntoGlobalValues(global_row,c_current,c_vals,c_cols);
    } else {
      C.InsertGlobalValues(global_row,c_current,c_vals,c_cols);
    }
#endif // 0

    // Remote matrix: Zero current counts again for matrix
    c_current=0;

    // Remote matrix: Do the "middle product"
    for(k=Arowptr[i]; k<Arowptr[i+1]; ++k) {
      int Ak      = Acolind[k];
      double Aval = Avals[k];
      // We're skipping local entries on this pass.
      if(!Bview.remote[Ak] || Aval==0) continue;

      int* Bcol_inds = Bview.indices[Ak];
      double* Bvals_k = Bview.values[Ak];

      for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
        int col=Bcol_inds[j];
        if(c_index[col]<0){
          c_cols[c_current]=col;
          c_vals[c_current]=Aval*Bvals_k[j];
          c_index[col]=c_current;
          c_current++;
        }
        else{
          c_vals[c_index[col]]+=Aval*Bvals_k[j];
        }
      }
    }
    // Remote matrix: Reset c_index and switch c_cols to GIDs
    for(k=0; k<c_current; k++){
      c_index[c_cols[k]]=-1;
      c_cols[k]=bcols_import[c_cols[k]];
    }
    // Remove matrix: Insert
    //
    // See above (on error handling).
#if 0
    err = C_filled ?
      C.SumIntoGlobalValues(global_row,c_current,c_vals,c_cols)
      :
      C.InsertGlobalValues(global_row,c_current,c_vals,c_cols);
#else
    if (C_filled) {
      C.SumIntoGlobalValues(global_row,c_current,c_vals,c_cols);
    } else {
      C.InsertGlobalValues(global_row,c_current,c_vals,c_cols);
    }
#endif // 0
  }

  // Since Multiply won't do this
  if(call_FillComplete_on_result)
    C.FillComplete(B.DomainMap(),A.RangeMap());

  delete [] dwork;
  delete [] iwork;
  delete [] c_index;
  return(0);
}





/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int MatrixMatrix::Tmult_A_B(const Epetra_CrsMatrix & A,
                           CrsMatrixStruct & Aview,
                           const Epetra_CrsMatrix & B,
                           CrsMatrixStruct& Bview,
                           Epetra_CrsMatrix& C,
			   bool call_FillComplete_on_result,
			   bool keep_all_hard_zeros){

  int i,rv;
  Epetra_Map* mapunion = 0;
  const Epetra_Map * colmap_B = &(B.ColMap());
  const Epetra_Map * colmap_C = &(C.ColMap());

  std::vector<int> Cremotepids;
  std::vector<int> Bcol2Ccol(B.ColMap().NumMyElements());
  std::vector<int> Bimportcol2Ccol;
  if(Bview.importMatrix) Bimportcol2Ccol.resize(Bview.importMatrix->ColMap_.NumMyElements());

#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM;
#endif

  // If the user doesn't want us to call FillComplete, use the general routine
  if(!call_FillComplete_on_result) {
#ifdef ENABLE_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM General Multiply")));
#endif
    rv=mult_A_B_general<int_type>(A,Aview,B,Bview,C,false,keep_all_hard_zeros);
    return rv;
  }

  // Is this a "clean" matrix
  bool NewFlag=!C.IndicesAreLocal() && !C.IndicesAreGlobal();

  // Does ExtractCrsDataPointers work?
  int *C_rowptr, *C_colind;
  double * C_vals;
  C.ExtractCrsDataPointers(C_rowptr,C_colind,C_vals);
  bool ExtractFailFlag=!C_rowptr || !C_colind || !C_vals;

  // It's a new matrix that hasn't been fill-completed, use the general routine
  if(!NewFlag && ExtractFailFlag){
#ifdef ENABLE_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM General Multiply")));
#endif
    rv=mult_A_B_general<int_type>(A,Aview,B,Bview,C,call_FillComplete_on_result,keep_all_hard_zeros);
    return rv;
  }


#ifdef ENABLE_MMM_TIMINGS
  if(NewFlag) MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Newmatrix CMap")));
  else MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Reuse CMap")));
#endif

  // If new, build & clobber a colmap for C
  if(NewFlag){
    if(Bview.importMatrix) {
      EPETRA_CHK_ERR( aztecoo_and_ml_compatible_map_union<int_type>(B,*Bview.importMatrix,mapunion,Cremotepids,Bcol2Ccol,Bimportcol2Ccol) );
      EPETRA_CHK_ERR( C.ReplaceColMap(*mapunion) );
    }
    else  {
      EPETRA_CHK_ERR( C.ReplaceColMap(B.ColMap()) );
      for(i=0;i<colmap_B->NumMyElements();i++) Bcol2Ccol[i]=i;

      // Copy B's remote list (if any)
      if(B.Importer())
        EPETRA_CHK_ERR( Epetra_Util::GetRemotePIDs(*B.Importer(),Cremotepids));
    }
  }

#ifdef ENABLE_MMM_TIMINGS
  if(NewFlag)  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Newmatrix Lookups")));
  else MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM Reuse Lookups")));
#endif

  // ********************************************
  // Setup Bcol2Ccol / Bimportcol2Ccol lookups
  // ********************************************
  // Note: If we ran the map_union, we have this information already

  if(!NewFlag) {
    if(colmap_B->SameAs(*colmap_C)){
      // Maps are the same: Use local IDs as the hash
      for(i=0;i<colmap_B->NumMyElements();i++)
        Bcol2Ccol[i]=i;
    }
    else {
      // Maps are not the same:  Use the map's hash
      for(i=0;i<colmap_B->NumMyElements();i++){
        Bcol2Ccol[i]=colmap_C->LID((int_type) colmap_B->GID64(i));
        if(Bcol2Ccol[i]==-1) EPETRA_CHK_ERR(-11);
      }
    }

    if(Bview.importMatrix){
      Bimportcol2Ccol.resize(Bview.importMatrix->ColMap_.NumMyElements());
      for(i=0;i<Bview.importMatrix->ColMap_.NumMyElements();i++){
      Bimportcol2Ccol[i]=colmap_C->LID((int_type) Bview.importMatrix->ColMap_.GID64(i));
      if(Bimportcol2Ccol[i]==-1) EPETRA_CHK_ERR(-12);
      }

    }
  }
#ifdef ENABLE_MMM_TIMINGS
  MM=Teuchos::null;
#endif

  // Call the appropriate core routine
  if(NewFlag) {
    EPETRA_CHK_ERR(mult_A_B_newmatrix<int_type>(A,B,Bview,Bcol2Ccol,Bimportcol2Ccol,Cremotepids,C,keep_all_hard_zeros));
  }
  else {
    // This always has a real map
    EPETRA_CHK_ERR(mult_A_B_reuse<int_type>(A,B,Bview,Bcol2Ccol,Bimportcol2Ccol,C,keep_all_hard_zeros));
  }

  // Cleanup
  delete mapunion;
  return 0;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int MatrixMatrix::mult_A_B(const Epetra_CrsMatrix & A,
                           CrsMatrixStruct & Aview,
                           const Epetra_CrsMatrix & B,
                           CrsMatrixStruct& Bview,
                           Epetra_CrsMatrix& C,
                           bool call_FillComplete_on_result,\
			   bool keep_all_hard_zeros){

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(A.RowMap().GlobalIndicesInt() && B.RowMap().GlobalIndicesInt()) {
    return Tmult_A_B<int>(A, Aview, B, Bview, C, call_FillComplete_on_result, keep_all_hard_zeros);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(A.RowMap().GlobalIndicesLongLong() && B.RowMap().GlobalIndicesLongLong()) {
    return Tmult_A_B<long long>(A, Aview, B, Bview, C, call_FillComplete_on_result, keep_all_hard_zeros);
  }
  else
#endif
    throw "EpetraExt::MatrixMatrix::mult_A_B: GlobalIndices type unknown";
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int jacobi_A_B_reuse(double omega,
                     const Epetra_Vector & Dinv,
                     const Epetra_CrsMatrix & A,
                     const Epetra_CrsMatrix & B,
                     CrsMatrixStruct& Bview,
                     std::vector<int> & Bcol2Ccol,
                     std::vector<int> & Bimportcol2Ccol,
                     Epetra_CrsMatrix& C){

#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi Reuse SerialCore")));
#endif

  // *****************************
  // Improved Parallel Gustavson in Local IDs
  // *****************************
  const Epetra_Map * colmap_C = &(C.ColMap());

  int m=A.NumMyRows();
  int n=colmap_C->NumMyElements();
  int i,j,k;

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));

  // DataPointers for B, Bimport
  int *Browptr, *Bcolind;
  double *Bvals;
  EPETRA_CHK_ERR(B.ExtractCrsDataPointers(Browptr,Bcolind,Bvals));

  int *Irowptr=0, *Icolind=0;
  double *Ivals=0;
  if(Bview.importMatrix){
    Irowptr = &Bview.importMatrix->rowptr_[0];
    Icolind = &Bview.importMatrix->colind_[0];
    Ivals   = &Bview.importMatrix->vals_[0];
  }

  // Data pointer for Dinv
  const double *Dvals = Dinv.Values();

  // DataPointers for C
  int *CSR_rowptr, *CSR_colind;
  double *CSR_vals;
  EPETRA_CHK_ERR(C.ExtractCrsDataPointers(CSR_rowptr,CSR_colind,CSR_vals));


  // The status array will contain the index into colind where this dude was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with -1's indicating that there are no entries yet.
  std::vector<int> c_status(n, -1);

  // Classic csr assembly
  int CSR_alloc=CSR_rowptr[m] - CSR_rowptr[0];
  int CSR_ip=0,OLD_ip=0;

  // For each row of C
  for(i=0; i<m; i++){
    double Dval = Dvals[i];

    // Entries of B
    for(k=Browptr[i]; k<Browptr[i+1]; k++){
      //      int Bk      = Bcolind[k];
      double Bval = Bvals[k];
      if(Bval==0) continue;
      int Ck=Bcol2Ccol[Bcolind[k]];

      // Assume no repeated entries in B
      c_status[Ck]=CSR_ip;
      CSR_colind[CSR_ip]=Ck;
      CSR_vals[CSR_ip]= Bvals[k];
      CSR_ip++;
    }

    // Entries of -omega * Dinv * A * B
    for(k=Arowptr[i]; k<Arowptr[i+1]; k++){
      int Ak=Acolind[k];
      double Aval = Avals[k];
      if(Aval==0) continue;

      if(Bview.targetMapToOrigRow[Ak] != -1){
        // Local matrix
        int Bk = Bview.targetMapToOrigRow[Ak];
        for(j=Browptr[Bk]; j<Browptr[Bk+1]; ++j) {
          int Cj=Bcol2Ccol[Bcolind[j]];

          if(c_status[Cj]<OLD_ip){
            // New entry
            if(CSR_ip >= CSR_alloc) EPETRA_CHK_ERR(-13);
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]= - omega * Dval * Aval * Bvals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]-= omega * Dval * Aval * Bvals[j];
        }
      }
      else{
        // Remote matrix
        int Ik = Bview.targetMapToImportRow[Ak];
        for(j=Irowptr[Ik]; j<Irowptr[Ik+1]; ++j) {
          int Cj=Bimportcol2Ccol[Icolind[j]];

          if(c_status[Cj]<OLD_ip){
            // New entry
            if(CSR_ip >= CSR_alloc) EPETRA_CHK_ERR(-14);
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]= - omega * Dval * Aval * Ivals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]-=omega * Dval * Aval * Ivals[j];
        }
      }
    }
    OLD_ip=CSR_ip;
  }

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi Reuse Final Sort")));
#endif
  // Sort the entries
  Epetra_Util::SortCrsEntries(m, &CSR_rowptr[0], &CSR_colind[0], &CSR_vals[0]);

  return 0;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int jacobi_A_B_newmatrix(double omega,
                         const Epetra_Vector & Dinv,
                         const Epetra_CrsMatrix & A,
                         const Epetra_CrsMatrix & B,
                         CrsMatrixStruct& Bview,
                         std::vector<int> & Bcol2Ccol,
                         std::vector<int> & Bimportcol2Ccol,
                         std::vector<int>& Cremotepids,
                         Epetra_CrsMatrix& C)
{
#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi NewMatrix SerialCore")));
#endif

  // *****************************
  // Improved Parallel Gustavson in Local IDs
  // *****************************
  const Epetra_Map * colmap_C = &(C.ColMap());
  int NumMyDiagonals=0; // Counter to speed up ESFC

  int m=A.NumMyRows();
  int n=colmap_C->NumMyElements();
  int i,j,k;

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));

  // DataPointers for B, Bimport
  int *Browptr, *Bcolind;
  double *Bvals;
  EPETRA_CHK_ERR(B.ExtractCrsDataPointers(Browptr,Bcolind,Bvals));

  // Data pointer for Dinv
  const double *Dvals = Dinv.Values();

  int *Irowptr=0, *Icolind=0;
  double *Ivals=0;
  if(Bview.importMatrix){
    Irowptr = &Bview.importMatrix->rowptr_[0];
    Icolind = (Bview.importMatrix->colind_.size()>0)?(&Bview.importMatrix->colind_[0]):0;
    Ivals   = (Bview.importMatrix->vals_.size()>0)?(&Bview.importMatrix->vals_[0]):0;
  }

  // MemorySetup: If somebody else is sharing this C's graphdata, make a new one.
  // This is needed because I'm about to walk all over the CrsGrapData...
  C.ExpertMakeUniqueCrsGraphData();

  // The status array will contain the index into colind where this entry was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with -1's indicating that there are no entries yet.
  std::vector<int> c_status(n, -1);

  // Classic csr assembly (low memory edition)
  int CSR_alloc=C_estimate_nnz(A,B);
  if(CSR_alloc < B.NumMyNonzeros()) CSR_alloc = B.NumMyNonzeros(); // update for Jacobi
  int CSR_ip=0,OLD_ip=0;
  Epetra_IntSerialDenseVector & CSR_rowptr = C.ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector & CSR_colind = C.ExpertExtractIndices();
  double *&                     CSR_vals   = C.ExpertExtractValues();

  CSR_rowptr.Resize(m+1);
  CSR_colind.Resize(CSR_alloc);
  resize_doubles(0,CSR_alloc,CSR_vals);

  // Static Profile stuff
  std::vector<int> NumEntriesPerRow(m);

  // For each row of C
  for(i=0; i<m; i++){
    bool found_diagonal=false;
    CSR_rowptr[i]=CSR_ip;
    double Dval = Dvals[i];

    // Entries of B
    for(k=Browptr[i]; k<Browptr[i+1]; k++){
      //int Bk      = Bcolind[k];
      double Bval = Bvals[k];
      if(Bval==0) continue;
      int Ck=Bcol2Ccol[Bcolind[k]];

      // Assume no repeated entries in B
      c_status[Ck]=CSR_ip;
      CSR_colind[CSR_ip]=Ck;
      CSR_vals[CSR_ip]= Bvals[k];
      CSR_ip++;
    }

    // Entries of -omega * Dinv * A * B
    for(k=Arowptr[i]; k<Arowptr[i+1]; k++){
      int Ak      = Acolind[k];
      double Aval = Avals[k];
      if(Aval==0) continue;

      if(Bview.targetMapToOrigRow[Ak] != -1){
        // Local matrix
        int Bk = Bview.targetMapToOrigRow[Ak];
        for(j=Browptr[Bk]; j<Browptr[Bk+1]; ++j) {
          int Cj=Bcol2Ccol[Bcolind[j]];

          if(Cj==i && !found_diagonal) {found_diagonal=true; NumMyDiagonals++;}

          if(c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]= - omega * Dval* Aval * Bvals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]-= omega * Dval * Aval * Bvals[j];
        }
      }
      else{
        // Remote matrix
        int Ik = Bview.targetMapToImportRow[Ak];
        for(j=Irowptr[Ik]; j<Irowptr[Ik+1]; ++j) {
          int Cj=Bimportcol2Ccol[Icolind[j]];

          if(Cj==i && !found_diagonal) {found_diagonal=true; NumMyDiagonals++;}

          if(c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]=CSR_ip;
            CSR_colind[CSR_ip]=Cj;
            CSR_vals[CSR_ip]= - omega * Dval * Aval * Ivals[j];
            CSR_ip++;
          }
          else
            CSR_vals[c_status[Cj]]-= omega * Dval * Aval * Ivals[j];
        }
      }
    }
    NumEntriesPerRow[i]=CSR_ip-CSR_rowptr[i];

    // Resize for next pass if needed
    if(CSR_ip + n > CSR_alloc){
      resize_doubles(CSR_alloc,2*CSR_alloc,CSR_vals);
      CSR_alloc*=2;
      CSR_colind.Resize(CSR_alloc);
    }
    OLD_ip=CSR_ip;
  }

  CSR_rowptr[m]=CSR_ip;

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi NewMatrix Final Sort")));
#endif

  // Sort the entries
  Epetra_Util::SortCrsEntries(m, &CSR_rowptr[0], &CSR_colind[0], &CSR_vals[0]);

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi NewMatrix Fast IE")));
#endif

  // Do a fast build of C's importer
  Epetra_Import * Cimport=0;
  int *RemotePIDs = Cremotepids.size()?&Cremotepids[0]:0;
  int NumExports=0;
  int *ExportLIDs=0, *ExportPIDs=0;
  if(Bview.importMatrix) {
    NumExports = Bview.importMatrix->ExportLIDs_.size();
    ExportLIDs = Bview.importMatrix->ExportLIDs_.size()?&Bview.importMatrix->ExportLIDs_[0]:0;
    ExportPIDs = Bview.importMatrix->ExportPIDs_.size()?&Bview.importMatrix->ExportPIDs_[0]:0;
  }
  else if(B.Importer()) {
    // Grab the exports from B proper
    NumExports = B.Importer()->NumExportIDs();
    ExportLIDs = B.Importer()->ExportLIDs();
    ExportPIDs = B.Importer()->ExportPIDs();
  }

  if(!C.ColMap().SameAs(B.DomainMap()))
    Cimport = new Epetra_Import(C.ColMap(),B.DomainMap(),Cremotepids.size(),RemotePIDs,NumExports,ExportLIDs,ExportPIDs);

#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi NewMatrix ESFC")));
#endif

  // Update the CrsGraphData
  C.ExpertStaticFillComplete(B.DomainMap(),A.RangeMap(),Cimport,0,NumMyDiagonals);

  return 0;
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int MatrixMatrix::Tjacobi_A_B(double omega,
                             const Epetra_Vector & Dinv,
                             const Epetra_CrsMatrix & A,
                             CrsMatrixStruct & Aview,
                             const Epetra_CrsMatrix & B,
                             CrsMatrixStruct& Bview,
                             Epetra_CrsMatrix& C,
                              bool call_FillComplete_on_result){
  int i,rv;
  Epetra_Map* mapunion = 0;
  const Epetra_Map * colmap_B = &(B.ColMap());
  const Epetra_Map * colmap_C = &(C.ColMap());

  std::vector<int> Cremotepids;
  std::vector<int> Bcol2Ccol(B.ColMap().NumMyElements());
  std::vector<int> Bimportcol2Ccol;
  if(Bview.importMatrix) Bimportcol2Ccol.resize(Bview.importMatrix->ColMap_.NumMyElements());

#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM;
#endif

  // If the user doesn't want us to call FillComplete, use the general routine
  if(!call_FillComplete_on_result) {
#ifdef ENABLE_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi General Multiply")));
#endif
    throw std::runtime_error("jacobi_A_B_general not implemented");
    //    rv=mult_A_B_general<int_type>(A,Aview,B,Bview,C,false);
    return rv;
  }

  // Is this a "clean" matrix
  bool NewFlag=!C.IndicesAreLocal() && !C.IndicesAreGlobal();

  // Does ExtractCrsDataPointers work?
  int *C_rowptr, *C_colind;
  double * C_vals;
  C.ExtractCrsDataPointers(C_rowptr,C_colind,C_vals);
  bool ExtractFailFlag=!C_rowptr || !C_colind || !C_vals;

  // It's a new matrix that hasn't been fill-completed, use the general routine
  if(!NewFlag && ExtractFailFlag){
#ifdef ENABLE_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi General Multiply")));
#endif
    throw std::runtime_error("jacobi_A_B_general not implemented");
    //    rv=mult_A_B_general<int_type>(A,Aview,B,Bview,C,call_FillComplete_on_result);
    return rv;
  }

#ifdef ENABLE_MMM_TIMINGS
  if(NewFlag) MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi Newmatrix CMap")));
  else MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi Reuse CMap")));
#endif

  // If new, build & clobber a colmap for C
  if(NewFlag){
    if(Bview.importMatrix) {
      EPETRA_CHK_ERR( aztecoo_and_ml_compatible_map_union<int_type>(B,*Bview.importMatrix,mapunion,Cremotepids,Bcol2Ccol,Bimportcol2Ccol) );
      EPETRA_CHK_ERR( C.ReplaceColMap(*mapunion) );
    }
    else  {
      EPETRA_CHK_ERR( C.ReplaceColMap(B.ColMap()) );
      for(i=0;i<colmap_B->NumMyElements();i++) Bcol2Ccol[i]=i;

      // Copy B's remote list (if any)
      if(B.Importer())
        EPETRA_CHK_ERR( Epetra_Util::GetRemotePIDs(*B.Importer(),Cremotepids));
    }
  }

#ifdef ENABLE_MMM_TIMINGS
  if(NewFlag) MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi Newmatrix Lookups")));
  else MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: Jacobi Reuse Lookups")));
#endif

  // ********************************************
  // Setup Bcol2Ccol / Bimportcol2Ccol lookups
  // ********************************************
  // Note: If we ran the map_union, we have this information already

  if(!NewFlag) {
    if(colmap_B->SameAs(*colmap_C)){
      // Maps are the same: Use local IDs as the hash
      for(i=0;i<colmap_B->NumMyElements();i++)
        Bcol2Ccol[i]=i;
    }
    else {
      // Maps are not the same:  Use the map's hash
      for(i=0;i<colmap_B->NumMyElements();i++){
        Bcol2Ccol[i]=colmap_C->LID((int_type) colmap_B->GID64(i));
        if(Bcol2Ccol[i]==-1) EPETRA_CHK_ERR(-11);
      }
    }

    if(Bview.importMatrix){
      Bimportcol2Ccol.resize(Bview.importMatrix->ColMap_.NumMyElements());
      for(i=0;i<Bview.importMatrix->ColMap_.NumMyElements();i++){
      Bimportcol2Ccol[i]=colmap_C->LID((int_type) Bview.importMatrix->ColMap_.GID64(i));
      if(Bimportcol2Ccol[i]==-1) EPETRA_CHK_ERR(-12);
      }

    }
  }


  // Call the appropriate core routine
  if(NewFlag) {
    EPETRA_CHK_ERR(jacobi_A_B_newmatrix<int_type>(omega,Dinv,A,B,Bview,Bcol2Ccol,Bimportcol2Ccol,Cremotepids,C));
  }
  else {
    // This always has a real map
    EPETRA_CHK_ERR(jacobi_A_B_reuse<int_type>(omega,Dinv,A,B,Bview,Bcol2Ccol,Bimportcol2Ccol,C));
  }

  // Cleanup
  delete mapunion;
  return 0;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int MatrixMatrix::jacobi_A_B(double omega,
                             const Epetra_Vector & Dinv,
                             const Epetra_CrsMatrix & A,
                             CrsMatrixStruct & Aview,
                             const Epetra_CrsMatrix & B,
                             CrsMatrixStruct& Bview,
                             Epetra_CrsMatrix& C,
                             bool call_FillComplete_on_result)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(A.RowMap().GlobalIndicesInt() && B.RowMap().GlobalIndicesInt()) {
    return Tjacobi_A_B<int>(omega,Dinv,A,Aview,B,Bview,C,call_FillComplete_on_result);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(A.RowMap().GlobalIndicesLongLong() && B.RowMap().GlobalIndicesLongLong()) {
      return Tjacobi_A_B<long long>(omega,Dinv,A,Aview,B,Bview,C,call_FillComplete_on_result);
    }
    else
#endif
    throw "EpetraExt::MatrixMatrix::jacobi_A_B: GlobalIndices type unknown";
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<typename int_type>
int MatrixMatrix::Tmult_AT_B_newmatrix(const CrsMatrixStruct & Atransview, const CrsMatrixStruct & Bview, Epetra_CrsMatrix & C,bool keep_all_hard_zeros) {
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef ENABLE_MMM_TIMINGS
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM-T Transpose")));
#endif


  /*************************************************************/
  /* 2/3) Call mult_A_B_newmatrix w/ fillComplete              */
  /*************************************************************/
#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM-T AB-core")));
#endif
  RCP<Epetra_CrsMatrix> Ctemp;

  // If Atrans has no Exporter, we can use C instead of having to create a temp matrix
  bool needs_final_export = Atransview.origMatrix->Exporter() != 0;
  if(needs_final_export)
    Ctemp = rcp(new Epetra_CrsMatrix(Copy,Atransview.origMatrix->RowMap(),Bview.origMatrix->ColMap(),0));
  else {
    EPETRA_CHK_ERR( C.ReplaceColMap(Bview.origMatrix->ColMap()) );
    Ctemp = rcp(&C,false);// don't allow deallocation
  }

  // Multiply
  std::vector<int> Bcol2Ccol(Bview.origMatrix->NumMyCols());
  for(int i=0; i<Bview.origMatrix->NumMyCols(); i++)
    Bcol2Ccol[i]=i;
  std::vector<int> Bimportcol2Ccol,Cremotepids;
  if(Bview.origMatrix->Importer())
    EPETRA_CHK_ERR( Epetra_Util::GetRemotePIDs(*Bview.origMatrix->Importer(),Cremotepids));

  EPETRA_CHK_ERR(mult_A_B_newmatrix<int_type>(*Atransview.origMatrix,*Bview.origMatrix,Bview,
                                              Bcol2Ccol,Bimportcol2Ccol,Cremotepids,
                                              *Ctemp,keep_all_hard_zeros));

  /*************************************************************/
  /* 4) ExportAndFillComplete matrix (if needed)               */
  /*************************************************************/
#ifdef ENABLE_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("EpetraExt: MMM-T ESFC")));
#endif

  if(needs_final_export)
    C.FusedExport(*Ctemp,*Ctemp->Exporter(),&Bview.origMatrix->DomainMap(),&Atransview.origMatrix->RangeMap(),false);

  return 0;
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int MatrixMatrix::mult_AT_B_newmatrix(const CrsMatrixStruct & Atransview, const CrsMatrixStruct & Bview, Epetra_CrsMatrix & C, bool keep_all_hard_zeros) {

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Atransview.origMatrix->RowMap().GlobalIndicesInt() && Bview.origMatrix->RowMap().GlobalIndicesInt()) {
    return Tmult_AT_B_newmatrix<int>(Atransview,Bview,C,keep_all_hard_zeros);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Atransview.origMatrix->RowMap().GlobalIndicesLongLong() && Bview.origMatrix->RowMap().GlobalIndicesLongLong()) {
    return Tmult_AT_B_newmatrix<long long>(Atransview,Bview,C,keep_all_hard_zeros);
  }
  else
#endif
    throw "EpetraExt::MatrixMatrix::mult_AT_B_newmatrix: GlobalIndices type unknown";
}



}//namespace EpetraExt
