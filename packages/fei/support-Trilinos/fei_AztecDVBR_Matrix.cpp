/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>

#ifdef HAVE_FEI_AZTECOO

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include <fei_mpi.h>

#ifndef FEI_SER

#define AZTEC_MPI AZTEC_MPI
#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif

#endif

#include <az_aztec.h>
#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_BlockMap.hpp>
#include <fei_Aztec_Vector.hpp>
#include <fei_AztecDVBR_Matrix.hpp>

//==============================================================================
AztecDVBR_Matrix::AztecDVBR_Matrix(Aztec_BlockMap& map, int* update)
  : amap_(map),
    Amat_(NULL),
    N_update_(0),
    update_(NULL),
    external_(NULL),
    extern_index_(NULL),
    update_index_(NULL),
    data_org_(NULL),
    orderingUpdate_(NULL),
    isLoaded_(false),
    isAllocated_(false),
    localNNZ_(0),
    nnzPerRow_(NULL),
    numRemoteBlocks_(0),
    remoteInds_(NULL),
    remoteBlockSizes_(NULL)
{

    N_update_ = amap_.getNumLocalBlocks();

    update_ = new int[N_update_];
    nnzPerRow_ = new int[N_update_];

    for(int i=0; i<N_update_; i++) {
       update_[i] = update[i];
       nnzPerRow_[i] = 0;
    }

    Amat_ = AZ_matrix_create(N_update_);
    Amat_->matrix_type = AZ_VBR_MATRIX;
    Amat_->matvec =
            (void(*)(double*,double*,AZ_MATRIX_STRUCT*,int*))AZ_VBR_matvec_mult;

    //now we can allocate and fill the rpntr array.
    Amat_->rpntr = NULL;
    calcRpntr();
}

//==============================================================================
AztecDVBR_Matrix::AztecDVBR_Matrix(const AztecDVBR_Matrix& src)
 : amap_(src.amap_),
   Amat_(NULL),
   N_update_(src.N_update_),
   update_(NULL),
   external_(NULL),
   extern_index_(NULL),
   update_index_(NULL),
   data_org_(NULL),
   orderingUpdate_(NULL),
   isLoaded_(src.isLoaded_),
   isAllocated_(src.isAllocated_),
   localNNZ_(src.localNNZ_),
   nnzPerRow_(NULL),
   numRemoteBlocks_(0),
   remoteInds_(NULL),
   remoteBlockSizes_(NULL)
{
//
//This copy constructor just takes a reference to src's amap_ (above), but
//'deep'-copies everything else. i.e., all arrays are allocated here and copied
//from src, etc.
//
//When this constructor completes, this matrix should be in the same state as
//the 'src' matrix. i.e., if src is already allocated, then this matrix will
//have its structure allocated. If src is already loaded (AZ_transform'd) then
//this matrix will have all arrays resulting from AZ_transform allocated too.
//The only thing this matrix won't get is the coefficient data from src.
//
   update_ = new int[N_update_];
   nnzPerRow_ = new int[N_update_];

   int i;
   for(i=0; i<N_update_; i++) {
      update_[i] = src.update_[i];
      nnzPerRow_[i] = src.nnzPerRow_[i];
   }

   Amat_ = AZ_matrix_create(N_update_);
   Amat_->matrix_type = AZ_VBR_MATRIX;
   Amat_->matvec =
            (void(*)(double*,double*,AZ_MATRIX_STRUCT*,int*))AZ_VBR_matvec_mult;

   Amat_->rpntr = NULL;
   calcRpntr();

   if (isAllocated_) {
      Amat_->bpntr = new int[N_update_+1];
      for(i=0; i<N_update_+1; i++) Amat_->bpntr[i] = src.Amat_->bpntr[i];

      int totalNNZBlks = Amat_->bpntr[N_update_];
      Amat_->bindx = new int[totalNNZBlks];
      for(i=0; i<totalNNZBlks; i++) Amat_->bindx[i] = src.Amat_->bindx[i];

      Amat_->indx = new int[totalNNZBlks+1];
      for(i=0; i<totalNNZBlks+1; i++) Amat_->indx[i] = src.Amat_->indx[i];

      Amat_->val = new double[localNNZ_];
      for(i=0; i<localNNZ_; i++) Amat_->val[i] = 0.0;
   }

   if (isLoaded_) {
      int dataOrgLength = src.data_org_[AZ_total_send] + AZ_send_list;
      data_org_ = (int*) AZ_allocate(dataOrgLength * sizeof(int));
      for(i=0; i<dataOrgLength; i++) data_org_[i] = src.data_org_[i];

      Amat_->data_org = data_org_;

      int extLength = src.data_org_[AZ_N_ext_blk];
      external_ = (int*) AZ_allocate(extLength * sizeof(int));
      extern_index_ = (int*) AZ_allocate(extLength * sizeof(int));
      for(i=0; i<extLength; i++) {
         external_[i] = src.external_[i];
         extern_index_[i] = src.extern_index_[i];
      }

      update_index_ = (int*) AZ_allocate(N_update_ * sizeof(int));
      orderingUpdate_ = new int[N_update_];
      for(i=0; i<N_update_; i++) {
         update_index_[i] = src.update_index_[i];
         orderingUpdate_[i] = src.orderingUpdate_[i];
      }

      int cpntrLength = N_update_ + src.data_org_[AZ_N_ext_blk] + 1;
      Amat_->cpntr = (int*) AZ_allocate(cpntrLength * sizeof(int));
      for(i=0; i<cpntrLength; i++) Amat_->cpntr[i] = src.Amat_->cpntr[i];
   }
}

//==============================================================================
AztecDVBR_Matrix::~AztecDVBR_Matrix(){

   if (isAllocated()) {
      delete [] Amat_->val;
      delete [] Amat_->bindx;
      delete [] Amat_->rpntr;
      delete [] Amat_->bpntr;
      delete [] Amat_->indx;

      delete [] remoteInds_;
      delete [] remoteBlockSizes_;

      setAllocated(false);
   }

   if (isLoaded()) {
      free(Amat_->cpntr);
      free(external_);
      free(extern_index_);
      free(update_index_);
      free(data_org_);
      delete [] orderingUpdate_;

      setLoaded(false);
   }

   delete [] update_;
   delete [] nnzPerRow_;
   localNNZ_ = 0;

   AZ_matrix_destroy(&Amat_);
   Amat_ = NULL;
}

//==============================================================================
int AztecDVBR_Matrix::getBlockMaps(Aztec_BlockMap** rowMap,
                                   Aztec_BlockMap** colMap) {
   *rowMap = & amap_;
   *colMap = & amap_;

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::getNumBlocksPerRow(int blkRow, int& nnzBlksPerRow) const {
//
//On return, nnzBlksPerRow will be the number of nonzero blocks in row blkRow.
//
   if (!isAllocated()) return(1);

   int index;

   if (!inUpdate(blkRow, index)) {
      FEI_CERR << "AztecDVBR_Matrix::getNumBlocksPerRow: ERROR: blkRow "
           << blkRow << " not in local update list." << FEI_ENDL;
      return(1);
   }

   nnzBlksPerRow = Amat_->bpntr[index+1] - Amat_->bpntr[index];

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::getNumBlocksPerRow(int* nnzBlksPerRow) const {
//
//nnzBlksPerRow must be allocated by the calling code.
//
//nnzBlksPerRow is a list of length number-of-local-block-rows, and
//nnzBlksPerRow[i] gives the number of nonzeros blocks in row i.
//
   if (!isAllocated()) return(1);

   for(int i=0; i<amap_.getNumLocalBlocks(); i++) {
      nnzBlksPerRow[i] = Amat_->bpntr[i+1] - Amat_->bpntr[i];
   }

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::getNumNonzerosPerRow(int blkRow, int& nnzPerRow) const {
//
//This function finds nnzPerRow, the number of nonzero *point* entries for
//row 'blkRow'.
//
   if (!isAllocated()) return(1);

   int index;

   if (!inUpdate(blkRow, index)) {
      FEI_CERR << "AztecDVBR_Matrix::getNumNonzerosPerRow: ERROR: blkRow "
           << blkRow << " not in local update list." << FEI_ENDL;
      return(1);
   }

   nnzPerRow = nnzPerRow_[index];

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::getNumNonzerosPerRow(int* nnzPerRow) const {
//
//nnzPerRow must be allocated by the calling code,
//length number-of-local-block-rows.
//
//This function fills nnzPerRow so that nnzPerRow[i] gives the
//number of nonzero *point* entries for row i.
//
   if (!isAllocated()) return(1);

   for(int i=0; i<amap_.getNumLocalBlocks(); i++) {
      nnzPerRow[i] = nnzPerRow_[i];
   }

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::getBlockSize(int blkRow, int blkCol,
                                   int& ptRows, int& ptCols) {
   int index;

   ptRows = 0;
   ptCols = 0;

   if (!inUpdate(blkRow, index)) {
      FEI_CERR << "AztecDVBR_Matrix::getBlockSize: ERROR: blkRow "
           << blkRow << " not in local update list." << FEI_ENDL;
      return(1);
   }

   ptRows = Amat_->rpntr[index+1] - Amat_->rpntr[index];

   int local = inUpdate(blkCol, index);

   if (local) {
      ptCols = Amat_->rpntr[index+1] - Amat_->rpntr[index];
   }
   else {
      index = AZ_find_index(blkCol, remoteInds_, numRemoteBlocks_);

      if (index < 0) return(1);
      else ptCols = remoteBlockSizes_[index];
   }

   return(0);
}

//==============================================================================
void AztecDVBR_Matrix::matvec(const Aztec_Vector& x, Aztec_Vector& y) const {
    	
// AztecDVBR_Matrix::matvec --- form y = Ax

    assert(isLoaded());

    int *proc_config = amap_.getProcConfig();
    double *b = (double*)x.startPointer();
    double *c = (double*)y.startPointer();

    AZ_VBR_matvec_mult(b, c, Amat_, proc_config);

    return;
}

//==============================================================================
void AztecDVBR_Matrix::put(double s){

   if (!isAllocated()) return;

   for(int i=0; i<localNNZ_; i++) {
      Amat_->val[i] = s;
   }

   return;
}

//==============================================================================
int AztecDVBR_Matrix::getBlockRow(int blkRow,
                                  double* val,
                                  int* blkColInds,
                                  int numNzBlks) const {

   if (!isAllocated()) return(1);

   int index;

   if (!inUpdate(blkRow, index)) {
      FEI_CERR << "AztecDVBR_Matrix::getBlockRow: ERROR: blkRow "
           << blkRow << " not in local update list." << FEI_ENDL;
      return(1);
   }

   //for each block, we need to find its block column index
   //in the bindx array, then go to that same position in the indx
   //array to find out how many point-entries are in that block.
   //We can then use the indx entry to go to the val array and get
   //the data.

   int nnzBlks = 0, nnzPts = 0;
   int err = getNumBlocksPerRow(blkRow, nnzBlks);
   if (err) return(err);
   err = getNumNonzerosPerRow(blkRow, nnzPts);
   if (err) return(err);
   
   if (numNzBlks != nnzBlks) return(1);

   int offset = 0;
   int blkCounter = 0;
   for(int indb = Amat_->bpntr[index]; indb<Amat_->bpntr[index+1]; indb++) {

      int numEntries = Amat_->indx[indb+1] - Amat_->indx[indb];
      int valOffset = Amat_->indx[indb];

      if (isLoaded()) {
         int ind = Amat_->bindx[indb];
         if (ind < N_update_) {
            blkColInds[blkCounter++] = update_[orderingUpdate_[ind]];
         }
         else {
            blkColInds[blkCounter++] = external_[ind-N_update_];
         }
      }
      else {
         blkColInds[blkCounter++] = Amat_->bindx[indb];
      }

      //ok, now we're ready to get the stuff.
      for(int i=0; i<numEntries; i++) {
         val[offset + i] = Amat_->val[valOffset + i];
      }

      offset += numEntries;
   }

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::putBlockRow(int blkRow,
                                   double* val,
                                   int* blkColInds,
                                   int numNzBlks) const {

   if (!isAllocated()) return(1);

   int index;

   if (!inUpdate(blkRow, index)) {
      FEI_CERR << "AztecDVBR_Matrix::putBlockRow: ERROR: blkRow "
           << blkRow << " not in local update list." << FEI_ENDL;
      return(1);
   }

   //for each incoming block, we need to find its block column index
   //in the bindx array, then go to that same position in the indx
   //array to find out how many (point) entries are in that block.
   //We can then use the indx entry to go to the val array and store
   //the data.

   int offset = 0;
   for(int blk = 0; blk<numNzBlks; blk++) {
      int indb = getBindxOffset(blkColInds[blk],
                                Amat_->bpntr[index], Amat_->bpntr[index+1]-1);

      if (indb < 0) messageAbort("putBlockRow: blk col not found in row.");

      int numEntries = Amat_->indx[indb+1] - Amat_->indx[indb];
      int valOffset = Amat_->indx[indb];

      //ok, now we're ready to store the stuff.
      for(int i=0; i<numEntries; i++) {
         Amat_->val[valOffset + i] = val[offset + i];
      }

      offset += numEntries;
   }

   return(0);
}

//==============================================================================
int AztecDVBR_Matrix::sumIntoBlockRow(int blkRow,
                                       double* val,
                                       int* blkColInds,
                                       int numNzBlks) const
{
  //
  //This function is the same as putBlockRow, except the values
  //are summed into any existing values rather than overwriting
  //them.
  //
  if (!isAllocated()) return(1);

  int index;

  if (!inUpdate(blkRow, index)) {
    FEI_CERR << "AztecDVBR_Matrix::sumIntoBlockRow: ERROR: blkRow "
	 << blkRow << " not in local update list." << FEI_ENDL;
    return(1);
  }

  //for each incoming block, we need to find its block column index
  //in the bindx array, then go to that same position in the indx
  //array to find out how many (point) entries are in that block.
  //We can then use the indx entry to go to the val array and store
  //the data.

  int offset = 0;

  for(int blk = 0; blk<numNzBlks; blk++) {
    int indb = getBindxOffset(blkColInds[blk],
			      Amat_->bpntr[index], Amat_->bpntr[index+1]-1);

    if (indb < 0) {
      FEI_CERR << "AztecDVBR_Matrix::sumIntoBlockRow: blk col "
	   << blkColInds[blk] << " not found in row " << blkRow << FEI_ENDL;
      abort();
    }

    int numEntries = Amat_->indx[indb+1] - Amat_->indx[indb];
    int valOffset = Amat_->indx[indb];

    //ok, now we're ready to store the stuff.
    for(int i=0; i<numEntries; i++) {
      Amat_->val[valOffset + i] += val[offset + i];
    }

    offset += numEntries;
  }

  return(0);
}
 
//==============================================================================
void AztecDVBR_Matrix::allocate(int* numNzBlks, int* blkColInds) {
//
// This function builds the structure of the matrix. i.e., does the
// memory allocation, etc.
//

   //calculate the bpntr array, which holds info about the number of
   //nonzero blocks per row.
   calcBpntr(numNzBlks);

   //we can now get the total number of nonzero blocks from the last
   //entry in bpntr.
   int totalNumNzBlks = Amat_->bpntr[N_update_];

   //now we can set the bindx array, which holds block column indices.
   setBindx(totalNumNzBlks, blkColInds);

   //and now we're ready to allocate and fill the indx array, which
   //holds info on the number of point entries in each nonzero block.
   calcIndx(totalNumNzBlks);

   //the last thing we need to do is allocate and initialize the val array.
   Amat_->val = new double[localNNZ_];

   for(int i=0; i<localNNZ_; i++) {
      Amat_->val[i] = 0.0;
   }

   setAllocated(true);
   return;
}

//==============================================================================
void AztecDVBR_Matrix::loadComplete() {
//
// This is where we call the Aztec function AZ_transform, which calculates
// communication parameters and re-orders the equations for use as a
// global distributed matrix.
//
   MPI_Comm thisComm = amap_.getCommunicator();

// Sync processors.

   MPI_Barrier(thisComm);

#ifndef FEI_SER
   int thisProc = 0;
   MPI_Comm_rank(thisComm, &thisProc);
#endif

   AZ_transform(amap_.getProcConfig(), &external_, Amat_->bindx, Amat_->val,
                update_, &update_index_, &extern_index_, &data_org_,
                N_update_, Amat_->indx, Amat_->bpntr, Amat_->rpntr,
                &(Amat_->cpntr), AZ_VBR_MATRIX);

   data_org_[AZ_internal_use] = 1;

   Amat_->data_org = data_org_;

//On return from AZ_transform, the array update_index_ contains a mapping
//to the local re-ordering of the indices of the update_ array. Now we will
//fill the orderingUpdate array with the reverse of that mapping. i.e., a
//record of how to get back to the original ordering of the update indices.

   orderingUpdate_ = new int[N_update_];
   for(int ii=0; ii<N_update_; ii++)
      orderingUpdate_[update_index_[ii]] = ii;

// Sync processors.
#ifndef FEI_SER
   MPI_Barrier(thisComm);
#endif

   setLoaded(true);

   return;
}

//==============================================================================
bool AztecDVBR_Matrix::readFromFile(const char *filename){
//
//readFromFile should be able to be called after the matrix is constructed,
//and before allocate has been called. i.e., example usage should include:
//
//  AztecDVBR_Matrix A(map, update);
//  A.readFromFile(fileName);
//  A.matvec(b, c);
//
//i.e., readFromFile can take the place of the allocate and loadComplete
//calls.
//
   FILE *infile = NULL;
   MPI_Comm thisComm = amap_.getCommunicator();

   MPI_Barrier(thisComm);

   infile = fopen(filename, "r");
   if (!infile) messageAbort("readFromFile: couldn't open file.");

   int* num_nz_blocks = NULL;
   int* blk_col_inds = NULL;

   readAllocateInfo(infile, num_nz_blocks, blk_col_inds);

   allocate(num_nz_blocks, blk_col_inds);

   delete [] num_nz_blocks;
   delete [] blk_col_inds;

   fclose(infile);
   infile = fopen(filename, "r");
    
   readMatrixData(infile);

   fclose(infile);

   loadComplete();

   return(true);
}

//==============================================================================
void AztecDVBR_Matrix::readAllocateInfo(FILE* infile,
                                     int*& num_nz_blocks,
                                     int*& blk_col_inds) {
//
//This function will read through infile and construct the lists
//num_nz_blocks (which is the number of nonzero blocks per row) and
//blk_col_inds (which is the block-column indices of those blocks).
//
//It is assumed that these two lists are empty when this function is
//called.

   int i;

   if (num_nz_blocks) delete [] num_nz_blocks;
   if (blk_col_inds) delete [] blk_col_inds;

   num_nz_blocks = new int[N_update_];

   //we'll use a 2-D array for constructing the set of block column indices,
   //because we need to keep them grouped by rows, and we aren't guaranteed
   //that they'll be grouped by rows in the file.
   int totalNumBlks = 0;
   int** blkColInds = new int*[N_update_];

   for(i=0; i<N_update_; i++) {
      num_nz_blocks[i] = 0;
      blkColInds[i] = NULL;
   }

   int blkRows, blkCols, rows, cols;
   char line[256];

   do {
      fgets(line,256,infile);
   } while(strchr(line,'%'));
   sscanf(line,"%d %d %d %d",&blkRows, &blkCols, &rows, &cols);

   if ((blkRows != blkCols) || (rows != cols))
      messageAbort("readAllocateInfo: non-square matrix not allowed.");

   int br, bc, pr, pc, index;

   while (!feof(infile)) {
      do {
         fgets(line,256,infile);
      } while(strchr(line,'%'));

      if(feof(infile))break;

      sscanf(line, "%d %d %d %d", &br, &bc, &pr, &pc);

      if (inUpdate(br, index)) {
         if ((bc < 0) || bc >= blkCols) {
            char mesg[80];
            sprintf(mesg,"readAllocateInfo: blkCols %d, 0-based col ind %d",
                    blkCols, bc);
            fclose(infile);
            messageAbort(mesg);
         }
         insertList(bc, blkColInds[index], num_nz_blocks[index]);
         totalNumBlks++;
      }
   }

   //so we've read the whole file, now flatten the 2-D list blkColInds
   //into the required 1-D list blk_col_inds.
   blk_col_inds = new int[totalNumBlks];

   int offset = 0;
   for(i=0; i<N_update_; i++) {
      for(int j=0; j<num_nz_blocks[i]; j++) {
         blk_col_inds[offset++] = blkColInds[i][j];
      }

      delete [] blkColInds[i];
   }

   delete [] blkColInds;
}

//==============================================================================
void AztecDVBR_Matrix::readMatrixData(FILE* infile) {

   int blkRows, blkCols, rows, cols;
   int br, bc, pr, pc, nnz, index;
   double* blockValues = NULL;
   char line[256];

   do {
      fgets(line,256,infile);
   } while(strchr(line,'%'));
   sscanf(line,"%d %d %d %d",&blkRows, &blkCols, &rows, &cols);

   while (!feof(infile))  {
      do {
         fgets(line,256,infile);
      } while(strchr(line,'%'));

      if(feof(infile))break;

      sscanf(line, "%d %d %d %d %d", &br, &bc, &pr, &pc, &nnz);

      if (inUpdate(br, index)){
         //br (block-row) is in the local row-space, so let's
         //plug this block into our matrix data structure.

         blockValues = new double[nnz];
         getValuesFromString(line, std::strlen(line)+1, blockValues, nnz);

         putBlockRow(br, blockValues, &bc, 1);

         delete [] blockValues;
      }
   }
}

/**=========================================================================**/
void AztecDVBR_Matrix::getValuesFromString(char *line, int len, double *values,
                                           int lenValues){
    (void)len;
    int i;

//first, we know that 'line' contains 5 integers at the beginning, separated
//by spaces, which we want to jump over. So we need the offset of the 5th
//space in 'line'.

    char *offset = &line[0];
    for(i=0; i<5; i++){
        offset = strchr(offset, ' ');
        offset++;
    }

//now we're ready to pick out the numbers to put into the values array.
    for(i=0; i<lenValues; i++){
        sscanf(offset,"%le",&(values[i]));
        offset = strchr(offset, ' ');
        offset++;
    }

    return;
}

/**=========================================================================**/
bool AztecDVBR_Matrix::writeToFile(const char *fileName) const {

   int thisProc = amap_.getProcConfig()[AZ_node];
   int numProcs = amap_.getProcConfig()[AZ_N_procs];
   MPI_Comm thisComm = amap_.getCommunicator();

   int numGlobalBlocks = amap_.getNumGlobalBlocks();
   int numGlobalPtEqns = amap_.globalSize();

   int numLocalBlocks = N_update_;

   FILE* file = NULL;

   for(int p=0; p<numProcs; p++) {
      MPI_Barrier(thisComm);

      if (p == thisProc) {
         if (thisProc == 0) {
            //open the file for writing and write the first line, which is
            //num-blk-rows num-block-cols num-pt-rows num-pt-cols

            file = fopen(fileName, "w");
            fprintf(file, "%d %d %d %d\n", numGlobalBlocks, numGlobalBlocks,
                                           numGlobalPtEqns, numGlobalPtEqns);
         }
         else {
            //open the file for appending
            file = fopen(fileName, "a");
         }

         //now loop over the local portion of the matrix, writing it to file,
         //one nonzer-block per line. Each line of the file will contain these
         //numbers separated by spaces:
         //blk-row-index
         //blk-col-index
         //blk-row-size (number of pt-rows in this block-row)
         //nnz          (number of nonzero point-entries in this block-entry)
         //nonzero1 nonzero2 ... nonzero<nnz> (the nonzeros for this block)

         for(int brow=0; brow<numLocalBlocks; brow++) {
            int bcolind1 = Amat_->bpntr[brow];
            int bcolind2 = Amat_->bpntr[brow+1];

            for(int ind=bcolind1; ind<bcolind2; ind++) {
               int nnzPts = Amat_->indx[ind+1] - Amat_->indx[ind];

               if (nnzPts <= 0) continue;

               int blkRowSize = Amat_->rpntr[brow+1]-Amat_->rpntr[brow];

               int globCol = -1;
               int lookup = Amat_->bindx[ind];
               if (isLoaded()) {
                  if (lookup < N_update_) {
                     globCol = update_[orderingUpdate_[lookup]];
                  }
                  else {
                     globCol = external_[lookup-N_update_];
                  }
               }
               else {
                  globCol = lookup;
               }

               int globalRow = update_[brow];
               if (isLoaded()) globalRow = update_[orderingUpdate_[brow]];

               fprintf(file, "%d %d %d %d ", globalRow, globCol,
                                          blkRowSize, nnzPts);

               int offset = Amat_->indx[ind];
               for(int i=0; i<nnzPts; i++) {
                  fprintf(file, "%20.13e ", Amat_->val[offset+i]);
               }
               fprintf(file, "\n");
            }
         }

         fclose(file);
      }
   }

   return(true);
}

//==============================================================================
int AztecDVBR_Matrix::inUpdate(int globalIndex, int& localIndex) const {
//
// This function determines whether globalIndex is in the local update set,
// and if it is, returns in localIndex the local index for it. If update_index_
// has already been allocated and set (by AZ_transform) then localIndex is
// taken from there.
// If globalIndex is not in the update set, inUpdate returns 0.
//
    localIndex = AZ_find_index(globalIndex, update_, N_update_);

    if(localIndex==-1)return(0);

    if(isLoaded_){
        localIndex = update_index_[localIndex];
    }

    return(1);
}

//==============================================================================
void AztecDVBR_Matrix::calcRpntr() {
//
//This function will use information from the Aztec_BlockMap 'amap_'
//to set the Amat_->rpntr array.
//
//rpntr[0..M] (where M = number-of-blocks)
//rpntr[0] = 0
//rpntr[k+1] - rpntr[k] = size of block k
//
   const int* blkSizes = amap_.getBlockSizes();

   Amat_->rpntr = new int[N_update_+1];

   Amat_->rpntr[0] = 0;

   for(int i=0; i<N_update_; i++) {
      Amat_->rpntr[i+1] = Amat_->rpntr[i] + blkSizes[i];
      if (blkSizes[i] < 0)
          messageAbort("allocate: negative block size.");
   }
}

//==============================================================================
void AztecDVBR_Matrix::calcBpntr(int* numNzBlks) {
//
//This function will set the Amat_->bpntr array.
//
//bpntr[0..M] (where M = number-of-blocks)
//bpntr[0] = 0
//bpntr[k+1]-bpntr[k] = number of nonzero blocks in row k
//
   Amat_->bpntr = new int[N_update_+1];

   Amat_->bpntr[0] = 0;

   for(int i=0; i<N_update_; i++) {
      Amat_->bpntr[i+1] = Amat_->bpntr[i] + numNzBlks[i];
   }
}

//==============================================================================
void AztecDVBR_Matrix::setBindx(int nnzBlks, int* blkColInds) {
//
//This function simply allocates and fills the Amat_->bindx array.
//
   Amat_->bindx = new int[nnzBlks];

   for(int i=0; i<nnzBlks; i++) {
      Amat_->bindx[i] = blkColInds[i];
      if (blkColInds[i] < 0)
         messageAbort("setBindx: negative block col index.");
   }
}

//==============================================================================
void AztecDVBR_Matrix::messageAbort(const char* mesg) const {
   FEI_CERR << "AztecDVBR_Matrix: ERROR: " << mesg << " Aborting." << FEI_ENDL;
   abort();
}

//==============================================================================
void AztecDVBR_Matrix::calcIndx(int nnzBlks) {
//
//This function allocates and fills the Amat_->indx array, which holds info
//on the number of entries in each nonzero block.
//
//indx[0..bpntr[M]], (where M = number of local block rows)
//indx[0] = 0
//indx[k+1]-indx[k] = number of entries in nonzero block k
//

   Amat_->indx = new int[nnzBlks+1];

   //we need to obtain block sizes for all local nonzero blocks. rpntr
   //gives us the sizes for the blocks with column indices in the local
   //update set, but we'll have to do some message passing to obtain the
   //sizes of blocks with column indices in other procs' update sets.

   int numProcs = amap_.getProcConfig()[AZ_N_procs];

   if (numProcs > 1) {
      //form a list of the column indices that are not local.
      calcRemoteInds(remoteInds_, numRemoteBlocks_);

      //now get sizes of blocks that correspond to remote rows.
      remoteBlockSizes_ = new int[numRemoteBlocks_];
      getRemoteBlkSizes(remoteBlockSizes_, remoteInds_, numRemoteBlocks_);
   }

   //now we're ready to set the block sizes in Amat_->indx.
   int index;

   Amat_->indx[0] = 0;

   for(int i=0; i<amap_.getNumLocalBlocks(); i++) {
      int rowBlkSize = Amat_->rpntr[i+1] - Amat_->rpntr[i];

      int colStart = Amat_->bpntr[i];
      int colEnd = Amat_->bpntr[i+1] - 1;

      for(int j=colStart; j<=colEnd; j++) {
         if (inUpdate(Amat_->bindx[j], index)) {
            int colBlkSize = Amat_->rpntr[index+1] - Amat_->rpntr[index];

            Amat_->indx[j+1] = Amat_->indx[j] + rowBlkSize*colBlkSize;
         }
         else { //it's a remoteIndex
            if (numProcs == 1) {
               char mesg[80];
               sprintf(mesg,"calcIndx: blk col index %d not in update set.",
                       Amat_->bindx[j]);
               messageAbort(mesg);
            }

            index = AZ_find_index(Amat_->bindx[j], remoteInds_,
                                  numRemoteBlocks_);
            if (index >= 0) {
               Amat_->indx[j+1] = Amat_->indx[j] +  
                                rowBlkSize*remoteBlockSizes_[index];
            }
            else { //if it wasn't in update or remoteInds, then panic!
               messageAbort("calcIndx: block column index not found.");
            }
         }
      } // end for j loop

      nnzPerRow_[i] = Amat_->indx[colEnd+1] - Amat_->indx[colStart];
   } // end for i loop

   localNNZ_ = Amat_->indx[nnzBlks];
}

//==============================================================================
int AztecDVBR_Matrix::getBindxOffset(int blkInd,
                                     int bpntrStart, int bpntrEnd) const {
//
//This function returns the index of blkInd in the bindx array,
//searching positions bindx[bpntrStart..bpntrEnd].
//
//If blkInd is not found, -1 is returned.
//
   for(int i=bpntrStart; i<=bpntrEnd; i++) {
      int ind = Amat_->bindx[i];
      int globalCol = -1;
      if (isLoaded()) {
         if (ind < N_update_) {
            globalCol = update_[orderingUpdate_[ind]];
         }
         else {
            globalCol = external_[ind-N_update_];
         }
      }
      else globalCol = ind;

      if (globalCol == blkInd) return(i);
   }

   return(-1);
}

//==============================================================================
void AztecDVBR_Matrix::calcRemoteInds(int*& remoteInds, int& len) {
//
//Form a list of the block column indices that are not in the local
//update set.
//
   int nnzBlks = Amat_->bpntr[amap_.getNumLocalBlocks()];
   int local;

   for(int i=0; i<nnzBlks; i++) {
      if (!inUpdate(Amat_->bindx[i], local)) {
         insertList(Amat_->bindx[i], remoteInds, len);
      }
   }
}

//==============================================================================
void AztecDVBR_Matrix::getRemoteBlkSizes(int* remoteBlkSizes,
					 int* remoteInds,
                                         int len)
{
  //
  //remoteInds is a sorted list of indices that correspond to rows
  //in remote processors' update lists. This function will spread the
  //indices to all processors so that they can provide the blk sizes,
  //then spread that information back to all processors.
  //
#ifdef FEI_SER
  return;
#else
  int numProcs = amap_.getProcConfig()[AZ_N_procs];
  int thisProc = amap_.getProcConfig()[AZ_node];
  MPI_Comm comm = amap_.getCommunicator();

   int* lengths = new int[numProcs];
   lengths[0] = 0;

   //gather up the lengths of the lists that each proc will be sending.
   MPI_Allgather(&len, 1, MPI_INT, lengths, 1, MPI_INT, comm);

   //now form a list of the offset at which each proc's contribution will
   //be placed in the all-gathered list.
   int* offsets = new int[numProcs];

   offsets[0] = 0;
   int totalLength = lengths[0];
   for(int i=1; i<numProcs; i++) {
      offsets[i] = offsets[i-1] + lengths[i-1];
      totalLength += lengths[i];
   }

   //now we can allocate the list to recv into.
   int* recvBuf = new int[totalLength];

   //now we're ready to do the gather.
   MPI_Allgatherv(remoteInds, len, MPI_INT, recvBuf, lengths, offsets,
                  MPI_INT, comm);

   //now we'll run through the list and put block sizes into a list of
   //the same length as the total recvBuf list.
   int* blkSizes = new int[totalLength];
   int index;

   for(int j=0; j<totalLength; j++) {
      if (inUpdate(recvBuf[j], index)) {
         blkSizes[j] = Amat_->rpntr[index+1]-Amat_->rpntr[index];
      }
      else blkSizes[j] = 0;
   }

   //now we'll reduce this info back onto all processors. We'll use MPI_SUM.
   //Since the sizes we did NOT supply hold a 0, and each spot in the list
   //should only have a nonzero size from 1 processor, the result will be
   //that each spot in the result list has the correct value.
   int* recvSizes = new int[totalLength];

   MPI_Allreduce(blkSizes, recvSizes, totalLength, MPI_INT, MPI_SUM, comm);

   //and finally, we just need to run our section of the list of recv'd sizes,
   //and transfer them into the remoteBlkSizes list.
   int offset = offsets[thisProc];
   for(int k=0; k<len; k++) {
      remoteBlkSizes[k] = recvSizes[offset + k];
      if (recvSizes[offset+k] <= 0)
         messageAbort("getRemoteBlkSizes: recvd a size <= 0.");
   }

   delete [] lengths;
   delete [] offsets;
   delete [] recvBuf;
   delete [] blkSizes;
   delete [] recvSizes;
#endif
}

//==============================================================================
void AztecDVBR_Matrix::insertList(int item, int*& list, int& len) {
//
//insert 'item' in 'list', if it's not already in there,
//and update the list's length, 'len'.
//
//We want to keep the list ordered, so we'll insert item in
//the list after the biggest existing entry that's smaller, and
//before the smallest existing entry that's bigger.
//

   if (len <= 0) {
      list = new int[1];
      list[0] = item;
      len = 1;
      return;
   }

   int index = AZ_find_index(item, list, len);

   if (index >= 0) return;

   int* newList = new int[len+1];

   //bring over the contents of the old list, putting in the new
   //one at the appropriate point.
   int inserted = 0;
   for(int i=0; i<len; i++) {
      if (!inserted) {
         if (list[i] < item) newList[i] = list[i];
         else {
            newList[i] = item;
            inserted = 1;
         }
      }
      else newList[i] = list[i-1];
   }

   //now put in the last list entry
   if (inserted) newList[len] = list[len-1];
   else newList[len] = item;

   //delete the old memory and reset the pointer.
   if (len > 0) delete [] list;
   list = newList;

   //update the length.
   len++;
}

#endif
//HAVE_FEI_AZTECOO

