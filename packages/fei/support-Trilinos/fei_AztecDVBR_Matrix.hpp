#ifndef _fei_AztecDVBR_Matrix_hpp_
#define _fei_AztecDVBR_Matrix_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This class provides a wrapper for the Aztec DVBR matrix data structure.
//
// Some functions return an int. This will be the ESI_error_code. It
// will be 0 if there were no errors, 1 if an error occurred.
//
namespace fei_trilinos {

class Aztec_BlockMap;

class AztecDVBR_Matrix {
 
  public:
    // Constructor.
    AztecDVBR_Matrix(Aztec_BlockMap& map, int* update);

    //Copy constructor.
    AztecDVBR_Matrix(const AztecDVBR_Matrix& src);

    virtual ~AztecDVBR_Matrix ();

    //query functions.

    int getBlockMaps(Aztec_BlockMap** rowMap, Aztec_BlockMap** colMap);
    int getNumBlocksPerRow(int blkRow, int& nnzBlksPerRow) const;
    int getNumNonzerosPerRow(int blkRow, int& nnzPerRow) const;
    int getNumBlocksPerRow(int* nnzBlksPerRow) const;
    int getNumNonzerosPerRow(int* nnzPerRow) const;
    int getBlockSize(int blkRow, int blkCol, int& ptRows, int& ptCols);

    // Mathematical functions.
    void matvec(const Aztec_Vector& x, Aztec_Vector& y) const;

    void put(double s);

    // ... to read matrix.

    int getBlockRow(int blk_row,
                    double* vals,
                    int* blk_col_inds,
                    int num_nz_blocks) const;
    
    // ... to write matrix.

    int putBlockRow(int blk_row,
                    double* vals,
                    int* blk_col_inds,
                    int num_nz_blocks) const;

    int sumIntoBlockRow(int blk_row,
                        double* vals,
                        int* blk_col_inds,
                        int num_nz_blocks) const;
 
    // configuration function.
    void allocate(int* num_nz_blocks, int* blk_col_inds);

    void loadComplete();

    bool isLoaded() const {return(isLoaded_);};
    void setLoaded(bool flag) {isLoaded_ = flag;};
    bool isAllocated() const {return(isAllocated_);};
    void setAllocated(bool flag) {isAllocated_ = flag;};

    AZ_MATRIX* getAZ_MATRIX_Ptr() const {return(Amat_);};

    bool readFromFile(const char *filename);
    bool writeToFile(const char *fileName) const;

    //numRemoteBlocks is the number of block-column-indices on this processor
    //that correspond to block-rows that reside on another processor.
    int getNumRemoteBlocks() {return(numRemoteBlocks_);};
    int* getRemoteBlockIndices() {return(remoteInds_);};
    int* getRemoteBlockSizes() {return(remoteBlockSizes_);};

    //Aztec-specific functions:
    int* getUpdate() {return(update_);};
    int* getUpdate_index() {return(update_index_);};
    int* getData_org() {return(data_org_);};

  private:
    int inUpdate(int globalIndex, int& localIndex) const;

    void readAllocateInfo(FILE* infile, int*& num_nz_blocks, int*& blk_col_inds);
    void readMatrixData(FILE* infile);

    void calcRpntr();
    void calcBpntr(int* nzBlksPerRow);
    void setBindx(int nnzBlks, int* blkColInds);
    void calcIndx(int nnzBlks);

    int getBindxOffset(int blkInd, int bpntrStart, int bpntrEnd) const;

    void calcRemoteInds(int*& remoteInds, int& len);
    void getRemoteBlkSizes(int* remoteBlkSizes, int* remoteInds, int len);
    void insertList(int item, int*& list, int& len);
    void getValuesFromString(char *line, int len, double *values,
                             int lenValues);
    void messageAbort(const char* mesg) const;

    Aztec_BlockMap& amap_;

    AZ_MATRIX *Amat_;

    int N_update_;
    int* update_;
    int* external_;
    int* extern_index_;
    int* update_index_;
    int* data_org_;
    int* orderingUpdate_;

    bool isLoaded_;
    bool isAllocated_;

    int localNNZ_;
    int* nnzPerRow_;

    int numRemoteBlocks_;
    int* remoteInds_;
    int* remoteBlockSizes_;
};

}//namespace fei_trilinos

#endif
