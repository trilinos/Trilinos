#ifndef _fei_Aztec_BlockMap_h_
#define _fei_Aztec_BlockMap_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This Aztec_BlockMap class is a wrapper that encapsulates the general
// information needed to describe the layout of an Aztec DVBR matrix.
// It is a companion/support class that goes with the data class wrappers
// Aztec_Vector and AztecDVBR_Matrix. Aztec_BlockMap inherits from
// Aztec_Map.
//
// Aztec_Map allows the storage and retrieval of information such as
// local and global sizes, the MPI communicator, and the proc_config array.
// Aztec_BlockMap allows the storage and retrieval of information that
// describes the partitioning and layout of a block matrix.
//

namespace fei_trilinos {

class Aztec_BlockMap : public Aztec_Map {
    
  public:
    Aztec_BlockMap(int globalSize, int localSize, int localOffset,
                   MPI_Comm comm,
                   int numGlobalBlocks, int numLocalBlocks,
                   int localBlockOffset, int* blockSizes);

    Aztec_BlockMap(const Aztec_BlockMap& map);       // copy constructor
    virtual ~Aztec_BlockMap(void);

    const int& getNumGlobalBlocks() const {return(numGlobalBlocks_);};
    const int& getNumLocalBlocks() const {return(numLocalBlocks_);};
    const int& getLocalBlockOffset() const {return(localBlockOffset_);};

    const int* getBlockSizes() const {return(blockSizes_);};

  private:

    void checkInput();

    int numGlobalBlocks_;
    int numLocalBlocks_;
    int localBlockOffset_;
    int* blockSizes_;
};

}// namespace fei_trilinos

#endif

