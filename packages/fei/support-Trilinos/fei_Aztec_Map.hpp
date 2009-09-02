#ifndef _fei_Aztec_Map_hpp_
#define _fei_Aztec_Map_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This Aztec_Map class is a wrapper that encapsulates the general
// information needed to describe the layout of an Aztec matrix or
// vector structure. It is a companion/support class that goes with
// the three data class wrappers Aztec_Vector, AztecDMSR_Matrix and
// AztecDVBR_Matrix (the Aztec_BlockMap specialization is also
// required for DVBR).
//
// Aztec_Map allows the storage and retrieval of information such as
// local and global sizes, the MPI communicator, and the proc_config array.
//
namespace fei_trilinos {

class Aztec_Map {
    
  public:
    Aztec_Map(int globalSize, int localSize, int localOffset,
              MPI_Comm comm);

    Aztec_Map(const Aztec_Map& map);            // copy constructor    
    virtual ~Aztec_Map(void);

    virtual const int& localSize() const {return(localSize_);};
    virtual const int& globalSize() const {return(globalSize_);};
    virtual const int& localOffset() const {return(localOffset_);};

    virtual MPI_Comm getCommunicator() const {return(comm_);};

    virtual int* getProcConfig() const {return(proc_config_);};

  private:
    void checkInput();

    int globalSize_;
    int localSize_;
    int localOffset_;

    MPI_Comm comm_;

    int* proc_config_;  //Aztec information container
};

}//namespace fei_trilinos

#endif

