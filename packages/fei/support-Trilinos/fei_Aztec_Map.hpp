#ifndef _fei_Aztec_Map_hpp_
#define _fei_Aztec_Map_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <vector>

//
// Aztec_Map is a wrapper that encapsulates the general
// information needed to describe the layout of an Aztec matrix or
// vector structure. It is a companion/support class that goes with
// the three data class wrappers Aztec_LSVector, AztecDMSR_Matrix and
// AztecDVBR_Matrix (the Aztec_BlockMap specialization is also
// required for DVBR).
//
// Aztec_Map allows the storage and retrieval of information such as
// local and global sizes, the MPI communicator, and the proc_config array.
//
namespace fei_trilinos {

class Aztec_Map {
    
  public:
    Aztec_Map(int globalSize, int N_update, const int* update, int localOffset,
              MPI_Comm comm);

    Aztec_Map(const Aztec_Map& map);            // copy constructor    
    virtual ~Aztec_Map(void);

    virtual const int& localSize() const {return(localSize_);}
    virtual const int& globalSize() const {return(globalSize_);}
    virtual const int& localOffset() const {return(localOffset_);}

    int* getUpdate()
    {
      return update.size()>0 ? &update[0] : NULL;
    }

    virtual MPI_Comm getCommunicator() const {return(comm_);}

    virtual int* getProcConfig()
    {
      return proc_config.size()>0 ? &proc_config[0] : NULL;
    }

    std::vector<int> proc_config;
    std::vector<int> update;
    int* external;
    int* update_index;
    int* extern_index;
    int* data_org;
    std::vector<int> orderingUpdate;

    bool az_transformed;

    int getTransformedEqn(int eqn) const {
      if (az_transformed == true) {
        return eqn<N_update_ ? update[orderingUpdate[eqn]] : external[eqn-N_update_];
      }
      return eqn;
    }

    bool inUpdate(int globalIndex, int& localIndex) const
    {
      localIndex = globalIndex - localOffset_;
      if (localIndex<0 || localIndex>=localSize_) {
        localIndex = -1;
        return false;
      }
      if (az_transformed == true) {
        localIndex = update_index[localIndex];
      }
      return true;
    }

  private:
    void checkInput();

    int globalSize_;
    int localSize_;
    int localOffset_;
    int N_update_;

    MPI_Comm comm_;
};

}//namespace fei_trilinos

#endif

