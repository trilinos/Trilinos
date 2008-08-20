/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_CommUtilsBase_hpp_
#define _fei_CommUtilsBase_hpp_

#include "fei_macros.hpp"
#include "fei_mpi.h"
#include "fei_CommCore.hpp"
#include "fei_ctg_set.hpp"
#include "snl_fei_RaggedTable.hpp"

#include <map>

namespace fei {

  typedef snl_fei::RaggedTable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> > comm_map;

  class CommUtilsBase {
  public:
    CommUtilsBase(MPI_Comm comm);
    virtual ~CommUtilsBase();

    int numProcs() const;
    int localProc() const;
    MPI_Comm getCommunicator() const;
    int Barrier() const;

    /** Scenario: The local processor has a list of processors to which data
        will be sent, but doesn't know which processors data will be received
        from. This method produces that list of processors to be received from.
        This is a collective method.
    */
    int mirrorProcs(std::vector<int>& toProcs,
                    std::vector<int>& fromProcs);

    /** Given a list of processors to send to, a scalar to send to each,
        and a list of processors to receive from, perform the exchange.

        @param sendProcs Input. List of processors to send to.
        @param sendData Input. List of data, same length as 'sendProcs', to be
        sent. (One item to be sent to each send proc.)
        @param recvProcs Input. List of processors to receive from.
        @param recvData Output. On exit, contains one item received from each
        recv proc.
        @return error-code 0 if successful
    */
    int exchangeIntData(std::vector<int>& sendProcs,
                        std::vector<int>& sendData,
                        std::vector<int>& recvProcs,
                        std::vector<int>& recvData);

    /** Given a comm-pattern, create and initialize its mirror...
     */
    int mirrorCommPattern(comm_map* inPattern,
                          comm_map*& outPattern);

    CommCore* commCore_;

  protected:
    MPI_Comm comm_;
    int numProcs_, localProc_;

  private:
    CommUtilsBase(const CommUtilsBase& src);
    CommUtilsBase& operator=(const CommUtilsBase& src);
  };//class CommUtilsBase

} //namespace fei

#endif // _fei_CommUtilsBase_hpp_

