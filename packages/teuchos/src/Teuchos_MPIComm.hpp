// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_MPICOMM_H
#define TEUCHOS_MPICOMM_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_RefCountPtr.hpp"



namespace Teuchos
{
  /**
   * Object representation of an MPI communicator.
   *
   * At present, groups are not implemented so the only communicator
   * is MPI_COMM_WORLD.
   */
  class MPIComm
    {
    public:
      /** empty ctor builds an object for MPI_COMM_WORLD */
      MPIComm();

#ifdef HAVE_MPI
      /** construct a MPIComm for a given MPI communicator */
      MPIComm(MPI_Comm comm);
#endif

      /** get an object representing MPI_COMM_WORLD */
      static MPIComm& world();


      /** return process rank */
      int getRank() const {return myRank_;}

      /** return number of processors in the communicator */
      int getNProc() const {return nProc_;}

      /** synchronize all the processors in the communicator */
      void synchronize() const ;

      /** {\bf Collective communications} */
      //@{
      /** all-to-all gather-scatter */
      void allToAll(void* sendBuf, int sendCount, int sendType,
                    void* recvBuf, int recvCount, int recvType) const ;
      /** variable-length gather-scatter */
      void allToAllv(void* sendBuf, int* sendCount, int* sendDisplacements,
                     int sendType,
                     void* recvBuf, int* recvCount,
                     int* recvDisplacements,
                     int recvType) const ;

      /** do a collective operation, scattering the results to all procs */
      void allReduce(void* input, void* result, int inputCount, int type,
                     int op) const ;

      /** gather to root */
      void gather(void* sendBuf, int sendCount, int sendType,
                  void* recvBuf, int recvCount, int recvType,
                  int root) const ;

      /** gather to all procs */
      void allGather(void* sendBuf, int sendCount, int sendType,
                     void* recvBuf, int recvCount, int recvType) const ;

      /** gather to all procs */
      void allGatherv(void* sendBuf, int sendCount, int sendType,
                      void* recvBuf, int* recvCount, int* recvDisplacements,
                      int recvType) const ;

      /** broadcast */
      void bcast(void* msg, int length, int type, int src) const ;
      //@}

#ifdef HAVE_MPI
      /** get the MPI_Comm communicator handle */
      MPI_Comm getComm() const {return comm_;}
#endif

      /** \name Constants */
      //@{
      /** \name Data types  */
      //@{
      /** */
      const static int INT;
      /** */
      const static int FLOAT;
      /** */
      const static int DOUBLE;
      /** */
      const static int CHAR;
      //@}

      /** \name Operations */
      //@{
      /** */
      const static int SUM;
      /** */
      const static int MIN;
      /** */
      const static int MAX;
      /** */
      const static int PROD;
      //@}

      // errCheck() checks the return value of an MPI call and throws
      // a ParallelException upon failure.
      static void errCheck(int errCode, const string& methodName);

#ifdef HAVE_MPI
      // getDataType converts a PMachine data type code to a MPI_Datatype
      static MPI_Datatype getDataType(int type);

      // getOp converts a PMachine operator code to a MPI_Op operator code.
      static MPI_Op getOp(int op);
#endif
    private:
#ifdef HAVE_MPI
      MPI_Comm comm_;
#endif

      int nProc_;
      int myRank_;

      /** common initialization function, called by all ctors */
      void init();
    };
}
#endif

