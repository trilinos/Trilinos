// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGMATRIX_HPP_
#define _ZOLTAN2_ALGMATRIX_HPP_

// #include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_MatrixPartitioningSolution.hpp>
// #include <Zoltan2_Algorithm.hpp>

#include <Zoltan2_PartitioningProblem.hpp>

// #include <sstream>
// #include <string>
// #include <bitset>

/*! \file Zoltan2_AlgMatrix.hpp
 *  \brief A set of algorithms for Matrix partitioning.  For now... 2D random Cartesian
 */

namespace Zoltan2{


/*! Block partitioning method.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param ids    an Identifier model
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  No special requirements on the
 *    identifiers.
 *
 */


template <typename Adapter>
class AlgMatrix : public Algorithm<Adapter>
{

public:
  typedef typename Adapter::lno_t lno_t;     // local ids
  typedef typename Adapter::gno_t gno_t;     // global ids
  typedef typename Adapter::scalar_t scalar_t;   // scalars
  typedef typename Adapter::part_t part_t;   // part numbers

  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

  typedef typename Adapter::base_adapter_t base_adapter_t;



  // Constructor
  AlgMatrix(const RCP<const Environment> &_env,
  	    const RCP<const Comm<int> > &_problemComm,
  	    const RCP<const MatrixAdapter<user_t, userCoord_t> > &_matrixAdapter)
    : mEnv(_env), mProblemComm(_problemComm), mMatrixAdapter(_matrixAdapter)
  {}


  // AlgMatrix(const RCP<const Environment> &_env,
  // 	    const RCP<const Comm<int> > &_problemComm,
  // 	    const RCP<const Adapter> &_matrixAdapter)
  //   : mEnv(_env), mProblemComm(_problemComm), mMatrixAdapter(_matrixAdapter)
  // {}

  // Partitioning method
  void partitionMatrix(const RCP<MatrixPartitioningSolution<Adapter> > &solution);


private:
  const RCP<const Environment> mEnv;
  const RCP<const Comm<int> > mProblemComm;

  const RCP<const MatrixAdapter<user_t, userCoord_t> > mMatrixAdapter;
  //const RCP<const Adapter > mMatrixAdapter;

  //typedef Tpetra::CrsMatrix<z2TestScalar, z2TestLO, z2TestGO> SparseMatrix;
  //  const RCP<const XpetraCrsMatrixAdapter<user_t, userCoord_t> > mMatrixAdapter;




};


////////////////////////////////////////////////////////////////////////////////
// Partitioning method
//
//   -- For now implementing 2D Random Cartesian
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgMatrix<Adapter>::partitionMatrix(const RCP<MatrixPartitioningSolution<Adapter> > &solution)
{
  mEnv->debug(DETAILED_STATUS, std::string("Entering AlgBlock"));

  int myrank = mEnv->myRank_;
  int nprocs = mEnv->numProcs_;

  ///////////////////////////////////////////////////////////////////////////
  // Determine processor dimension d for 2D block partitioning d = sqrt(p) 
  // Determine processor row and processor column for this rank for 2D partitioning
  ///////////////////////////////////////////////////////////////////////////
  int procDim = sqrt(nprocs);
  assert(procDim * procDim == nprocs); // TODO: Should check this earlier and produce more useful error
  
  int myProcRow = myrank / procDim;
  int myProcCol = myrank % procDim;
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Create 1D Random partitioning problem to create partitioning for the 2D partitioning
  ///////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Create parameters for 1D partitioning
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params1D("Params for 1D partitioning");
  //  params1D.set("algorithm", "random"); //TODO add support for random
  params1D.set("algorithm", "block");
  params1D.set("imbalance_tolerance", 1.1);
  params1D.set("num_global_parts", procDim);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Create Zoltan2 partitioning problem 
  //    -- Alternatively we could simply call algorithm directly
  //////////////////////////////////////////////////////////////////////
  MatrixAdapter<user_t, userCoord_t> *adapterPtr = 
    const_cast<MatrixAdapter<user_t, userCoord_t>*>(mMatrixAdapter.getRawPtr());
  

  Zoltan2::PartitioningProblem<base_adapter_t> *problem = 
    new Zoltan2::PartitioningProblem<base_adapter_t>(adapterPtr, &params1D);

  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Solve the problem, get the solution
  //////////////////////////////////////////////////////////////////////
  problem->solve();

  //  Zoltan2::PartitioningSolution<SparseMatrixAdapter_t> solution = problem.getSolution();
  const Zoltan2::PartitioningSolution<base_adapter_t> &solution1D = problem->getSolution();

  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Get part assignments for each local_id
  //////////////////////////////////////////////////////////////////////
  //  size_t numGlobalParts = solution1D->getTargetGlobalNumberOfParts();
  
  const part_t *parts = solution1D.getPartListView();

  //////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Create column Ids ArrayRCP colIDs to store in solution
  //    This will define which column IDs are allowed in a particular processor
  //    column.
  ///////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Group gids corresponding to local ids based on process column 
  //////////////////////////////////////////////////////////////////////
  const gno_t *rowIds;
  mMatrixAdapter->getRowIDsView(rowIds);

  size_t nLocIDs = mMatrixAdapter->getLocalNumRows();

  std::vector<std::vector<gno_t> > idsInProcColI(procDim);
  Teuchos::ArrayRCP<gno_t> colIDs;

  for(size_t i=0; i<nLocIDs; i++)
  {
    // Nonzeros with columns of index rowIds[i] belong to some processor
    // in processor column parts[i]
    idsInProcColI[ parts[i] ].push_back(rowIds[i]);
  }
  //////////////////////////////////////////////////////////////////////

  delete problem; // delete 1D partitioning problem


  //////////////////////////////////////////////////////////////////////
  // Communicate gids to process col roots 
  //////////////////////////////////////////////////////////////////////

  // Loop over each of the processor columns
  for(int procColI=0; procColI<procDim; procColI++)
  {

    /////////////////////////////////////////////////////////////////
    // This rank is root of procColI, need to receive
    /////////////////////////////////////////////////////////////////
    if(myProcCol==procColI && myProcRow==0) 
    {
      assert(myrank==procColI); // Could make this above conditional?


      std::set<gno_t> gidSet; // to get sorted list of gids

      ////////////////////////////////////////////////////////////
      // Copy local gids to sorted gid set
      ////////////////////////////////////////////////////////////
      for(int i=0;i<idsInProcColI[procColI].size();i++)
      {
        gidSet.insert(idsInProcColI[procColI][i]);
      }
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // Receive gids from remote processors, insert into set
      ////////////////////////////////////////////////////////////
      std::vector<gno_t> recvBuf;
      for(int src=0;src<nprocs;src++)
      {
        if(src!=myrank)
  	{
          int buffSize;

  	  Teuchos::receive<int,int>(*mProblemComm, src, 1, &buffSize);

          if(buffSize>0)
  	  {
            recvBuf.resize(buffSize);

  	    Teuchos::receive<int,gno_t>(*mProblemComm, src, buffSize, recvBuf.data());

            for(int i=0;i<buffSize;i++)
  	    {
              gidSet.insert(recvBuf[i]);
  	    }
  	  }
  	}
      }
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      //Copy data to std::vector 
      ////////////////////////////////////////////////////////////
      colIDs.resize(gidSet.size());

      typename std::set<gno_t>::const_iterator iter;
      int indx=0;
      for (iter=gidSet.begin(); iter!=gidSet.end(); ++iter)
      {
        colIDs[indx] = *iter;
  	indx++;
      }
      ////////////////////////////////////////////////////////////
    }
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    // This rank is not root of procColI, need to senddata block to root
    /////////////////////////////////////////////////////////////////
    else
    {
      int sizeToSend = idsInProcColI[procColI].size();

      //is the dst proc info correct here?

      // Root of procColI is rank procColI
      Teuchos::send<int,int>(*mProblemComm,1,&sizeToSend,procColI);

      Teuchos::send<int,gno_t>(*mProblemComm,sizeToSend,idsInProcColI[procColI].data(),procColI);
    }
    /////////////////////////////////////////////////////////////////

  }
  //////////////////////////////////////////////////////////////////////

  // Free memory if possible


  //////////////////////////////////////////////////////////////////////
  // Communicate gids from process col root down process column to other procs 
  //////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////
  // This rank is root of its processor column, need to send
  /////////////////////////////////////////////////////////////////
  if(myProcRow==0) 
  {
    ////////////////////////////////////////////////////////////
    // Send data to remote processors in processor column
    ////////////////////////////////////////////////////////////
    for(int procColRank=0; procColRank<procDim; procColRank++)
    {
      // Convert from proc column rank to global rank
      int dstRank = procColRank*procDim + myProcCol;

      if(dstRank!=myrank)
      {
        int sizeToSend = colIDs.size();

  	Teuchos::send<int,int>(*mProblemComm,1,&sizeToSend,dstRank);
        Teuchos::send<int,gno_t>(*mProblemComm,sizeToSend,colIDs.get(),dstRank);

      }
    }
    ////////////////////////////////////////////////////////////

  }
  /////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////
  // This rank is not root of processor, need to recv data from root
  /////////////////////////////////////////////////////////////////
  else
  {
    // Src rank is root of processor column = myProcCol
    int srcRank = myProcCol;

    int buffSize;
    Teuchos::receive<int,int>(*mProblemComm, srcRank, 1, &buffSize);

    colIDs.resize(buffSize);
    Teuchos::receive<int,gno_t>(*mProblemComm, srcRank, buffSize, colIDs.get());
  }
  /////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Create domain/range IDs (same for both now) Array RCP to store in solution
  //   Created from equal division of column IDs
  ///////////////////////////////////////////////////////////////////////////

  // Split processor column ids into nDim parts
  // Processor column i ids split between ranks 0+i, nDim+i, 2*nDim+i, ...

  // 

  ArrayRCP<gno_t> domRangeIDs;     // Domain/Range vector Ids assigned to this process

  size_t nCols = colIDs.size();
  lno_t nColsDivDim = nCols / procDim;
  lno_t nColsModDim = nCols % procDim;

  size_t sColIndx;
  size_t domRangeSize;

  // This proc will have nColsDivDim+1 domain/range ids
  if(myProcRow < nColsModDim)
  {
    sColIndx = myProcRow * (nColsDivDim+1);
    domRangeSize = nColsDivDim+1;
  }
  // This proc will have nColsDivDim domain/range ids
  else
  {
    sColIndx =  nColsModDim*(nColsDivDim+1) + (myProcRow-nColsModDim) * nColsDivDim;
    domRangeSize = nColsDivDim;
  }

  domRangeIDs.resize(domRangeSize);

  for(size_t i=0;i<domRangeSize;i++)
  {
    domRangeIDs[i] = colIDs[sColIndx+i];
  }  
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Create row IDs Array RCP to store in solution
  //   Created from union of domRangeIDs
  //
  // Procs 0, 1, ... nDim-1 share the same set of rowIDs
  // Procs nDim, nDim+1, ..., 2*nDim-1 share the same set of rowIDs
  ///////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Create subcommunicator for processor row
  ////////////////////////////////////////////////////////////////////// 
  std::vector<int> subcommRanks(procDim);

  for(unsigned int i=0; i<procDim; i++)
  {
    subcommRanks[i] = myProcRow*procDim + i;
  }

  ArrayView<const int> subcommRanksView = Teuchos::arrayViewFromVector (subcommRanks);

  RCP<Teuchos::Comm<int> > rowcomm = mProblemComm->createSubcommunicator(subcommRanksView);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Determine max number of columns in this process row
  //////////////////////////////////////////////////////////////////////
  size_t maxProcRowNCols=0;

  Teuchos::reduceAll<int, size_t>(*rowcomm,Teuchos::REDUCE_MAX,1,&domRangeSize,&maxProcRowNCols);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Communicate all domRangeIDs to all processes in procRow
  //////////////////////////////////////////////////////////////////////
  gno_t MAXVAL = std::numeric_limits<gno_t>::max();

  std::vector<gno_t> locRowIDs(maxProcRowNCols,MAXVAL);

  Teuchos::ArrayRCP<gno_t> rowIDs(maxProcRowNCols * procDim);

  // Copy local domRangeIDs into local row IDs, "extra elements" will have
  // value MAXVAL
  for(size_t i=0;i<domRangeIDs.size();i++)
  {
    locRowIDs[i] = domRangeIDs[i];
  }    

  // Gather all ids onto all processes in procRow
  Teuchos::gatherAll<int,gno_t>(*rowcomm,maxProcRowNCols, locRowIDs.data(),
                                maxProcRowNCols*procDim, rowIDs.get());
  //////////////////////////////////////////////////////////////////////

  // Free local data
  std::vector<int>().swap(locRowIDs);


  //////////////////////////////////////////////////////////////////////
  // Insert local domRangeIDs into row IDs set, filter out extra items
  //////////////////////////////////////////////////////////////////////
  std::set<gno_t> setRowIDs;

  for(size_t i=0;i<rowIDs.size();i++)
  {
    if(rowIDs[i]!=MAXVAL)
    {
      setRowIDs.insert(rowIDs[i]);
    }
  }  
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Resize rowIDs array and copy data to array (now sorted)
  //////////////////////////////////////////////////////////////////////
  rowIDs.resize(setRowIDs.size());

  typename std::set<gno_t>::const_iterator iter;
  size_t indx=0;

  for(iter=setRowIDs.begin(); iter!=setRowIDs.end(); ++iter)
  {
    rowIDs[indx] = *iter;
    indx++;
  }  
  //////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Finished. Store partition in solution.
  ///////////////////////////////////////////////////////////////////////////
  solution->setIDLists(rowIDs,colIDs,domRangeIDs,domRangeIDs);
  ///////////////////////////////////////////////////////////////////////////

  mEnv->debug(DETAILED_STATUS, std::string("Exiting AlgMatrix"));
}
////////////////////////////////////////////////////////////////////////////////



}   // namespace Zoltan2

#endif
