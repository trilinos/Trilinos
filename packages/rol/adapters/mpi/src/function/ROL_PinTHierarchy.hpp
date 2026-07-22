// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTHIERARCHY_HPP
#define ROL_PINTHIERARCHY_HPP

#include <unordered_map>

#include "ROL_TimeStamp.hpp"
#include "ROL_PinTVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_PinTCommunicationUtilities.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

namespace ROL {

template<typename Real> 
class PinTHierarchy {

  using V  = ROL::Vector<Real>;
  using PV = ROL::PartitionedVector<Real>;

  using size_type = typename std::vector<Real>::size_type;
  template<typename T> using Ptr = ROL::Ptr<T>;

private:

  // internal state members
  bool isInitialized_;

  Ptr<const std::vector<TimeStamp<Real>>>         userTimeStamps_;    // these are the untouched ones the user passes in
  Ptr<std::vector<TimeStamp<Real>>>               timeStamps_;        // these are used internally
  std::vector<Ptr<std::vector<TimeStamp<Real>>>>  stamps_;            // these are used internally on each level

  // preconditioner settings
  int maxLevels_;               // must turn on multigrid 

  std::vector<Ptr<const PinTCommunicators>> communicators_;
  Ptr<const PinTVectorCommunication<Real>> vectorComm_; // object for sending whole vectors between processors

  bool rebalance_;

  // This class encapsulates the processor redistribution required
  // for both restriction and prolongation of a vector
  struct TimeStepCommMap {
    // This vector is a map from the locally computed vector
    // to its place on another processor. For instance if there are
    // three time steps locally on a processor
    // hen this vector will have 3 processor destinations, one for
    // each coarse vector.
    std::vector<int> sendToProc_;

    // This vector is a map from the coarse vector on the reciving processor
    // to where it is computed on another processor. For instance if there are
    // three time steps in the coarse vector on this processor, 
    // then this vector will have 3 processor destinations, one for
    // each coarse vector.
    std::vector<int> recvFromProc_;
  };

  std::vector<TimeStepCommMap> restrictOptMap_;
  std::vector<TimeStepCommMap> restrictSimMap_;
  std::vector<TimeStepCommMap> prolongOptMap_;
  std::vector<TimeStepCommMap> prolongSimMap_;

  /**
   * \brief Get the time stamps for by level
   *
   * The current implementation is recursive and expensive, but its the easiest
   * way to do this right now. 
   */
  Ptr<std::vector<TimeStamp<Real>>> buildTimeStampsByLevel(int level,
                                                           const ROL::Ptr<const PinTCommunicators> & comm) const
  {
    assert(level>=0); // precondition

    // base case
    if(level==0)
      return timeStamps_;

    int timeRank = comm->getTimeRank();

    Ptr<std::vector<TimeStamp<Real>>> higherLevel = getTimeStampsByLevel(level-1);

    int myFirstIndex = 0;
    int myLength = int(higherLevel->size());
    MPI_Exscan(&myLength,&myFirstIndex,1,MPI_INT,MPI_SUM,comm->getTimeCommunicator());

    // send your end point to the right
    Real coarseTimeEnd = (higherLevel->end()-1)->t[0];
    if(timeRank!=comm->getTimeSize()-1 && ((myFirstIndex+higherLevel->size()) % 2)==0) {
      int tag = 0;
      MPI_Send(&coarseTimeEnd,1,MPI_DOUBLE,comm->getTimeRank()+1,tag,comm->getTimeCommunicator());
    }

    // recieve the left time node on the coarse grid
    Real coarseTimeStart = 0.0;
    if(timeRank!=0 && (myFirstIndex % 2)==0) {
      int tag = 0;
      MPI_Recv(&coarseTimeStart,1,MPI_DOUBLE,comm->getTimeRank()-1,tag,comm->getTimeCommunicator(),MPI_STATUS_IGNORE);
    }
    else {
      coarseTimeStart = higherLevel->begin()->t[0];
    }

    // get the size of the array
    int currentStamps = -1; 
    if(higherLevel->size() % 2 ==0) {
      currentStamps = higherLevel->size()/2;
    }
    else {
      currentStamps = (higherLevel->size()-1)/2+(myFirstIndex+1)%2;
    }

    Ptr<std::vector<TimeStamp<Real>>> currentLevel 
        = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>(currentStamps);

    // you offset based on your starting index (and if you are on rank 0
    int offset = -1;
    if(myFirstIndex % 2 ==0)
      offset = 0;
    else 
      offset = 1;

    for(size_t k=0;k<currentLevel->size();k++) {
      currentLevel->at(k).t.resize(2);

      if(k==0) {
        currentLevel->at(k).t.at(0) = coarseTimeStart;
        currentLevel->at(k).t.at(1) = higherLevel->at(2*k+offset).t.at(1);
      }
      else {
        currentLevel->at(k).t.at(0) = higherLevel->at(2*k+offset-1).t.at(0);
        currentLevel->at(k).t.at(1) = higherLevel->at(2*k+offset).t.at(1);
      }
    }

    return currentLevel;
  }

  /** This method is used to hide the details of rebalancing or not. It's
    * goal is to provide the user an abstraction that seperates the pointer
    * from the restriction algorithm. This is usually paired with "send" and
    * "recv" methods which when not rebalanced, end up doing nothing.
    */
  ROL::Ptr<ROL::Vector<Real>>
  getLocalVector(int k,
                       int myRank,
                       const ROL::Ptr<ROL::PinTVector<Real>> & pint_output,
                       const ROL::Ptr<ROL::Vector<Real>> & scratch,
                       const TimeStepCommMap & commMap) const
  {
    if(commMap.sendToProc_[k]!=myRank)
      return scratch->clone();
    else 
      return pint_output->getVectorPtr(k);
  }

  /** This method sends a vector to a processor it should belong to.
    *
    * \param[in] vec Vector to send to the coarse processor
    * \param[in] index Local coarse index if a straightforward coarsening had been used
    * \param[in] communicators Communicators for the time domain
    * \param[in] vectorComm Object to communicate entire vector between processors
    */
  void sendToProcs(ROL::Vector<Real> & vec,
                         int index,
                         const ROL::PinTCommunicators & communicators,
                         const ROL::PinTVectorCommunication<Real> & vectorComm,
                         const TimeStepCommMap & commMap) const
  {
    int myRank    = communicators.getTimeRank();

    // nothing to do, all in place
    if(commMap.sendToProc_[index]==myRank)
      return;
    
    // nothing to do if reblance isn't on (a no-op
    MPI_Comm comm = communicators.getTimeCommunicator();

    int targetRank = commMap.sendToProc_[index];

    vectorComm.send(comm,targetRank,vec,index);
  }

  /** This method sends a vector to a processor it should belong to.
    *
    * \param[in] vec Vector to send to the coarse processor
    * \param[in] index Local coarse index if a straightforward coarsening had been used
    * \param[in] communicators Communicators for the time domain
    * \param[in] vectorComm Object to communicate entire vector between processors
    */
  void sendToAllProcs(const std::vector<ROL::Ptr<ROL::Vector<Real>>> & input,
                      const ROL::PinTCommunicators & communicators,
                      const ROL::PinTVectorCommunication<Real> & vectorComm,
                      const TimeStepCommMap & commMap) const
  {
    MPI_Comm comm = communicators.getTimeCommunicator();
    int myRank    = communicators.getTimeRank();

    // this little bit of trickery is used to make sure the default value is 0 (e.g. index 0)
    struct DefaultInt { int value=0; };
    std::unordered_map<int,DefaultInt> procToIndex;

    // loop over all vectors you need to recieve, and recieve them
    for(size_t i=0;i<commMap.sendToProc_.size();i++) {
      int proc = commMap.sendToProc_[i];
 
      if(proc==myRank) continue;

      int index = procToIndex[proc].value;

      ROL::Ptr<Vector<Real>> vec = input[i];

      vectorComm.send(comm,proc,*vec,index /* the tag */);

      // increment the index value for this processor (this is used as the tag
      procToIndex[proc].value++;
    }
  }

  /** Recieve a set of vectors that have already been sent.
    *
    * \param[in] pint_output Vector to fill
    * \param[in] communicators Communicators for the time domain
    * \param[in] vectorComm Object to communicate entire vector between processors
    * \param[in] commMap Object describing how the coarse vectors are distributed.
    */
  void recvAllFromProcs(ROL::PinTVector<Real> & pint_output,
                        const ROL::PinTCommunicators & communicators,
                        const ROL::PinTVectorCommunication<Real> & vectorComm,
                        const TimeStepCommMap & commMap) const
  {
    MPI_Comm comm = communicators.getTimeCommunicator();
    int myRank    = communicators.getTimeRank();

    // this little bit of trickery is used to make sure the default value is 0 (e.g. index 0)
    struct DefaultInt { int value=0; };
    std::unordered_map<int,DefaultInt> procToIndex;

    // loop over all vectors you need to recieve, and recieve them
    for(size_t i=0;i<commMap.recvFromProc_.size();i++) {
      int proc = commMap.recvFromProc_[i];
 
      if(proc==myRank) continue;

      int index = procToIndex[proc].value;

      ROL::Ptr<Vector<Real>> coarseVec = pint_output.getVectorPtr(i);

      vectorComm.recv(comm,proc,*coarseVec,index /* the tag */);

      // increment the index value for this processor (this is used as the tag
      procToIndex[proc].value++;
    }
  }

public:

  //! Default constructor
  PinTHierarchy()
    : isInitialized_(false)
    , maxLevels_(-1)
    , rebalance_(false)
  { }

  /**
   * \brief Constructor
   *
   * Build a parallel-in-time constraint with a specified step constraint. This specifies
   * any communication you might need between processors and steps using "stencils".
   * Multigrid in time preconditioning is disabled by default.
   * 
   * \param[in] stepConstraint Constraint for a single step.
   * \param[in] initialCond Initial condition
   */
  PinTHierarchy(const Ptr<std::vector<TimeStamp<Real>>> & timeStamps)
    : isInitialized_(false)
    , maxLevels_(-1)
    , rebalance_(false)
  { 
    initialize(timeStamps);
  }

  /**
   * Turn on multigrid preconditioning in time with a specified number of levels.
   *
   * \param[in] maxLevels Largest number of levels
   */
  void setMaxLevels(int maxLevels)
  {
    maxLevels_ = maxLevels;
  }

  /**
   * \brief Get the time stamps for by level
   *
   * The current implementation is recursive and expensive, but its the easiest
   * way to do this right now. 
   */
  Ptr<std::vector<TimeStamp<Real>>> getTimeStampsByLevel(int level) const
  {
    return stamps_[level];
  }


  /** \brief Initialize this class, setting up parallel distribution.
   
   */
  void initialize(const Ptr<std::vector<TimeStamp<Real>>> & timeStamps)
  {
    // initialize user member variables
    userTimeStamps_ = timeStamps;

    // build up the internally used time stamps
    ////////////////////////////////////////////////////////////////
   
    timeStamps_ = makePtr<std::vector<TimeStamp<Real>>>(timeStamps->size());

    for(size_type k=0;k<timeStamps->size();k++) {
      timeStamps_->at(k).t = timeStamps->at(k).t;
    }

    isInitialized_ = true;
  }

   // restriction and prolongation functions
   ///////////////////////////////////////////////////////////////////////////////////////////
   
   ROL::Ptr<const PinTCommunicators> getLevelCommunicators(int level) const
   { 
     if(not rebalance_) return communicators_[0]; 
     else               return communicators_[level]; 
   }

   /** Indicator if a specified level is active on this rank.
     */
   bool levelIsActiveOnMyRank(int level) const
   { 
     // all levels are active when not rebalanced
     if(not rebalance_) 
       return true;

     return communicators_[level]!=ROL::nullPtr; 
   }

   /** 
    * \brief Builds up communicators and time stamps for each level
    *
    * This currently doubles the size of the time step as it coarsens. This calls
    * the recursive helper function buildLevelDataStructures to do most of the work
    *
    * \param[in] level_0_ref Reference vector that is used to pull out all the communication
    *                        devices.
    */
   void buildLevels(const ROL::Ptr<const ROL::PinTCommunicators> & pintComm,
                    const ROL::Ptr<const ROL::PinTVectorCommunication<Real>> & vectorComm,
                    bool rebalance=false)
   {
     rebalance_ = rebalance;

     // set the vector communicator
     vectorComm_ = vectorComm;

     // we allocate the communicators with null
     communicators_.resize(maxLevels_);
     communicators_[0] = pintComm;

     // we allocate the time stamps vectors with null
     stamps_.resize(maxLevels_);
     stamps_[0] = timeStamps_;

     // allocate the restriction/prolognation maps (on the zeroth level they are empty)
     restrictOptMap_.resize(maxLevels_);
     restrictSimMap_.resize(maxLevels_);
    
     // build prolong optimization variable map
     prolongOptMap_.resize(maxLevels_);
     if(rebalance_) {
       int procSize = pintComm->getTimeSize();
       int rank = pintComm->getTimeRank();
       MPI_Comm mpiComm = pintComm->getTimeCommunicator();
       TimeStepCommMap & commMap = prolongOptMap_[0];

       int recvSize = stamps_[0]->size();
       int recvTarget = -1;
       if(rank % 2 == 0)
         recvTarget = rank/2;
       else
         recvTarget = (rank-1)/2;
       MPI_Send(&recvSize,1,MPI_INT,recvTarget,0,mpiComm);

       commMap.recvFromProc_.resize(stamps_[0]->size(),recvTarget);

       if(2*rank+1<procSize) {
         int sizeA = -1;
         int sizeB = -1;
         MPI_Recv(&sizeA,1,MPI_INT,2*rank  ,0,mpiComm,MPI_STATUS_IGNORE);
         MPI_Recv(&sizeB,1,MPI_INT,2*rank+1,0,mpiComm,MPI_STATUS_IGNORE);

         commMap.sendToProc_.resize(sizeA+sizeB,-1);

         for(int i=0;i<sizeA;i++) 
           commMap.sendToProc_[i      ] = 2*rank;
         for(int i=0;i<sizeB;i++) 
           commMap.sendToProc_[i+sizeA] = 2*rank+1;
       }
     }
     else {
       int rank = pintComm->getTimeRank();
       TimeStepCommMap & commMap = prolongOptMap_[0];
       commMap.recvFromProc_.resize(stamps_[0]->size(),rank);
       commMap.sendToProc_.resize(stamps_[0]->size(),rank);
     }

     // build prolong simulation variable map
     prolongSimMap_.resize(maxLevels_);
     if(rebalance_) {
       int procSize = pintComm->getTimeSize();
       int rank = pintComm->getTimeRank();
       MPI_Comm mpiComm = pintComm->getTimeCommunicator();
       TimeStepCommMap & commMap = prolongSimMap_[0];

       int recvSize = 2*stamps_[0]->size();
       int recvTarget = -1;
       if(rank % 2 == 0)
         recvTarget = rank/2;
       else
         recvTarget = (rank-1)/2;

       MPI_Send(&recvSize,1,MPI_INT,recvTarget,0,mpiComm);

       commMap.recvFromProc_.resize(2*stamps_[0]->size(),recvTarget);

       if(2*rank+1<procSize) {
         int sizeA = -1;
         int sizeB = -1;
         MPI_Recv(&sizeA,1,MPI_INT,2*rank  ,0,mpiComm,MPI_STATUS_IGNORE);
         MPI_Recv(&sizeB,1,MPI_INT,2*rank+1,0,mpiComm,MPI_STATUS_IGNORE);

         commMap.sendToProc_.resize(sizeA+sizeB,-1);

         for(int i=0;i<sizeA;i++) 
           commMap.sendToProc_[i      ] = 2*rank;
         for(int i=0;i<sizeB;i++) 
           commMap.sendToProc_[i+sizeA] = 2*rank+1;
       }
     }
     else {
       int rank = pintComm->getTimeRank();
       TimeStepCommMap & commMap = prolongSimMap_[0];
       commMap.recvFromProc_.resize(2*stamps_[0]->size(),rank);
       commMap.sendToProc_.resize(2*stamps_[0]->size(),rank);
     }

     buildLevelDataStructures(1);
   }

   //! Recursive function that builds up the communicators and time stamps on coarse grids
   void buildLevelDataStructures(int level)
   {
     // protect sanity
     assert(level>0);

     // don't go too deep (base case)
     if(level>=maxLevels_) {
       return;
     }

     auto comm = getLevelCommunicators(level-1);
 
     if(comm==ROL::nullPtr)
       return;

     int rank = comm->getTimeRank();

     stamps_[level] = buildTimeStampsByLevel(level,comm);
     int origSize = stamps_[level]->size();

     if(rebalance_) {
       // this processor is no longer participating (base case)
       if(ROL::is_nullPtr(comm)) {
         return;
       }

       // this will subdivide the communicators by two
       communicators_[level] = comm->buildCoarseCommunicators();

       // rebalance the coarse stamps
       ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> coarseStamps 
         = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>();
       ROL::PinT::exportToCoarseDistribution_TimeStamps(*stamps_[level],*coarseStamps,*comm,0);

       stamps_[level] = coarseStamps;
     }

     // build restriction optimization maps
     if(rebalance_) {
       MPI_Comm mpiComm = comm->getTimeCommunicator();
       TimeStepCommMap & commMap = restrictOptMap_[level];

       int sendTarget = -1;
       if(rank % 2 == 0)
         sendTarget = rank/2;
       else
         sendTarget = (rank-1)/2;

       commMap.sendToProc_.resize(origSize,sendTarget);

       MPI_Send(&origSize,1,MPI_INT,sendTarget,level,mpiComm);

       commMap.recvFromProc_.resize(stamps_[level]->size(),-1);

       if(commMap.recvFromProc_.size()>0) {
         int sizeA = -1;
         int sizeB = -1;
         MPI_Recv(&sizeA,1,MPI_INT,2*rank  ,level,mpiComm,MPI_STATUS_IGNORE);
         MPI_Recv(&sizeB,1,MPI_INT,2*rank+1,level,mpiComm,MPI_STATUS_IGNORE);

         TEUCHOS_ASSERT(static_cast<long long>(sizeA) + static_cast<long long>(sizeB) == static_cast<long long>(stamps_[level]->size()));

         for(int i=0;i<sizeA;i++) 
           commMap.recvFromProc_[i      ] = 2*rank;
         for(int i=0;i<sizeB;i++) 
           commMap.recvFromProc_[i+sizeA] = 2*rank+1;
       }
     }
     else {
       TimeStepCommMap & commMap = restrictOptMap_[level];
 
       commMap.sendToProc_.resize(stamps_[level]->size(),rank);
       commMap.recvFromProc_.resize(stamps_[level]->size(),rank);
     }

     // build restriction simulation maps
     if(rebalance_) {
       MPI_Comm mpiComm = comm->getTimeCommunicator();
       TimeStepCommMap & commMap = restrictSimMap_[level];

       int sendTarget = -1;
       if(rank % 2 == 0)
         sendTarget = rank/2;
       else
         sendTarget = (rank-1)/2;

       int simSize = 2*origSize;
       commMap.sendToProc_.resize(simSize,sendTarget);

       MPI_Send(&simSize,1,MPI_INT,sendTarget,level,mpiComm);

       commMap.recvFromProc_.resize(2*stamps_[level]->size(),-1);

       if(commMap.recvFromProc_.size()>0) {
         int sizeA = -1;
         int sizeB = -1;
         MPI_Recv(&sizeA,1,MPI_INT,2*rank  ,level,mpiComm,MPI_STATUS_IGNORE);
         MPI_Recv(&sizeB,1,MPI_INT,2*rank+1,level,mpiComm,MPI_STATUS_IGNORE);

         TEUCHOS_ASSERT(static_cast<long long>(sizeA) + static_cast<long long>(sizeB) == 2 * static_cast<long long>(stamps_[level]->size()));

         for(int i=0;i<sizeA;i++) 
           commMap.recvFromProc_[i      ] = 2*rank;
         for(int i=0;i<sizeB;i++) 
           commMap.recvFromProc_[i+sizeA] = 2*rank+1;
       }
     }
     else {
       TimeStepCommMap & commMap = restrictSimMap_[level];
 
       commMap.sendToProc_.resize(2*stamps_[level]->size(),rank);
       commMap.recvFromProc_.resize(2*stamps_[level]->size(),rank);
     }

     // build prolongation optimization maps
     if(rebalance_) {
       int procSize = comm->getTimeSize();
       MPI_Comm mpiComm = comm->getTimeCommunicator();
       TimeStepCommMap & commMap = prolongOptMap_[level];

       int recvTarget = -1;
       if(rank % 2 == 0)
         recvTarget = rank/2;
       else
         recvTarget = (rank-1)/2;

       commMap.recvFromProc_.resize(stamps_[level]->size(),recvTarget);

       int recvSize = stamps_[level]->size();
       MPI_Send(&recvSize,1,MPI_INT,recvTarget,level,mpiComm);

       if(2*rank+1<procSize) {
         int sizeA = -1;
         int sizeB = -1;
         MPI_Recv(&sizeA,1,MPI_INT,2*rank  ,level,mpiComm,MPI_STATUS_IGNORE);
         MPI_Recv(&sizeB,1,MPI_INT,2*rank+1,level,mpiComm,MPI_STATUS_IGNORE);

         commMap.sendToProc_.resize(sizeA+sizeB,-1);

         for(int i=0;i<sizeA;i++) 
           commMap.sendToProc_[i      ] = 2*rank;
         for(int i=0;i<sizeB;i++) 
           commMap.sendToProc_[i+sizeA] = 2*rank+1;
       }
     }
     else {
       TimeStepCommMap & commMap = prolongOptMap_[level];
 
       commMap.sendToProc_.resize(stamps_[level]->size(),rank);
       commMap.recvFromProc_.resize(stamps_[level]->size(),rank);
     }

     // build prolongation simulation maps
     if(rebalance_) {
       int procSize = comm->getTimeSize();
       MPI_Comm mpiComm = comm->getTimeCommunicator();
       TimeStepCommMap & commMap = prolongSimMap_[level];

       int recvTarget = -1;
       if(rank % 2 == 0)
         recvTarget = rank/2;
       else
         recvTarget = (rank-1)/2;

       commMap.recvFromProc_.resize(2*stamps_[level]->size(),recvTarget);

       int recvSize = 2*stamps_[level]->size();
       MPI_Send(&recvSize,1,MPI_INT,recvTarget,level,mpiComm);

       if(2*rank+1<procSize) {
         int sizeA = -1;
         int sizeB = -1;
         MPI_Recv(&sizeA,1,MPI_INT,2*rank  ,level,mpiComm,MPI_STATUS_IGNORE);
         MPI_Recv(&sizeB,1,MPI_INT,2*rank+1,level,mpiComm,MPI_STATUS_IGNORE);

         commMap.sendToProc_.resize(sizeA+sizeB,-1);

         for(int i=0;i<sizeA;i++) 
           commMap.sendToProc_[i      ] = 2*rank;
         for(int i=0;i<sizeB;i++) 
           commMap.sendToProc_[i+sizeA] = 2*rank+1;
       }
     }
     else {
       TimeStepCommMap & commMap = prolongSimMap_[level];
 
       commMap.sendToProc_.resize(2*stamps_[level]->size(),rank);
       commMap.recvFromProc_.resize(2*stamps_[level]->size(),rank);
     }

     buildLevelDataStructures(level+1);
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<PinTVector<Real>> allocateSimVector(const Vector<Real> & level_0_ref,int level) const
   {
     // when rebalnacing these are null as we coarsen
     if(not levelIsActiveOnMyRank(level))
       return ROL::nullPtr;

     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     ROL::Ptr<const PinTCommunicators> comm = getLevelCommunicators(level);
     auto vectorComm = pint_ref.vectorCommunicationPtr();
   
     int totalSteps = 0;
     Ptr<std::vector<TimeStamp<Real>>> stamps = getTimeStampsByLevel(level);
     int mySteps = int(stamps->size());
     MPI_Allreduce(&mySteps,&totalSteps,1,MPI_INT,MPI_SUM,comm->getTimeCommunicator());

     return makePtr<PinTVector<Real>>(comm,vectorComm,pint_ref.getVectorPtr(0)->clone(),totalSteps,mySteps,1,2);
   }

   /**
    * \brief Allocate a optimization space vector at a particular multigrid level.
    */
   Ptr<PinTVector<Real>> allocateOptVector(const Vector<Real> & level_0_ref,int level) const
   {
     // when rebalanacing these are null as we coarsen
     if(not levelIsActiveOnMyRank(level))
       return ROL::nullPtr;

     ROL::Ptr<const PinTCommunicators> comm = getLevelCommunicators(level);
  
     TEUCHOS_ASSERT(comm!=Teuchos::null);

     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto vectorComm = pint_ref.vectorCommunicationPtr();
     
     int totalSteps = 0;
     Ptr<std::vector<TimeStamp<Real>>> stamps = getTimeStampsByLevel(level);
     int mySteps = int(stamps->size());
     MPI_Allreduce(&mySteps,&totalSteps,1,MPI_INT,MPI_SUM,comm->getTimeCommunicator());

     return  makePtr<PinTVector<Real>>(comm,vectorComm,pint_ref.getVectorPtr(0)->clone(),totalSteps,mySteps,1,1);
   }
   
   /**
    * \brief Restrict a simulation space vector
    *
    * Currently doing the following. Let the coarse distribution be defined locally
    * to contain Np steps. The fine discretization virtual variable is denoted by a V, while the coarse
    * discretization is denoted by a v (similarly U and u for the primary variable). A superscript indicates the relative processor index.
    * If no index is included then it is assumed to be 0 (denoting no relative change in
    * processor index). The restriction for the sim vector thus computes
    *
    *   v_i = V_{2*i} + V_{2*i+1} + V_{2*i+2}      0 <= i < Np - 1
    *   u_i = U_{2*i} + U_{2*i+1} + U_{2*i+2}      0 <= i < Np - 1
    *
    * This ensures that for all interior points the coarse grid is the sum of the three neighbors. These
    * averages are dividecd by three.
    * At the boundaries of processors the same thing holds, but with the appropriate sharing. Finally,
    * at the origin for v and the terminus for u only the relavant two points are taken and the sum is divided
    * by two.
    *
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    */
   void restrictSimVector(const Vector<Real> & input,Vector<Real> & output,int inputLevel)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     restrictSimVector(ROL::makePtrFromRef(pint_input),ROL::makePtrFromRef(pint_output),restrictSimMap_[inputLevel+1],inputLevel); 
   }

   void restrictSimVector(const ROL::Ptr<const ROL::PinTVector<Real>> & pint_input,
                          const ROL::Ptr<      ROL::PinTVector<Real>> & pint_output,
                          int inputLevel)
   {
     restrictSimVector(pint_input,pint_output,restrictSimMap_[inputLevel+1],inputLevel); 
   }

   void restrictSimVector(const ROL::Ptr<const ROL::PinTVector<Real>> & pint_input,
                          const ROL::Ptr<      ROL::PinTVector<Real>> & pint_output,
                          const TimeStepCommMap & commMap,
                          int inputLevel)
   {
     // nothing to do in this case
     if(not levelIsActiveOnMyRank(inputLevel)) {
       return;
     }

     int timeRank = pint_input->communicators().getTimeRank();
     int numSteps = pint_input->numOwnedSteps();
     std::pair<int,int> fneRange = pint_input->ownedStepRange();

     auto scratch = pint_input->getRemoteBufferPtr(0)->clone();

     int crs_i = 0;
     for(int i=0;i<numSteps;i++) {
       // only do evens
       if((i+fneRange.first) % 2 != 0) 
         continue;

       ROL::Ptr<Vector<Real>> output_u = getLocalVector(2*crs_i,timeRank,pint_output,scratch,commMap);
       output_u->axpy(1.0,*pint_input->getVectorPtr(2*i)); 

       // send to coarse
       sendToProcs(*output_u,                   
                    2*crs_i,                  
                    pint_input->communicators(),
                    *pint_input->vectorCommunicationPtr(),
                    commMap); 

       ROL::Ptr<Vector<Real>> output_v = getLocalVector(2*crs_i+1,timeRank,pint_output,scratch,commMap);
       output_v->axpy(1.0,*pint_input->getVectorPtr(2*i+1)); 

       // send to coarse
       sendToProcs(*output_v,                   
                    2*crs_i+1,                  
                    pint_input->communicators(),
                    *pint_input->vectorCommunicationPtr(),
                    commMap); 

       crs_i++;
     }

     // if the coarse level is inactive on this proccessor skip the recieve step
     if(levelIsActiveOnMyRank(inputLevel+1)) {
       recvAllFromProcs(*pint_output,
                         pint_input->communicators(),          
                         *pint_input->vectorCommunicationPtr(),
                         commMap);
     }
   }

   void restrictOptVector(const ROL::Ptr<const ROL::PinTVector<Real>> & pint_input,
                          const ROL::Ptr<      ROL::PinTVector<Real>> & pint_output,
                          int inputLevel)
   {
     restrictOptVector(pint_input,pint_output,restrictOptMap_[inputLevel+1],inputLevel); 
   }

   void restrictOptVector(const ROL::Ptr<const ROL::PinTVector<Real>> & pint_input,
                          const ROL::Ptr<      ROL::PinTVector<Real>> & pint_output,
                          const TimeStepCommMap & commMap,
                          int inputLevel)
   {
     // nothing to do in this case
     if(not levelIsActiveOnMyRank(inputLevel)) {
       return;
     }

     int fineRank = pint_input->communicators().getTimeRank();

     // communicate points on the left of this interval
     pint_input->boundaryExchangeLeftToRight();
     auto leftStart = pint_input->getRemoteBufferPtr(0)->clone();
     auto scratch = leftStart->clone();
     leftStart->set(*pint_input->getRemoteBufferPtr(0));

     std::pair<int,int> fneRange = pint_input->ownedStepRange();

     // handle interior
     int crsLocal = 0;
     for(int fineLocal=0;fineLocal<pint_input->numOwnedSteps();fineLocal++) {
       int fineIndex = fineLocal+fneRange.first;
 
       ROL::Ptr<Vector<Real>> output = getLocalVector(crsLocal,fineRank,pint_output,scratch,commMap);

       if(fineIndex == 0) {
         output->set(*pint_input->getVectorPtr(0));
       }
       else {
         if(fineLocal>=0)
           output->set(     *pint_input->getVectorPtr(fineLocal));
         else
           output->set(     *leftStart);

         output->axpy(1.0,*pint_input->getVectorPtr(fineLocal+1));
         output->scale(1.0/2.0);

         fineLocal++;
       }

       // send to coarse
       sendToProcs(*output,                   
                    crsLocal,                  
                    pint_input->communicators(),
                    *pint_input->vectorCommunicationPtr(),
                    commMap); 
       crsLocal++;
     }

     // if the coarse level is inactive on this proccessor skip the recieve step
     if(levelIsActiveOnMyRank(inputLevel+1)) {
       recvAllFromProcs(*pint_output,
                         pint_input->communicators(),          
                         *pint_input->vectorCommunicationPtr(),
                         commMap);
     }
   }

   /**
    * \brief Restrict a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void restrictOptVector(const Vector<Real> & input,Vector<Real> & output,int inputLevel)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     restrictOptVector(ROL::makePtrFromRef(pint_input),ROL::makePtrFromRef(pint_output),restrictOptMap_[inputLevel+1],inputLevel); 
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector 
    * \param[in] input The input vector to be "prolonged" 
    * \param[in] output The output vector resulting from prolongation, must be at level inputLevel-1
    */
   void prolongSimVector(const ROL::Ptr<const PinTVector<Real>> & pint_input,
                         const ROL::Ptr<PinTVector<Real>> & pint_output,
                         const TimeStepCommMap & commMap,
                         int inputLevel)
   {
     // nothing to do in this case
     if(levelIsActiveOnMyRank(inputLevel)) {
       int timeRank = pint_output->communicators().getTimeRank();

       // communicate points on the left of this interval
       pint_input->boundaryExchangeLeftToRight();
       auto u_fineStart = pint_input->getRemoteBufferPtr(0)->clone();
       auto v_fineStart = pint_input->getRemoteBufferPtr(1)->clone();
       u_fineStart->set(*pint_input->getRemoteBufferPtr(0));             // this is a little unpalatable since we are
       v_fineStart->set(*pint_input->getRemoteBufferPtr(1));             // allocating more memory. it might be better for
       // the user to pass in a buffer
       //
       auto scratch = u_fineStart->clone();

       // communicate points on the right of this interval
       pint_input->boundaryExchangeRightToLeft();
       auto u_fineEnd = pint_input->getRemoteBufferPtr(0);
       auto v_fineEnd = pint_input->getRemoteBufferPtr(1);

       std::vector<ROL::Ptr<ROL::Vector<Real>>> local_store;

       // for(int fne_i=0;fne_i<fineSteps;fne_i++) {
       int fne_i = 0;
       std::pair<int,int> crsRange = pint_input->ownedStepRange();
       for(int crs_i=0;crs_i<pint_input->numOwnedSteps();crs_i++) {
         int crsIndex = crs_i+crsRange.first;

         // fine indices
         ROL::Ptr<Vector<Real>> output_u;
         ROL::Ptr<Vector<Real>> output_v;

         if(crsIndex>0) {
           int u_fne_index = 2*fne_i;   
           int v_fne_index = 2*fne_i+1;   

           int u_crs_index = -1; 
           int v_crs_index = -1; 

           int coarse_si_0 = crs_i-1;     // coarse step index

           u_crs_index = 2*coarse_si_0;
           v_crs_index = 2*coarse_si_0+1;

           output_u = getLocalVector(u_fne_index,timeRank,pint_output,scratch,commMap);
           output_v = getLocalVector(v_fne_index,timeRank,pint_output,scratch,commMap);

           local_store.push_back(output_u);
           local_store.push_back(output_v);

           // this special case is used when the first fine time step does not belong to the 
           // coarse grid
           if(coarse_si_0 < 0) {
             output_u->set(*u_fineStart);
             output_v->set(*v_fineStart);
           }
           else {
             output_u->set(*pint_input->getVectorPtr(u_crs_index));
             output_v->set(*pint_input->getVectorPtr(v_crs_index));
           }

           u_crs_index += 2;
           v_crs_index += 2;

           output_u->axpy(1.0,*pint_input->getVectorPtr(u_crs_index));
           output_v->axpy(1.0,*pint_input->getVectorPtr(v_crs_index));

           // average the two points
           output_u->scale(0.5);
           output_v->scale(0.5);
 
           fne_i++;
         }

         // inject into the fine index
         {
           int u_fne_index = 2*fne_i;   
           int v_fne_index = 2*fne_i+1;   

           int u_crs_index = 2*crs_i;
           int v_crs_index = 2*crs_i+1;

           output_u = getLocalVector(u_fne_index,timeRank,pint_output,scratch,commMap);
           output_v = getLocalVector(v_fne_index,timeRank,pint_output,scratch,commMap);

           local_store.push_back(output_u);
           local_store.push_back(output_v);

           output_u->set(*pint_input->getVectorPtr(u_crs_index));
           output_v->set(*pint_input->getVectorPtr(v_crs_index));

           fne_i++;
         }
       }

       sendToAllProcs(local_store,
           pint_output->communicators(),          
           *pint_output->vectorCommunicationPtr(),
           commMap);
     }

     // if the fine level is inactive on this proccessor skip the recieve step
     if(levelIsActiveOnMyRank(inputLevel-1)) {
       recvAllFromProcs(*pint_output,
                         pint_output->communicators(),          
                         *pint_output->vectorCommunicationPtr(),
                         commMap);
     }
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] input The input vector to be "prolonged" 
    * \param[in] output The output vector resulting from prolongation, must be at level inputLevel-1
    * \param[in] inputLevel Multigrid level of the input vector 
    */
   void prolongSimVector(const ROL::Ptr<const PinTVector<Real>> & pint_input,
                         const ROL::Ptr<PinTVector<Real>> & pint_output,
                         int inputLevel)
   {
     prolongSimVector(pint_input,pint_output,prolongSimMap_[inputLevel-1],inputLevel); 
   }

   /**
    * \brief Prolong a simulation space vector
    *
    *   X_{2*i+1} = x_i                           -1 <= i < Np      (injection variables, including virtual variable)
    *   X_{2*i+2} = (x_i + x_{i+1})/2              0 <= i < Np - 1  (averaged variables, created at fine level)
    *
    * Note: Currently if the timeRank==0, the initial condition is assumed to be zero. 
    *
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    */
   void prolongSimVector(const Vector<Real> & input,Vector<Real> & output,int inputLevel)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     prolongSimVector(ROL::makePtrFromRef(pint_input),ROL::makePtrFromRef(pint_output),inputLevel); 
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector 
    * \param[in] input The input vector to be "prolonged" 
    * \param[in] output The output vector resulting from prolongation, must be at level inputLevel-1
    */
   void prolongOptVector(const ROL::Ptr<const PinTVector<Real>> & pint_input,
                         const ROL::Ptr<PinTVector<Real>> & pint_output,
                         const TimeStepCommMap & commMap,
                         int inputLevel)
   {
     // nothing to do in this case
     if(levelIsActiveOnMyRank(inputLevel)) {
       int fineRank = pint_output->communicators().getTimeRank();

       // communicate points on the right of this interval
       pint_input->boundaryExchangeRightToLeft();
       auto rightStart = pint_input->getRemoteBufferPtr(0)->clone();
       rightStart->set(*pint_input->getRemoteBufferPtr(0));
       auto scratch = rightStart->clone();

       std::pair<int,int> crsRange = pint_input->ownedStepRange();

       // handle interior
       int fineLocal = 0;

       std::vector<ROL::Ptr<ROL::Vector<Real>>> local_store;

       for(int crsLocal=0;crsLocal<pint_input->numOwnedSteps();crsLocal++) {
         int crsIndex = crsLocal+crsRange.first;

         ROL::Ptr<Vector<Real>> output;

         if(crsIndex==0) {
           output = getLocalVector(fineLocal,fineRank,pint_output,scratch,commMap);
           local_store.push_back(output);

           output->set(*pint_input->getVectorPtr(0)); 

           fineLocal++;
         }
         else {
           // send the first time step
           output = getLocalVector(fineLocal,fineRank,pint_output,scratch,commMap);
           local_store.push_back(output);

           if(crsLocal<pint_input->numOwnedSteps())
             output->set(*pint_input->getVectorPtr(crsLocal)); 
           else
             output->set(*rightStart);

           fineLocal++;

           /////////////////////////////////////////////////////////

           // send the second time step
           output = getLocalVector(fineLocal,fineRank,pint_output,scratch,commMap);
           local_store.push_back(output);

           if(crsLocal<pint_input->numOwnedSteps())
             output->set(*pint_input->getVectorPtr(crsLocal)); 
           else
             output->set(*rightStart);

           fineLocal++;
         }
       }

       sendToAllProcs(local_store,
           pint_output->communicators(),          
           *pint_output->vectorCommunicationPtr(),
           commMap);
     }

     // if the fine level is inactive on this proccessor skip the recieve step
     if(levelIsActiveOnMyRank(inputLevel-1)) {
       recvAllFromProcs(*pint_output,
                         pint_output->communicators(),          
                         *pint_output->vectorCommunicationPtr(),
                         commMap);
     }
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] input The input vector to be "prolonged" 
    * \param[in] output The output vector resulting from prolongation, must be at level inputLevel-1
    * \param[in] inputLevel Multigrid level of the input vector 
    */
   void prolongOptVector(const ROL::Ptr<const PinTVector<Real>> & pint_input,
                         const ROL::Ptr<PinTVector<Real>> & pint_output,
                         int inputLevel)
   {
     prolongOptVector(pint_input,pint_output,prolongOptMap_[inputLevel-1],inputLevel); 
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector 
    * \param[in] input The input vector to be "prolonged" 
    * \param[in] output The output vector resulting from prolongation, must be at level inputLevel-1
    */
   void prolongOptVector(const Vector<Real> & input,Vector<Real> & output,int inputLevel)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     prolongOptVector(ROL::makePtrFromRef(pint_input),ROL::makePtrFromRef(pint_output),inputLevel); 
   }

}; // ROL::PinTConstraint


} // namespace ROL 

#endif // ROL_PINTCONSTRAINT_HPP
