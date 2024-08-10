
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTVECTOR_H
#define ROL_PINTVECTOR_H

#include <ostream>

#include "mpi.h"

#include "ROL_Vector.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_StdVector.hpp"

#include "ROL_PinTCommunicators.hpp"
#include "ROL_PinTVectorCommunication.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

/** @ingroup la_group
    \class ROL::PinTVector
    \brief Defines a parallel in time vector.
*/

namespace ROL {

template <class Real>
class PinTVector
 : public Vector<Real>
{

protected:
  // Class to build all the communciators

  // members
  bool isInitialized_;

  Ptr<const PinTCommunicators> communicators_;

  Ptr<Vector<Real>> localVector_;
  int globalSteps_;
  int localSteps_;
  std::vector<int> stencil_;

  // parallel distribution information
  int stepStart_;
  int stepEnd_;
  int replicate_;
  int bufferSize_; // boundary exchange buffer size

  mutable Ptr<PinTVector<Real>> dual_;

  Ptr<PartitionedVector<Real>> stepVectors_;
  std::vector<Ptr<Vector<Real>>> bufferVectors_;

  std::vector<Ptr<Vector<Real>>> leftVectors_;
  std::vector<Ptr<Vector<Real>>> rightVectors_;

  Ptr<const PinTVectorCommunication<Real>> vectorComm_;


  // Using parallel communication and a linear decomposition
  // determine where this processor lives in the global
  // scheme of things
  void computeStepStartEnd()
  {
    int numRanks = communicators_->getTimeSize();
    int myRank   = communicators_->getTimeRank();

    bool stepOff = false;

    if(localSteps_<0 && globalSteps_>=0) {

      // determine which steps are owned by this processor
      int stepsPerRank = globalSteps_ / numRanks;
      int remainder    = globalSteps_ % numRanks; 

      stepStart_ = 0;

      if(myRank<remainder) {
        stepStart_ = myRank*(stepsPerRank+1);
        stepEnd_   = (myRank+1)*(stepsPerRank+1);
      }
      else if(myRank==remainder) {
        stepStart_ = myRank*(stepsPerRank+1);
        stepEnd_   = (myRank+1)*stepsPerRank + myRank;
      }
      else if(myRank>remainder) {
        stepStart_ = myRank*stepsPerRank + remainder;
        stepEnd_   = (myRank+1)*stepsPerRank + remainder;
      }

      localSteps_ = stepEnd_-stepStart_;
    }

   //  else if(localSteps_>=0) 
    {
      stepOff = true;

      int foundTotalSteps = 0;
      MPI_Allreduce(&localSteps_,&foundTotalSteps,1,MPI_INT,MPI_SUM,communicators_->getTimeCommunicator());

      if(globalSteps_>=0) 
        assert(foundTotalSteps==globalSteps_);
      else 
        globalSteps_ = foundTotalSteps;

      stepStart_ = 0;
      MPI_Exscan(&localSteps_,&stepStart_,1,MPI_INT,MPI_SUM,communicators_->getTimeCommunicator());
      if(myRank==0) stepStart_ = 0;

      stepEnd_ = stepStart_+localSteps_;
    }

    // else {
    //   std::stringstream ss;
    //   ss << "ROL::PinTVector::must assign positive number to one of local or global steps: found local = " 
    //      << localSteps_ << ", global = " << globalSteps_;
    //   throw std::logic_error(ss.str());
    // }
 
    int input = stepEnd_>globalSteps_;
    int test = 0;
    MPI_Allreduce(&input,&test,1,MPI_INT,MPI_SUM,communicators_->getTimeCommunicator());

    if(test) {
      std::cout << "P" << myRank << ". step range = " << stepStart_ << " " << stepEnd_ << " of " << globalSteps_ << " " << stepOff << std::endl;
    }
  }

  void allocateBoundaryExchangeVectors(int sz)
  {
    bufferVectors_.resize(sz);

    for(int i=0;i<sz;i++) {
      bufferVectors_[i]  = localVector_->clone();
      bufferVectors_[i]->set(*localVector_);      // make sure each subvector is initialized
    }
  }

public:
  
  typedef enum {SEND_AND_RECV, SEND_ONLY, RECV_ONLY} ESendRecv;

  PinTVector()
    : isInitialized_(false)
    , globalSteps_(-1)
    , localSteps_(-1)
    , stepStart_(-1)
    , stepEnd_(-1)
    , replicate_(-1)
    , bufferSize_(-1)
  {}

  PinTVector(const PinTVector & v)
    : isInitialized_(false)
    , globalSteps_(-1)
    , localSteps_(-1)
    , stepStart_(-1)
    , stepEnd_(-1)
    , replicate_(-1)
    , bufferSize_(-1)
  {
    initialize(v.communicators_,v.vectorComm_,v.localVector_,v.globalSteps_,v.localSteps_,v.bufferSize_/v.replicate_,v.replicate_,v.stepStart_,v.stepEnd_);
  }

  PinTVector(const Ptr<const PinTCommunicators> & comm,
             const Ptr<const PinTVectorCommunication<Real>> & vectorComm,
             const Ptr<Vector<Real>> & localVector,
             int steps,
             const std::vector<int> & stencil,
             int replicate=1)
    : isInitialized_(false)
    , globalSteps_(-1)
    , localSteps_(-1)
    , stepStart_(-1)
    , stepEnd_(-1)
    , replicate_(-1)
    , bufferSize_(-1)
  {
    initialize(comm,vectorComm,localVector,steps,-1,1,replicate,-1,-1);
    TEUCHOS_ASSERT(false);
  }

  PinTVector(const Ptr<const PinTCommunicators> & comm,
             const Ptr<const PinTVectorCommunication<Real>> & vectorComm,
             const Ptr<Vector<Real>> & localVector,
             int globalSteps,
             int localSteps,
             int bufferSize,
             int replicate=1)
    : isInitialized_(false)
    , globalSteps_(-1)
    , localSteps_(-1)
    , stepStart_(-1)
    , stepEnd_(-1)
    , replicate_(-1)
    , bufferSize_(-1)
  {
    initialize(comm,vectorComm,localVector,globalSteps,localSteps,bufferSize,replicate,-1,-1);
  }

  virtual ~PinTVector() {}

  void initialize(const Ptr<const PinTCommunicators> & comm,
                  const Ptr<const PinTVectorCommunication<Real>> & vectorComm,
                  const Ptr<Vector<Real>> & localVector,
                  int globalSteps,
                  int localSteps,
                  int bufferSize,
                  int replicate,
                  int stepStart,
                  int stepEnd)
  {
    replicate_     = replicate;
    bufferSize_    = bufferSize*replicate;
    communicators_ = comm;
    localVector_   = localVector;
    globalSteps_   = globalSteps;
    localSteps_    = localSteps;
    vectorComm_    = vectorComm; // makePtr<PinTVectorCommunication_StdVector<Real>>();

    if(stepStart<0)
      computeStepStartEnd();
    else {
      stepStart_ = stepStart;
      stepEnd_ = stepEnd;
    }

    assert(stepStart_>=0);
    assert(stepEnd_>stepStart_);
    assert(stepEnd_<=globalSteps_);

    allocateBoundaryExchangeVectors(bufferSize_);

    std::vector<Ptr<Vector<Real>>> stepVectors;
    stepVectors.resize(replicate*(stepEnd_ - stepStart_));

    // build up local vectors
    for(int i=0;i<(int) stepVectors.size();i++) {
      stepVectors[i]  = localVector_->clone();
      stepVectors[i]->set(*localVector_);      // make sure each subvector is initialized
    }

    stepVectors_ = makePtr<PartitionedVector<Real>>(stepVectors);

    isInitialized_ = true;
  }

  /** How many vectors are "owned" by this processor. This also includes
      any duplicate vectors.
      
      It may be the case that
      this processor belongs to a sub group of processors. All processors
      in that sub group will return the same number of owned vectors.
    */
  int numOwnedVectors() const 
  { 
    return stepVectors_->numVectors(); 
  }

  /** How many steps are "owned" by this processor. 
      
      It may be the case that
      this processor belongs to a sub group of processors. All processors
      in that sub group will return the same number of owned steps.
    */
  int numOwnedSteps() const 
  { 
    return numOwnedVectors()/replicate_;
  }

  std::pair<int,int> ownedStepRange() const
  { return std::make_pair(stepStart_,stepEnd_); }

  /** What is the stencil used to build this vector?
 
      This accessor is directly based on what the user intiializes the
      vector with. Its used by ROL algorithms and constraints to ensure
      the communication patterns are understood.
    */
  const std::vector<int> & stencil() const { return stencil_; }

  /** What is the communicators object used to build this vector?
    */
  const PinTCommunicators & communicators() const { return *communicators_; }

  /** What is the communicators object used to build this vector?
    */
  Ptr<const PinTCommunicators> communicatorsPtr() const { return communicators_; }

  /** What is the communicators object used to build this vector?
    */
  Ptr<const PinTVectorCommunication<Real>> vectorCommunicationPtr() const { return vectorComm_; }

  /** \brief Determine if an index is valid.
   */
  bool isValidIndex(int i) const
  {
    if(i<0 || i>=numOwnedVectors()) 
      return false;

    return true;
  }

  /** Get a vector pointer. 
   */
  Ptr<Vector<Real>> getVectorPtr(int i) const
  {
    if(not isValidIndex(i))
      return nullPtr;

    return stepVectors_->get(i);
  }

  /** Get a vector pointer to the local buffer. 
   */
  Ptr<Vector<Real>> getRemoteBufferPtr(int i) const
  {
    assert( 0<=i || i < static_cast<int>(bufferVectors_.size()) );

    return bufferVectors_[i];
  }

  /** \brief Exchange unknowns with neighboring processors.
      
      Send the last vector on the sending processor to the remote buffer on the receiving
      processor. In this case the time rank of the sending process is less than the recieving.

      \note This method is const because it doesn't change the behavior
            of the vector. It does ensure the components are correct on processor.
  */
  void boundaryExchangeLeftToRight(ESendRecv   send_recv=SEND_AND_RECV) const
  {
    // build the send buffer 
    std::vector<Ptr<Vector<Real>>> sendBuffer;
    for(int i=0;i<bufferSize_;i++) {
        sendBuffer.push_back(getVectorPtr(numOwnedVectors()-bufferSize_+i));
    }

    boundaryExchangeLeftToRight(sendBuffer,send_recv);
  }

  /** \brief Exchange unknowns with neighboring processors.
      
      Send the last vector on the sending processor to the remote buffer on the receiving
      processor. In this case the time rank of the sending process is less than the recieving.

      \note This method is const because it doesn't change the behavior
            of the vector. It does ensure the components are correct on processor.
  */
  void boundaryExchangeLeftToRight(const std::vector<Ptr<Vector<Real>>> & sendBuffer,
                                   ESendRecv   send_recv=SEND_AND_RECV) const
  {
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("PinTVector::boundaryExchangeLeftToRight");

    assert( static_cast<int>(sendBuffer.size()) <= bufferSize_ );

    MPI_Comm timeComm = communicators_->getTimeCommunicator();
    int      myRank   = communicators_->getTimeRank();

    bool sendToRight   = stepEnd_ < globalSteps_;
    bool recvFromLeft = stepStart_ > 0;

    // this allows finer granularity of control of send recieve
    // and will permit some blocking communication
    if(send_recv==SEND_ONLY) {
     recvFromLeft = false;
    }
    if(send_recv==RECV_ONLY) {
     sendToRight = false;
    }
    // do nothing if(send_recv==SEND_AND_RECV) 
    
    // send from left to right
    int bufferSize = int(sendBuffer.size());
    for(int i=0;i<bufferSize;i++) {

      if(sendToRight) {
        vectorComm_->send(timeComm,myRank+1,*sendBuffer[i],i); // this is "owned"
      }
      
      if(recvFromLeft) {
        vectorComm_->recv(timeComm,myRank-1,*getRemoteBufferPtr(i),false,i);                  
      }
    }

    timer->stop("PinTVector::boundaryExchangeLeftToRight");
  }

  /** \brief Exchange unknowns with neighboring processors.
      
      Send the first vector on the sending processor to the remote buffer on the receiving
      processor. In this case the time rank of the sending process is greater than the recieving.

      \note This method is const because it doesn't change the behavior
            of the vector. It does ensure the components are correct on processor.
  */
  void boundaryExchangeRightToLeft(ESendRecv   send_recv=SEND_AND_RECV) const
  {
    // build the send buffer 
    std::vector<Ptr<Vector<Real>>> sendBuffer;
    for(int i=0;i<bufferSize_;i++) {
        sendBuffer.push_back(getVectorPtr(i));
    }

    boundaryExchangeRightToLeft(sendBuffer,send_recv);
  }

  /** \brief Exchange unknowns with neighboring processors.
      
      Send the first vector on the sending processor to the remote buffer on the receiving
      processor. In this case the time rank of the sending process is greater than the recieving.

      \note This method is const because it doesn't change the behavior
            of the vector. It does ensure the components are correct on processor.
  */
  void boundaryExchangeRightToLeft(const std::vector<Ptr<Vector<Real>>> & sendBuffer,
                                   ESendRecv   send_recv=SEND_AND_RECV) const
  {
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("PinTVector::boundaryExchangeRightToLeft");

    MPI_Comm timeComm = communicators_->getTimeCommunicator();
    int      myRank   = communicators_->getTimeRank();

    bool recvFromRight = stepEnd_ < globalSteps_;
    bool sendToLeft   = stepStart_ > 0;

    // this allows finer granularity of control of send recieve
    // and will permit some blocking communication
    if(send_recv==SEND_ONLY) {
      recvFromRight = false;
    }
    if(send_recv==RECV_ONLY) {
      sendToLeft = false;
    }
    // do nothing if(send_recv==SEND_AND_RECV) 

    // send from right to left
    int bufferSize = int(sendBuffer.size());
    for(int i=0;i<bufferSize;i++) {
      if(sendToLeft)
        vectorComm_->send(timeComm,myRank-1,*sendBuffer[i],i); // this is "owned"

      if(recvFromRight)
        vectorComm_->recv(timeComm,myRank+1,*getRemoteBufferPtr(i),false,i);
    }

    timer->stop("PinTVector::boundaryExchangeRightToLeft");
  }

  /** \brief Exchange unknowns with neighboring processors.
      
      Using the stencil information, do global communication with time neighbor
      processors.

      \note This method is const because it doesn't change the behavior
            of the vector. It does ensure the components are correct on processor.
  */
  void boundaryExchange(ESendRecv   send_recv=SEND_AND_RECV) const
  {
    MPI_Comm timeComm = communicators_->getTimeCommunicator();
    int      myRank   = communicators_->getTimeRank();

    bool sendToRight   = stepEnd_ < globalSteps_;
    bool recvFromRight = stepEnd_ < globalSteps_;

    bool sendToLeft   = stepStart_ > 0;
    bool recvFromLeft = stepStart_ > 0;

    // this allows finer granularity of control of send recieve
    // and will permit some blocking communication
    if(send_recv==SEND_ONLY) {
     recvFromRight = false;
     recvFromLeft = false;
    }
    if(send_recv==RECV_ONLY) {
     sendToRight = false;
     sendToLeft = false;
    }
    // do nothing if(send_recv==SEND_AND_RECV) 
    
    // send from left to right
    for(std::size_t i=0;i<stencil_.size();i++) {
      int offset = stencil_[i];
      if(offset >= 0)
        continue;

      if(sendToRight)
        vectorComm_->send(timeComm,myRank+1,*getVectorPtr(numOwnedVectors()+offset)); // this is "owned"
      
      if(recvFromLeft)
        vectorComm_->recv(timeComm,myRank-1,*getRemoteBufferPtr(offset),false);                  
    }

    // send from right to left
    for(std::size_t i=0;i<stencil_.size();i++) {
      int offset = stencil_[i];
      if(offset <= 0)
        continue;

      if(sendToLeft)
        vectorComm_->send(timeComm,myRank-1,*getVectorPtr(offset-1)); // this is "owned"
      
      if(recvFromRight)
        vectorComm_->recv(timeComm,myRank+1,*getRemoteBufferPtr(numOwnedSteps()+offset-1),false);
    }
  }

  /** \brief Exchange unknowns with neighboring processors using summation.
      
      Using the stencil information, do global communication with time neighbor
      processors. This sums in the the reverse direction.

      \note This method is not const because it does change the behavior
            of the vector.
  */
  void boundaryExchangeSumInto(ESendRecv   send_recv=SEND_AND_RECV)
  {
    MPI_Comm timeComm = communicators_->getTimeCommunicator();
    int      myRank   = communicators_->getTimeRank();

    bool sendToRight   = stepEnd_ < globalSteps_;
    bool recvFromRight = stepEnd_ < globalSteps_;

    bool sendToLeft   = stepStart_ > 0;
    bool recvFromLeft = stepStart_ > 0;

    // this allows finer granularity of control of send recieve
    // and will permit some blocking communication
    if(send_recv==SEND_ONLY) {
     recvFromRight = false;
     recvFromLeft = false;
    }
    if(send_recv==RECV_ONLY) {
     sendToRight = false;
     sendToLeft = false;
    }
    // do nothing if(send_recv==SEND_AND_RECV) 

    // send from left to right
    for(std::size_t i=0;i<stencil_.size();i++) {
      int offset = stencil_[i];
      if(offset >= 0)
        continue;

      if(recvFromRight) {
        vectorComm_->recv(timeComm,myRank+1,*getVectorPtr(numOwnedSteps()+offset),true); // this is "owned"
      }
      
      if(sendToLeft) {
        vectorComm_->send(timeComm,myRank-1,*getRemoteBufferPtr(offset));                  
      }
    }

    // send from right to left
    for(std::size_t i=0;i<stencil_.size();i++) {
      int offset = stencil_[i];
      if(offset <= 0)
        continue;

      if(recvFromLeft)
        vectorComm_->recv(timeComm,myRank-1,*getVectorPtr(offset-1),true); // this is "owned"
      
      if(sendToRight)
        vectorComm_->send(timeComm,myRank+1,*getRemoteBufferPtr(numOwnedSteps()+offset-1));
    }
  }

  // Vector virtual methods

  /** \brief Compute \f$y \leftarrow y + x\f$, where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector to be added to \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \mathtt{*this} + x\f$.

             ---
  */
  virtual void plus( const Vector<Real> &x ) override
  {
    typedef PinTVector<Real> PinTVector; 
    const PinTVector &xs = dynamic_cast<const PinTVector&>(x);

    stepVectors_->plus(*xs.stepVectors_);
  }


  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \alpha (\mathtt{*this}) \f$.

             ---
  */
  virtual void scale( const Real alpha ) override
  {
    stepVectors_->scale(alpha);
  }


  /** \brief Compute \f$ \langle y,x \rangle \f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector that forms the dot product with \f$\mathtt{*this}\f$.
             @return         The number equal to \f$\langle \mathtt{*this}, x \rangle\f$.

             ---
  */
  virtual Real dot( const Vector<Real> &x ) const override
  {
    // this is probably very inefficient way to do this... oh well!
    typedef PinTVector<Real> PinTVector; 
    const PinTVector &xs = dynamic_cast<const PinTVector&>(x);

    // this won't work for Real!=double...oh well!
    Real subdot = stepVectors_->dot(*xs.stepVectors_);
    if(communicators_->getSpaceRank()!=0)  // the space vectors compute the right 'dot', but we need to get rid of them
      subdot = 0.0;

    // NOTE: it might be cheaper to use the TimeCommunicator on the all reduce and then
    //       broadcast to the space communicator
    Real dot = 0;
    MPI_Allreduce(&subdot,&dot,1,MPI_DOUBLE,MPI_SUM,communicators_->getParentCommunicator());

    return dot;
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mathtt{*this}\f$.

             @return         A nonnegative number equal to the norm of \f$\mathtt{*this}\f$.

             ---
  */
  virtual Real norm() const override
  {
    // this is probably very inefficient way to do this... oh well!
    return std::sqrt(this->dot(*this));
  }

  /** \brief Clone to make a new (uninitialized) vector.

             @return         A reference-counted pointer to the cloned vector.

             Provides the means of allocating temporary memory in ROL.

             ---             
  */
  virtual ROL::Ptr<Vector<Real>> clone() const override
  {
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("PinTVector::clone");

    auto v = clonePinT();
    
    timer->stop("PinTVector::clone");

    return v;
  }

  /** \brief Clone that provides direct access to a PinT vector.

             @return         A reference-counted pointer to the cloned vector.
  */
  ROL::Ptr<PinTVector<Real>> clonePinT() const
  {
    return makePtr<PinTVector<Real>>(*this);
  }

  virtual void print( std::ostream &outStream) const override
  {
    stepVectors_->print(outStream);
  }

  virtual void applyUnary( const Elementwise::UnaryFunction<Real> &f ) override {
    stepVectors_->applyUnary(f);
  }

  /** \brief Set \f$y \leftarrow x\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = x\f$.
             Uses #zero and #plus methods for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void set( const Vector<Real> &x ) override {
    typedef PinTVector<Real> PinTVector; 
    const PinTVector &xs = dynamic_cast<const PinTVector&>(x);

    stepVectors_->set(*xs.stepVectors_);
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of @b x.
             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = \mathtt{*this} + \alpha x \f$.
             Uses #clone, #set, #scale and #plus for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void axpy( const Real alpha, const Vector<Real> &x ) override {
    typedef PinTVector<Real> PinTVector; 
    const PinTVector &xs = dynamic_cast<const PinTVector&>(x);

    stepVectors_->axpy(alpha,*xs.stepVectors_);
  }

  /** \brief Zeros all entries of the vector, including boundary exchange fields.
   */
  virtual void zeroAll() 
  {
    stepVectors_->zero();

    for(std::size_t i=0;i<bufferVectors_.size();i++)
      bufferVectors_[i]->zero();
  }

#if 0
  /** \brief Set to zero vector.

             Uses #scale by zero for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void zero();

  /** \brief Return dimension of the vector space.

             @return The dimension of the vector space, i.e., the total number of basis vectors.

             Overload if the basis is overloaded.

             ---
  */
  virtual int dimension() const;
#endif

}; // class Vector

} // namespace ROL

#endif
