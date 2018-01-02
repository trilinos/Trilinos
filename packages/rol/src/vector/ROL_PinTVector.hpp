
// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_PINTVECTOR_H
#define ROL_PINTVECTOR_H

#include <ostream>

#include "mpi.h"

#include "ROL_Vector.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup la_group
    \class ROL::PinTVector
    \brief Defines a parallel in time vector.
*/

namespace ROL {

/** 
 * This class handles PinT communication. Its
 * based on the approach taken in X-Braid.
 */
class PinTCommunicators {
public:
  PinTCommunicators(MPI_Comm parent,int spatialProcs)
  { 
    parentComm_ = parent; 

    int myGlobalRank = -1;
    int globalProcs = -1;
    MPI_Comm_size(parentComm_, &globalProcs);
    MPI_Comm_rank(parentComm_, &myGlobalRank);

    // make sure they divide evenly (for sanity!)
    assert(globalProcs % spatialProcs == 0);

    int space = myGlobalRank / spatialProcs;
    int time  = myGlobalRank % spatialProcs;

    // this decomposition comes from X-Braid 
    MPI_Comm_split(parentComm_,space,myGlobalRank,&spaceComm_); 
    MPI_Comm_split(parentComm_, time,myGlobalRank, &timeComm_); 

    MPI_Comm_size(timeComm_, &timeSize_);
    MPI_Comm_rank(timeComm_, &timeRank_);
  }

  // cleanup
  ~PinTCommunicators()
  { MPI_Comm_free(&spaceComm_);  
    MPI_Comm_free(&timeComm_); }

   MPI_Comm getParentCommunicator() const { return parentComm_; }
   MPI_Comm getSpaceCommunicator() const { return spaceComm_; }
   MPI_Comm getTimeCommunicator() const { return timeComm_; }

   int getTimeRank() const { return timeRank_; }
   int getTimeSize() const { return timeSize_; }

private:

  MPI_Comm parentComm_;

  MPI_Comm spaceComm_;
  MPI_Comm timeComm_;
  int timeSize_;
  int timeRank_;
};

template <class Real> 
class PinTVectorCommunication {
public:
  void send(MPI_Comm comm,int rank,Vector<Real> & source) const
  {
    const std::vector<Real> & std_source = *dynamic_cast<const StdVector<Real>&>(source).getVector();

    MPI_Send(const_cast<Real*>(&std_source[0]),int(std_source.size()),MPI_DOUBLE,rank,0,comm);
  }

  void recv(MPI_Comm comm,int rank,Vector<Real> & dest) const
  {
    std::vector<Real> & std_dest = *dynamic_cast<StdVector<Real>&>(dest).getVector();

    MPI_Recv(&std_dest[0],int(std_dest.size()),MPI_DOUBLE,rank,0,comm,MPI_STATUS_IGNORE);
  }
};

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
  int steps_;
  std::vector<int> stencil_;

  // parallel distribution information
  int stepStart_;
  int stepEnd_;

  Ptr<PartitionedVector<Real>> stepVectors_;
  std::vector<Ptr<Vector<Real>>> leftVectors_;
  std::vector<Ptr<Vector<Real>>> rightVectors_;

  Ptr<PinTVectorCommunication<Real>> vectorComm_;

  // Using parallel communication and a linear decomposition
  // determine where this processor lives in the global
  // scheme of things
  void computeStepStartEnd(int steps)
  {
    int numRanks = communicators_->getTimeSize();
    int myRank   = communicators_->getTimeRank();

    // determine which steps are owned by this processor
    {
      int stepsPerRank = steps / numRanks;
      int remainder    = steps % numRanks; 

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
    }

    assert(stepStart_>=0);
    assert(stepEnd_>stepStart_);
  }

  void allocateBoundaryExchangeVectors()
  {
    int numLeft = 0;
    int numRight = 0;
 
    for(int i=0;i<stencil_.size();i++) {
      if(stencil_[i]<0) 
        numLeft = std::max(numLeft,std::abs(stencil_[i]));
      else if(stencil_[i]>0) 
        numRight = std::max(numRight,stencil_[i]);
    }

    // there is a slight over allocation here if the stencil is sparse
    leftVectors_.resize(numLeft);
    for(int i=0;i<numLeft;i++) {
      leftVectors_[i]  = localVector_->clone();
      leftVectors_[i]->scale(0.0);      // make sure each subvector is initialized
    }

    // there is a slight over allocation here if the stencil is sparse
    rightVectors_.resize(numRight);
    for(int i=0;i<numRight;i++) {
      rightVectors_[i]  = localVector_->clone();
      rightVectors_[i]->scale(0.0);      // make sure each subvector is initialized
    }
  }

public:

  PinTVector()
    : isInitialized_(false)
  {}

  PinTVector(const PinTVector & v)
  {
    isInitialized_ = v.isInitialized_;
    communicators_ = v.communicators_;
    localVector_   = v.localVector_;
    steps_         = v.steps_;
    stencil_       = v.stencil_;
    vectorComm_    = v.vectorComm_;
    stepVectors_   = dynamicPtrCast<PartitionedVector<Real>>(v.stepVectors_->clone());

    computeStepStartEnd(steps_);
  }

  PinTVector(const Ptr<const PinTCommunicators> & comm,
             const Ptr<Vector<Real>> & localVector,
             int steps,
             const std::vector<int> & stencil)
    : isInitialized_(false)
  {
    initialize(comm,localVector,steps,stencil);
  }

  virtual ~PinTVector() {}

  void initialize(const Ptr<const PinTCommunicators> & comm,
                  const Ptr<Vector<Real>> & localVector,
                  int steps,
                  const std::vector<int> & stencil)
  {
    communicators_ = comm;
    localVector_   = localVector;
    steps_         = steps;
    stencil_       = stencil;
    vectorComm_  = makePtr<PinTVectorCommunication<Real>>();

    computeStepStartEnd(steps_);
    allocateBoundaryExchangeVectors();
    
    // build up local vectors
    std::vector<Ptr<Vector<Real>>> stepVectors(stepEnd_-stepStart_);
    for(int i=stepStart_;i<stepEnd_;i++) {
      stepVectors[i-stepStart_]  = localVector_->clone();
      stepVectors[i-stepStart_]->set(*localVector_);      // make sure each subvector is initialized
    }

    stepVectors_ = makePtr<PartitionedVector<Real>>(stepVectors);

    isInitialized_ = true;
  }

  /** How many steps are "owned" by this processor. 
      
      It may be the case that
      this processor belongs to a sub group of processors. All processors
      in that sub group will return the same number of owned steps.
    */
  int numOwnedSteps() const { return stepVectors_->numVectors(); }

  /** What is the stencil used to build this vector?
 
      This accessor is directly based on what the user intiializes the
      vector with. Its used by ROL algorithms and constraints to ensure
      the communication patterns are understood.
    */
  const std::vector<int> & stencil() { return stencil_; }


  /** Determine if an index is valid including the stencil.
      An index is valid if is in [min(stencil),max(stencil)+numOwnedSteps()-1]
   */
  bool isValidIndex(int i) const
  {
    if(0<=i && i <numOwnedSteps()) 
      return true;

    // these are "neighbor" unowned vectors (left)
    if(-leftVectors_.size() <= i && i < 0)
      return true;

    // these are "neighbor" unowned vectors (right)
    if(numOwnedSteps() <= i && i<numOwnedSteps()+rightVectors_.size())
      return true;

    return false;
  }

  /** Get a vector pointer. This range is valid from i in [min(stencil),max(stencil)+numOwnedSteps()-1]
   */
  Ptr<Vector<Real>> getVectorPtr(int i) const
  {
    assert(isValidIndex(i));

    // these are owned vectors
    if(0<=i && i<numOwnedSteps()) 
      return stepVectors_->get(i);

    // these are "neighbor" unowned vectors (left)
    if(-leftVectors_.size() <= i && i < 0) 
      return leftVectors_[leftVectors_.size()+i];

    // these are "neighbor" unowned vectors (right)
    if(numOwnedSteps() <= i && i<numOwnedSteps()+rightVectors_.size())
      return rightVectors_[i-numOwnedSteps()];

    return nullPtr;
  }

  /** \brief Exchange unknowns with neighboring processors.
      
      Using the stencil information, do global communication with time neighbor
      processors.
  */
  void boundaryExchange()
  {
    MPI_Comm timeComm = communicators_->getTimeCommunicator();
    int      myRank   = communicators_->getTimeRank();

    bool sendToRight   = stepEnd_ < steps_-1;
    bool recvFromRight = stepEnd_ < steps_-1;

    bool sendToLeft   = stepStart_ > 0;
    bool recvFromLeft = stepStart_ > 0;
    
    // send from left to right
    for(std::size_t i=0;i<stencil_.size();i++) {
      int offset = stencil_[i];
      if(offset >= 0)
        continue;

      if(sendToRight)
        vectorComm_->send(timeComm,myRank+1,*getVectorPtr(numOwnedSteps()+offset)); // this is "owned"
      
      if(recvFromLeft)
        vectorComm_->recv(timeComm,myRank-1,*getVectorPtr(offset));                  
    }

    // send from right to left
    for(std::size_t i=0;i<stencil_.size();i++) {
      int offset = stencil_[i];
      if(offset <= 0)
        continue;

      if(sendToLeft)
        vectorComm_->send(timeComm,myRank-1,*getVectorPtr(offset-1)); // this is "owned"
      
      if(recvFromRight)
        vectorComm_->recv(timeComm,myRank+1,*getVectorPtr(numOwnedSteps()+offset-1));
    }
  }

  // Vector virtual methods

  /** \brief Compute \f$y \leftarrow y + x\f$, where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector to be added to \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \mathtt{*this} + x\f$.

             ---
  */
  virtual void plus( const Vector<Real> &x )
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
  virtual void scale( const Real alpha ) 
  {
    stepVectors_->scale(alpha);
  }


  /** \brief Compute \f$ \langle y,x \rangle \f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector that forms the dot product with \f$\mathtt{*this}\f$.
             @return         The number equal to \f$\langle \mathtt{*this}, x \rangle\f$.

             ---
  */
  virtual Real dot( const Vector<Real> &x ) const
  {
    // this is probably very inefficient way to do this... oh well!
    typedef PinTVector<Real> PinTVector; 
    const PinTVector &xs = dynamic_cast<const PinTVector&>(x);

    // this won't work for Real!=double...oh well!
    Real subdot = stepVectors_->dot(*xs.stepVectors_);
    if(communicators_->getTimeRank()!=0) 
      subdot = 0.0;

    Real dot = 0;
    MPI_Allreduce(&subdot,&dot,1,MPI_DOUBLE,MPI_SUM,communicators_->getParentCommunicator());

    return dot;
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mathtt{*this}\f$.

             @return         A nonnegative number equal to the norm of \f$\mathtt{*this}\f$.

             ---
  */
  virtual Real norm() const
  {
    // this is probably very inefficient way to do this... oh well!
    return std::sqrt(this->dot(*this));
  }

  /** \brief Clone to make a new (uninitialized) vector.

             @return         A reference-counted pointer to the cloned vector.

             Provides the means of allocating temporary memory in ROL.

             ---             
  */
  virtual ROL::Ptr<Vector<Real>> clone() const
  {
    return makePtr<PinTVector<Real>>(*this);
  }

  virtual void print( std::ostream &outStream) const override
  {
    stepVectors_->print(outStream);
  }

#if 0
  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of @b x.
             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = \mathtt{*this} + \alpha x \f$.
             Uses #clone, #set, #scale and #plus for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void axpy( const Real alpha, const Vector<Real> &x );

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


  /** \brief Set \f$y \leftarrow x\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = x\f$.
             Uses #zero and #plus methods for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void set( const Vector<Real> &x );
#endif

}; // class Vector

} // namespace ROL

#endif
