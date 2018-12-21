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


#ifndef ROL_PINTHIERARCHY_HPP
#define ROL_PINTHIERARCHY_HPP

#include "ROL_TimeStamp.hpp"
#include "ROL_PinTVector.hpp"
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

  bool useRepart_;

public:

  //! Default constructor
  PinTHierarchy()
    : isInitialized_(false)
    , maxLevels_(-1)
    , useRepart_(false)
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
    , useRepart_(false)
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
    assert(level>=0); // precondition

    // base case
    if(level==0)
      return timeStamps_;

    Ptr<std::vector<TimeStamp<Real>>> higherLevel = getTimeStampsByLevel(level-1);

    // let's start easy!
    assert((higherLevel->size()-1) % 2 == 0);

    Ptr<std::vector<TimeStamp<Real>>> currentLevel 
        = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>((higherLevel->size()-1)/2+1);

    // build up the current level
    for(size_t k=0;k<currentLevel->size();k++) {
      currentLevel->at(k).t.resize(2);

      if(k==0) {
        currentLevel->at(k).t.at(0) = higherLevel->at(2*k).t.at(0);
        currentLevel->at(k).t.at(1) = higherLevel->at(2*k).t.at(1);
      }
      else {
        currentLevel->at(k).t.at(0) = higherLevel->at(2*k-1).t.at(0);
        currentLevel->at(k).t.at(1) = higherLevel->at(2*k+0).t.at(1);
      }
    }

    return currentLevel;
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
   { return communicators_[level]; }

   /** 
    * \brief Builds up communicators and time stamps for each level
    *
    * This currently doubles the size of the time step as it coarsens. This calls
    * the recursive helper function buildLevelDataStructures to do most of the work
    *
    * \param[in] level_0_ref Reference vector that is used to pull out all the communication
    *                        devices.
    */
   void buildLevels(const ROL::Ptr<ROL::PinTCommunicators> & pintComm,
                    const ROL::Ptr<const ROL::PinTVectorCommunication<Real>> & vectorComm)
   {
     // set the vector communicator
     vectorComm_ = vectorComm;

     // we allocate the communicators with null
     communicators_.resize(maxLevels_);
     communicators_[0] = pintComm;

     // we allocate the time stamps vectors with null
     stamps_.resize(maxLevels_);
     stamps_[0] = timeStamps_;

     buildLevelDataStructures(1);

     useRepart_ = true;
   }

   //! Recursive function that builds up the communicators and time stamps on coarse grids
   void buildLevelDataStructures(int level)
   {
     // protect sanity
     assert(level>0);

     // don't go too deep (base case)
     if(level>=maxLevels_)
       return;

     // this processor is no longer participating (base case)
     if(ROL::is_nullPtr(communicators_[level-1]))
       return;

     // this will subdivide the communicators by two
     communicators_[level] = communicators_[level-1]->buildCoarseCommunicators();

     // we now build the satmps to be distributed to other processors, we are halfing the size
     auto & fineStamps = *stamps_[level-1];
     std::vector<ROL::TimeStamp<Real>> sourceStamps((fineStamps.size()-1)/2+1); 

     for(size_type i=0;i<fineStamps.size();i+=2) {
       auto & stamp = sourceStamps[i/2];

       // buld a stamp skipping one step
       stamp.t.resize(2);

       if(i==0) {
         stamp.t[0] = fineStamps[i].t[0];
         stamp.t[1] = fineStamps[i].t[1];
       }
       else {
         stamp.t[0] = fineStamps[i-1].t[0];
         stamp.t[1] = fineStamps[i].t[1];
       }
     }

     // export to coarse distribution always requires a reference to the stamps, however that is only filled
     // if its on the right processor. The builtStamps states if that vector was filled.
     // ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>();
     // bool builtStamps = PinT::exportToCoarseDistribution_TimeStamps(sourceStamps,*stamps,*communicators_[level-1],0);

     // setup false step (this is to ensure the control is properly handled
     if(false) 
     {
       stamps_[level] = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>(sourceStamps);

/*
       Real ta = stamps_[level]->at(0).t.at(0);
       Real tb = stamps_[level]->at(0).t.at(1);
       Real dt =  tb-ta;
 
       // the serial constraint should never see this!
       ROL::TimeStamp<Real> stamp; 
       stamp.t.resize(2);
       stamp.t.at(0) = ta-dt;
       stamp.t.at(1) = ta;

       stamps_[level]->insert(stamps_[level]->begin(),stamp);
*/
     }

     buildLevelDataStructures(level+1);
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<Vector<Real>> allocateSimVector(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = pint_ref.communicatorsPtr();
     auto vectorComm = pint_ref.vectorCommunicationPtr();
     
     Ptr<std::vector<TimeStamp<Real>>> stamps = getTimeStampsByLevel(level);

     return buildStatePinTVector(comm,vectorComm,comm->getTimeSize()*stamps->size(),pint_ref.getVectorPtr(0)->clone());
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<Vector<Real>> allocateOptVector(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = pint_ref.communicatorsPtr();
     auto vectorComm = pint_ref.vectorCommunicationPtr();
     
     Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel(level);
     return buildControlPinTVector(comm,vectorComm,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
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
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void restrictSimVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     int numSteps = pint_input.numOwnedSteps();
     std::pair<int,int> fneRange = pint_input.ownedStepRange();
     std::pair<int,int> crsRange = pint_input.ownedStepRange();

     // make sure they line up, we are assuming we are only coarsening even indexed points
     assert(fneRange.first == 2*crsRange.first);

     int crs_i = 0;
     for(int i=0;i<numSteps;i++) {
       // only do evens
       if((i+fneRange.first) % 2 != 0) 
         continue;

       int coarseStep = crs_i+crsRange.second;

       ROL::Vector<Real> & out_u = *pint_output.getVectorPtr(2*crs_i);
       out_u.axpy(1.0,*pint_input.getVectorPtr(2*i)); 

       ROL::Vector<Real> & out_v = *pint_output.getVectorPtr(2*crs_i+1);
       out_v.axpy(1.0,*pint_input.getVectorPtr(2*i+1)); 

       crs_i++;
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
   void restrictOptVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     // handle interior
     for(int k=0;k<pint_output.numOwnedSteps();k++) {
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(k);

       out.axpy(1.0/2.0,*pint_input.getVectorPtr(2*k+0)); 
       out.axpy(1.0/2.0,*pint_input.getVectorPtr(2*k+1)); 
     }
   }

   /**
    * \brief Prolong a simulation space vector
    *
    *   X_{2*i+1} = x_i                           -1 <= i < Np      (injection variables, including virtual variable)
    *   X_{2*i+2} = (x_i + x_{i+1})/2              0 <= i < Np - 1  (averaged variables, created at fine level)
    *
    * Note: Currently if the timeRank==0, the initial condition is assumed to be zero. 
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void prolongSimVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input); // coarse vector
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     int fineSteps = pint_output.numOwnedSteps();

     // handle left exterior points
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     pint_input.boundaryExchangeLeftToRight();
     
     for(int fne_i=0;fne_i<fineSteps;fne_i++) {

       // fine indices
       int u_fne_index = 2*fne_i;   
       int v_fne_index = 2*fne_i+1;   

       if(fne_i % 2 == 0) {
         int crs_i = fne_i/2;
         int u_crs_index = 2*crs_i;
         int v_crs_index = 2*crs_i+1;

         pint_output.getVectorPtr(u_fne_index)->set(*pint_input.getVectorPtr(u_crs_index));
         pint_output.getVectorPtr(v_fne_index)->set(*pint_input.getVectorPtr(v_crs_index));
       }
       else {
         int coarse_si_0 = (fne_i-1)/2;
         int coarse_si_1 = coarse_si_0+1;
         int u_crs_index = -1; 
         int v_crs_index = -1;
 
         u_crs_index = 2*coarse_si_0;
         v_crs_index = 2*coarse_si_0+1;

         pint_output.getVectorPtr(u_fne_index)->set(*pint_input.getVectorPtr(u_crs_index));
         pint_output.getVectorPtr(v_fne_index)->set(*pint_input.getVectorPtr(v_crs_index));

         u_crs_index = 2*coarse_si_1;
         v_crs_index = 2*coarse_si_1+1;

         pint_output.getVectorPtr(u_fne_index)->axpy(1.0,*pint_input.getVectorPtr(u_crs_index));
         pint_output.getVectorPtr(v_fne_index)->axpy(1.0,*pint_input.getVectorPtr(v_crs_index));

         pint_output.getVectorPtr(u_fne_index)->scale(0.5);
         pint_output.getVectorPtr(v_fne_index)->scale(0.5);
       }
     }
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void prolongOptVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     // handle interior
     for(int k=0;k<pint_input.numOwnedSteps();k++) {

       pint_output.getVectorPtr(2*k+0)->set(*pint_input.getVectorPtr(k)); 
       pint_output.getVectorPtr(2*k+1)->set(*pint_input.getVectorPtr(k)); 
     }
   }

}; // ROL::PinTConstraint


} // namespace ROL 

#endif // ROL_PINTCONSTRAINT_HPP

