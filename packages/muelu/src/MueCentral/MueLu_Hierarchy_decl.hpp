// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_HIERARCHY_DECL_HPP
#define MUELU_HIERARCHY_DECL_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Ptr.hpp>

#include <Xpetra_ConfigDefs.hpp>  // global_size_t
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>

#include <MueLu_SmootherCloner.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Hierarchy_fwd.hpp"

#include "MueLu_Types.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManager.hpp"  // no fwd declaration because constructor of FactoryManager is used as a default parameter of Setup()
#include "MueLu_KeepType.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

enum class ConvergenceStatus {
  Converged,
  Unconverged,
  Undefined
};

/*!
  @class Hierarchy
  @brief Provides methods to build a multigrid hierarchy and apply multigrid cycles.

  Allows users to manually populate operators at different levels within
  a multigrid method and push them into the hierarchy via SetLevel()
  and/or to supply factories for automatically generating prolongators,
  restrictors, and coarse level discretizations.  Additionally, this class contains
  an apply method that supports V and W cycles.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class Hierarchy : public BaseClass {
#undef MUELU_HIERARCHY_SHORT
#include "MueLu_UseShortNames.hpp"

  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType MagnitudeType;

  //! Data struct for defining stopping criteria of multigrid iteration
  struct ConvData {
    ConvData()
      : maxIts_(1)
      , tol_(-STS::magnitude(STS::one())) {}
    ConvData(LO maxIts)
      : maxIts_(maxIts)
      , tol_(-STS::magnitude(STS::one())) {}
    ConvData(MagnitudeType tol)
      : maxIts_(10000)
      , tol_(tol) {}
    ConvData(std::pair<LO, MagnitudeType> p)
      : maxIts_(p.first)
      , tol_(p.second) {}

    LO maxIts_;
    MagnitudeType tol_;
  };

 public:
  //! @name Constructors/Destructors
  //@{

  //! Default constructor.
  Hierarchy();
  //! Constructor that labels the hierarchy.
  Hierarchy(const std::string& label);

  //! Constructor
  Hierarchy(const RCP<Matrix>& A);

  //! Constructor
  Hierarchy(const RCP<Matrix>& A, const std::string& label);

  //! Destructor.
  virtual ~Hierarchy() {}

  //@}

  //! @name Set/Get Methods.
  //@{

  //!
  static CycleType GetDefaultCycle() { return MasterList::getDefault<std::string>("cycle type") == "V" ? VCYCLE : WCYCLE; }
  static int GetDefaultCycleStartLevel() { return MasterList::getDefault<int>("W cycle start level"); }
  static bool GetDefaultImplicitTranspose() { return MasterList::getDefault<bool>("transpose: use implicit"); }
  static bool GetDefaultFuseProlongationAndUpdate() { return MasterList::getDefault<bool>("fuse prolongation and update"); }
  static Xpetra::global_size_t GetDefaultMaxCoarseSize() { return MasterList::getDefault<int>("coarse: max size"); }
  static int GetDefaultMaxLevels() { return MasterList::getDefault<int>("max levels"); }
  static bool GetDefaultPRrebalance() { return MasterList::getDefault<bool>("repartition: rebalance P and R"); }

  Xpetra::global_size_t GetMaxCoarseSize() const { return maxCoarseSize_; }
  bool GetImplicitTranspose() const { return implicitTranspose_; }
  bool GetFuseProlongationAndUpdate() const { return fuseProlongationAndUpdate_; }

  void SetMaxCoarseSize(Xpetra::global_size_t maxCoarseSize) { maxCoarseSize_ = maxCoarseSize; }
  void SetPRrebalance(bool doPRrebalance) { doPRrebalance_ = doPRrebalance; }
  void SetPRViaCopyrebalance(bool doPRViaCopyrebalance) { doPRViaCopyrebalance_ = doPRViaCopyrebalance; }
  void SetImplicitTranspose(const bool& implicit) { implicitTranspose_ = implicit; }
  void SetFuseProlongationAndUpdate(const bool& fuse) { fuseProlongationAndUpdate_ = fuse; }

  //@}

  //!

  template <class S2, class LO2, class GO2, class N2>
  friend class Hierarchy;

 private:
  int LastLevelID() const { return Levels_.size() - 1; }
  void DumpCurrentGraph(int level) const;

 public:
  //! Add a level at the end of the hierarchy
  void AddLevel(const RCP<Level>& level);

  //! Add a new level at the end of the hierarchy
  void AddNewLevel();

  //! Retrieve a certain level from hierarchy.
  RCP<Level>& GetLevel(const int levelID = 0);

  int GetNumLevels() const;
  int GetGlobalNumLevels() const;

  MagnitudeType GetRate() const { return rate_; }

  // This function is global
  double GetOperatorComplexity() const;

  // This function is global
  double GetSmootherComplexity() const;

  //! Helper function
  void CheckLevel(Level& level, int levelID);

  void SetMatvecParams(RCP<ParameterList> matvecParams);

  //! Multi-level setup phase: build a new level of the hierarchy.
  /*!  This method is aimed to be used in a loop building the hierarchy level by level. See Hierarchy::Setup(manager, startLevel, numDesiredLevels) for an example of usage.
   *
   *   @param coarseLevelID ID of the level to be built.
   *   @param fineLevelManager defines how to build missing data of the fineLevel (example: aggregates)
   *   @param coarseLevelManager defines how to build the level
   *   @param nextLevelManager defines how the next coarse level will be built. This is used to post corresponding request before building the coarse level to keep useful data.

   CoarseLevel is considered to be the last level if:
    - input parameter isLastLevel == true
    or
    - Ac->getRowMap()->getGlobalNumElements() <= maxCoarseSize_
   Method return true if CoarseLevel is the last level.

   Pre-condition:
    * FineLevel:
       - must have kept useful data (TODO: not tested yet)
       - must be Teuchos::null when Setup is called for finest level (Setup then automatically calls Request for "Smoother" and "CoarseSolver")
    * CoarseLevel:
       - already allocated (using Hierarchy::AddLevel())
       - requests already posted
         (exception: for finest level (=fineLevelManager==null) requests are called within setup routine)
    * NextLevel:
       - do not need to be allocate but could (FIXME: will be deleted if lastlevel...).
       - should be null when Setup is called for last level

   Post-condition:
    * FineLevel:
       - temporary data have been used and released (this condition is not tested)
    * CoarseLevel:
       - built, requests have been used
       - if it is the last level (due to input parameter isLastLevel or getGlobalNumElements() <= maxCoarseSize_),
         then the coarse solver factory of the factory manager have been used instead of the smoother factory.
    * NextLevel:
      If input parameter isLastLevel == false:
       - have been allocated
       - requests already posted.
  */
  bool Setup(int coarseLevelID, const RCP<const FactoryManagerBase> fineLevelManager /* = Teuchos::null */, const RCP<const FactoryManagerBase> coarseLevelManager,
             const RCP<const FactoryManagerBase> nextLevelManager = Teuchos::null);

  //!
  void Setup(const FactoryManagerBase& manager = FactoryManager(), int startLevel = 0, int numDesiredLevels = GetDefaultMaxLevels());

  void SetupRe();

  //! Clear impermanent data from previous setup
  void Clear(int startLevel = 0);
  void ExpertClear();

  //! Returns multigrid cycle type (supports VCYCLE and WCYCLE)
  CycleType GetCycle() const { return Cycle_; }

  //! Supports VCYCLE and WCYCLE types.
  void SetCycle(CycleType Cycle) { Cycle_ = Cycle; }

  void SetCycleStartLevel(int cycleStart) { WCycleStartLevel_ = cycleStart; }

  //! Specify damping factor alpha such that x = x + alpha*P*c, where c is the coarse grid correction.
  void SetProlongatorScalingFactor(double scalingFactor) { scalingFactor_ = scalingFactor; }

  /*!
    @brief Apply the multigrid preconditioner.

    In theory, more general cycle types than just V- and W-cycles are possible.  However,
    the enumerated type CycleType would have to be extended.

    @param B right-hand side of linear problem
    @param X initial and final (approximate) solution of linear problem
    @param ConvData struct which stores convergence criteria (maximum number of multigrid iterations or stopping tolerance)
    @param InitialGuessIsZero Indicates whether the initial guess is zero
    @param startLevel index of starting level to build multigrid hierarchy (default = 0)
  */
  ConvergenceStatus Iterate(const MultiVector& B, MultiVector& X, ConvData conv = ConvData(),
                            bool InitialGuessIsZero = false, LO startLevel = 0);

  /*!
    @brief Print matrices in the multigrid hierarchy to file.

    @param[in] start start level
    @param[in] end   end level

    Default behavior is to print system and transfer matrices from the entire hierarchy.
    Files are named "A_0.m", "P_1.m", "R_1.m", etc, and are in matrix market coordinate format.
  */
  void Write(const LO& start = -1, const LO& end = -1, const std::string& suffix = "");

  //@}

  //! @name Permanent storage
  //@{

  //! Call Level::Keep(ename, factory) for each level of the Hierarchy.
  void Keep(const std::string& ename, const FactoryBase* factory = NoFactory::get());

  //! Call Level::Delete(ename, factory) for each level of the Hierarchy.
  void Delete(const std::string& ename, const FactoryBase* factory = NoFactory::get());

  //! Call Level::AddKeepFlag for each level of the Hierarchy.
  void AddKeepFlag(const std::string& ename, const FactoryBase* factory = NoFactory::get(), KeepType keep = MueLu::Keep);

  //! Call Level::RemoveKeepFlag for each level of the Hierarchy
  void RemoveKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keep = MueLu::All);

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  /*! @brief Print the Hierarchy with some verbosity level to a FancyOStream object.

      @param[in] out The Teuchos::FancyOstream.
      @param[in] verbLevel Controls amount of output.
  */
  void describe(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_HIGH) const;

  //! Hierarchy::print is local hierarchy function, thus the statistics can be different from global ones
  void print(std::ostream& out = std::cout, const VerbLevel verbLevel = (MueLu::Parameters | MueLu::Statistics0)) const;

  /*! Indicate whether the multigrid method is a preconditioner or a solver.

    This is used in conjunction with the verbosity level to determine whether the residuals can be printed.
  */
  void IsPreconditioner(const bool flag);

  //@}

  void EnableGraphDumping(const std::string& filename, int levelID = 1) {
    isDumpingEnabled_ = true;
    dumpLevel_        = levelID;
    dumpFile_         = filename;
  }

  void setlib(Xpetra::UnderlyingLib inlib) { lib_ = inlib; }
  Xpetra::UnderlyingLib lib() { return lib_; }

  //! force recreation of cached description_ next time description() is called:
  void ResetDescription() {
    description_ = "";
  }

  void AllocateLevelMultiVectors(int numvecs, bool forceMapCheck = false);
  void DeleteLevelMultiVectors();

 protected:
  const RCP<const FactoryManagerBase>& GetLevelManager(const int levelID) const {
    return levelManagers_[levelID];
  }

 private:
  //! Copy constructor is not implemented.
  Hierarchy(const Hierarchy& h);

  //! Decide if the residual needs to be computed
  bool IsCalculationOfResidualRequired(const LO startLevel, const ConvData& conv) const;

  /*!
  \brief Decide if the multigrid iteration is converged

  We judge convergence by comparing the current \c residualNorm
  to the user given \c convergenceTolerance and then return the
  appropriate \c ConvergenceStatus
  */
  ConvergenceStatus IsConverged(const Teuchos::Array<MagnitudeType>& residualNorm,
                                const MagnitudeType convergenceTolerance) const;

  //! Print \c residualNorm for this \c iteration to the screen
  void PrintResidualHistory(const LO iteration,
                            const Teuchos::Array<MagnitudeType>& residualNorm) const;

  //! Compute the residual norm and print it depending on the verbosity level
  ConvergenceStatus ComputeResidualAndPrintHistory(const Operator& A, const MultiVector& X,
                                                   const MultiVector& B, const LO iteration,
                                                   const LO startLevel, const ConvData& conv, MagnitudeType& previousResidualNorm);

  //! Container for Level objects
  Array<RCP<Level> > Levels_;

  //! We replace coordinates GIDs to make them consistent with matrix GIDs,
  //! even if user does not do that.  Ideally, though, we should completely
  //! remove any notion of coordinate GIDs, and deal only with LIDs, assuming
  //! that they are consistent with matrix block IDs
  void ReplaceCoordinateMap(Level& level);

  //! Minimum size of a matrix on any level. If we fall below that, we stop
  //! the coarsening
  Xpetra::global_size_t maxCoarseSize_;

  //! Potential speed up of the setup by skipping R construction, and using
  //! transpose matrix-matrix product for RAP
  bool implicitTranspose_;

  //! Potential speed up of the solve by fusing prolongation and update steps.
  //! This can lead to more iterations to round-off error accumulation.
  bool fuseProlongationAndUpdate_;

  //! Potential speed up of the setup by skipping rebalancing of P and R, and
  //! doing extra import during solve
  bool doPRrebalance_;
  bool doPRViaCopyrebalance_;  // fully explicit, needed for CombinePFactory

  //! Hierarchy may be used in a standalone mode, or as a preconditioner
  bool isPreconditioner_;

  //! V- or W-cycle
  CycleType Cycle_;

  //! Level at which to start W-cycle
  int WCycleStartLevel_;

  //! Scaling factor to be applied to coarse grid correction.
  double scalingFactor_;

  //! Epetra/Tpetra mode
  Xpetra::UnderlyingLib lib_;

  //! cache description to avoid recreating in each call to description() - use ResetDescription() to force recreation in Setup, SetupRe, etc.
  mutable std::string description_ = "";  // mutable so that we can lazily initialize in description(), which is declared const

  /*!
  @brief Graph dumping

  If enabled, we dump the graph on a specified level into a specified file
  */
  bool isDumpingEnabled_;
  // -1 = dump all levels, -2 = dump nothing
  int dumpLevel_;
  std::string dumpFile_;

  //! Convergece rate
  MagnitudeType rate_;

  //! Level managers used during the Setup
  Array<RCP<const FactoryManagerBase> > levelManagers_;

  //! Caching (Multi)Vectors used in Hierarchy::Iterate()
  int sizeOfAllocatedLevelMultiVectors_;
  Array<RCP<MultiVector> > residual_, coarseRhs_, coarseX_, coarseImport_, coarseExport_, correction_;

};  // class Hierarchy

}  // namespace MueLu

#define MUELU_HIERARCHY_SHORT
#endif  // MUELU_HIERARCHY_DECL_HPP
