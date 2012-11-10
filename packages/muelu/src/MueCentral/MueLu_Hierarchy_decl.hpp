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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
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
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Hierarchy_fwd.hpp"

#include "MueLu_Types.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryManager.hpp" // no fwd declaration because constructor of FactoryManager is used as a default parameter of Setup()
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase_fwd.hpp"
#include "MueLu_PFactory_fwd.hpp"
#include "MueLu_RFactory_fwd.hpp"
#include "MueLu_SmootherFactoryBase_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_HierarchyHelpers_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_KeepType.hpp"

namespace MueLu {
  /*!
    @class Hierarchy
    @brief Provides methods to build a multigrid hierarchy and apply multigrid cycles.

    Allows users to manually populate operators at different levels within
    a multigrid method and push them into the hierarchy via SetLevel()
    and/or to supply factories for automatically generating prolongators,
    restrictors, and coarse level discretizations.  Additionally, this class contains
    an apply method that supports V and W cycles.
  */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps> //TODO: or BlockSparseOp ?
  class Hierarchy : public BaseClass {
#undef MUELU_HIERARCHY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    //! Default constructor.
    Hierarchy();

    //! Constructor
    Hierarchy(const RCP<Matrix> & A);

    //! Destructor.
    virtual ~Hierarchy();

    //@}

    //! @name Set/Get Methods.
    //@{

    //!
    void SetMaxCoarseSize(Xpetra::global_size_t const &maxCoarseSize);

    //!
    Xpetra::global_size_t GetMaxCoarseSize() const;

  private:
    int LastLevelID() const;

  public:

    //! Add a level at the end of the hierarchy
    void AddLevel(const RCP<Level> & level);

    //! Add a new level at the end of the hierarchy
    void AddNewLevel();


    //! Retrieve a certain level from hierarchy.
    RCP<Level> & GetLevel(const int levelID = 0);

    LO GetNumLevels() const;

    //! Indicate that Iterate should use tranpose of prolongator for restriction operations.
    void SetImplicitTranspose(const bool &implicit);

    //! If true is returned, iterate will use tranpose of prolongator for restriction operations.
    bool GetImplicitTranspose() const;

    //@}

    //! @name Populate Methods.
    //@{

    /*!
      @brief Constructs components of the hierarchy.

      Invoke a set of factories to populate (construct prolongation,
      restriction, coarse level discretizations, and smoothers in this
      order) a multigrid Hierarchy starting with information on 'startLevel'
      and continuing for at most 'numDesiredLevels'.
    */
    Teuchos::ParameterList FullPopulate(const FactoryBase & PFact,
                                        const FactoryBase & RFact,
                                        const TwoLevelFactoryBase & AcFact,
                                        const SmootherFactory & SmooFact,
                                        const int &startLevel = 0, const int &numDesiredLevels = 10);

    /*! @brief Populate hierarchy with A's, R's, and P's.

    Invoke a set of factories to populate (construct prolongation,
    restriction, and coarse level discretizations in this
    order) a multigrid Hierarchy starting with information on 'startLevel'
    and continuing for at most 'numDesiredLevels'.

    @return  List containing starting and ending level numbers, operator complexity, \#nonzeros in the fine
    matrix, and the sum of nonzeros all matrices (including the fine).
    */
    Teuchos::ParameterList FillHierarchy(const PFactory & PFact, const RFactory & RFact,
                                         const TwoLevelFactoryBase & AcFact,
                                         const int startLevel = 0, const int numDesiredLevels = 10);
    // FillHierarchy

    //! Helper function
    void CheckLevel(Level& level, int levelID);

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
    bool Setup(int coarseLevelID, const Teuchos::Ptr<const FactoryManagerBase> fineLevelManager /* = Teuchos::null */, const Teuchos::Ptr<const FactoryManagerBase> coarseLevelManager,
               const Teuchos::Ptr<const FactoryManagerBase> nextLevelManager = Teuchos::null);

    //!
    Teuchos::ParameterList Setup(const FactoryManagerBase & manager = FactoryManager(), const int &startLevel = 0, const int &numDesiredLevels = 10); // Setup()

    /*! @brief Set solve method for coarsest level.

    @param smooFact  fully constructed SmootherFactory
    @param pop       whether to use pre, post, or both pre and post smoothing

    Note: Whether the SmootherFactory builds both a pre- and post-smoother can be also be
    controlled by SmootherFactory::SetSmootherPrototypes. This approach is a bit cumbersome,
    however.
    */
    //TODO: remove PRE/POST

    void SetCoarsestSolver(SmootherFactoryBase const &smooFact, PreOrPost const &pop = BOTH);

    /*! @brief Construct smoothers on all levels but the coarsest.

    Invoke a set of factories to construct smoothers within
    a multigrid Hierarchy starting with information on 'startLevel'
    and continuing for at most 'numDesiredLevels'.

    Note: last level smoother will not be set here. Use SetCoarsestSolver()
    to define a smoother for the last level. Otherwise, a direct solve is
    assumed
    */
    void SetSmoothers(SmootherFactory const & smooFact, LO const & startLevel = 0, LO numDesiredLevels = -1); //SetSmoothers()

    // #define GimmeNorm(someVec, someLabel);

    /*!
      @brief Apply the multigrid preconditioner.

      In theory, more general cycle types than just V- and W-cycles are possible.  However,
      the enumerated type CycleType would have to be extended.

      @param B right-hand side of linear problem
      @param nIts number of multigrid cycles to perform
      @param X initial and final (approximate) solution of linear problem
      @param InitialGuessIsZero Indicates whether the initial guess is zero
      @param Cycle Supports VCYCLE and WCYCLE types.
    */
    void Iterate(MultiVector const &B, LO nIts, MultiVector &X, //TODO: move parameter nIts and default value = 1
                 const bool &InitialGuessIsZero = false, const CycleType &Cycle = VCYCLE, const LO &startLevel = 0);

    //@}

    //! @name Permanent storage
    //@{

    //! Call Level::Keep(ename, factory) for each level of the Hierarchy.
    void Keep(const std::string & ename, const FactoryBase* factory = NoFactory::get());

    //! Call Level::Delete(ename, factory) for each level of the Hierarchy.
    void Delete(const std::string& ename, const FactoryBase* factory = NoFactory::get());

    //! Call Level::AddKeepFlag for each level of the Hierarchy.
    void AddKeepFlag(const std::string & ename, const FactoryBase* factory = NoFactory::get(), KeepType keep = MueLu::Keep);

    //! Call Level::RemoveKeepFlag for each level of the Hierarchy
    void RemoveKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep = MueLu::All);

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    /*! Indicate whether the multigrid method is a preconditioner or a solver.

      This is used in conjunction with the verbosity level to determine whether the residuals can be printed.
    */
    void IsPreconditioner(const bool flag);

    //@}

  private:
    //! Copy constructor is not implemented.
    Hierarchy(const Hierarchy &h);

    //! vector of Level objects
    Array<RCP<Level> > Levels_;

    Xpetra::global_size_t maxCoarseSize_;
    bool implicitTranspose_;
    bool isPreconditioner_;

  }; //class Hierarchy

} //namespace MueLu

#define MUELU_HIERARCHY_SHORT
#endif // MUELU_HIERARCHY_DECL_HPP
