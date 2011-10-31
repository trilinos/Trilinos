#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_FactoryManager.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_HierarchyHelpers.hpp"

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

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    //! Default constructor.
    Hierarchy()
      : maxCoarseSize_(50), implicitTranspose_(false)
    ;

    //! Constructor
    Hierarchy(const RCP<Operator> & A)
      :  maxCoarseSize_(50), implicitTranspose_(false)
    ;

    //! Destructor.
    virtual ~Hierarchy() ;

    //@}

    //! @name Set/Get Methods.
    //@{

    //!
    void SetMaxCoarseSize(Xpetra::global_size_t const &maxCoarseSize) ;

    //!
    Xpetra::global_size_t GetMaxCoarseSize() const ;

  private:
    int LastLevelID() const ;

  public:

    //! Add a level at the end of the hierarchy
    void AddLevel(const RCP<Level> & level) ;

    //! Add a new level at the end of the hierarchy
    void AddNewLevel() ;


    //! Retrieve a certain level from hierarchy.
    RCP<Level> & GetLevel(const int levelID = 0) ;

    LO GetNumLevels() const ;

    //! Indicate that Iterate should use tranpose of prolongator for restriction operations.
    void SetImplicitTranspose(const bool &implicit) ;

    //! If true is returned, iterate will use tranpose of prolongator for restriction operations.
    bool GetImplicitTranspose() const ;

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
                                        const int &startLevel = 0, const int &numDesiredLevels = 10) ;

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
                                         const int startLevel = 0, const int numDesiredLevels = 10) ;
    // FillHierarchy

    Teuchos::ParameterList Setup(const FactoryManager & manager = FactoryManager(), const int &startLevel = 0, const int &numDesiredLevels = 10) ; // Setup()

    /*! @brief Set solve method for coarsest level.

    @param smooFact  fully constructed SmootherFactory 
    @param pop       whether to use pre, post, or both pre and post smoothing 

    Note: Whether the SmootherFactory builds both a pre- and post-smoother can be also be
    controlled by SmootherFactory::SetSmootherPrototypes. This approach is a bit cumbersome, 
    however.
    */
    //TODO: remove PRE/POST

    void SetCoarsestSolver(SmootherFactoryBase const &smooFact, PreOrPost const &pop = BOTH) ;

    /*! @brief Construct smoothers on all levels but the coarsest.

    Invoke a set of factories to construct smoothers within 
    a multigrid Hierarchy starting with information on 'startLevel' 
    and continuing for at most 'numDesiredLevels'. 

    Note: last level smoother will not be set here. Use SetCoarsestSolver()
    to define a smoother for the last level. Otherwise, a direct solve is
    assumed
    */
    void SetSmoothers(SmootherFactory const & smooFact, LO const & startLevel = 0, LO numDesiredLevels = -1) ; //SetSmoothers()

    // #define GimmeNorm(someVec, someLabel) ;

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
                 const bool &InitialGuessIsZero = false, const CycleType &Cycle = VCYCLE, const LO &startLevel = 0)
    ; //Iterate()

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const ;
    
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;; //class Hierarchy

} //namespace MueLu

// TODO: We need a Set/Get function to change the CycleType (for when Iterate() calls are embedded in a Belos Preconditionner for instance).

#define MUELU_HIERARCHY_SHORT
#endif //ifndef MUELU_HIERARCHY_HPP
