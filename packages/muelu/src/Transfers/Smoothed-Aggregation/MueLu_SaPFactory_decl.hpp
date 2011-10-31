#ifndef MUELU_SAPFACTORY_DECL_HPP
#define MUELU_SAPFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <iostream>

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  /*!
    @class SaPFactory class.
    @brief Factory for building Smoothed Aggregation prolongators.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SaPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{
  
    /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
    */
    SaPFactory(RCP<PFactory> InitialPFact = Teuchos::null, RCP<SingleLevelFactoryBase> AFact = Teuchos::null)
      : initialPFact_(InitialPFact), AFact_(AFact),
        dampingFactor_(4./3), diagonalView_("current") ;
  
    //! Destructor.
    virtual ~SaPFactory() ;
  
    //@}

    //! @name Set methods.
    //@{

    //! Set prolongator smoother damping factor.
    void SetDampingFactor(Scalar dampingFactor) ;

    //! Change view of diagonal.
    void SetDiagonalView(std::string const& diagView) ;
    //@}

    //! @name Get methods.
    //@{

    //! Returns prolongator smoother damping factor.
    Scalar GetDampingFactor() ;

    //! Returns current view of diagonal.
    std::string GetDiagonalView() ;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const ;

    //@}

    //! @name Build methods.
    //@{
  
    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
      //FIXME what does the return code mean (unclear in MueMat)?
      */
    void Build(Level& fineLevel, Level &coarseLevel) const ;

    void BuildP(Level &fineLevel, Level &coarseLevel) const ; //Build()

    //@}

    /*
    //TODO
    function [this] = SaPFactory(CoalesceFact,AggFact, diagonalView) //copy ctor
    function SetDiagonalView(this, diagonalView)
    */


  private:

    //! Input factories
    RCP<PFactory> initialPFact_;        //! Ptentative Factory
    RCP<SingleLevelFactoryBase> AFact_; //! A Factory
    
    //! Factory parameters
    Scalar dampingFactor_;
    std::string diagonalView_;

  }; //class SaPFactory

} //namespace MueLu

#define MUELU_SAPFACTORY_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_SAPFACTORY_DECL_HPP
