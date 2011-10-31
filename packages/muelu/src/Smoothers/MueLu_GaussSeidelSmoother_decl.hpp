#ifndef MUELU_GAUSSSEIDELSMOOTHER_DECL_HPP
#define MUELU_GAUSSSEIDELSMOOTHER_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <Xpetra_Operator.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels < void, LocalOrdinal, Node>::SparseOps>
  class GaussSeidelSmoother : public SmootherPrototype <Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    GaussSeidelSmoother(LO sweeps = 1, SC omega = 1.0) : nSweeps_(sweeps), omega_(omega) ;

    virtual ~GaussSeidelSmoother() ;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    //! @name Setup and apply methods.
    //@{

    //! Set up the smoother. (Right now, just grab A from the Level.)
    void Setup(Level & level) ;

    /*! Solve A*x = b approximately with the smoother.
      @param x  unknown vector
      @param b  right-hand side
      @param InitialGuessIsZero if true, indicates that x is zero, and that some flops might be avoided
    */
    void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero = false) const ; // Apply ()
    
    //@}

    //! @name Utilities.
    //@{

    RCP <SmootherPrototype> Copy() const ;

    //@}

  private:

    LO nSweeps_;       // < ! sweeps
    SC omega_;         // < ! relaxation parameter
    RCP <Operator> A_;

  }; //class GaussSeidelSmoother

} //namespace MueLu

#define MUELU_GAUSSSEIDELSMOOTHER_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_GAUSSSEIDELSMOOTHER_DECL_HPP
