#ifndef MUELU_MERGEDSMOOTHER_DECL_HPP
#define MUELU_MERGEDSMOOTHER_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  class Level;

  /*!
    @class MergedSmoother
  */
template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class MergedSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  public:
    // UNUSED: typedef vector::size_type size_type;
    //         ArrayRCP<RCP<SmootherPrototype>>::size_type

    //! @name Constructors / destructors
    //@{

    //! Constructor
    MergedSmoother(ArrayRCP<RCP<SmootherPrototype> > & smootherList, bool verbose=false)
    ;
    
    //! Copy constructor (performs a deep copy of input object)
    MergedSmoother(const MergedSmoother& src) 
    ;

    //! Copy method (performs a deep copy of input object)
    // TODO: Copy() should be virtual (if a subclass of MergedSmoother is created later) ?
    RCP<SmootherPrototype> Copy() const ;

    //! Destructor
    virtual ~MergedSmoother() ;
    //@}

    //! @name Set/Get methods
    //@{

    void StandardOrder() ;
    void ReverseOrder()  ;
 
    bool GetReverseOrder() const ; // TODO: GetOrder() is a better name (+ enum to define order)

    // UNUSED // TODO: GetSmoother() do not take into account the reverseOrder option. Might be confusing... To be changed in MueMat too
    // const SmootherPrototype & GetSmoother(size_type Smoother) const ;

    //  TODO  const ArrayRCP<const RCP<const SmootherPrototype> > & GetSmootherList() const ;
    const ArrayRCP<const RCP<SmootherPrototype> > /* & */ GetSmootherList() const ;

    // UNUSED size_type GetNumSmoothers() const ;
    //@}

    void DeclareInput(Level &currentLevel) const;

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up. */
    void Setup(Level &level) ;

    /*! @brief Apply

    Solves the linear system <tt>AX=B</tt> using the smoothers of the list.

    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const
    ;
      
    //@}
    
    //! @name Utilities.
    //@{
    void Print(std::string prefix) const ;

    void CopyParameters(RCP<SmootherPrototype> src) // TODO: wrong prototype. We do not need an RCP here.
    ;
  
    ArrayRCP<RCP<SmootherPrototype> > SmootherListDeepCopy(const ArrayRCP<const RCP<SmootherPrototype> >& srcSmootherList) ;

    //@}

  private:
    // List of smoothers. It is an ArrayRCP of RCP because:
    //  1) I need a vector of pointers (to avoid slicing problems) 
    //  2) I can use an std::vector insead of an ArrayRCP but then the constructor will do a copy of user input
    ArrayRCP<RCP<SmootherPrototype> > smootherList_;

    //
    bool reverseOrder_;

    // tmp, for debug
    bool verbose_;

  }; //class MergedSmoother

} //namespace MueLu

#define MUELU_MERGED_SMOOTHER_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_MERGEDSMOOTHER_DECL_HPP
