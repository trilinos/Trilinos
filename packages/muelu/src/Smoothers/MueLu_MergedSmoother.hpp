#ifndef MUELU_MERGED_SMOOTHER_HPP
#define MUELU_MERGED_SMOOTHER_HPP

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

  private:
    // List of smoothers. It is an ArrayRCP of RCP because:
    //  1) I need a vector of pointers (to avoid slicing problems) 
    //  2) I can use an std::vector insead of an ArrayRCP but then the constructor will do a copy of user input
    ArrayRCP<RCP<SmootherPrototype> > smootherList_;

    //
    bool reverseOrder_;

    // tmp, for debug
    bool verbose_;

  protected:
    RCP<Teuchos::FancyOStream> out_;

  public:
    // UNUSED: typedef vector::size_type size_type;
    //         ArrayRCP<RCP<SmootherPrototype>>::size_type

    //! @name Constructors / destructors
    //@{

    //! Constructor
    MergedSmoother(ArrayRCP<RCP<SmootherPrototype> > & smootherList, bool verbose=false)
      : smootherList_(smootherList), reverseOrder_(false), verbose_(verbose), out_(this->getOStream())
    {
      // TODO: check that on each method TEST_FOR_EXCEPTION(smootherList == Teuchos::null, MueLu::Exceptions::RuntimeError, "");

      SmootherPrototype::SetType("MergedSmoother");
      SmootherPrototype::IsSetup(false);
    }
    
    //! Copy constructor (performs a deep copy of input object)
    MergedSmoother(const MergedSmoother& src) 
      : reverseOrder_(src.reverseOrder_), verbose_(src.verbose_), out_(this->getOStream())
    {
      // Deep copy of src.smootherList_
      smootherList_ = SmootherListDeepCopy(src.GetSmootherList());
    }

    //! Copy method (performs a deep copy of input object)
    // TODO: Copy() should be virtual (if a subclass of MergedSmoother is created later) ?
    RCP<SmootherPrototype> Copy() const { 
      if (verbose_)
        *out_ << "MueLu::MergedSmoother<>::Copy()" << std::endl;

      return rcp(new MergedSmoother(*this));
    }

    //! Destructor
    virtual ~MergedSmoother() {}
    //@}

    //! @name Set/Get methods
    //@{

    void StandardOrder() { reverseOrder_ = true;  }
    void ReverseOrder()  { reverseOrder_ = false; }
 
    bool GetReverseOrder() const { return reverseOrder_; } // TODO: GetOrder() is a better name (+ enum to define order)

    // UNUSED // TODO: GetSmoother() do not take into account the reverseOrder option. Might be confusing... To be changed in MueMat too
    // const SmootherPrototype & GetSmoother(size_type Smoother) const {  return smootherList_[i]; }

    //  TODO  const ArrayRCP<const RCP<const SmootherPrototype> > & GetSmootherList() const { return smootherList_; }
    const ArrayRCP<const RCP<SmootherPrototype> > /* & */ GetSmootherList() const { return smootherList_; }

    // UNUSED size_type GetNumSmoothers() const { return smootherList_.size(); }
    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up. */
    void Setup(Level &level) {
      if (verbose_)
        *out_ << "MueLu::MergedSmoother<>::Setup() - isSetup==|" << SmootherPrototype::IsSetup() << "|" << std::endl;

      if (SmootherPrototype::IsSetup()) return; // skip setup
      
      for (typename ArrayView<RCP<SmootherPrototype> >::iterator it = smootherList_.begin(); it != smootherList_.end(); ++it) {
        try {
          (*it)->Setup(level);          
        } catch(MueLu::Exceptions::RuntimeError e) { 
          std::string msg = "MueLu::MergedSmoother<>::Setup(): Runtime Error.\n One of the underlying smoother throwed the following exception: \n"; msg += e.what();
          throw MueLu::Exceptions::RuntimeError(msg); 
        }
      }
      
      SmootherPrototype::IsSetup(true);
    }

    /*! @brief Apply

    Solves the linear system <tt>AX=B</tt> using the smoothers of the list.

    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false)
    {
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:Apply(): Setup() has not been called");

      if (verbose_)
        *out_ << "MueLu::MergedSmoother<>::Apply(): reverseOrder_==" << reverseOrder_ << std::endl;

      if (reverseOrder_ == false) {
        
        //for (vector<SmootherPrototype>::iterator it = smootherList_.begin(); it != smootherList_.end(); ++it)
        for (typename ArrayView<RCP<SmootherPrototype> >::iterator it = smootherList_.begin(); it != smootherList_.end(); ++it) {
          try {
            if (verbose_)
              *out_ << "MueLu::MergedSmoother<>::Apply() ." << std::endl;
            (*it)->Apply(X,B,InitialGuessIsZero);
          } catch(MueLu::Exceptions::RuntimeError e) { 
            std::string msg = "MueLu::MergedSmoother<>::Apply(): Runtime Error.\n One of the underlying smoothers throwed the following exception: \n"; msg += e.what();
            throw MueLu::Exceptions::RuntimeError(msg); 
          }
        }
        
      } else {
   
        // NOTE: reverse_iterators not provided by ArrayRCP :-(
        // for (...::reverse_iterator rit = smootherList_.rbegin(); rit != smootherList_.rend(); ++rit) {
        for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = smootherList_.size() - 1; i >= 0; i--)
          try {
            smootherList_[i]->Apply(X,B,InitialGuessIsZero);
          } catch(MueLu::Exceptions::RuntimeError e) { 
            std::string msg = "MueLu::MergedSmoother<>::Apply(): Runtime Error. One of the underlying smoothers throws the following exception: \n"; msg += e.what();
            throw MueLu::Exceptions::RuntimeError(msg); 
          }
      }
        
    }
      
    //@}
    
    //! @name Utilities.
    //@{
    void Print(std::string prefix) const {
      throw(Exceptions::NotImplemented("MueLu::MergedSmoother<>::Print() is not implemented"));
    }

    void CopyParameters(RCP<SmootherPrototype> src) // TODO: wrong prototype. We do not need an RCP here.
    {
      *out_ << "MueLu::MergedSmoother<>::CopyParameters(): Warning, do not work as advertised for the moment (complex reuse desactivated)" << std::endl;

      if (verbose_)
        *out_ << "MueLu::MergedSmoother<>::CopyParameters()" << std::endl;
      
      RCP<MergedSmoother> srcMergedSmoother = rcp_dynamic_cast<MergedSmoother>(src); // TODO: check if dynamic cast fails

      reverseOrder_ = srcMergedSmoother->GetReverseOrder();

      {
        //TODO        const ArrayRCP<const RCP<const SmootherPrototype> > & srcSmootherList  = srcMergedSmoother->GetSmootherList();
        const ArrayRCP<const RCP< SmootherPrototype> > & srcSmootherList  = srcMergedSmoother->GetSmootherList();
        const ArrayRCP<const RCP<SmootherPrototype> >       & thisSmootherList = smootherList_;
        TEST_FOR_EXCEPTION(srcSmootherList == Teuchos::null, MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:CopyParameters(): thisSmootherList == Teuchos::null"); // might be allowed later if useful
        TEST_FOR_EXCEPTION(srcSmootherList == Teuchos::null, MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:CopyParameters(): srcSmootherList == Teuchos::null");

        // If the smootherList of 'this' and 'src' contains the same type of smoothers, we can transfert parameters from src to 'this' in order to tentatively reuse the current setup information of each smoothers.
        // Note that the reuse of the setup phase of the MergedSmoother 'src' can be implemented for a larger set of cases (and more complicated cases), but it does not seems useful for now.

        bool reuse = true; // true == can we transfert parameters of smoothers one by one or do we have to copy the whole list of src?

        // test 1: same list size
        reuse = reuse && (thisSmootherList.size() == srcSmootherList.size());
        
        if (reuse) {

          //  test 2: one-by-one comparison of smoother types
          for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = 0; i < srcSmootherList.size(); i++) {

            // Following test should never throw in our use cases because 'src' is a prototype and 'this' is a real smoother so they don't share any data. We may allow such case later if useful.
            TEST_FOR_EXCEPTION((thisSmootherList[i] == srcSmootherList[i]) && (thisSmootherList[i] != Teuchos::null) , MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:CopyParameters(): internal logic error");

            reuse = reuse && ((thisSmootherList[i] == Teuchos::null && srcSmootherList[i] == Teuchos::null) ||
                              thisSmootherList[i]->GetType() == srcSmootherList[i]->GetType());
          }
        }

        reuse=false;//TODO: temporary disactivated.

        if (reuse) {
          bool isSetup = true;

          if (verbose_)
            *out_ << "MueLu::MergedSmoother<>::CopyParameters() -> reuse" << std::endl;
      
          // Call CopyParameters for each smoothers and update IsSetup status of the MergedSmoother
          for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = 0; i < srcSmootherList.size(); i++) {
            if (srcSmootherList[i] != Teuchos::null) {
              TEST_FOR_EXCEPTION(srcSmootherList[i] == Teuchos::null, MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:CopyParameters(): internal logic error");

              thisSmootherList[i]->CopyParameters(srcSmootherList[i]);
              isSetup = isSetup && thisSmootherList[i]->IsSetup();
            }
            SmootherPrototype::IsSetup(isSetup);
          }

        } else {

          if (verbose_)
            *out_ << "MueLu::MergedSmoother<>::CopyParameters() -> no reuse" << std::endl;

          // No reuse: copy srcSmootherList.
          smootherList_ = Teuchos::null;
          smootherList_ = SmootherListDeepCopy(srcSmootherList);
          SmootherPrototype::IsSetup(false);
        }

      }
    }
  
    ArrayRCP<RCP<SmootherPrototype> > SmootherListDeepCopy(const ArrayRCP<const RCP<SmootherPrototype> >& srcSmootherList) {
      ArrayRCP<RCP<SmootherPrototype> > newSmootherList(srcSmootherList.size());
      for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = 0; i < srcSmootherList.size(); i++) {
        newSmootherList[i] = srcSmootherList[i]->Copy();
      }
    
      return newSmootherList;
    }

    //@}

  }; //class MergedSmoother

} //namespace MueLu

#define MUELU_MERGED_SMOOTHER_SHORT

#endif //ifndef MUELU_MERGED_SMOOTHER_HPP
