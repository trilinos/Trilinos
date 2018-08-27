#ifndef FROSCH_XPETRAOPERATOR_DECL_HPP
#define FROSCH_XPETRAOPERATOR_DECL_HPP


#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVector.hpp>

using namespace Teuchos;

namespace FROSch {
    
    /*!  @brief Wraps an existing MueLu::Hierarchy as a Xpetra::Operator.
     */
    template <class Scalar = Xpetra::Operator<>::scalar_type,
    class LocalOrdinal = typename Xpetra::Operator<Scalar>::local_ordinal_type,
    class GlobalOrdinal = typename Xpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
    class Node = typename Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
    class FROSch_XpetraOperator : public Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
        typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>                     Matrix;
        typedef FROSch::TwoLevelPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>     TwoLevelPrec;
        typedef RCP<Teuchos::ParameterList>                     ParameterListPtr;
        //protected:
        //XpetraOperator() { }
        public:
        
        //! @name Constructor/Destructor
        //@{
        
        //! Constructor
        FROSch_XpetraOperator(const RCP<TwoLevelPrec>& T) : TwoLevel_(T) { }
        
        //! Destructor.
        virtual ~FROSch_XpetraOperator() { }
        
        //@}
        
        //! Returns the Tpetra::Map object associated with the domain of this operator.
        Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
            
            return TwoLevel_->getDomainMap();
        }
        
        //! Returns the Tpetra::Map object associated with the range of this operator.
        Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
            return TwoLevel_->getRangeMap();
        }
        
        //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
        /*!
         \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
         \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
         */
        void apply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                   Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                   Scalar beta  = Teuchos::ScalarTraits<Scalar>::one()) const{
            try {
                
                
                // X is supposed to live in the range map of the operator (const rhs = B)
                RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> Xop = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(TwoLevel_->getRangeMap(),X.getNumVectors());
                RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Yop = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(TwoLevel_->getDomainMap(),Y.getNumVectors());
                TEUCHOS_TEST_FOR_EXCEPTION(TwoLevel_->getRangeMap()->isSameAs(*(Xop->getMap())) == false, std::logic_error,
                                           "FROSch::XpetraOperator::apply: map of X is incompatible with range map of TwoLevelPreconditioner");
                TEUCHOS_TEST_FOR_EXCEPTION(TwoLevel_->getDomainMap()->isSameAs(*(Yop->getMap())) == false, std::logic_error,
                                           "FROSch::XpetraOperator::apply: map of Y is incompatible with domain map of TwoLevelPreconditioner");

                
                Y.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
                TwoLevel_->apply(X,Y,NO_TRANS,Teuchos::ScalarTraits<Scalar>::one(),Teuchos::ScalarTraits<Scalar>::zero());
            } catch (std::exception& e) {
                
                //FIXME add message and rethrow
                std::cerr << "Caught an exception in FROSch::XpetraOperator::apply():" << std::endl
                << e.what() << std::endl;
            }
        }
        
        //! Indicates whether this operator supports applying the adjoint operator.
        bool hasTransposeApply() const { return false; }
        
        template <class NewNode>
        Teuchos::RCP< FROSch_XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, NewNode> >
        clone(const RCP<NewNode>& new_node) const {
            return Teuchos::rcp (new FROSch_XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, NewNode> (TwoLevel_->template clone<NewNode> (new_node)));
        }
        
        //! @name MueLu specific
        //@{
        
        //! Direct access to the underlying MueLu::Hierarchy.
        RCP<TwoLevelPrec> GetTwoLevelPreconditioner() const { return TwoLevel_; }
        
        //@}
        
        private:
        RCP<TwoLevelPrec> TwoLevel_;
    };
    
} // namespace

#endif // MUELU_XPETRAOPERATOR_DECL_HPP

