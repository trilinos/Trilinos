#ifndef PIKE_RXN_DATA_TRANSFER_1_TO_3_HPP
#define PIKE_RXN_DATA_TRANSFER_1_TO_3_HPP

#include "Pike_DataTransfer.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos { template<typename T> class Comm; }
namespace pike { class BlackBoxModelEvaluator; }

namespace pike_test {

  class RxnDT1To3 : public pike::DataTransfer {

  public:

    RxnDT1To3(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
	      const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
	      const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me3);
    
    //@{ DataTransfer derived methods
    
    virtual std::string name() const;

    virtual bool doTransfer(const pike::Solver& solver);

    virtual bool transferSucceeded() const;

    virtual const std::vector<std::string>& getSourceModelNames() const;
    
    virtual const std::vector<std::string>& getTargetModelNames() const;
    //@}

  private:
    Teuchos::RCP<pike::MultiphysicsDistributor> mpd_;
    Teuchos::RCP<pike::BlackBoxModelEvaluator> me1_;
    Teuchos::RCP<pike::BlackBoxModelEvaluator> me3_;
    std::vector<std::string> sourceModelNames_;
    std::vector<std::string> targetModelNames_;
    pike::MultiphysicsDistributor::ApplicationIndex eq1_;
    pike::MultiphysicsDistributor::ApplicationIndex eq3_;
    pike::MultiphysicsDistributor::TransferIndex eq1ToEq3_;
  };

  /** \brief non-member ctor
      \relates RxnDT1To3
  */
  Teuchos::RCP<pike_test::RxnDT1To3> 
  rxnDT1To3(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
	    const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
	    const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me3);
  

}

#endif
