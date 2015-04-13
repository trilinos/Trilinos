#ifndef PIKE_VANDERPOL_DATA_TRANSFER_1_TO_2_HPP
#define PIKE_VANDERPOL_DATA_TRANSFER_1_TO_2_HPP

#include "Pike_DataTransfer.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos { template<typename T> class Comm; }
namespace pike { class BlackBoxModelEvaluator; }

namespace pike_test {

  class VanderPolDT1To2 : public pike::DataTransfer {

  public:

    VanderPolDT1To2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		    const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		    const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2);
    
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
    Teuchos::RCP<pike::BlackBoxModelEvaluator> me2_;
    std::vector<std::string> sourceModelNames_;
    std::vector<std::string> targetModelNames_;
    pike::MultiphysicsDistributor::ApplicationIndex eq1_;
    pike::MultiphysicsDistributor::ApplicationIndex eq2_;
    pike::MultiphysicsDistributor::TransferIndex eq1ToEq2_;
  };

  /** \brief non-member ctor
      \relates VanderPolDT1To2
  */
  Teuchos::RCP<pike_test::VanderPolDT1To2> 
  vanderPolDT1To2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2);
  

}

#endif
