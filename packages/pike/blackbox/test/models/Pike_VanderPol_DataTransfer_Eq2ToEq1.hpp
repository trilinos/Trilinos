#ifndef PIKE_VANDERPOL_DATA_TRANSFER_2_TO_1_HPP
#define PIKE_VANDERPOL_DATA_TRANSFER_2_TO_1_HPP

#include "Pike_DataTransfer.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos { template<typename T> class Comm; }
namespace pike { class BlackBoxModelEvaluator; }

namespace pike_test {

  class VanderPolDT2To1 : public pike::DataTransfer {

  public:

    VanderPolDT2To1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
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
    pike::MultiphysicsDistributor::TransferIndex eq2ToEq1_;
  };

  /** \brief non-member ctor
      \relates VanderPolDT2To1
  */
  Teuchos::RCP<pike_test::VanderPolDT2To1> 
  vanderPolDT2To1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2);
  

}

#endif
