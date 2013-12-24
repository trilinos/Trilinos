
#ifndef PIKE_STATUS_TESTS_COMPOSITE_HPP
#define PIKE_STATUS_TESTS_COMPOSITE_HPP

#include "Pike_StatusTest.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <vector>

namespace pike {

  class Composite : 
    public pike::StatusTest,
    public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    enum CompositeType {
      AND,
      OR
    };

    //! Default ctor.
    Composite(const pike::Composite::CompositeType type);
    
    //! Ctor that takes 2 tests.  If you have more that 2 tests, they can be added via the addTest() method.
    Composite(const pike::Composite::CompositeType type,
	      const Teuchos::RCP<pike::StatusTest>& t1,
	      const Teuchos::RCP<pike::StatusTest>& t2);
    
    void addTest(const Teuchos::RCP<pike::StatusTest>& t);
    
    pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE);
    
    pike::SolveStatus getStatus() const;
    
    void reset();

    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=verbLevel_default) const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:
    void checkAnd(const pike::Solver& solver, const CheckType checkType);
    void checkOr(const pike::Solver& solver, const CheckType checkType);

  private:
    CompositeType type_;
    std::vector<Teuchos::RCP<pike::StatusTest> > tests_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    pike::SolveStatus status_;

    typedef std::vector<Teuchos::RCP<pike::StatusTest> >::iterator TestIterator;
    typedef std::vector<Teuchos::RCP<pike::StatusTest> >::const_iterator TestConstIterator;
  };

  //! Non-member ctor.
  Teuchos::RCP<pike::Composite> composite(const pike::Composite::CompositeType type);

}

#endif
