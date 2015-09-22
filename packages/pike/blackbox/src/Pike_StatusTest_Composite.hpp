
#ifndef PIKE_STATUS_TESTS_COMPOSITE_HPP
#define PIKE_STATUS_TESTS_COMPOSITE_HPP

#include "Pike_StatusTest.hpp"
#include <vector>

namespace Teuchos { class ParameterList; }

namespace pike {

  class StatusTestAbstractFactory;

  class Composite : 
    public pike::StatusTest {

  public:

    enum CompositeType {
      AND,
      OR
    };

    //! Default ctor.
    Composite(const pike::Composite::CompositeType type = OR);
    
    void addTest(const Teuchos::RCP<pike::StatusTest>& t);
    
    pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE);
    
    pike::SolveStatus getStatus() const;
    
    void reset();

    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=verbLevel_default) const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList,
			  const pike::StatusTestAbstractFactory& factory);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:
    void checkAnd(const pike::Solver& solver, const CheckType checkType);
    void checkOr(const pike::Solver& solver, const CheckType checkType);

  private:
    CompositeType type_;
    std::vector<Teuchos::RCP<pike::StatusTest> > tests_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    Teuchos::RCP<Teuchos::ParameterList> myParameters_;
    pike::SolveStatus status_;
    
    
    typedef std::vector<Teuchos::RCP<pike::StatusTest> >::iterator TestIterator;
    typedef std::vector<Teuchos::RCP<pike::StatusTest> >::const_iterator TestConstIterator;
  };

  //! Non-member ctor.
  Teuchos::RCP<pike::Composite> composite(const pike::Composite::CompositeType type = pike::Composite::OR);

}

#endif
