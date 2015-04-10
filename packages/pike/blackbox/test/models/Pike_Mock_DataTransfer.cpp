#include "Pike_Mock_DataTransfer.hpp"
#include "Pike_Solver.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

namespace pike_test {
  
  MockDataTransfer::MockDataTransfer(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
					 const std::string& myName,
					 const std::vector<std::string>& sourceModelNames,
					 const std::vector<std::string>& targetModelNames)
    : comm_(comm),
      name_(myName),
      sourceModelNames_(sourceModelNames),
      targetModelNames_(targetModelNames)
  { }

  std::string MockDataTransfer::name() const
  { return name_; }
  
  bool MockDataTransfer::doTransfer(const pike::Solver& solver)
  { return true; }

  bool MockDataTransfer::transferSucceeded() const
  { return true; }

  const std::vector<std::string>& MockDataTransfer::getSourceModelNames() const
  { return sourceModelNames_; }

  const std::vector<std::string>& MockDataTransfer::getTargetModelNames() const
  { return targetModelNames_; }

  // non-member ctor
  Teuchos::RCP<pike_test::MockDataTransfer> 
  mockDataTransfer(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		   const std::string& name,
		   const std::vector<std::string>& sourceModelNames,
		   const std::vector<std::string>& targetModelNames)
  {
    return Teuchos::rcp(new pike_test::MockDataTransfer(comm,name,sourceModelNames,targetModelNames));
  }

}
