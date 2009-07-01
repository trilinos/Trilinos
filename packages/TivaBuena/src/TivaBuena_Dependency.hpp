#ifndef TIVABUENA_DEPENDENCY_HPP_
#define TIVABUENA_DEPENDENCY_HPP_
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterEntry.hpp"
namespace TivaBuena{


/**
 * This class represents a depndency between to elements in a Parameter List.
 * 
 * @see DependencyList 
 */
class Dependency{
public:
	/**
	 * Allows two dependecies to be compared.
	 */
	class DepComp{
	public:
		bool operator () (const Teuchos::RCP<Dependency> dep1, const Teuchos::RCP<Dependency> dep2) const{
			return dep1->getDependent() >= dep2->getDependent();
		}
	};

	/**
	 * Enum classifying various types of dependencies.
	 */
	enum Type{VisualDep, ValidatorDep, NumberValidatorAspectDep, NumberArrayLengthDep};

	/**
	 * Constructs a Dependency
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param type The type of dependency.
	 */
	Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, Type type);

	/**
	 * Gets the dependee of the dependency.
	 *
	 *  @return The dependee of the dependency.
	 */
	const Teuchos::ParameterEntry* getDependee() const;

	/**
	 * Gets the dependet of the dependency.
	 *
	 * @return The dependent of the dependency.
	 */
	const Teuchos::ParameterEntry* getDependent() const;

	/**
	 * Determines whether or not the potentialParentList constains the dependee's parent
	 * ParameterList.
	 *
	 * @param The ParameterList that potentially could be the parent of the
	 * dependee's parent ParameterList.
	 * @return True of the potentialParentList contians the dependee's parent
	 * ParameterList, false otherwise.
	 */
	bool isDependeeParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList);

	/**
	 * Determines whether or not the potentialParentList constains the dependent's parent
	 * ParameterList.
	 *
	 * @param The ParameterList that potentially could be the parent of the
	 * dependent's parent ParameterList.
	 * @return True of the potentialParentList contians the dependent's parent
	 * ParameterList, false otherwise.
	 */
	bool isDependentParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList);

	/**
	 * Gets the name of the dependee parameter.
	 *
	 * @return The name of the dependee parameter.
	 */
	const std::string& getDependeeName() const;

	/**
	 * Gets the name of the dependent parameter.
	 *
	 * @return The name of the dependent parameter.
	 */
	const std::string& getDependentName() const;

	/**
	 * Gets the type of the dependency.
	 *
	 * @return The type of dependency.
	 */
	Type getType() const;

	/**
	 * Evaluates the dependency and makes any appropriate changes to the
	 * dependee based on the dependent.
	 */
	virtual void evaluate() =0;

protected:
	/**
	 * The dependee is the parameter being depended upon.
	 */
	Teuchos::ParameterEntry* dependee;

	/**
	 * The dependent is the parameter that dependes on another parameter.
	 */
	Teuchos::ParameterEntry* dependent;

	/**
	 * The ParameterList containing the dependent parameter.
	 */
	Teuchos::RCP<Teuchos::ParameterList> dependentParentList;

	/**
	 * The ParameterList containing the dependee parameter.
	 */
	Teuchos::RCP<Teuchos::ParameterList> dependeeParentList;

	/**
	 * The name of the dependent and dependee parameters.
	 */
	std::string dependentName, dependeeName;

private:
	/**
	 * The type of dependency.
	 */
	Type type;

	/**
	 * Validates the dependency to make sure it's valid. Should be called in a subclasses
	 * constructor.
	 */
	virtual void validateDep() = 0;

	/**
	 * Determines whether or not a ParameterList or any of it's children lists contain a specific
	 * ParameterList.
	 *
	 * @param parentList The ParameterList to search.
	 * @param listname The name of the ParameterList for which is being searched.
	 * @return True if the parentList or and of it's children ParameterLists contains the list
	 * specified by the listname parameter.
	 */
	bool doesListContainList(Teuchos::RCP<Teuchos::ParameterList> parentList, std::string listname);
};


}
#endif //TIVABUENA_DEPENDENCY_HPP_
