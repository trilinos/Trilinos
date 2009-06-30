#ifndef TIVABUENA_DEPENDENCYSHEET_HPP_
#define TIVABUENA_DEPENDENCYSHEET_HPP_
#include "TivaBuena_Dependency.hpp"

namespace TivaBuena{

/**
 * A Parameter List that keeps track of dependencies between it's various elements.
 * Note that a DependencySheet never acts on these dependencies. It mearly keeps
 * track of them.
 */
class DependencySheet{
public:
	/**
	 * Convience typedef
	 */
	typedef std::set<Teuchos::RCP<TivaBuena::Dependency>, TivaBuena::Dependency::DepComp > DepSet;

	/**
	 * Convience typedef
	 */
	typedef Teuchos::map<const Teuchos::ParameterEntry*, DepSet > DepMap;

	/**
	 * Constructs an empty DependencySheet with no name
	 */
	DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList);

	/**
	 * Constructs a DependencySheet.
	 *
	 * @param name Name of the Parameter List
	 */
	DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList, const std::string &name);

	/**
	 * Adds a dependency to the list.
	 * 
	 * @param dependency The dependency to be added.
	 * @return True if the addition was sucessful, false otherwise.
	 */
	bool addDependency(const Teuchos::RCP<TivaBuena::Dependency> dependency);

	/**
	 * Removes a particular dependency between two parameters.
	 *
	 * @param dependency The dependency to be removed.
	 * @return True if the removal was sucessfull, false otherwise.
	 */
	bool removeDependency(Teuchos::RCP<TivaBuena::Dependency> dependency);

	/**
	 * Determines whether or not a parameter is depended upon by another parameter or
	 * parameter list.
	 *
	 * @parameter name The paramteter to be checked for dependents.
	 * @return True if the parameter you're checking has other dependents, false otherwise.
	 */
	bool hasDependents(const Teuchos::ParameterEntry *dependee) const;

	/**
	 * Returns a map of all the parameters and parameter lists that depend on the parameter specified.
	 *
	 * @param parameterName The parameter whose dependencies are in question.
	 * @return A map of all the parameters and parameter lists that depend on the specified parameter,
	 * and the associated dependencies.
	 */
	const DepSet& getDependenciesForParameter(const Teuchos::ParameterEntry *dependee) const;

	/**
	 * Prints out a list of the dependencies in the DependencySheet
	 */
	void printDeps();

	/**
	 * Returns an iterator to the beginning of all the dependees in the Dependent Parameter List.
	 *
	 * @return An iterator to the beginning of all the dependees in the Dependent Parameter List.
	 */
	DepMap::iterator depBegin();

	/**
	 * Returns an iterator to the end of all of the dependees in the Dependent Parameter List.
	 *
	 * @return An iterator to the end of all of the dependees in the Dependent Parameter List.
	 */
	DepMap::iterator depEnd();

	/**
	 * Returns a const iterator to the beginning of all the dependees in the Dependent Parameter List.
	 *
	 * @return A const iterator to the beginning of all the dependees in the Dependent Parameter List.
	 */
	DepMap::const_iterator depBegin() const;

	/**
	 * Returns a const iterator to the end of all of the dependees in the Dependent Parameter List.
	 *
	 * @return A const iterator to the end of all of the dependees in the Dependent Parameter List.
	 */
	DepMap::const_iterator depEnd() const;

private:
	/**
	 * A map containing all the depenecies for a list.
	 */
	DepMap dependencies;

	/**
	 * The Name of the dependency sheet.
	 */
	std::string name;

	/**
	 * The root parameterlist that this dependency sheet is associated with.
	 */
	Teuchos::RCP<Teuchos::ParameterList> rootList;

	/**
	 * Validates whether or not the dependee and dependet of a dependency exist
	 * within the root ParameterList.
	 */
	void validateExistanceInRoot(Teuchos::RCP<TivaBuena::Dependency> dependency);
};


}
#endif //TIVABUENA_DEPENDENCYSHEET_HPP_
