#ifndef TIVABUENA_DEPENDENCYSHEET_HPP_
#define TIVABUENA_DEPENDENCYSHEET_HPP_
#include "TivaBuena_Dependency.hpp"

namespace TivaBuena{

/**
 * A Dependency sheet keeps track of dependencies between various elements located
 * somewhere withing a "Root List". All dependencies added to a DependencySheet
 * must have dependents and dependees who are either in the Root List or one of
 * its sublists.
 *
 * Note that a DependencySheet never acts on these dependencies. It mearly keeps
 * track of them.
 */
class DependencySheet{
public:
	/**
	 * Convience typedef representing a set of dependencies.
	 */
	typedef std::set<Teuchos::RCP<TivaBuena::Dependency>, TivaBuena::Dependency::DepComp > DepSet;

	/**
	 * Convience typedef. Maps dependee parameter entries to a set of their corresponding
	 * dependencies.
	 */
	typedef Teuchos::map<const Teuchos::ParameterEntry*, DepSet > DepMap;

	/**
	 * Constructs an empty DependencySheet with no name.
	 */
	DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList);

	/**
	 * Constructs a DependencySheet.
	 *
	 * @param name Name of the Dependency Sheet.
	 */
	DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList, const std::string &name);

	/**
	 * Adds a dependency to the sheet.
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
	 * Determines whether or not a parameter is depended upon by any another
	 * parameters or parameter lists.
	 *
	 * @parameter name The paramteter to be checked for dependents.
	 * @return True if the parameter you're checking has other dependents, false otherwise.
	 */
	bool hasDependents(const Teuchos::ParameterEntry *dependee) const;

	/**
	 * Returns a set of all the dependencies associated with a particular dependee.
	 *
	 * @param dependee The parameter whose dependencies are sought. 
	 * @return A set of all dependencies associated with the dependee parameter.
	 * */
	const DepSet& getDependenciesForParameter(const Teuchos::ParameterEntry *dependee) const;


	/**
	 * Returns an iterator to the beginning of all the dependees in the sheet.
	 *
	 * @return An iterator to the beginning of all the dependees in the sheet.
	 */
	DepMap::iterator depBegin();

	/**
	 * Returns an iterator to the end of all of the dependees in the sheet.
	 *
	 * @return An iterator to the end of all of the dependees in the sheet.
	 */
	DepMap::iterator depEnd();

	/**
	 * Returns a const iterator to the beginning of all the dependees in the sheet.
	 *
	 * @return A const iterator to the beginning of all the dependees in the sheet.
	 */
	DepMap::const_iterator depBegin() const;

	/**
	 * Returns a const iterator to the end of all of the dependees in the sheet.
	 *
	 * @return A const iterator to the end of all of the dependees in the sheet.
	 */
	DepMap::const_iterator depEnd() const;

	/**
	 * Prints out a list of the dependencies in the DependencySheet
	 */
	void printDeps();

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
