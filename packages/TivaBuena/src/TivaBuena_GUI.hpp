#ifndef TIVABUENA_GUI_HPP_
#define TIVABUENA_GUI_HPP_
#include "TivaBuena_DependencySheet.hpp"
#include "TivaBuena_metawindow.hpp"
namespace TivaBuena{
	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters);

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet);
}


#endif //TIVABUENA_GUI_HPP_
