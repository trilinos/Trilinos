#include "TivaBuena_GUI.hpp"
#include "TivaBuena_SpecificParameterEntryValidators.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "TivaBuena_StandardDependencies.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

int main(){
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * !!!!!!!!!!!!!!!!!!!!!! NOTICE: CURRENTLY BAD INFORMATION, DON'T READ THIS EXAMPLE!!!!!!!!!!!!
   */
  /*
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              ATTENTION              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * !!!!   PLEASE VIEW THE BASIC EXAMPLE FIRST BEFORE READING THIS EXAMPLE. IT PROVIDES FUNDAMENTAL    !!!! 
   * !!!!   KNOWLEDGE THAT WILL BE VERY HELPFUL IN UNDERSTANDING THIS EXAMPLE.                          !!!!
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */ 


  /*
   * One of the most powerful features of the TivaBuena package is it's ability to understand Dependencies
   * between various parameters. In order to take advantage of this capability we cannot use
   * a standard Teuchos::ParameterList. In this case we use the DependentParameterList class from the
   * TivaBuena package. Everything is still the same as before, and we do initial setup just like we 
   * did in the basic example. Except instead of using a Teuchos::ParameterList we now use a
   * TivaBuena::DependentParameterList.
   */
  Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);


  /*
   * We'll start off by adding a few parameters.
   */
  Teuchos::Array<int> cheeseArray(10, 0);
  My_deplist->set("Cheese array stuff:", cheeseArray, "Array stuff");
  My_deplist->set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_deplist->set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  My_deplist->set("ShowPrecs", false, "Whether or not to show the preconditioners");
  Teuchos::ParameterList&
    Prec_List0 = My_deplist->sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  Prec_List0.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List0.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\nfactorization are left out of the factors.");

  /*
   * Here we create a Dependency. There are several different types of dependencies but they all follow the same
   * general outline: A dependent (which can be a parameter or parameter list) is dependent upon a dependee
   * (which is always a non-array parameter). Different dependencies usually have different requirements of 
   * their dependees and dependts. So be sure to check a dependecies documentation if you're at all unsure
   * whether or not you're using one correctly.
   *
   * Below is a bool visual dependency. What this means is the dependent's visibility to the user
   * is determined by the dependee's boolean value. Here the dependent is the ParameterList "Preconditioner".
   * The dependee is a boolean parameter called ShowPrec. We only want the Preconditioner ParameterList to show if the ShowPrecs
   * parameter is set to true, so we give the boolean value of "true" as the showIf argument.
   *
   * If we were to write out this dependency as a sentance, it would read like this:
   * Only show the Preconditioner list if the ShowPrecs parameter is set to true.
   */
  /*Teuchos::RCP<TivaBuena::BoolVisualDependency> precDep1 =
    Teuchos::RCP<TivaBuena::BoolVisualDependency>(
      new TivaBuena::BoolVisualDependency(
        "Preconditioner",
	"ShowPrecs",
	true));*/

  /*
   * Once we have created the depenency we add it to our dependent parameter list using the addDependency function.
   * Dependencies can also be removed using the removeDependency function.
   */
/*  My_deplist->addDependency(precDep1);

  My_deplist->set(
   	"Favorite Cheese",
	"American",
	"Your favorite type of cheese");
  My_deplist->set("Swiss rating", 0, "How you rate swiss on a scail of 1 to 10");*/

  /*
   * Here we are creating a StringVisualDependency. The dependent this time is a parameter called Swiss rating. We only want
   * the swiss rating parameter to show if the Favorite Cheese parameter has the value Swiss. So we make Favorite Cheese the
   * dependee, set the desired value to Swiss, and the showIf argument to true. If we were to state this dependency as a
   * sentence it would look something like this:
   * "Show the "Swiss rating" parameter when the "Favorite Cheese" parameter has the value "Swiss".
   */
   /*Teuchos::RCP<TivaBuena::StringVisualDependency> swissDep1 = 
      Teuchos::RCP<TivaBuena::StringVisualDependency>(
        new TivaBuena::StringVisualDependency(
	  "Swiss rating",
	  "Favorite Cheese",
	  "Swiss",
	  true));*/

  /*
   * We then add the dependency.
   */
   /*My_deplist->addDependency(swissDep1);

   My_deplist->set("No awesome param", true, "Whether or not the awesome parameter should be shown");
   My_deplist->set("Awesomeness", 10, "How awesome you think dep lists are.");*/

  /*
   * Some times you only want a dependent to be shown when another parameter is NOT equal to a value.
   * Here the dependent is the Awesomeness parameter and we only want it to be shown when the No awesome param
   * is equal to false. So we pass false as the showIf argument. As a sentance, the dependnecy would read
   * like this:
   *  Only show the "Awesomeness" parameter when the "No awesome param" is false.
   */
  /*Teuchos::RCP<TivaBuena::BoolVisualDependency> awesomeDep1 =
    Teuchos::RCP<TivaBuena::BoolVisualDependency>( new TivaBuena::BoolVisualDependency(
      "No awesome param",
      "Awesomeness",
      false));*/

  /*
   * Here we call the GUI. Before the user is able to input anything into the GUI, all dependencies are evaulated.
   * That way, any default dependee values will be expressed properly on their dependents.
   */
  TivaBuena::getInput(My_deplist);

  /*
   * Create the output stream.
   */
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  /*
   * Printing out the Dep List that was used to construct the GUI.
   */
  std::cout << "Dep List: \n";
  Teuchos::writeParameterListToXmlOStream(*My_deplist, *out);


  /*
   * Final notes:
   *
   * Remember: When making dependencies you must use the exact names of the the parameter and parameter lists
   * when specifying the dependent and the dependee. If this is not done properly, then the dependency will
   * have no effect on the GUI representation of the parameters.
   *
   * Remmeber: Different depenencies have different requirements. Be sure to check the documentation of the
   * dependency your using to make sure you're getting it right. Otherwise your (or worse, your users) might
   * have some nasty problems down the road.
   *
   * Remember: All depenencies are evaluted before the user can even interact with the GUI. Make sure that your
   * dependee's default values will result in desireable behavior.
   *
   * Remember: If you want to use Dependencies, you have to use the DependentParameterList from the TivaBuena
   * package, not the standard ParameterList from the Teuchos package. Also, you can't use a standard
   * ParameterList as your valid parameter list and a DependentParameterList as uer userInputList or vice
   * versa.
   *
   * Remember: Arrays can't be dependees. If you would like this functionality please contact the author
   * of this package (Kurtis Nusbaum: klnusbaum@gmail.com), because he's thinking about including it in 
   * a future release of this package.
   */
  return 0;
}
