#include "TivaBuena_GUI.hpp"
#include "TivaBuena_SpecificParameterEntryValidators.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "TivaBuena_StandardDependencies.hpp"
#include "TivaBuena_DependencySheet.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

int main(){
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
   * between various parameters. In order to take advantage of this capability need to use the new class
   * "DependencySheet". For the most part, everything is still the same as before, and we do initial setup 
   * just like we did in the basic example. But now, we will also create our own DependencySheet.
   *
   */
  Teuchos::RCP<Teuchos::ParameterList> My_deplist = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
  /**
   * Notice how we specify My_deplist as the "Root List". When a dependency sheet is constructed it assumes
   * all dependencies that will be added to it will be located in the Root List or one of the sublists in
   * the Root List. Make sure when you're adding a dependency to the Dependency sheet that your dependent
   * and dependee (more on those later) are located somewhere in the "Root List".
   */
  Teuchos::RCP<TivaBuena::DependencySheet> depSheet1 =
  	Teuchos::RCP<TivaBuena::DependencySheet>(new TivaBuena::DependencySheet(My_deplist, "My dep sheet"));


  /*
   * We'll start off by adding a few parameters. Pretty basic stuff.
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
   * Here we create a Dependency. There are several different types of dependencies but they all follow the 
   * same general outline: A dependent (which can be a parameter or parameter list) is dependent upon a 
   * dependee (which is always a non-array parameter). Different dependencies usually have different 
   * requirements of their dependees and dependent. So be sure to check a dependecies documentation if you're 
   * at all unsure whether or not you're using one correctly. Also, if you use a dependency incorrectly, 
   * you'll be notified by an error and the GUI will never execute. It's alwasy important to make sure the 
   * GUI can at leaste run. Most errors that result from improperly formed dependencies will be caught before 
   * the user even sees the GUI.
   *
   * Below is a bool visual dependency. What this means is the dependent's visibility to the user
   * is determined by the dependee's boolean value. Here the dependent is the ParameterList "Preconditioner".
   * The dependee is a boolean parameter called ShowPrec. We only want the Preconditioner ParameterList to show 
   * if the ShowPrecs parameter is set to true, so we give the boolean value of "true" as the showIf argument.
   * We also specify the parent parameter lists of both the dependent and the dependee.
   *
   * Note how both the Preconditioner sublist and ShowPrecs parameter are located in our "Root List" (which also
   * happens to be their parent ParameterList).
   *
   * If we were to write out this dependency as a sentance, it would read like this:
   * Only show the Preconditioner list if the ShowPrecs parameter is set to true.
   */
  Teuchos::RCP<TivaBuena::BoolVisualDependency> precDep1 =
    Teuchos::RCP<TivaBuena::BoolVisualDependency>(
      new TivaBuena::BoolVisualDependency(
	"ShowPrecs",
	My_deplist,
        "Preconditioner",
	My_deplist,
	true));

  /*
   * Once we have created the depenency we add it to our dependent sheet using the addDependency function.
   * Dependencies can also be removed using the removeDependency function.
   */
  depSheet1->addDependency(precDep1);


  /*
   * Next we'll add two more parameters
   */
  My_deplist->set(
   	"Favorite Cheese",
	"American",
	"Your favorite type of cheese");
  My_deplist->set("Swiss rating", 0, "How you rate swiss on a scail of 1 to 10");

  /*
   * Here we are creating a StringVisualDependency. The dependent this time is a parameter called Swiss rating. We only want
   * the swiss rating parameter to show if the Favorite Cheese parameter has the value Swiss. So we make Favorite Cheese the
   * dependee, set the desired value to Swiss, and the showIf argument to true. If we were to state this dependency as a
   * sentence it would look something like this:
   * Show the "Swiss rating" parameter when the "Favorite Cheese" parameter has the value "Swiss".
   */
   Teuchos::RCP<TivaBuena::StringVisualDependency> swissDep1 = 
      Teuchos::RCP<TivaBuena::StringVisualDependency>(
        new TivaBuena::StringVisualDependency(
	  "Swiss rating",
	  My_deplist,
	  "Favorite Cheese",
	  My_deplist,
	  "Swiss",
	  true));

  /*
   * We then add the dependency.
   */
   depSheet1->addDependency(swissDep1);

  /*
   * Here we add two more parameters.
   */
   My_deplist->set("No awesome param", true, "Whether or not the awesome parameter should be shown");
   My_deplist->set("Awesomeness", 10, "How awesome you think dep lists are.");

  /*
   * Some times you only want a dependent to be shown when another parameter is NOT equal to a value.
   * Here the dependent is the Awesomeness parameter and we only want it to be shown when the No awesome param
   * is equal to false. So we pass false as the showIf argument. As a sentance, the dependnecy would read
   * like this:
   * Only show the "Awesomeness" parameter when the "No awesome param" is false.
   */
  Teuchos::RCP<TivaBuena::BoolVisualDependency> awesomeDep1 =
    Teuchos::RCP<TivaBuena::BoolVisualDependency>( new TivaBuena::BoolVisualDependency(
      "No awesome param",
      My_deplist,
      "Awesomeness",
      My_deplist,
      false));

  /*
   * Lets make a sublist to put in to our "Root List".
   */
  Teuchos::ParameterList& waterList = My_deplist->sublist("Water", false, "A sublist about a lovely liquid.");

  /*
   * And let's put a few parameters in the sublist.
   */
  waterList.set("Number Of Buckets", 3, "How many buckets we have");
  waterList.set("Amount in Buckets", Teuchos::Array<double>(3, 2.5), "How much water is in each bucket");

  /*
   * Obviously the number of buckets we have is going to affect the length of our array. If the number of
   * buckets gets changed to 5, we'll need an array length of 5. To solve this problem we'll use a
   * NumberArrayLengthDependency. If written in a sentence, the dependency says this:
   * The number of entry's in the "Amount in Buckets" array is dependent upon the "Number Of Buckets" parameter.
   *
   * Note how this time the parent list of the dependent and dependee are not My_deplist. They are the actuall list
   * that the parameters a located in.
   */
  Teuchos::RCP<TivaBuena::NumberArrayLengthDependency> arrayLengthDep = Teuchos::RCP<TivaBuena::NumberArrayLengthDependency>(
  	new TivaBuena::NumberArrayLengthDependency(
	"Number Of Buckets", 
	Teuchos::sublist(My_deplist,"Water"),
	"Amount in Buckets", 
	Teuchos::sublist(My_deplist,"Water")));

  /*
   * Then add the dependency to the dependency sheet.
   */
  depSheet1->addDependency(arrayLengthDep);

  /*
   * Here we call the GUI. Before the user is able to input anything into the GUI, all dependencies are evaulated.
   * That way, any default dependee values will be expressed properly on their dependents.
   */
  TivaBuena::getInput(My_deplist, depSheet1);

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
   * Remember: It's always a good idea to make sure you're GUI works. Basic issues will not allow the GUI
   * to be even displayed to the user. So if your GUI can simply launch, that means you can trust everything
   * to workout pretty well.
   *
   * Remember: Dependents and Dependees don't have to have the same parent list. They just have to be located
   * some where in the "Root List" of the dependency sheet.
   *
   * Remember: When making dependencies you must use the exact names of the the parameter and parameter lists
   * when specifying the dependent and the dependee. If you mispell the names, and there is no parameter or
   * parameter list by the name in the Root List, the GUI won't even start up. It'll just throw and error.
   * Worse, if you mispelled a name and there is a parameter or parameter list by that name in the root list,
   * then your GUI will probably behave erratically.
   *
   * Remmeber: Different depenencies have different requirements. Be sure to check the documentation of the
   * dependency your using to make sure you're getting it right. Otherwise your (or worse, your users) might
   * have some nasty problems down the road. Most of the time, if you try to preform a dependency incorrectly
   * your program will compile but the GUI will throw an error before it ever opens.
   *
   * Remember: All depenencies are evaluted before the user can ever interact with the GUI. Make sure that your
   * dependee's default values will result in desireable behavior.
   *
   * Remember: Arrays can't be dependees. If you would like this functionality please contact the author
   * of this package (Kurtis Nusbaum: klnusbaum@gmail.com), because he's thinking about including it in 
   * a future release of this package.
   */
  return 0;
}
