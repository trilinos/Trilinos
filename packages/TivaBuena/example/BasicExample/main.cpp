#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "TivaBuena_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[])
{
 /*
  * Welcome to the TivaBuena Package!
  *
  * This package was designed to help rapid development of GUIs for existing and new projects
  * using the Trilinos Framework. Using the ParameterList class found in the Teuchos package
  * and the new DependentParameterList class provided by the TivaBuena Package, this package 
  * will allow you to use ParameterLists to define a set of values you whish to obtain from
  * the user. You may then pass this ParameterList along with another empty ParameterList
  * to the function getInput. This function will dynamically generate a GUI based on your
  * ParameterList, display the GUI to the user, obtain input from the user, and then store
  * the users input back into the empty ParameterList. Let's take a look at an example to see
  * how this all works.
  *
  * Before you Start:
  * We recommend you have at least a basic understanding of what a Teuchos::RCP is. While not
  * crucial to the understanding of these examples, undestanding RCPs allow you to more
  * easily understand what is going on in the examples.
  */

  /* 
   * First we create an empty parameter list. We will use this to define
   * all of the parameters we wish to obtain from the user. This type of ParameterList is
   * commonly known as the "Valid Parameter List".
   */
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);

  /* 
   * Creating parameters in this list can be easily done using the set function.
   * The first argument is the name of the Parameter. The second is the default
   * value for the parameter. The third is a short description of what the parameter is for.
   */
  My_List->set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_List->set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  /* 
   * Validators are useful for making sure that an input has only a particular set of values.
   * For the "Solver" option, we will create a validator that will automatically
   * create documentation for this parameter but will also help in validation.
   * Here we use a StringToIntegralParameterEntryValidator and a tuple
   * to specify which string values are valid for the "Solver" option.
   */
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    solverValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )
        ,"Solver"
        )
      );
  My_List->set(
    "Solver"
    ,"GMRES"
    ,"The type of solver to use."
    ,solverValidator
    );

  /* 
   * The TivaBuean Package can also handle arrays.
   * Here we create a Teuchos::Array class of 10 doubles
   * representing an initial guess for a linear solver .
   */
  Teuchos::Array<double> doubleArray( 10, 0.0 );
  
  My_List->set("Initial Guess", doubleArray, "The initial guess as a RCP to an array object.");

  /* 
   * A hierarchy of parameter lists can be constructed using {\tt Teuchos::ParameterList}.  This 
   * means another parameter list is a valid {\it value} in any parameter list.  To create a sublist
   * in a parameter list and obtain a reference to it:
   */
  Teuchos::ParameterList&
    Prec_List = My_List->sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  /*
   * Now this parameter list can be filled with values:
   */
  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.");

  /*
   * The getInput function starts up a TivaBuena GUI and lets the user start to input parameter values. When the user
   * has completed their data entry, the function will finish right after all of the input values are stored in My_List
   * userInputList.
   */
  TivaBuena::getInput(My_List);

  /*
   * Here we can print out what the user entered in nice XML format.
   */
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::writeParameterListToXmlOStream(*My_List, *out);
  

  /*
   * A Few Final Notes
   *
   * After calling the getInput function, the userInputList will be an exact copy of the valid parameter list that you
   * specified in the first arguement (with the new values the users specified of course). This includes and validators or
   * dependencies you man have used. What this means is that you can reuse
   * the userInputList again, by feeding IT in as the valid list the next time you call the getInput function. If you don't
   * care about maintaining the integerity of your valid parameter list, you may use it as both the arguments for the getInput
   * function.
   *
   * The GUI can only handle certain types of parameters. They are:
   * int
   * short
   * double
   * float
   * bool
   * std::string
   * Teuchos::Array<int>
   * Teuchos::Array<short>
   * Teuchos::Array<double>
   * Teuchos::Array<float>
   * Teuchos::Array<string>
   * If you give it a ParameterList containing a parameter that is not of one of the types specified above,
   * the parameter will still be displayed in the GUI. However, the user will not be able to modify it's value
   * and the userInputList will contain the default value specified by the valid parameter list.
   *
   * That's it for now. Be sure Check out the other examples to see some of the more advanced 
   * features of the TivaBuena package. If you have any suggestions or feature requests, please send them to
   * klnusbaum@gmail.com.
   */
  return 0;
}

