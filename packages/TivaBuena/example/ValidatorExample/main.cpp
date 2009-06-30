#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "TivaBuena_GUI.hpp"
#include "TivaBuena_SpecificParameterEntryValidators.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "TivaBuena_StandardDependencies.hpp"
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
  * Every supported data type (with the exception of bool and Teuchos::Array<bool>) has 
  * at least one validator. Lets take a look at them.
  */
 
 /*
  * The validator for number types has been templated. The validator has two constuctors.
  * The first takes three arguments:
  * 1. Minimum allowed value (inclusive).
  * 2. Maximum allowed value (exclusive).
  * 3. The step. This is the how much the value of the parameter should be changed
  * when it is told to increase or decreas in the GUI. Play around with this value a bit
  * to get a good idea of what it does. It really only has meaning when used with the GUI.
  *
  * The second constructor only takes one argument:
  * 1. The step (see above for explaination).
  * If you use this second constructor, no minimum or maximum will be set. If you would 
  * like your validator to have a minimum and no maximum you may call the setMin function. 
  * The same can be done with the setMax if you wish to have a maximum and no minimum.
  */

 /*
  * First we create an empty parameter list. We will use this to define
  * all of the parameters we wish to obtain from the user.
  */
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::RCP<Teuchos::ParameterList>(new Teuchos::ParameterList);
  

 /*
  * Here we create a validator for an int parameter. It's minimum value is 0, and it's maximum value is 10.
  * The step is 1 (the default value for the argument).
  */
  Teuchos::RCP<TivaBuena::EnhancedNumberValidator<int> > intVali = 
  	Teuchos::rcp(new TivaBuena::EnhancedNumberValidator<int>(0,10));

  /*
   * We then create an int parameter and use intVali as the
   * validator argument.
   */
  My_List->set("Int", 8, "Int tester", intVali);

  /*
   * Here we create an int validator with a minimum of 0, a maximum of 100
   * and a step value of 10. Try running the example program and press the up
   * and down buttons that appear in the edit box of this parameter. That's
   * the best way to explain what the step parameter specifies.
   */
  Teuchos::RCP<TivaBuena::EnhancedNumberValidator<int> > intStepVali =
     Teuchos::rcp(new TivaBuena::EnhancedNumberValidator<int>(0,100,10));
  My_List->set("Step int", 10, "Step int tester", intStepVali);


  /*
   * Now suppose we wanted to make a validator for a short parameter that only
   * had a minimum, and no maximum. First we creat it.
   */
  Teuchos::RCP<TivaBuena::EnhancedNumberValidator<short> > shortVali = 
  	Teuchos::rcp(new TivaBuena::EnhancedNumberValidator<short>());

  /*
   * We then call the setMin function with a value of 0.
   */
  shortVali->setMin(0);

  /*
   * We then apply the validator to a short parameter.
   */
  My_List->set("Short", (short)4, "short tester", shortVali);

  /*
   * Floats and Doubles have an extra argument that can be tacked on to their constructor,
   * the precision argument. This controls how many decimals are displayed to the user in 
   * the GUI, NOT THE ACTUALL PRECISION USED TO STORE THE VALUE! Here we set the step of the
   * double validator to 1e-6 and the precision to 6 decimal places. Try running the program
   * to see it in action.
   */
  Teuchos::RCP<TivaBuena::EnhancedNumberValidator<double> > doubleVali = 
  	Teuchos::rcp(new TivaBuena::EnhancedNumberValidator<double>(0,20,1e-6, 6));
  My_List->set("Double", (double)4.5, "double tester", doubleVali);


  /*
   * If you create a string parameter and specify no validator the GUI will interpert this parameter
   * as a "free string". This means the user is allowed to type in anything they want. However, you may
   * use one of the two accepted string validators to controll what it entered by the user.
   */

  /*
   * The first validator comes from the teuchos packaged. It is called the StringToIntegralParameterEntryValidator.
   * Amoung other things, you may specify a list of acceptable values for a specific string parameter. The template 
   * arugment doesn't really matter as far as the GUI is concerened, so go ahead and use 
   * what ever your particular program calls for. Here we use int as the template argument.
   * This validator is for the Solver parameter and allows only the values GMRES, CG, and TFQMR to be supplied by the user.
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
    ,solverValidator  // This will be validated by solverValidator right here!
    );
  
  /*
   * The other validator that you may use is the FileNameValidator. This makes sure that the user supplies
   * a valid filename for this parameter.
   */
  Teuchos::RCP<TivaBuena::FileNameValidator> filnameVali = 
  	Teuchos::rcp(new TivaBuena::FileNameValidator);
  My_List->set("filename", "", "filename tester", filnameVali);

  /*
   * Array validators may also be used as well. For arrays containing numbers, simply use the ArrayNumberValidator wrapper class.
   * The ArrayNumberValidator takes an ordinary EnhancedNumberValidator as an argument for its constructor, and then uses that
   * validator to validate each entry in the array.
   *
   * If you would like to emulate the functionality of the StringToIntegralParameterEntryValidator for an array, use the
   * ArrayStringValidator wrapper class.
   *
   * If you would like to emulate the functionality of the FileNameValidator for an array, use the ArrayFileNameValidator
   * wrapper class.
   *
   * Examples of all these are shown below.
   */
  Teuchos::Array<int> intArray(10,0);
  Teuchos::Array<std::string> stringArray(10,"CG");
  Teuchos::Array<std::string> filenameArray(3,"~/");

  My_List->set("IntArray",
  	intArray,
	"intarray tester", 
	Teuchos::RCP<TivaBuena::ArrayNumberValidator<int> >(new TivaBuena::ArrayNumberValidator<int>(
	  Teuchos::RCP<TivaBuena::EnhancedNumberValidator<int> >(
	  	new TivaBuena::EnhancedNumberValidator<int>(0,20,5)	
	))));


  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    optionsValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "Option 1", "Option 2", "Option 3", "Option 4" )
        ,"Options"
        )
      );

  My_List->set("StringArray", 
  	stringArray,
	"string tester", 
  	Teuchos::RCP<TivaBuena::ArrayStringValidator<int> >(new TivaBuena::ArrayStringValidator<int>(optionsValidator))); 

  Teuchos::RCP<TivaBuena::FileNameValidator> arrayFilnameVali = 
  	Teuchos::rcp(new TivaBuena::FileNameValidator);
  
  My_List->set("Filename Array", 
  	filenameArray,
	"filename array tester",
  	Teuchos::RCP<TivaBuena::ArrayFileNameValidator>(new TivaBuena::ArrayFileNameValidator(arrayFilnameVali)));

  /*
   * Here we print ouf the user entered values in XML format.
   */
  TivaBuena::getInput(My_List);
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::writeParameterListToXmlOStream(*My_List, *out);

  /*
   * Final Notes:
   *
   * Remember: You default value for a parameter should be within the valid range
   * for the validator you are using on it. I this is not the case, an error will
   * be thrown before the GUI can even start up.
   *
   * Remember: Make sure you're using the right type of validator with the right
   * type of parameter, especially when it comes to arrays. If you don't use the
   * correct validator on the parameter, an error will be thrown before the GUI
   * can even startup.
   *
   * Careful! Reusing validators on multiple parameters is perfectly legal, but
   * be careful. Things like the NumberValidatorDependency can cause validators
   * min's and max's to change during the running of the GUI. Because we're
   * using pointers this means any parameter using the same validator will have
   * it's min's and max's changed too!
   */
  return 0;
}

