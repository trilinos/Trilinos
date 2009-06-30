#ifndef TIVABUENA_STANDARDDEPENDCIES_HPP_
#define TIVABUENA_STANDARDDEPENDCIES_HPP_
#include "TivaBuena_InvalidDependencyException.hpp"
#include "TivaBuena_Dependency.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "TivaBuena_SpecificParameterEntryValidators.hpp"
namespace TivaBuena{

class VisualDependency : public Dependency{
public:
	VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList);

	virtual void evaluate() = 0;

	/**
	 * Determines whether or not the dependent is currently visible.
	 */
	bool isDependentVisible();

protected:

	/**
	 * Whether or not the dependent is currently visible.
	 */
	bool dependentVisible;
private:

	virtual void validateDep() = 0;
};


class ValidatorDependency : public Dependency{
public:
	ValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList);

	virtual void evaluate() = 0;

private:
	virtual void validateDep() = 0;

};

/**
 * A string visual depdencies says the following about the relationship between two elements in a Dependent Parameter List:
 * Depending on wether or not the dependee has a particular value, the dependent may or may not be displayed to the user in a GUI.
 * 
 * The dependee of a StringVisualDependency must be of type string and can't be an array. The dependent may be any type of
 * parameter or parameter list.
 */
class StringVisualDependency : public VisualDependency{
public:
	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param value The value of the depndee that affects the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to the same value as specified by the value parameter.
	 * If false, the dependent will be shown only when the depndee is set to a value other than the one specified by the value parameter.
	 */
	StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, std::string value, bool showIf);

	void evaluate();

private:
	/**
	 * The value used to deteremine the visiblity of the dependent.
	 */
	std::string value;

	/**
	 * Whether or not to show the dependent if the dependee is set to the value.
	 */
	bool showIf;

	void validateDep();
};

/**
 * A bool visual depdencies says the following about the relationship between to elements in a Dependent Parameter List:
 * Depending on wether or not the dependee is true or false, the dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a BoolVisualDependency must be of type bool and can't be an array. The dependent may be any type of parameter
 * or parameter list.
 */
class BoolVisualDependency : public VisualDependency{
public:
	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 * If false, the dependent will be shown only when the dependee is false.
	 */
	BoolVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf);

	void evaluate();

	/**
	 * Whether or not to show the dependent if the dependee is set to the value.
	 */
	bool showIf;

	void validateDep();
};

/**
 * A number visual depdencies says the following about the relationship between to elements in a Dependent Parameter List:
 * Depending on wether or not the dependee has a certain value, the dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a NumberVisualDependency must be a number type and can't be an array. The dependent may be any type of parameter
 * or parameter list.
 */
template <class S>
class NumberVisualDependency : public VisualDependency{
public:
	/**
	 * Constructs a NumberVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param func A function that takes the dependees value, does some calculations on it, and then
	 * returns a value. If this value is greater than 0, the dependent is show. If the value returned is
	 * less than or equal to zero, the dependent is not shown. If no fuction is specified, the direct
	 * value of the dependee will be used to determine the dependents visibility in a similar fashion (postive
	 * numbers causing the dependent to be displayed and 0 or negative numbers cuasing the dependent to be
	 * hidden).
	 */
	NumberVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, S (*func)(S) =0)
		:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList)
	{
		this->func = func;
		validateDep();
	}

	void evaluate(){
		int functionVal = runFunction(dependeeParentList->get<S>(dependeeName));
		functionVal > 0 ? dependentVisible = true : dependentVisible = false;
	}

private:
	/**
	 * the function used to determine the
	 * visibility of the dependent.
	 */
	S (*func)(S);

	/**
	 * run the function on the argument and returns the value of the fucntion. if no function is specified,
	 * the argument is simple returned.
	 *
	 * @param argument the value to use as an argument for the function.
	 * @return the result of running the function with the value. if no function is specified,
	 * the argument is simple returned.
	 */
	S runFunction(S argument) const{
		if(func !=0)
			return (*func)(argument);
		else
			return argument;
	}	

	void validateDep(){
		if(dependee->isType<int>()
		&& dependee->isType<short>()
		&& dependee->isType<double>()
		&& dependee->isType<float>())
			throw InvalidDependencyException("The dependee of a"
			"Number Visual Dependency must be of a supported number type!\n"
			"Problem dependee: " + dependeeName + "\n"
			"Actual type: " + dependee->getAny().typeName() + "\n"
			"Dependent: " + dependentName);
	}
};

/**
 * A NumberValidatorAspectDependency says the following about the relationship between two parameters in a dependent parameter list:
 * depending of the value of the dependee a particular aspect of the dependents validator will have a certain value.
 *
 * A NumberValidatorAspectDependency must have the following characteristics:
 * <ul>
 * 	<li>The dependee must be of a number type and can't be an array.</li>
 * 	<li>The dependent must have an enhancednumberparameter validator or an arrayNumberValidator.</li>
 * 	<li>The validator specified in the constructor must be the same validator being used on the dependent parameter.</li>
 * 	<li>The template type, dependee type, and dependent type must all be the same.</li>
 * <ul>
 *
 * This dependency can have an interesting effect on your program. The class will modifies
 * dependents validator. If that same validator is used more than once, every other 
 * parameter that used that validator will also see the change. Make sure that you are
 * aware of this when making your validators.
 */
template <class S>
class NumberValidatorAspectDependency : public Dependency{
public:
	/**
	 * An enum specifying the aspect of the
	 * validator that should be modified
	 */
	enum ValidatorAspect{
		Min,
		Max,
		Step
	};

	/**
	 * Constructs a NumberValidatorDependency
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, 
	Teuchos::RCP<TivaBuena::EnhancedNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberValidatorAspectDep)
	{
		this->validator = validator;
		this->aspect = aspect;
		this->func = func;
		validateDep();
	}

	/**
	 * Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, 
	Teuchos::RCP<TivaBuena::ArrayNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberValidatorAspectDep)
	{
		this->validator = validator->getPrototype();
		this->aspect = aspect;
		this->func = func;
		validateDep();
	}

	void evaluate(){
		S newAspectValue = runFunction(dependeeParentList->get<S>(dependeeName));
		switch(aspect){
			case NumberValidatorAspectDependency<S>::Min:
				validator->setMin(newAspectValue);
				break;
			case NumberValidatorAspectDependency<S>::Max:
				validator->setMax(newAspectValue);
				break;
			case NumberValidatorAspectDependency<S>::Step:
				validator->setStep(newAspectValue);
				break;
		}
	}

private:
	/**
	 * The aspect of the validator to be modified.
	 */
	ValidatorAspect aspect;

	/**
	 * The validator to be modified.
	 */
	Teuchos::RCP<TivaBuena::EnhancedNumberValidator<S> > validator;

	/**
	 * The function used to calculate the new value of the
	 * aspect of the validator.
	 */
	S (*func)(S);
	
	/**
	 * Runs the dependencies function on the given argument and returns
	 * the value that function returns.
	 *
	 * @param The value to run the function on.
	 * @return The value the function returned.
	 */
	S runFunction(S argument) const{
		if(func !=0)
			return (*func)(argument);
		else
			return argument;
	}	

	void validateDep(){
		if(dependee->isType<int>()
		&& dependee->isType<short>()
		&& dependee->isType<double>()
		&& dependee->isType<float>())
			throw InvalidDependencyException("The dependee of a"
			"Number Validator Aspect Dependency must be of a supported number type!\n"
			"Problem dependee: " + dependeeName + "\n"
			"Actual type: " + dependee->getAny().typeName() + "\n"
			"Dependent: " + dependentName);
		if(dependent->validator() == Teuchos::null)
			throw InvalidDependencyException("The dependent of an"
			"Number Validator Aspect Dependency must have an EnhancedNumberValidator "
			"or an ArrayNumberValidator\n"
			"Problem dependent: " + dependentName + "\n" 
			"Dependee: " + dependeeName);
		if(validator != dependent->validator())
			throw InvalidDependencyException("The dependent's validator and the validator specified "
			"in the constructor must be the same for a Number Validator Aspect Dependency!"
			"Problem dependent: " + dependentName + "\n" 
			"Dependee: " + dependeeName);
		if(typeid(S) != dependee->getAny().type() || typeid(S) != dependent->getAny().type())
			throw InvalidDependencyException("The dependent type, dependee type, and EnhancedNumberValidator "
			"template type must all be the same for a Number Validator Aspect Dependency.\n"
			"Dependee: " + dependeeName + "\n"
			"Dependee Type: " + dependee->getAny().typeName() + "\n"
			"Dependent: " + dependentName + "\n" 
			"Dependent Type: " + dependent->getAny().typeName() + "\n"
			"Validator Template Type: " + typeid(S).name());
	}	
};

/**
 * An ArrayLengthDependency says the following about the relationship between two parameters:
 * The length of the dependent's array depends on the value of the dependee.
 *
 * An ArrayLengthDependency must have the following characteristics:
 * <ul>
 * 	<li>The dependee must be either of type int or type short.</li>
 * 	<li>The dependent must be an array</li>
 * </ul>
 */
class NumberArrayLengthDependency : public Dependency{
public:
	/**
	 * Constructs an ArrayLengthDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param func A function specifying how the arrays length 
	 * should be calculated from the dependees value.
	 */
	NumberArrayLengthDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, unsigned int (*func)(int) = 0);
	
	void evaluate();

private:
	/**
	 * The function used to calculate the new value of the
	 * arrays length.
	 */
	unsigned int (*func)(int);
	
	/**
	 * Runs the dependencies function on the given argument and returns
	 * the value that function returns.
	 *
	 * @param The value to run the function on.
	 * @return The value the function returned.
	 */
	unsigned int runFunction(int argument) const;

	/**
	 * Modifies the length of an array.
	 *
	 * @param newLength The new length the array should be.
	 * @param dependentValue The index of the dependent array that is going to be changed.
	 */
	template <class S>
	void modifyArrayLength(unsigned int newLength);

	void validateDep();
};

/**
 * A StringValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A StringValidatorDependency must have the following characterisitics:
 * <ul>
 * 	<li>The dependee must be of type string</li>
 * 	<li>The dependent must already have a validator assigned to it.</li>
 * 	<li>All of the validators that may be applied to the dependent must match the type of the
 * 	validator that is currently on the dependent.</li>
 * </ul>
 */
class StringValidatorDependency : public ValidatorDependency{
public:
	/**
	 * Conveniece typedef
	 */
	typedef std::map<std::string, Teuchos::RCP<const Teuchos::ParameterEntryValidator> > ValueToValidatorMap;

	/**
	 * Constructs a StringValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param valuesAndValidators A map associating string values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that is not in the valuesAndValidators map,
	 * this is the validator that will be assigned to the dependent.
	 */
	StringValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,  ValueToValidatorMap valuesAndValidators, 
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator);

	void evaluate();

private:
	/**
	 * The default validator to be used if a request is made for a value that does not
	 * appear in the valuesAndValidators map.
	 */
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator;

	/**
	 * A map associating particular dependee values with validators that
	 * could be placed on the dependent.
	 */
	ValueToValidatorMap valuesAndValidators;
	
	void validateDep();
};

/**
 * A BoolValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A StringValidatorDependency must have the following characterisitics:
 * <ul>
 * 	<li>The dependee must be of type bool</li>
 * 	<li>The dependent must already have a validator assigned to it.</li>
 * 	<li>The "true" and "false" validators that may be applied to the dependent must match the type of the
 * 	validator that is currently on the dependent.</li>
 * </ul>
 */
class BoolValidatorDependency : public ValidatorDependency{
public:
	/**
	 * Constructs a BoolValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
	 * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
	 */
	BoolValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator);

	void evaluate();

private:
	/**
	 * The validators to be used when the dependee is either
	 * true or false.
	 */
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator, falseValidator;

	void validateDep();
};


/**
 * A RangeValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A RangeValidatorDependency achieves this by associating ranges of numerical values with validators.
 * If the dependees value falls within the one of the ranges, the validator associated with the range is
 * used on the dependent. All ranges are inclusive.
 *
 * A RangeValidatorDependency must have the following characterisitics:
 * <ul>
 * 	<li>The dependee must be a number type</li>
 * 	<li>The dependent must already have a validator assigned to it.</li>
 * 	<li>All of the validators that may be applied to the dependent must match the type of the
 * 	validator that is currently on the dependent.</li>
 * </ul>
 */
template<class S>
class RangeValidatorDependency : public ValidatorDependency{
public:
	/**
	 * Convenience typedef
	 */
	typedef std::pair<S,S> Range;
	typedef std::map<Range, Teuchos::RCP<const Teuchos::ParameterEntryValidator> > RangeToValidatorMap;

	/**
	 * Constructs a RangeValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param rangesAndValidators A map associating ranges of values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that does not fall within any of the ranges in
	 * the rangesAndValidators map, this is the validator that will be assigned to the dependent.
	 */
	RangeValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, RangeToValidatorMap rangesAndValidators,
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
		:ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList)
	{
		this->defaultValidator = defaultValidator;
		this->rangesAndValidators = rangesAndValidators;
		validateDep();
	}

	void evaluate(){
		typename RangeToValidatorMap::const_iterator it;
		S dependeeValue = dependeeParentList->get<S>(dependeeName);
		for(it = rangesAndValidators.begin(); it != rangesAndValidators.end(); it++){
			S min = it->first.first;
			S max = it->first.second;
			if(dependeeValue >= min && dependeeValue <=max){
				dependent->setValidator(it->second);
				return;
			}
		}

	}
private:
	/**
	 * The default validator
	 */
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> defaultValidator;

	/**
	 * A map associating ranges with validators.
	 */
	RangeToValidatorMap rangesAndValidators;
	
	void validateDep(){
		if(dependee->isType<int>()
		&& dependee->isType<short>()
		&& dependee->isType<double>()
		&& dependee->isType<float>())
			throw InvalidDependencyException("The dependee of a"
			"Range Validator Dependency must be of a supported number type!\n"
			"Problem dependee: " + dependeeName + "\n"
			"Actual type: " + dependee->getAny().typeName() + "\n"
			"Dependent: " + dependentName);
		typename RangeToValidatorMap::const_iterator it;
		for(it = rangesAndValidators.begin(); it != rangesAndValidators.end(); it++){
			if(typeid(dependent->validator()) != typeid(it->second))
			throw InvalidDependencyException("The validator of a dependent of a"
			"Range Validator Dependency must be the same type as all of the validators\n"
			"in the rangesAndValidators map.\n"
			"Note this means that the dependent must have an initial validator\n"
			"Problem dependent: " + dependentName + "\n"
			"Validator Type: " + typeid(dependent->validator()).name() + "\n"
			"One of the validators in the rangesAndValidators map is of type: " + typeid(it->second).name());
		}
				
	}
};


}
#endif //TIVABUENA_STANDARDDEPENDCIES_HPP_
