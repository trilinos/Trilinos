#include "TivaBuena_StandardDependencies.hpp"
#include "TivaBuena_ArrayHelperFunctions.hpp"
#include <QDir>

namespace TivaBuena{

VisualDependency::VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList)
	:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList,  Dependency::VisualDep){}

bool VisualDependency::isDependentVisible(){
	return dependentVisible;
}

ValidatorDependency::ValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList)
	:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::ValidatorDep){}

StringVisualDependency::StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, std::string value, bool showIf)
	:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList)
{
	this->showIf = showIf;
	this->value = value;
	validateDep();
}

void StringVisualDependency::validateDep(){
	if(!dependee->isType<std::string>())
		throw InvalidDependencyException("The dependee of a "
		"String Visual Dependency must be of type string!\n"
		"Problem dependee: " + dependeeName + "\n"
		"Actual type: " + dependee->getAny().typeName() + "\n"
		"Dependent: " + dependentName);
}

void StringVisualDependency::evaluate(){
	std::string dependeeValue = dependeeParentList->get<std::string>(dependeeName);
	if((dependeeValue == value && showIf) || (dependeeValue != value && !showIf))
		dependentVisible = true;
	else 
		dependentVisible = false;
}


BoolVisualDependency::BoolVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf)
	:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList)
{
	this->showIf = showIf;
	validateDep();
}

void BoolVisualDependency::evaluate(){
	bool dependeeValue = dependeeParentList->get<bool>(dependeeName);
	if((dependeeValue && showIf) || (!dependeeValue && !showIf))
		dependentVisible = true;
	else 
		dependentVisible = false;
}

void BoolVisualDependency::validateDep(){
	if(!dependee->isType<bool>())
		throw InvalidDependencyException("The dependee of a "
		"Bool Visual Dependency must be of type bool!\n"
		"Problem dependee: " + dependeeName + "\n"
		"Actual type: " + dependee->getAny().typeName() + "\n"
		"Dependent: " + dependentName);
}

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, unsigned int (*func)(int))
	:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberArrayLengthDep)
{
	this->func = func;
	validateDep();
}

template <class S>
void NumberArrayLengthDependency::modifyArrayLength(unsigned int newLength){
	const Teuchos::Array<S> originalArray = Teuchos::any_cast<Teuchos::Array<S> >(dependent->getAny()); 
	Teuchos::RCP<const EnhancedNumberValidator<S> > potentialValidator = Teuchos::null;
	if(!dependent->validator().is_null()){
		potentialValidator = Teuchos::rcp_dynamic_cast<const ArrayNumberValidator<S> >(dependent->validator(),true)->getPrototype();
	}
	Teuchos::Array<S> newArray;
	unsigned int i;
	for(i=0; i<originalArray.size() && i<newLength; i++){
		newArray.append(originalArray[i]);
	}
	for(;i<newLength;i++){
			if(potentialValidator.is_null()){
				newArray.append(0);
			}
			else{
				newArray.append(potentialValidator->min());
			}
	}
	dependent->setValue(
		newArray,
		false,
		dependent->docString(),
		dependent->validator()
		);
}

template<>
void NumberArrayLengthDependency::modifyArrayLength<std::string>(unsigned int newLength){
	const Teuchos::Array<std::string> originalArray = Teuchos::any_cast<Teuchos::Array<std::string> >(dependent->getAny()); 
	Teuchos::Array<std::string> newArray;
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = dependent->validator();
	unsigned int i;
	for(i=0; i<originalArray.size() && i<newLength; i++)
		newArray.append(originalArray[i]);
	for(;i<newLength;i++){
		if(Teuchos::is_null(validator))
			newArray.append(" ");
		else if(!Teuchos::is_null(Teuchos::rcp_dynamic_cast<const ArrayFileNameValidator>(validator)))
			newArray.append(QDir::homePath().toStdString());
		else
			newArray.append(validator->validStringValues()->at(0));
	}
	dependent->setValue(
		newArray,
		false,
		dependent->docString(),
		dependent->validator()
		);
}

void NumberArrayLengthDependency::evaluate(){
	int newLength = runFunction(Teuchos::getValue<int>(*dependee));
	QString dependentType = determineArrayType(dependent);
	if(dependentType.contains(intId))
		modifyArrayLength<int>(newLength);
	else if(dependentType.contains(shortId))
		modifyArrayLength<short>(newLength);
	else if(dependentType.contains(doubleId))
		modifyArrayLength<double>(newLength);
	else if(dependentType.contains(floatId))
		modifyArrayLength<float>(newLength);
	else if(dependentType.contains(stringId))
		modifyArrayLength<std::string>(newLength);
		
}

unsigned int NumberArrayLengthDependency::runFunction(int argument) const{
	if(func !=0)
		return (*func)(argument);
	else
		return argument;
}

void NumberArrayLengthDependency::validateDep(){
	if(!dependee->isType<int>()
	&& !dependee->isType<short>())
		throw InvalidDependencyException("The dependee for an "
		"Array Length Dependency must be either of type short or of type int\n" 
		"Problem dependee: " + dependeeName + "\n" 
		"Actual type: " + dependee->getAny().typeName() + "\n"
		"Dependent: " + dependentName);
	if(!doesParameterContainArray(dependent))
		throw InvalidDependencyException("The dependent of an "
		"Array Length Dependency must be an array\n"
		"Problem dependent: " + dependentName + "\n" 
		"Actual type: " + dependent->getAny().typeName() + "\n"
		"Dependee: " + dependeeName);
}

StringValidatorDependency::StringValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,
	ValueToValidatorMap valuesAndValidators, Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
	:ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList)
{
	this->valuesAndValidators = valuesAndValidators;
	this->defaultValidator = defaultValidator;
	validateDep();
}

void StringValidatorDependency::evaluate(){
	std::string dependeeValue = dependeeParentList->get<std::string>(dependeeName);
	if(valuesAndValidators.find(dependeeValue) == valuesAndValidators.end())
		dependent->setValidator(defaultValidator);
	else
		dependent->setValidator(valuesAndValidators[dependeeValue]);
}

void StringValidatorDependency::validateDep(){
	if(!dependee->isType<std::string>())
		throw InvalidDependencyException("The dependee of a "
		"String Validator Dependency must be of type string!\n"
		"Problem dependee: " + dependeeName + "\n"
		"Actual type: " + dependee->getAny().typeName() + "\n"
		"Dependent: " + dependentName);
	for(ValueToValidatorMap::const_iterator it = valuesAndValidators.begin(); it != valuesAndValidators.end(); it++){
		if(typeid(dependent->validator()) != typeid(it->second))
		throw InvalidDependencyException("The validator of a dependent of a "
		"String Validator Dependency must be the same type as all of the validators\n"
		"in the valuesAndValidators map.\n"
		"Note this means that the dependent must have an initial validator\n"
		"Problem dependent: " + dependentName + "\n"
		"Validator Type: " + typeid(dependent->validator()).name() + "\n"
		"Problem Map Value: " + it->first + "\n"
		"The Validator type associated with the Map Value: " + typeid(it->second).name());
	}
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator, Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator)
	:ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList)
{
	this->trueValidator = trueValidator;
	this->falseValidator = falseValidator;
	validateDep();
}

void BoolValidatorDependency::evaluate(){
	bool dependeeValue = dependeeParentList->get<bool>(dependeeName);
	dependeeValue ? dependent->setValidator(trueValidator) : dependent->setValidator(falseValidator);
}

void BoolValidatorDependency::validateDep(){
	if(!dependee->isType<bool>())
		throw InvalidDependencyException("The dependee of a "
		"Bool Validator Dependency must be of type string!\n"
		"Problem dependee: " + dependeeName + "\n"
		"Actual type: " + dependee->getAny().typeName() + "\n"
		"Dependent: " + dependentName);
	if(typeid(dependent->validator()) != typeid(trueValidator))
		throw InvalidDependencyException("The validator of a dependent of a"
		"Bool Validator Dependency must be the same type as the \"true\" validator.\n"
		"Note this means that the dependent must have an initial validator\n"
		"Problem dependent: " + dependentName + "\n"
		"Validator Type: " + typeid(dependent->validator()).name() + "\n"
		"Type of the \"true\" validator: " + typeid(trueValidator).name());
	if(typeid(dependent->validator()) != typeid(falseValidator))
		throw InvalidDependencyException("The validator of a dependent of a"
		"Bool Validator Dependency must be the same type as the \"false\" validator.\n"
		"Note this means that the dependent must have an initial validator\n"
		"Problem dependent: " + dependentName + "\n"
		"Validator Type: " + typeid(dependent->validator()).name() + "\n"
		"Type of the \"false\" validator: " + typeid(falseValidator).name());
}


}

