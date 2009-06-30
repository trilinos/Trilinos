#include "TivaBuena_ArrayHelperFunctions.hpp"
#include "TivaBuena_Types.hpp"

namespace TivaBuena{
bool doesParameterContainArray(const Teuchos::ParameterEntry *parameter){
	QString typeName = QString::fromStdString(parameter->getAny(false).typeName());
	return typeName.contains("Teuchos") && typeName.contains("Array");	
}

QStringList getValues(QString& values){
	values = values.remove("{");
	values = values.remove("}");
	values = values.remove(" ");
	return values.split(",");
}

QString determineArrayType(Teuchos::ParameterEntry *parameter){
	Teuchos::any anyArray = parameter->getAny();
	if(anyArray.type() == typeid(Teuchos::Array<int>))
		return intId;
	else if(anyArray.type() == typeid(Teuchos::Array<short>))
		return shortId;
	else if(anyArray.type() == typeid(Teuchos::Array<double>))
		return doubleId;
	else if(anyArray.type() == typeid(Teuchos::Array<float>))
		return floatId;
	else if(anyArray.type() == typeid(Teuchos::Array<std::string>))
		return stringId;
	else
		return unrecognizedId;		
}

}

