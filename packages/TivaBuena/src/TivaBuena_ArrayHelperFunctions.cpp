#include "TivaBuena_ArrayHelperFunctions.hpp"

namespace TivaBuena{
bool doesParameterContainArray(const Teuchos::ParameterEntry *parameter){
	QString typeName = QString::fromStdString(parameter->getAny(false).typeName());
	return typeName.contains("Teuchos") && typeName.contains("Array");	
}

QStringList getValues(QString& values){
	values = values.remove("{");
	values = values.remove("}");
	QStringList toReturn = values.split(",");
	for(int i = 0; i < toReturn.size(); i++){
		if(toReturn[i].at(0) == QChar(' ')){
			toReturn[i] = toReturn[i].remove(0,1);
		}
	}
	return toReturn;
}

QString determineArrayType(Teuchos::ParameterEntry *parameter){
	Teuchos::any anyArray = parameter->getAny();
	if(anyArray.type() == typeid(Teuchos::Array<int>)){
		return intId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<short>)){
		return shortId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<double>)){
		return doubleId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<float>)){
		return floatId;
	}
	else if(anyArray.type() == typeid(Teuchos::Array<std::string>)){
		return stringId;
	}
	else{
		return unrecognizedId;		
	}
}

template <>
Teuchos::Array<std::string> fromStringToArray<std::string>(QString arrayString){
	arrayString = arrayString.remove("{");
	arrayString = arrayString.remove("}");
	QStringList tempValues = arrayString.split(",");
	for(int i = 0; i < tempValues.size(); i++){
		if(tempValues[i].at(0) == QChar(' ')){
			tempValues[i] = tempValues[i].remove(0,1);
		}
	}
	QList<QVariant> values;
	for(int i = 0; i<tempValues.size(); i++){
		values.append(tempValues[i]);
	}
	Teuchos::Array<std::string> toReturn;
	for(int i = 0; i<values.size(); i++){
		toReturn.append(values[i].value<QString>().toStdString());	
	}
	return toReturn;

}


}

