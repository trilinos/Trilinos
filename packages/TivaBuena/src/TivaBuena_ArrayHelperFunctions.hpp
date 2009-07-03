#ifndef TIVABUENA_DOESPARAMETERCONTAINARRAY_HPP_
#define TIVABUENA_DOESPARAMETERCONTAINARRAY_HPP_
#include <QStringList>
#include <QVariant>
#include "TivaBuena_Types.hpp"
#include "Teuchos_ParameterEntry.hpp"
namespace TivaBuena{

/**
 * Determines whether or not a ParameterEntry contains an array.
 *
 * @return True if the ParameterEntry contains an array, false otherwise.
 */
bool doesParameterContainArray(const Teuchos::ParameterEntry *parameter);

/**
 * Takes a string representing an array, formats it, and returns
 * a QStringList containing each value in the array.
 *
 * @param values A QString containing the values in the array.
 * @return A QStringList containing the values in the array.
 */
QStringList getValues(QString& values);

/**
 * Determines the type of array stored in a parameter.
 *
 * @param parameter The parameters whose array type is in question.
 * @return A QString containing the type of array in the parameter.
 */
QString determineArrayType(Teuchos::ParameterEntry *parameter);

/*template <class S>
Teuchos::Array<S> fromStringToArray(QString arrayString);

template <>
Teuchos::Array<std::string> fromStringToArray(QString arrayString);*/

template <class S>
Teuchos::Array<S> fromStringToArray(QString arrayString){
	arrayString = arrayString.remove("{");
	arrayString = arrayString.remove("}");
	arrayString = arrayString.remove(" ");
	QStringList tempValues = arrayString.split(",");
	QList<QVariant> values;
	for(int i = 0; i<tempValues.size(); i++){
		values.append(tempValues[i]);
	}
	Teuchos::Array<S> toReturn;
	for(int i = 0; i<values.size(); i++){
		toReturn.append(values[i].value<S>());	
	}
	return toReturn;

}

template <>
Teuchos::Array<std::string> fromStringToArray<std::string>(QString arrayString);

}
#endif /* TIVABUENA_DOESPARAMETERCONTAINARRAY_HPP_ */
