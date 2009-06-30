#ifndef TIVABUENA_DOESPARAMETERCONTAINARRAY_HPP_
#define TIVABUENA_DOESPARAMETERCONTAINARRAY_HPP_
#include "Teuchos_ParameterEntry.hpp"
#include <QStringList>
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

}
#endif /* TIVABUENA_DOESPARAMETERCONTAINARRAY_HPP_ */
