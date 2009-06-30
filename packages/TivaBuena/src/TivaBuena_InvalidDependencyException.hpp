#ifndef TIVABUENA_INVALIDDEPENDENCYEXCEPTION_HPP_
#define TIVABUENA_INVALIDDEPENDENCYEXCEPTION_HPP_
#include <stdexcept>

namespace TivaBuena {


/**
 * Thrown when some aspect of a Dependency has been determined to be invalid.
 */
class InvalidDependencyException : public std::logic_error{
public: 
	/**
	 * Constructs an InvalidDependencyException
	 */
	InvalidDependencyException(const std::string& what_arg) 
	: std::logic_error(what_arg){}
};



}
#endif //TIVABUENA_INVALIDDEPENDENCYEXCEPTION_HPP_
