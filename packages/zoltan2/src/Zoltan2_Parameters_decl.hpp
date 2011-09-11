// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_PARAMETERS_DECL_HPP_
#define _ZOLTAN2_PARAMETERS_DECL_HPP_

#include <Zoltan2_config.h>
#include <string>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>

/*! \file Zoltan2_Parameters_decl.hpp

  This file contains parameter-related declarations.
*/

// Had to redefine this type from Teuchos_ParameterEntryValidator.hpp.
// Compiler stumbled on it.
#include <Teuchos_RCP.hpp>
typedef Teuchos::RCP<const Teuchos::Array<std::string> > ValidStringsList;


namespace Zoltan2{

template <typename Integral>
  class IntegerRangeListValidator : public Teuchos::ParameterEntryValidator
{
private:
  Integral _min;
  Integral _max;

  static const std::string _listDelim;
  static const std::string _rangeDelim;
  static const std::string _allText;

  static void checkValid(char c); 
  static bool listSaysAll(std::string &l);
  static int breakRange(std::string &range, std::string &from, std::string &to);

public:
  // Constructor: any Integral is valid
  IntegerRangeListValidator();

  // Constructor: only Integrals in [validMin,validMax] are valid
  IntegerRangeListValidator(Integral validMin, Integral validMax); 

  // Implementation of Teuchos::ParameterEntryValidator interface

  const std::string getXMLTypeName() const; 

  void printDoc(std::string const& docString, std::ostream &out) const;

  //Teuchos::ValidStringsList validStringValues() const ;
  ValidStringsList validStringValues() const ;

  void validate( Teuchos::ParameterEntry  const& entry,
    std::string const& paramName, std::string const& sublistName
    ) const;

  void validateAndModify( std::string const& paramName,
    std::string const& sublistName, Teuchos::ParameterEntry * entry
    ) const;
}; // end class

// Helpers for IntegralRangeList parameter type

template <typename Integral>
  bool validIntegralRangeList(const Teuchos::Array<Integral> &vals);

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::Array<Integral> &vals);

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::ParameterEntry &e);

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::Array<Integral> &vals);

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::ParameterEntry &e);

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e);

template <typename Integral>
  void printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl);

// Function declarations

void createValidParameterList(Teuchos::ParameterList &pl);

}

#endif
