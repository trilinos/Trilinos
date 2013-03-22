// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#ifndef TEUCHOS_DEPENDENCYXMLCONVERTERDB_HPP
#define TEUCHOS_DEPENDENCYXMLCONVERTERDB_HPP

/*! \file Teuchos_DependencyXMLConverterDB.hpp
 * \brief A database for DependencyXMLConverters.
*/

//Must include this and not just DependencyXMLConverter.hpp
//because of convience macros
#include "Teuchos_StandardDependencyXMLConverters.hpp"
//Once again, done for the macros.
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_XMLParameterListReader.hpp"


namespace Teuchos {

class Dependency;

/** \brief Provides ability to lookup DependencyXMLConverterDB
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT DependencyXMLConverterDB {
public:

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param A dummy dependency representing the type of dependency the converter
   * will convert.
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(RCP<const Dependency> dependency,
    RCP<DependencyXMLConverter> converterToAdd);
  
  //@}

  /** \name Converter Functions */
  //@{
  
  /** \brief Get an appropriate DependencyXMLConverter given a 
   *  ParameterEntry.
   *
   * \param dependency The ParameterEntryDependency for which a 
   * converter is desired.
   */
  static RCP<const DependencyXMLConverter> getConverter(
    const Dependency& dependency);

  /** \brief Get an appropriate DependencyXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const DependencyXMLConverter> 
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a dependency converts the
   * dependency to XML.
   *
   * \param dependency Dependency to Convert.
   * \param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   * \param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   *
   * \return XML representation of the dependency.
   */
  static XMLObject convertDependency(
    RCP<const Dependency> dependency,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap); 

  /**
   * \brief Given an XMLObject converts the XMLObject 
   * to a Dependency.
   *
   * \param xmlObject The XMLObject to convert into a depdendency.
   * \param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   * \param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   *
   * \return A Dependency that was represented by the XML.
   */
  static RCP<Dependency> convertXML(
    const XMLObject& xmlObject, 
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap); 
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out){
    out << "Known DependencyXMLConverters: " << std::endl;
    for(
      ConverterMap::const_iterator it = getConverterMap().begin();
      it != getConverterMap().end();
      ++it)
    {
      out << "\t" << it->first <<std::endl;
    }
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief convience class. */
  typedef std::map<std::string, RCP<DependencyXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<DependencyXMLConverter> > 
    ConverterPair;

  /** \brief Gets the default converter to be used to convert
   * Dependencies.
   */
  static ConverterMap& getConverterMap();
  
  //@}

};

} // end namespace Teuchos


//
// Helper Macros
//

/** \brief Adds converter to the list of DependencyXMLConverters
 * so that all dependencies of DEP_TYPE will be converted using
 * CONVERTER.
 */
#define TEUCHOS_ADD_DEP_CONVERTER(DEP_TYPE, CONVERTER) \
   Teuchos::DependencyXMLConverterDB::addConverter( \
        Teuchos::DummyObjectGetter< DEP_TYPE >:: \
          getDummyObject(), \
        Teuchos::rcp(new CONVERTER)); 

/**
 * \brief Adds converters for NumberVisualDepednency,
 * RangeValidatorDepencny, and NumberArrayLengthDependency 
 * which are templated on type T to the list of available 
 * converters.
 */
#define TEUCHOS_ADD_TEMPLATED_NUMBER_DEPS(T) \
  TEUCHOS_ADD_NUMBER_VISUAL_DEP(T); \
  TEUCHOS_ADD_RANGE_VALIDATOR_DEP(T); \
  TEUCHOS_ADD_ARRAY_MODIFIER_DEP_GROUP(T);

/**
 * \brief Adds a NumberVisualDependencyXMLConverter temeplated on
 * type T to the list of available converters.
 */
#define TEUCHOS_ADD_NUMBER_VISUAL_DEP(T) \
   Teuchos::DependencyXMLConverterDB::addConverter( \
      Teuchos::DummyObjectGetter<Teuchos::NumberVisualDependency< T > >:: \
      getDummyObject(), \
      Teuchos::rcp(new Teuchos::NumberVisualDependencyXMLConverter< T >));

/**
 * \brief Adds a RangeValidatorDependencyXMLConverter temeplated on
 * type T to the list of available converters.
 */
#define TEUCHOS_ADD_RANGE_VALIDATOR_DEP(T) \
   Teuchos::DependencyXMLConverterDB::addConverter( \
      Teuchos::DummyObjectGetter<Teuchos::RangeValidatorDependency< T > >:: \
        getDummyObject(), \
      Teuchos::rcp(new Teuchos::RangeValidatorDependencyXMLConverter< T >));
/**
 * \brief Adds a NumberArrayLengthDependencyXMLConverter tmeplated on
 * type DEPENDEE_TYPE and DEPENDENT_TYPE to the list of available converters.
 */
#define TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP(DEPENDEE_TYPE , DEPENDENT_TYPE) \
   Teuchos::DependencyXMLConverterDB::addConverter( \
      Teuchos::DummyObjectGetter<Teuchos::NumberArrayLengthDependency< \
        DEPENDEE_TYPE , DEPENDENT_TYPE > >::getDummyObject(), \
        Teuchos::rcp(new Teuchos::NumberArrayLengthDependencyXMLConverter< \
        DEPENDEE_TYPE , DEPENDENT_TYPE >));

/**
 * \brief Adds a TwoDRowDependencyXMLConverter tmeplated on
 * type DEPENDEE_TYPE and DEPENDENT_TYPE to the list of available converters.
 */
#define TEUCHOS_ADD_TWODROW_DEP(DEPENDEE_TYPE , DEPENDENT_TYPE) \
   Teuchos::DependencyXMLConverterDB::addConverter( \
      Teuchos::DummyObjectGetter<Teuchos::TwoDRowDependency< \
        DEPENDEE_TYPE , DEPENDENT_TYPE > >::getDummyObject(), \
        Teuchos::rcp(new Teuchos::TwoDRowDependencyXMLConverter< \
        DEPENDEE_TYPE , DEPENDENT_TYPE >));
/**
 * \brief Adds a TwoDColDependencyXMLConverter tmeplated on
 * type DEPENDEE_TYPE and DEPENDENT_TYPE to the list of available converters.
 */
#define TEUCHOS_ADD_TWODCOL_DEP(DEPENDEE_TYPE , DEPENDENT_TYPE) \
   Teuchos::DependencyXMLConverterDB::addConverter( \
      Teuchos::DummyObjectGetter<Teuchos::TwoDColDependency< \
        DEPENDEE_TYPE , DEPENDENT_TYPE > >::getDummyObject(), \
        Teuchos::rcp(new Teuchos::TwoDColDependencyXMLConverter< \
        DEPENDEE_TYPE , DEPENDENT_TYPE >));

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
/**
 * \brief Adds several ArrayModifierDependencies templated on
 * DEPENDEE_TYPE and several standard dependent types.
 */
#define TEUCHOS_ADD_ARRAY_MODIFIER_DEP_GROUP(DEPENDEE_TYPE) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , std::string) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , int) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , long long int) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , double) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , float) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , std::string) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , int) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , long long int) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , double) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , float)  \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , std::string) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , int) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , long long int) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , double) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , float) 
#else
/**
 * \brief Adds several ArrayModifierDependencies templated on
 * DEPENDEE_TYPE and several standard dependent types.
 */
#define TEUCHOS_ADD_ARRAY_MODIFIER_DEP_GROUP(DEPENDEE_TYPE) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , std::string) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , int) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , double) \
  TEUCHOS_ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , float) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , std::string) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , int) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , double) \
  TEUCHOS_ADD_TWODROW_DEP( DEPENDEE_TYPE , float)  \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , std::string) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , int) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , double) \
  TEUCHOS_ADD_TWODCOL_DEP( DEPENDEE_TYPE , float) 
#endif


#endif // TEUCHOS_DEPENDENCYXMLCONVERTERDB_HPP
