#!/bin/bash




# create _decl and _def file for facade class with kenner stored in $1
function createFiles {
  facade_xml=$1
  echo "Facade xml " $facade_xml
  filekenner=`xsltproc gen_name.xsl ${facade_xml}`
  templatestr=`xsltproc gen_template.xsl ${facade_xml}` 
  echo "The filekenner is " $filekenner

  # Facade class files
  facade_decl=MueLu_Facade${filekenner}_decl.hpp
  facade_def=MueLu_Facade${filekenner}_def.hpp

  # delete old facade class files
  rm -f ${facade_decl}
  rm -f ${facade_def}

  # create new (empty) template files and fix class names throughout all files
  # Inplace use of sed (sed -i) is not supported with BSD sed.
  sedcommandstr=s/XYZNAMEXYZ/${filekenner}/g
  sed $sedcommandstr MueLu_Facade_decl.tmpl > ${facade_decl}
  sed $sedcommandstr MueLu_Facade_def.tmpl > ${facade_def}

  # Default input parameters
  echo '    // obtain ParameterList with default input parameters for this facade class' >> ${facade_def}
  echo '    // Note all parameters are of type string (we use it for string replacement)' >> ${facade_def}
  echo '    std::string defaultString = ' >> ${facade_def}
  echo '"<ParameterList name=\"Input\">"' >> ${facade_def}
  echo '"<Parameter name=\"MueLu preconditioner\" type=\"string\" value=\"undefined\"/>"' >> ${facade_def}
  xsltproc gen_defaults.xsl ${facade_xml} >> ${facade_def}
  echo '"</ParameterList>"' >> ${facade_def}
  echo ';' >> ${facade_def}
  echo '    Teuchos::RCP<ParameterList> defaultList = Teuchos::getParametersFromXmlString(defaultString);' >> ${facade_def}
  echo '    // validate user input parameters (and set defaults if necessary)' >> ${facade_def}
  echo '    Teuchos::ParameterList inputParameters = paramList;' >> ${facade_def}
  echo '    inputParameters.validateParametersAndSetDefaults(*defaultList);' >> ${facade_def}
  echo '    TEUCHOS_TEST_FOR_EXCEPTION(inputParameters.get<std::string>("MueLu preconditioner") == "undefined", MueLu::Exceptions::RuntimeError, "Facade'${filekenner}': undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");' >> ${facade_def}
  echo '' >> ${facade_def}

  # Template string
  echo '    // create copy of template string which is updated with in-place string replacements' >> ${facade_def}
  echo '    // template string for preconditioner layout (factory based parameters)' >> ${facade_def}
  echo '    std::string finalString =' >> ${facade_def}
  xsltproc gen_template.xsl ${facade_xml} >> ${facade_def}
  echo ';' >> ${facade_def}
  echo '' >> ${facade_def}

  # Logical code
  echo '    // logical code for more complicated distinctions' >> ${facade_def}
  xsltproc gen_logic.xsl ${facade_xml} >> ${facade_def}
  echo '' >> ${facade_def}
  echo '    // end logical code' >> ${facade_def}
  echo '' >> ${facade_def}

  # Loop over input parameters
  echo '    // loop over all input parameters' >> ${facade_def}
  echo '    for(Teuchos::ParameterList::ConstIterator it = inputParameters.begin(); it != inputParameters.end(); it++) {' >> ${facade_def}
  echo '      // form replacement string' >> ${facade_def}
  echo '      std::string par_name = inputParameters.name(it);' >> ${facade_def}
  echo '      std::stringstream ss;' >> ${facade_def}
  echo '      ss << "XXX" << par_name << "YYY";' >> ${facade_def}
  echo '' >> ${facade_def}

  # Update file string
  echo '      // update final string with parameters' >> ${facade_def}
  echo '      Teuchos::ParameterEntry par_entry = inputParameters.entry(it);' >> ${facade_def}
  echo '      this->ReplaceString(finalString,' >> ${facade_def}
  echo '              ss.str(), Teuchos::toString(par_entry.getAny()));' >> ${facade_def}
  echo '    }' >> ${facade_def}
  echo '' >> ${facade_def}

  echo '    Teuchos::RCP<ParameterList> ret = Teuchos::getParametersFromXmlString(finalString);' >> ${facade_def}
  echo '    return ret;' >> ${facade_def}
  echo '  }' >> ${facade_def}
  echo '' >> ${facade_def}

  echo '} // end namespace MueLu' >> ${facade_def}
  echo '#endif // PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_'${filekenner}'_DEF_HPP_' >> ${facade_def}
}

# create facade files
createFiles def_facadeSimple2x2.xml
createFiles def_facadeBGS2x2.xml
