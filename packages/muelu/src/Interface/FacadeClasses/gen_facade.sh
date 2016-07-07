#!/bin/bash


facade_xml=def_facade_Simple2x2.xml

filekenner=`xsltproc gen_name.xsl ${facade_xml}`
templatestr=`xsltproc gen_template.xsl ${facade_xml}` 
echo "The filekenner is " $filekenner

# delete old facade class files
rm -f MueLu_Facade_${filekenner}_decl.hpp
rm -f MueLu_Facade_${filekenner}_def.hpp

# create new (empty) template files
cp MueLu_Facade_decl.tmpl MueLu_Facade_${filekenner}_decl.hpp
cp MueLu_Facade_def.tmpl MueLu_Facade_${filekenner}_def.hpp

# fix class names throughout all files
sedcommandstr=s/XYZNAMEXYZ/${filekenner}/g
sed -i $sedcommandstr MueLu_Facade_${filekenner}_decl.hpp
sed -i $sedcommandstr MueLu_Facade_${filekenner}_def.hpp

xsltproc gen_logic.xsl ${facade_xml} >> MueLu_Facade_${filekenner}_def.hpp
echo '    // end logical code' >> MueLu_Facade_${filekenner}_def.hpp
echo '' >> MueLu_Facade_${filekenner}_def.hpp    
echo '    // loop over all input parameters' >> MueLu_Facade_${filekenner}_def.hpp
echo '    for(Teuchos::ParameterList::ConstIterator it = inputParameters.begin(); it != inputParameters.end(); it++) {' >> MueLu_Facade_${filekenner}_def.hpp
echo '      // form replacement string' >> MueLu_Facade_${filekenner}_def.hpp
echo '      std::string par_name = inputParameters.name(it);' >> MueLu_Facade_${filekenner}_def.hpp
echo '      std::stringstream ss;' >> MueLu_Facade_${filekenner}_def.hpp
echo '      ss << "XXX" << par_name << "YYY";' >> MueLu_Facade_${filekenner}_def.hpp
echo '' >> MueLu_Facade_${filekenner}_def.hpp
echo '      // update final string with parameters' >> MueLu_Facade_${filekenner}_def.hpp
echo '      Teuchos::ParameterEntry par_entry = inputParameters.entry(it);' >> MueLu_Facade_${filekenner}_def.hpp
echo '      this->ReplaceString(finalString,' >> MueLu_Facade_${filekenner}_def.hpp
echo '              ss.str(), Teuchos::toString(par_entry.getAny()));' >> MueLu_Facade_${filekenner}_def.hpp
echo '    }' >> MueLu_Facade_${filekenner}_def.hpp
echo '' >> MueLu_Facade_${filekenner}_def.hpp
echo '    Teuchos::RCP<ParameterList> ret = Teuchos::getParametersFromXmlString(finalString);' >> MueLu_Facade_${filekenner}_def.hpp
echo '    return ret;' >> MueLu_Facade_${filekenner}_def.hpp
echo '  }' >> MueLu_Facade_${filekenner}_def.hpp
echo '' >> MueLu_Facade_${filekenner}_def.hpp
echo '  // Note all parameters are of type string (we use it for string replacement)' >> MueLu_Facade_${filekenner}_def.hpp
echo '  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>' >> MueLu_Facade_${filekenner}_def.hpp
echo '  const std::string FacadeClassBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::defaultParams_ =' >> MueLu_Facade_${filekenner}_def.hpp
echo '"<ParameterList name=\"Input\">"' >> MueLu_Facade_${filekenner}_def.hpp
echo '"<Parameter name=\"MueLu preconditioner\" type=\"string\" value=\"undefined\"/>"' >> MueLu_Facade_${filekenner}_def.hpp
xsltproc gen_defaults.xsl ${facade_xml} >> MueLu_Facade_${filekenner}_def.hpp
echo '"</ParameterList>"' >> MueLu_Facade_${filekenner}_def.hpp
echo ';' >> MueLu_Facade_${filekenner}_def.hpp

echo '  // template string for preconditioner layout (factory based parameters)' >> MueLu_Facade_${filekenner}_def.hpp
echo '  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>' >> MueLu_Facade_${filekenner}_def.hpp
echo '  const std::string FacadeClassBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stringTemplate_ =' >> MueLu_Facade_${filekenner}_def.hpp
xsltproc gen_template.xsl ${facade_xml} >> MueLu_Facade_${filekenner}_def.hpp
echo ';' >> MueLu_Facade_${filekenner}_def.hpp


echo '} // end namespace MueLu' >> MueLu_Facade_${filekenner}_def.hpp
echo '#endif'                   >> MueLu_Facade_${filekenner}_def.hpp

