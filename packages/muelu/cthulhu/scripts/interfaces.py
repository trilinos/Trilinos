import sys
sys.path.append('lib/')

import os
from string import Template
from ConfigParser import SafeConfigParser
from XpetraLib import *


def buildFuncLineInterface( functionNode ):

    # <name> = function name
    name = functionNode.xpath('name')[0].text

    # <type> = return type of the function
    type = functionNode.xpath('type')[0].xpath("string()")
    
    # <argsstring>
    argsstring = functionNode.xpath('argsstring')[0].text
    
    # briefdescription
    briefdescription = functionNode.xpath("briefdescription")[0].xpath("string()")

    #hack for Vector
    if 'magnitudeType' in type: type = 'typename ' + type
    if functionNode.xpath('//compoundname')[0].text == 'Tpetra::Vector':
        if name in ['dot','norm1','norm2','normInf','normWeighted','meanValue'] and 'ArrayView' in argsstring: return ''

    #
    if len(type) > 0 :
        declStr = type + " " + name + argsstring
    else:
        declStr = name + argsstring
#    declStr = declStr.rstrip()

    if 'TPETRA_DEPRECATED' in type: return ''
    if "const =0" in argsstring: return '' #hack for CrsMatrix

    # hack for MultiVector
    if name == "scale" and "Teuchos::ArrayView< const Scalar > alpha" in argsstring: return ''
    if name == "scale" and "const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A" in argsstring: return ''

    if name in conf_RemoveRefFunctionList: declStr = declStr.replace('&', '')
        
    descStr = "    //! " + briefdescription.lstrip().rstrip() + "\n"
    declStr = "    virtual " + declStr + "= 0;"

    return descStr + declStr + "\n" + "\n"
####

xml_dir = '../../../../packages/tpetra/doc/xml/'
conf_dir = 'interfaces/conf/'
tmpl_dir = 'interfaces/tmpl/'
out_dir = '../src/'

for file in os.listdir(conf_dir):
    basename, extension = os.path.splitext(file)
    if extension == ".conf":

#### READ CONFIG ####
        parser = SafeConfigParser()
        parser.read(conf_dir + file)
        
        conf_XMLheaders = xml_dir + parser.get('io', 'XMLheaders')
        conf_XMLclass   = xml_dir + parser.get('io', 'XMLclass')
        conf_template   = tmpl_dir + parser.get('io', 'template')
        conf_output     = parser.get('io', 'output')
        
        conf_SkipFunctionList = set(parser.get('function', 'skip').split(';'))
        conf_RemoveRefFunctionList = set(parser.get('function', 'removeref').split(';'))
        conf_SkipHeaderList = set(parser.get('header', 'skip').split(';'))
#
        
        template = open(conf_template, 'r').read()
        out = Template(template)
        
        className = buildClassDefinition(conf_XMLclass)
        out = out.substitute(
            TMPL_HEADERS=buildHeader(className, 'interfaces.py'),
            TMPL_INCLUDES=buildInclude(conf_XMLheaders, conf_SkipHeaderList),
            TMPL_TEMPLATE_PARAM=buildTemplateParam(conf_XMLclass),
            TMPL_CLASS=className,
            TMPL_INHERITANCE='  ' + parser.get('inheritance', 'parent').rstrip(),
            TMPL_DESTRUCTOR=buildDestructor(className),
            TMPL_PUBLIC_FUNCTIONS=buildClassFunctions(conf_XMLclass, conf_SkipFunctionList, buildFuncLineInterface),
            TMPL_FOOTERS=buildFooter(className)
            )
        f = open(out_dir + conf_output, 'w')
        f.write(out)
        f.close()
