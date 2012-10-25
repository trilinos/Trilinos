import sys
sys.path.append('lib/')

import os
from string import Template
from ConfigParser import SafeConfigParser
from XpetraLib import *
from XpetraLibConfig import *


def buildFuncLineTpetra( functionNode ):

#TODO: clean up
    tree = etree.parse(conf_XMLclass)
    root = tree.getroot() # root == <doxygen>
    classNode = root[0]   # classNode == <compounddef>
    
    fullClassName = classNode.xpath('compoundname')[0].text # Tpetra::Map
    baseClassName = fullClassName.lstrip('Tpetra::')        # Map
    className = 'Tpetra'+baseClassName                      # TpetraMap
##

    # <name> = function name
    name = functionNode.xpath('name')[0].text
    if name == baseClassName: name = className
    if name == '~'+baseClassName: name = '~'+className
        
    # <type> = return type of the function
    type = functionNode.xpath('type')[0].xpath("string()")
    
    # <argsstring>
    argsstring = functionNode.xpath('argsstring')[0].text

    # skip deprecated functions
    if 'TPETRA_DEPRECATED' in type: return ''
    
    # hack for Vector:
    # - add missing 'typename'
    # - do not add MultiVector inherited methods
    if 'magnitudeType' in type: type = 'typename ' + type
    if functionNode.xpath('//compoundname')[0].text == 'Tpetra::Vector':
        if name in ['dot','norm1','norm2','normInf','normWeighted','meanValue'] and 'ArrayView' in argsstring: return ''
    if functionNode.xpath('//compoundname')[0].text == 'Tpetra::Vector':
        if name in ['replaceGlobalValue','sumIntoGlobalValue','replaceLocalValue','sumIntoLocalValue'] and 'size_t vectorIndex' in argsstring: return ''

    # hack for MultiVector
    #  if name == "scale" and "Teuchos::ArrayView< const Scalar > alpha" in argsstring: return ''
    #  if name == "scale" and "const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A" in argsstring: return ''

    #hack for CrsMatrix
    if name == "getLocalRowCopy" and "const =0" in argsstring: return '' 
    if name == "TpetraCrsMatrix" and "const RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph" in argsstring: return ''
    if className == "TpetraCrsMatrix" and "const =0" in argsstring: argsstring = argsstring.replace("const =0", "const")
    
    #hack for RowMatrix
    if className == "TpetraRowMatrix" and "const =0" in argsstring: argsstring = argsstring.replace("const =0", "const")

    # <param> -> get list of arg name as a string 'GIDList, nodeIDList, LIDList'
    # Simple version
    #        paramList = functionNode.xpath('param/declname/text()')
    #        paramStr  = ', '.join(param for param in paramList)
    
    # More complete version
    paramStr  = ''
    paramNodes = functionNode.xpath('param')
    #print name
    for paramNode in paramNodes:
        n = paramNode.xpath('declname')[0].xpath("string()")
        if paramNode.xpath('type')[0].xpath("string()") in conf_TypeWrapped:
            paramStr += "toTpetra(" + n + ")"
        else:
            paramStr += n
            
        paramStr += ", "
            
    paramStr = paramStr.rstrip(', ')
                    
    # briefdescription
    briefdescription = functionNode.xpath("briefdescription")[0].xpath("string()")
    
    if len(type) > 0:
        declStr = type + " " + name + argsstring
    else:
        declStr = name + argsstring
    declStr = declStr.rstrip()

    if name in conf_RemoveRefFunctionList: declStr = declStr.replace('&', '')
    
    descStr = "    //! " + briefdescription.lstrip().rstrip() + "\n"
    defStr  = "    " + declStr

    if name != className and name != "~"+className:
        defStr += " { "
        defStr += "XPETRA_MONITOR(\"" + className + "::" + name + "\"); "
        if len(type) > 0 and type != 'void': defStr += 'return '
        if type in conf_TypeWrapped: defStr += "toXpetra("
        defStr += conf_memberName + "->" + name
        defStr += "(" + paramStr
        if type in conf_TypeWrapped: defStr += ")"
        defStr += "); }"

    # constructor
    if name == className:
        defStr += "\n      " + ": " + conf_memberName + "(Teuchos::rcp(new " + fullClassName + "< "+templateParam+" >"
        defStr += "(" + paramStr + "))) { "
        defStr += " }"
      
    # destructor
    if name == '~'+className:
        defStr += " { "
        defStr += " }"
        
    return descStr + defStr + "\n" + "\n";

####

xml_dir = trilinosRoot_dir + '/packages/tpetra/doc/xml/'
conf_dir = 'tpetra/conf/'
tmpl_dir = 'tpetra/tmpl/'
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
        conf_memberName = parser.get('member', 'name')
        conf_TypeWrapped = set(parser.get('type', 'wrapped').split(';'))
#
        
        template = open(conf_template, 'r').read()
        out = Template(template)
        
        className = buildClassDefinition(conf_XMLclass, 'Tpetra')
        templateParam = buildTemplateParam2(conf_XMLclass)
        
        out = out.substitute(
            TMPL_HEADERS=buildHeader(className, 'tpetra.py'),
            TMPL_INCLUDES=buildInclude(conf_XMLheaders, conf_SkipHeaderList),
            TMPL_TEMPLATE_PARAM=buildTemplateParam(conf_XMLclass),
            TMPL_CLASS=className,
            TMPL_INHERITANCE='  ' + parser.get('inheritance', 'parent').rstrip(),
            TMPL_DESTRUCTOR=buildDestructor(className),
            TMPL_PUBLIC_FUNCTIONS=buildClassFunctions(conf_XMLclass, conf_SkipFunctionList, buildFuncLineTpetra),
            TMPL_FOOTERS=buildFooter(className)
            )
        f = open(out_dir + conf_output, 'w')
        f.write(out)
        f.close()

