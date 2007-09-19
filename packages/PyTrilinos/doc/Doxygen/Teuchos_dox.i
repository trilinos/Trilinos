
// File: index.xml

// File: classTeuchos_1_1EmptyXMLError.xml
%feature("docstring") Teuchos::EmptyXMLError "

Thrown when attempting to parse an empty XML std::string.

C++ includes: Teuchos_XMLObject.hpp ";

%feature("docstring")  Teuchos::EmptyXMLError::EmptyXMLError "Teuchos::EmptyXMLError::EmptyXMLError(const std::string &what_arg) ";


// File: classTeuchos_1_1FileInputSource.xml
%feature("docstring") Teuchos::FileInputSource "

Instantiation of XMLInputSource class for reading XML from a file.

C++ includes: Teuchos_FileInputSource.hpp ";

%feature("docstring")  Teuchos::FileInputSource::FileInputSource "FileInputSource::FileInputSource(const std::string &filename)

Constructor. ";

%feature("docstring")  Teuchos::FileInputSource::~FileInputSource "virtual Teuchos::FileInputSource::~FileInputSource()

Destructor. ";

%feature("docstring")  Teuchos::FileInputSource::stream "RCP<
XMLInputStream > FileInputSource::stream() const

Create a FileInputStream. ";


// File: classTeuchos_1_1Exceptions_1_1InvalidParameter.xml
%feature("docstring") Teuchos::Exceptions::InvalidParameter "C++
includes: Teuchos_ParameterListExceptions.hpp ";

%feature("docstring")
Teuchos::Exceptions::InvalidParameter::InvalidParameter "Teuchos::Exceptions::InvalidParameter::InvalidParameter(const
std::string &what_arg) ";


// File: classTeuchos_1_1Exceptions_1_1InvalidParameterName.xml
%feature("docstring") Teuchos::Exceptions::InvalidParameterName "C++
includes: Teuchos_ParameterListExceptions.hpp ";

%feature("docstring")
Teuchos::Exceptions::InvalidParameterName::InvalidParameterName "Teuchos::Exceptions::InvalidParameterName::InvalidParameterName(const
std::string &what_arg) ";


// File: classTeuchos_1_1Exceptions_1_1InvalidParameterType.xml
%feature("docstring") Teuchos::Exceptions::InvalidParameterType "C++
includes: Teuchos_ParameterListExceptions.hpp ";

%feature("docstring")
Teuchos::Exceptions::InvalidParameterType::InvalidParameterType "Teuchos::Exceptions::InvalidParameterType::InvalidParameterType(const
std::string &what_arg) ";


// File: classTeuchos_1_1Exceptions_1_1InvalidParameterValue.xml
%feature("docstring") Teuchos::Exceptions::InvalidParameterValue "C++
includes: Teuchos_ParameterListExceptions.hpp ";

%feature("docstring")
Teuchos::Exceptions::InvalidParameterValue::InvalidParameterValue "Teuchos::Exceptions::InvalidParameterValue::InvalidParameterValue(const
std::string &what_arg) ";


// File: classTeuchos_1_1ParameterList.xml
%feature("docstring") Teuchos::ParameterList "

Templated parameter list.

Parameters can be added and retreived with the templated \"get\" and
\"set\" functions. These parameters can any data type which uses value
sementics (e.g. double, float, int, *double, *float, *int, ...) which
includes other parameter lists, allowing for a hierarchy of parameter
lists. These parameters can also be pointers to vectors or functions.

Use static_cast<T>() when the type is ambiguous.

Both char* and std::string std::map to are stored as strings
internally.

C++ includes: Teuchos_ParameterList.hpp ";


// File: classTeuchos_1_1ParameterList_1_1PrintOptions.xml
%feature("docstring") Teuchos::ParameterList::PrintOptions "

Utility class for setting and passing in print options.

C++ includes: Teuchos_ParameterList.hpp ";

%feature("docstring")
Teuchos::ParameterList::PrintOptions::PrintOptions "Teuchos::ParameterList::PrintOptions::PrintOptions() ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::indent "PrintOptions& Teuchos::ParameterList::PrintOptions::indent(int
_indent) ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::showTypes
"PrintOptions& Teuchos::ParameterList::PrintOptions::showTypes(bool
_showTypes) ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::showFlags
"PrintOptions& Teuchos::ParameterList::PrintOptions::showFlags(bool
_showFlags) ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::showDoc "PrintOptions& Teuchos::ParameterList::PrintOptions::showDoc(bool
_showDoc) ";

%feature("docstring")
Teuchos::ParameterList::PrintOptions::incrIndent "PrintOptions&
Teuchos::ParameterList::PrintOptions::incrIndent(int indents) ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::indent "int Teuchos::ParameterList::PrintOptions::indent() const ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::showTypes
"bool Teuchos::ParameterList::PrintOptions::showTypes() const ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::showFlags
"bool Teuchos::ParameterList::PrintOptions::showFlags() const ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::showDoc "bool Teuchos::ParameterList::PrintOptions::showDoc() const ";

%feature("docstring")  Teuchos::ParameterList::PrintOptions::copy "PrintOptions Teuchos::ParameterList::PrintOptions::copy() const ";


// File: classTeuchos_1_1ParameterListAcceptor.xml
%feature("docstring") Teuchos::ParameterListAcceptor "

Base class objects that can accept a parameter list.

ToDo: Finish Documentation!

C++ includes: Teuchos_ParameterListAcceptor.hpp ";

%feature("docstring")
Teuchos::ParameterListAcceptor::~ParameterListAcceptor "Teuchos::ParameterListAcceptor::~ParameterListAcceptor() ";


// File: structTeuchos_1_1ScalarTraits.xml
%feature("docstring") Teuchos::ScalarTraits "

This structure defines some basic traits for a scalar field type.

Scalar traits are an essential part of templated codes. This structure
offers the basic traits of the templated scalar type, like defining
zero and one, and basic functions on the templated scalar type, like
performing a square root.

The functions in the templated base unspecialized struct are designed
not to compile (giving a nice compile-time error message) and
therefore specializations must be written for Scalar types actually
used.

The default defined specializations are provided for int, float, and
double.

ScalarTraits can be used with the Arbitrary Precision Library (
http://crd.lbl.gov/~dhbailey/mpdist/ ) by configuring Teuchos with
--enable-teuchos-arprec and giving the appropriate paths to ARPREC.
Then ScalarTraits has the specialization: mp_real.

If Teuchos is configured with --enable-teuchos-stdcomplex then
ScalarTraits also has a parital specialization for all std::complex
numbers of the form std::complex<T>.

C++ includes: Teuchos_ScalarTraits.hpp ";


// File: structTeuchos_1_1ScalarTraits_3_01char_01_4.xml
%feature("docstring") Teuchos::ScalarTraits< char > " ";


// File: structTeuchos_1_1ScalarTraits_3_01double_01_4.xml
%feature("docstring") Teuchos::ScalarTraits< double > " ";


// File: structTeuchos_1_1ScalarTraits_3_01float_01_4.xml
%feature("docstring") Teuchos::ScalarTraits< float > " ";


// File: structTeuchos_1_1ScalarTraits_3_01int_01_4.xml
%feature("docstring") Teuchos::ScalarTraits< int > " ";


// File: classTeuchos_1_1StringInputSource.xml
%feature("docstring") Teuchos::StringInputSource "

Instantiation of XMLInputSource class for reading XML from a
std::string.

C++ includes: Teuchos_StringInputSource.hpp ";

%feature("docstring")  Teuchos::StringInputSource::StringInputSource "StringInputSource::StringInputSource(const std::string &text)

Constructor. ";

%feature("docstring")  Teuchos::StringInputSource::~StringInputSource
"virtual Teuchos::StringInputSource::~StringInputSource()

Destructor. ";

%feature("docstring")  Teuchos::StringInputSource::stream "RCP<
XMLInputStream > StringInputSource::stream() const

Create a StringInputStream. ";


// File: classTeuchos_1_1Time.xml
%feature("docstring") Teuchos::Time "

Basic wall-clock timer class.

To time a section of code, place it in between calls to start() and
stop().

For std::exception safety and correct behavior in reentrant code, this
class should generally be used only through the Teuchos::TimeMonitor
mechanism.

C++ includes: Teuchos_Time.hpp ";

%feature("docstring")  Teuchos::Time::Time "Teuchos::Time::Time(const
std::string &name, bool start=false)

Construct with a descriptive name. ";

%feature("docstring")  Teuchos::Time::start "void
Teuchos::Time::start(bool reset=false)

Starts the timer. ";

%feature("docstring")  Teuchos::Time::stop "double
Teuchos::Time::stop()

Stops the timer. ";

%feature("docstring")  Teuchos::Time::totalElapsedTime "double
Teuchos::Time::totalElapsedTime(bool readCurrentTime=false) const

Returns the total time accumulated by this timer. This should be
called only when the clock is stopped.. ";

%feature("docstring")  Teuchos::Time::reset "void
Teuchos::Time::reset()

Resets the cummulative time and number of times this timer has been
called. Does not affect any other state. ";

%feature("docstring")  Teuchos::Time::isRunning "bool
Teuchos::Time::isRunning() const

Indicates if this timer is currently running, i.e., if it has been
started but not yet stopped.

It is necessary to know if a timer is running to avoid incorrectly
starting or stopping in reentrant code. ";

%feature("docstring")  Teuchos::Time::name "const std::string&
Teuchos::Time::name() const

Return the name of this timer. ";

%feature("docstring")  Teuchos::Time::incrementNumCalls "void
Teuchos::Time::incrementNumCalls()

Increment the number of times this timer has been called. ";

%feature("docstring")  Teuchos::Time::numCalls "int
Teuchos::Time::numCalls() const

Return the number of times this timer has been called. ";


// File: classTeuchos_1_1TypeNameTraits_3_01ParameterList_01_4.xml
%feature("docstring") Teuchos::TypeNameTraits< ParameterList > "

Traits specialization.

C++ includes: Teuchos_ParameterList.hpp ";


// File: structTeuchos_1_1UndefinedScalarTraits.xml
%feature("docstring") Teuchos::UndefinedScalarTraits "";


// File: classTeuchos_1_1XMLInputSource.xml
%feature("docstring") Teuchos::XMLInputSource "

XMLInputSource represents a source of XML input that can be parsed to
produce an XMLObject.

The source might be a file, a socket, a std::string. The XMLObject is
created with a call to the getObject() method.

The source gets its data from a XMLInputStream object that is created
(internally) to work with this source.

getObject() is implemented with EXPAT if Teuchos is configured with
--enable-expat.

C++ includes: Teuchos_XMLInputSource.hpp ";

%feature("docstring")  Teuchos::XMLInputSource::XMLInputSource "Teuchos::XMLInputSource::XMLInputSource()

Empty constructor. ";

%feature("docstring")  Teuchos::XMLInputSource::~XMLInputSource "virtual Teuchos::XMLInputSource::~XMLInputSource()

Destructor. ";

%feature("docstring")  Teuchos::XMLInputSource::stream "virtual
RCP<XMLInputStream> Teuchos::XMLInputSource::stream() const =0

Virtual input source interface. ";

%feature("docstring")  Teuchos::XMLInputSource::getObject "XMLObject
XMLInputSource::getObject() const

Get an object by invoking the TreeBuildingXMLHandler on the input
data. ";


// File: classTeuchos_1_1XMLInputStream.xml
%feature("docstring") Teuchos::XMLInputStream "

XMLInputStream represents an XML input stream that can be used by a
XMLInputSource.

C++ includes: Teuchos_XMLInputStream.hpp ";

%feature("docstring")  Teuchos::XMLInputStream::XMLInputStream "Teuchos::XMLInputStream::XMLInputStream()

Constructor. ";

%feature("docstring")  Teuchos::XMLInputStream::~XMLInputStream "virtual Teuchos::XMLInputStream::~XMLInputStream()

Destructor. ";

%feature("docstring")  Teuchos::XMLInputStream::readBytes "virtual
unsigned int Teuchos::XMLInputStream::readBytes(unsigned char *const
toFill, const unsigned int maxToRead)=0

Read up to maxToRead bytes from the stream. ";

%feature("docstring")  Teuchos::XMLInputStream::curPos "unsigned int
XMLInputStream::curPos() const

Identify current position. ";


// File: classTeuchos_1_1XMLObject.xml
%feature("docstring") Teuchos::XMLObject "

Representation of an XML data tree. XMLObject is a ref-counted handle
to a XMLObjectImplem object, allowing storage by reference.

C++ includes: Teuchos_XMLObject.hpp ";


// File: classTeuchos_1_1XMLObjectImplem.xml
%feature("docstring") Teuchos::XMLObjectImplem "

The XMLObjectImplem class takes care of the low-level implementation
details of XMLObject.

C++ includes: Teuchos_XMLObjectImplem.hpp ";

%feature("docstring")  Teuchos::XMLObjectImplem::XMLObjectImplem "XMLObjectImplem::XMLObjectImplem(const std::string &tag)

Construct with a 'tag'. ";

%feature("docstring")  Teuchos::XMLObjectImplem::deepCopy "XMLObjectImplem * XMLObjectImplem::deepCopy() const

Deep copy. ";

%feature("docstring")  Teuchos::XMLObjectImplem::addAttribute "void
XMLObjectImplem::addAttribute(const std::string &name, const
std::string &value)

Add a [name, value] attribute. ";

%feature("docstring")  Teuchos::XMLObjectImplem::addChild "void
XMLObjectImplem::addChild(const XMLObject &child)

Add a child XMLObject. ";

%feature("docstring")  Teuchos::XMLObjectImplem::addContent "void
XMLObjectImplem::addContent(const std::string &contentLine)

Add a content line. ";

%feature("docstring")  Teuchos::XMLObjectImplem::getTag "const
std::string& Teuchos::XMLObjectImplem::getTag() const

Return the tag std::string. ";

%feature("docstring")  Teuchos::XMLObjectImplem::hasAttribute "bool
Teuchos::XMLObjectImplem::hasAttribute(const std::string &name) const

Determine whether an attribute exists. ";

%feature("docstring")  Teuchos::XMLObjectImplem::getAttribute "const
std::string& Teuchos::XMLObjectImplem::getAttribute(const std::string
&name) const

Look up an attribute by name. ";

%feature("docstring")  Teuchos::XMLObjectImplem::numChildren "int
XMLObjectImplem::numChildren() const

Return the number of children. ";

%feature("docstring")  Teuchos::XMLObjectImplem::getChild "const
XMLObject & XMLObjectImplem::getChild(int i) const

Look up a child by its index. ";

%feature("docstring")  Teuchos::XMLObjectImplem::numContentLines "int
Teuchos::XMLObjectImplem::numContentLines() const

Get the number of content lines. ";

%feature("docstring")  Teuchos::XMLObjectImplem::getContentLine "const std::string& Teuchos::XMLObjectImplem::getContentLine(int i)
const

Look up a content line by index. ";

%feature("docstring")  Teuchos::XMLObjectImplem::print "void
XMLObjectImplem::print(std::ostream &os, int indent) const

Print to stream with the given indentation level. Output will be well-
formed XML. ";

%feature("docstring")  Teuchos::XMLObjectImplem::toString "std::string XMLObjectImplem::toString() const

Write as a std::string. Output may be ill-formed XML. ";

%feature("docstring")  Teuchos::XMLObjectImplem::header "std::string
XMLObjectImplem::header(bool strictXML=false) const

Write the header. ";

%feature("docstring")  Teuchos::XMLObjectImplem::terminatedHeader "std::string XMLObjectImplem::terminatedHeader(bool strictXML=false)
const

Write the header terminated as <Header>. ";

%feature("docstring")  Teuchos::XMLObjectImplem::footer "std::string
Teuchos::XMLObjectImplem::footer() const

Write the footer. ";


// File: classTeuchos_1_1XMLParameterListReader.xml
%feature("docstring") Teuchos::XMLParameterListReader "

Writes an XML object to a parameter list.

C++ includes: Teuchos_XMLParameterListReader.hpp ";

%feature("docstring")
Teuchos::XMLParameterListReader::toParameterList "ParameterList
XMLParameterListReader::toParameterList(const XMLObject &xml) const

Write the given XML object to a parameter list ";


// File: classTeuchos_1_1XMLParameterListWriter.xml
%feature("docstring") Teuchos::XMLParameterListWriter "

Writes a ParameterList to an XML object.

C++ includes: Teuchos_XMLParameterListWriter.hpp ";

%feature("docstring")  Teuchos::XMLParameterListWriter::toXML "XMLObject XMLParameterListWriter::toXML(const ParameterList &p) const

Write the given list to an XML object ";


// File: classTeuchos_1_1XMLParser.xml
%feature("docstring") Teuchos::XMLParser "

XMLParser consumes characters from an XMLInputStream object, parsing
the XML and using a TreeBuildingXMLHandler to construct an XMLObject.

C++ includes: Teuchos_XMLParser.hpp ";

%feature("docstring")  Teuchos::XMLParser::XMLParser "Teuchos::XMLParser::XMLParser(RCP< XMLInputStream > is)

Constructor. ";

%feature("docstring")  Teuchos::XMLParser::~XMLParser "Teuchos::XMLParser::~XMLParser()

Destructor. ";

%feature("docstring")  Teuchos::XMLParser::parse "XMLObject
XMLParser::parse()

Consume the XMLInputStream to build an XMLObject. ";


// File: namespace@0.xml


// File: namespaceTeuchos.xml
%feature("docstring")  Teuchos::Exceptions::haveSameValues "bool
Teuchos::haveSameValues(const ParameterList &list1, const
ParameterList &list2) ";

%feature("docstring")  Teuchos::Exceptions::parameterList "RCP<ParameterList> Teuchos::parameterList()

Nonmember constructor. ";

%feature("docstring")  Teuchos::Exceptions::parameterList "RCP<ParameterList> Teuchos::parameterList(const std::string &name)

Nonmember constructor. ";

%feature("docstring")  Teuchos::Exceptions::parameterList "RCP<ParameterList> Teuchos::parameterList(const ParameterList &source)

Nonmember constructor. ";

%feature("docstring")
Teuchos::Exceptions::throwScalarTraitsNanInfError "void
Teuchos::throwScalarTraitsNanInfError(const std::string &errMsg) ";

%feature("docstring")
Teuchos::Exceptions::updateParametersFromXmlFile "void
Teuchos::updateParametersFromXmlFile(const std::string
&xmlFileName,Teuchos::ParameterList *paramList)

Reads XML parameters from a file and updates those already in the
given parameter list.

Parameters:
-----------

xmlFileName:  [in] The file name containing XML parameter list
specification.

paramList:  [in/out] On input, *paramList may be empty or contain some
parameters and sublists. On output, parameters and sublist from the
file xmlFileName will be set or overide those in *paramList. ";

%feature("docstring")
Teuchos::Exceptions::updateParametersFromXmlString "void
Teuchos::updateParametersFromXmlString(const std::string
&xmlStr,Teuchos::ParameterList *paramList)

Reads XML parameters from a std::string and updates those already in
the given parameter list.

Parameters:
-----------

xmlStr:  [in] String containing XML parameter list specification.

paramList:  [in/out] On input, *paramList may be empty or contain some
parameters and sublists. On output, parameters and sublist from the
file xmlStr will be set or overide those in *paramList. ";

%feature("docstring")
Teuchos::Exceptions::writeParameterListToXmlOStream "void
Teuchos::writeParameterListToXmlOStream(const Teuchos::ParameterList
&paramList,std::ostream &xmlOut)

Write parameters and sublists in XML format to an std::ostream.

Parameters:
-----------

paramList:  [in] Contains the parameters and sublists that will be
written to file.

xmlOut:  [in] The stream that will get the XML output. ";

%feature("docstring")
Teuchos::Exceptions::writeParameterListToXmlFile "void
Teuchos::writeParameterListToXmlFile(const Teuchos::ParameterList
&paramList,const std::string &xmlFileName)

Write parameters and sublist to an XML file.

Parameters:
-----------

paramList:  [in] Contains the parameters and sublists that will be
written to file.

xmlFileName:  [in] The file name that will be create to contain the
XML version of the parameter list specification. ";


// File: namespaceTeuchos_1_1Exceptions.xml


// File: Teuchos__FileInputSource_8cpp.xml


// File: Teuchos__FileInputSource_8hpp.xml


// File: Teuchos__ParameterList_8cpp.xml
%feature("docstring")  Teuchos::filterValueToString "std::string
@0::filterValueToString(const Teuchos::ParameterEntry &entry) ";


// File: Teuchos__ParameterList_8hpp.xml


// File: Teuchos__ParameterListAcceptor_8cpp.xml


// File: Teuchos__ParameterListAcceptor_8hpp.xml


// File: Teuchos__ParameterListExceptions_8hpp.xml


// File: Teuchos__ScalarTraits_8cpp.xml
%feature("docstring")  returnFloatZero "float @0::returnFloatZero()
";

%feature("docstring")  returnDoubleZero "double
@0::returnDoubleZero() ";


// File: Teuchos__ScalarTraits_8hpp.xml


// File: Teuchos__StringInputSource_8cpp.xml


// File: Teuchos__StringInputSource_8hpp.xml


// File: Teuchos__Time_8cpp.xml


// File: Teuchos__Time_8hpp.xml


// File: Teuchos__XMLInputSource_8cpp.xml


// File: Teuchos__XMLInputSource_8hpp.xml


// File: Teuchos__XMLInputStream_8cpp.xml


// File: Teuchos__XMLInputStream_8hpp.xml


// File: Teuchos__XMLObject_8cpp.xml


// File: Teuchos__XMLObject_8hpp.xml


// File: Teuchos__XMLObjectImplem_8cpp.xml


// File: Teuchos__XMLObjectImplem_8hpp.xml


// File: Teuchos__XMLParameterListHelpers_8cpp.xml


// File: Teuchos__XMLParameterListHelpers_8hpp.xml


// File: Teuchos__XMLParameterListReader_8cpp.xml


// File: Teuchos__XMLParameterListReader_8hpp.xml


// File: Teuchos__XMLParameterListWriter_8cpp.xml


// File: Teuchos__XMLParameterListWriter_8hpp.xml


// File: Teuchos__XMLParser_8cpp.xml


// File: Teuchos__XMLParser_8hpp.xml


// File: dir_daa2443682d0f547b84c7fa838636502.xml


// File: dir_2c9d975476051ec1d3cf6e8ee401cf57.xml


// File: ParameterList_2cxx__main_8cpp-example.xml

