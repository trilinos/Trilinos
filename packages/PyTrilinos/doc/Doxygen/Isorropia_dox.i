
// File: index.xml

// File: classIsorropia_1_1Colorer.xml
%feature("docstring") Isorropia::Colorer "

Interface (abstract base class) for computing a new coloring and
describing the result.

The colors returned have values between 1 and C, where C is the number
of colors used.

C++ includes: Isorropia_Colorer.hpp ";

%feature("docstring")  Isorropia::Colorer::~Colorer "virtual
Isorropia::Colorer::~Colorer()

Destructor ";

%feature("docstring")  Isorropia::Colorer::color "virtual void
Isorropia::Colorer::color(bool forceColoring=false)=0

Method which does the work of computing a new coloring.

Parameters:
-----------

forceColoring:  Optional argument defaults to false. Depending on the
implementation, color() should only perform a coloring the first time
it is called, and subsequent repeated calls are no-ops. If the user's
intent is to re- compute the coloring (e.g., if parameters or other
inputs have been changed), then setting this flag to true will force a
new coloring to be computed. ";

%feature("docstring")  Isorropia::Colorer::numColors "virtual int
Isorropia::Colorer::numColors() const

Method which returns the number (global) of colors used.

The overall number of colors used. All colors used for all vertices
are between 1 and this value (included).

See:   Isorropia::Operator::numProperties() ";

%feature("docstring")  Isorropia::Colorer::numElemsWithColor "virtual
int Isorropia::Colorer::numElemsWithColor(int color) const

Return the number of local elements of a given color.

Parameters:
-----------

color:  The wanted color.

The number of local of the asked color.

See:   Isorropia::Operator::numElemsWithProperty() ";

%feature("docstring")  Isorropia::Colorer::elemsWithColor "virtual
void Isorropia::Colorer::elemsWithColor(int color, int *elementList,
int len) const

Fill user-allocated list (of length len) with the local element ids
for LOCAL elements of the given color.

Parameters:
-----------

color:  the wanted color

elementList:  an array to receive local elements of the given color

len:  the number of elements wanted

See:   Isorropia::Operator::elemsWithProperty() ";

%feature("docstring")  Isorropia::Colorer::extractColorsView "virtual
int Isorropia::Colorer::extractColorsView(int &size, const int
*&array) const

Give access of the color assignments array that is owned by the
current processor.

Parameters:
-----------

size:  Number of elements in the array.

array:  Pointer to the color assignements array inside the object.

This pointer is only significant if the object still exists.
Otherwise, you must use

See:  Isorropia::Operator::extractPartsCopy()

Isorropia::Operator::extractPropertiesView() ";

%feature("docstring")  Isorropia::Colorer::extractColorsCopy "virtual
int Isorropia::Colorer::extractColorsCopy(int len, int &size, int
*array) const

Copy a part of the color assignments array.

Parameters:
-----------

len:  of the array given by the user.

size:  Number of elements in the array.

array:  Array of color assignments. Allocated by the user with a size
of at least len elements.

Memory space which is not useful in the array is not initialized or
used in this method.

See:   Isorropia::Operator::extractPropertiesCopy() ";


// File: classIsorropia_1_1CostDescriber.xml
%feature("docstring") Isorropia::CostDescriber "

Interface (abstract base class) for describing the weights or costs
associated with the vertices and/or edges or hyperedges of the object
to be partitioned, ordered or colored.

A CostDescriber object is created by the application. If no
CostDescriber is supplied by the application, sensible default weights
should be used.

C++ includes: Isorropia_CostDescriber.hpp ";

%feature("docstring")  Isorropia::CostDescriber::~CostDescriber "virtual Isorropia::CostDescriber::~CostDescriber()

Destructor ";


// File: classIsorropia_1_1Exception.xml
%feature("docstring") Isorropia::Exception "

A simple extension of std::exception with constructors that accept a
message in the form of a 'const char*' or a 'stdstring'.

C++ includes: Isorropia_Exception.hpp ";

%feature("docstring")  Isorropia::Exception::Exception "Isorropia::Exception::Exception(const char *msg)  throw () constructor
that accepts a const char-ptr ";

%feature("docstring")  Isorropia::Exception::Exception "Isorropia::Exception::Exception(std::string msg)  throw () constructor
that accepts a std::string ";

%feature("docstring")  Isorropia::Exception::~Exception "Isorropia::Exception::~Exception()  throw () destructor ";

%feature("docstring")  Isorropia::Exception::what "const char *
Isorropia::Exception::what() const  throw () return const char-ptr of
exception message ";


// File: classIsorropia_1_1LevelScheduler.xml
%feature("docstring") Isorropia::LevelScheduler "

Interface (abstract base class) for an operator that computes a
partitioning of local elements into levels. On each process the levels
begin at 1 and increase by 1. All elements in the same level can be
used in calculation concurrently. All elements in level i are
dependent upon the completion of the calculations involving the
elements in level i-1.

C++ includes: Isorropia_LevelScheduler.hpp ";

%feature("docstring")  Isorropia::LevelScheduler::~LevelScheduler "virtual Isorropia::LevelScheduler::~LevelScheduler()

Destructor ";

%feature("docstring")  Isorropia::LevelScheduler::schedule "virtual
void Isorropia::LevelScheduler::schedule(bool forceScheduling=false)=0

Method which does the work of computing a new level schedule.

Parameters:
-----------

forceScheduling:  Optional argument defaults to false. Depending on
the implementation, schedule() should only perform a scheduling the
first time it is called, and subsequent repeated calls are no-ops. If
the user's intent is to re- compute the scheduling (e.g., if
parameters or other inputs have been changed), then setting this flag
to true will force a new scheduling to be computed. ";

%feature("docstring")  Isorropia::LevelScheduler::numLevels "virtual
int Isorropia::LevelScheduler::numLevels() const

Method which returns the number of levels.

The number of levels on the local process. Levels begin at 1 and
increase by 1.

See:   Isorropia::Operator::numProperties() ";

%feature("docstring")  Isorropia::LevelScheduler::numElemsWithLevel "virtual int Isorropia::LevelScheduler::numElemsWithLevel(int level)
const

Return the number of elements in a given level.

Parameters:
-----------

level:  The wanted level.

The number of elements in the level.

See:   Isorropia::Operator::numElemsWithProperty() ";

%feature("docstring")  Isorropia::LevelScheduler::elemsWithLevel "virtual void Isorropia::LevelScheduler::elemsWithLevel(int level, int
*elementList, int len) const

Fill user-allocated list (of length len) with the local ID for each
element in the given level.

Parameters:
-----------

level:  the wanted level

elementList:  an array to receive local elements of the given level

len:  the number of elements wanted

See:   Isorropia::Operator::elemsWithProperty() ";


// File: classIsorropia_1_1Operator.xml
%feature("docstring") Isorropia::Operator "

Interface (abstract base class) for computing a new
partitioning/coloring/ ordering and exploiting their results.

If the accessor methods are called before the computation of the
result (by a method like compute()) has been called, behavior is not
well defined. Implementations will either return empty/erroneous data,
or throw an exception. In most cases, implementations will probably
call compute_partitioning() internally in a constructor or factory
method, so this won't usually be an issue.

C++ includes: Isorropia_Operator.hpp ";

%feature("docstring")  Isorropia::Operator::alreadyComputed "virtual
bool Isorropia::Operator::alreadyComputed() const =0

Query whether the computation has already been called.

True if the computation has already been done, False otherwise. ";

%feature("docstring")  Isorropia::Operator::numProperties "virtual
int Isorropia::Operator::numProperties() const =0

Return the number of different values used for \"properties\".

For example, the number of colors or the number of parts used for the
overall graph/matrix.

Global number of values for properties

Infact, it returns the upper bound of the interval of taken values.
For example, for the colors \"1,2,4\" , it will return \"4\" ";

%feature("docstring")  Isorropia::Operator::numLocalProperties "virtual int Isorropia::Operator::numLocalProperties() const =0

Return the number of different values used for \"properties\" for this
process only.

Local number of values for properties ";

%feature("docstring")  Isorropia::Operator::numElemsWithProperty "virtual int Isorropia::Operator::numElemsWithProperty(int property)
const =0

Return the number of LOCAL elements with the given property.

Parameters:
-----------

property:  Value of the property to consider.

Number of local elems which have this property. ";

%feature("docstring")  Isorropia::Operator::elemsWithProperty "virtual void Isorropia::Operator::elemsWithProperty(int property, int
*elementList, int len) const =0

Fill user-allocated list (of length len) with the local element ids of
the LOCAL elements with the given property.

Parameters:
-----------

property:  Value of the property to consider.

elementList:  User allocated array (of size at least len) of local ID
that have the asked property.

len:  Maximum lenght for the array. If len is greater than the result
of numElemsWithProperty() for property, only the first and relevant
elements are filled.

Memory space which is not useful in the array is not initialized or
used in this method. ";

%feature("docstring")  Isorropia::Operator::extractPropertiesView "virtual int Isorropia::Operator::extractPropertiesView(int &size,
const int *&array) const =0

Give access of the property array that is owned by the current
processor.

Parameters:
-----------

size:  Number of elements in the array.

array:  Pointer to the the properties array inside the object.

This pointer is only significant if the object still exists.
Otherwise, you must use

See:   Isorropia::Operator::extractPropertiesCopy().

Isorropia::Operator::extractPropertiesCopy() ";

%feature("docstring")  Isorropia::Operator::extractPropertiesCopy "virtual int Isorropia::Operator::extractPropertiesCopy(int len, int
&size, int *array) const =0

Copy a part of the property array.

Parameters:
-----------

len:  of the array given by the user.

size:  Number of elements in the array.

array:  Array of properties. Allocated by the user with a size of at
least len elements.

Memory space which is not useful in the array is not initialized or
used in this method.

See:   Isorropia::Operator::extractPropertiesView() ";

%feature("docstring")  Isorropia::Operator::~Operator "virtual
Isorropia::Operator::~Operator()

Destructor ";

%feature("docstring")  Isorropia::Operator::setParameters "virtual
void Isorropia::Operator::setParameters(const Teuchos::ParameterList
&paramlist)=0

Set parameters for the Operator instance. The contents of the input
paramlist object are copied into an internal ParameterList attribute.
Instances of this interface should not retain a reference to the input
ParameterList after this method returns.

Parameters:
-----------

paramlist:  List of parameters that the user wants to use. ";

%feature("docstring")  Isorropia::Operator::compute "virtual void
Isorropia::Operator::compute(bool forceRecomputing=false)=0

Method which does the work of computing a new
partitioning/coloring/ordering, depending on the child class used.

Parameters:
-----------

forceRecomputing:  Optional argument defaults to false. Depending on
the implementation, compute() should only perform a computation the
first time it is called, and subsequent repeated calls are no-ops. If
the user's intent is to re- compute the results (e.g., if parameters
or other inputs have been changed), then setting this flag to true
will force a new result to be computed. ";


// File: classIsorropia_1_1Orderer.xml
%feature("docstring") Isorropia::Orderer "

Interface (abstract base class) for computing a new ordering and
describing the layout of elements in the new order.

If the methods which describe the new ordering (e.g., operator[],
etc.) are called before order() has been called, behavior is not well
defined. Implementations will either return empty/erroneous data, or
throw an exception. In most cases, implementations will probably call
order() internally in a constructor or factory method, so this won't
usually be an issue.

C++ includes: Isorropia_Orderer.hpp ";

%feature("docstring")  Isorropia::Orderer::~Orderer "virtual
Isorropia::Orderer::~Orderer()

Destructor ";

%feature("docstring")  Isorropia::Orderer::order "virtual void
Isorropia::Orderer::order(bool forceOrdering=false)=0

Method which does the work of computing a new ordering.

Parameters:
-----------

forceOrdering:  Optional argument defaults to false. Depending on the
implementation, compute_partitioning() should only perform a
repartitioning the first time it is called, and subsequent repeated
calls are no-ops. If the user's intent is to re-compute the
partitioning (e.g., if parameters or other inputs have been changed),
then setting this flag to true will force a new partitioning to be
computed. ";

%feature("docstring")  Isorropia::Orderer::extractPermutationView "virtual int Isorropia::Orderer::extractPermutationView(int &size,
const int *&array) const

Give access of the \"direct\" permutation vector that is owned by the
current processor.

Parameters:
-----------

size:  Number of elements in the array.

array:  Pointer to the the part assignements array inside the object.

This pointer is only significant if the object still exists.
Otherwise, you must use

See:  Isorropia::Operator::extractPartsCopy()

Isorropia::Operator::extractPropertiesView() ";

%feature("docstring")  Isorropia::Orderer::extractPermutationCopy "virtual int Isorropia::Orderer::extractPermutationCopy(int len, int
&size, int *array) const

Copy a part of the \"direct\" permutation vector.

Parameters:
-----------

len:  of the array given by the user.

size:  Number of elements in the array.

array:  Direct permutation vector. Allocated by the user with a size
of at least len elements.

Memory space which is not useful in the array is not initialized or
used in this method.

See:   Isorropia::Operator::extractPropertiesCopy() ";


// File: classIsorropia_1_1Partitioner.xml
%feature("docstring") Isorropia::Partitioner "

Interface (abstract base class) for computing a new partitioning and
describing the layout of elements in the new partition (the parts).

If the methods which describe the new partitioning (e.g., operator [],
elemsInPart()) are called before compute_partitioning() has been
called, behavior is not well defined. Implementations will either
return empty/erroneous data, or throw an exception. In most cases,
implementations will probably call compute_partitioning() internally
in a constructor or factory method, so this won't usually be an issue.

C++ includes: Isorropia_Partitioner.hpp ";

%feature("docstring")  Isorropia::Partitioner::~Partitioner "virtual
Isorropia::Partitioner::~Partitioner()

Destructor ";

%feature("docstring")  Isorropia::Partitioner::partition "virtual
void Isorropia::Partitioner::partition(bool
forceRepartitioning=false)=0

Method which does the work of computing a new partitioning.
Implementations of this interface will typically be constructed with
an object or information describing the existing ('old') partitioning.
This method computes a 'new' rebalanced partitioning for that input
data.

Parameters:
-----------

forceRepartitioning:  Optional argument defaults to false. Depending
on the implementation, partitioning() should only perform a
repartitioning the first time it is called, and subsequent repeated
calls are no-ops. If the user's intent is to re-compute the
partitioning (e.g., if parameters or other inputs have been changed),
then setting this flag to true will force a new partitioning to be
computed.

See:   Isorropia::Operator::compute() ";

%feature("docstring")  Isorropia::Partitioner::numElemsInPart "virtual int Isorropia::Partitioner::numElemsInPart(int part) const =0

Return the number of LOCAL elements in a given part.

Parameters:
-----------

part:  the part ID we want to know the number of local elements.

number of local elements that belongs to the given part.

See:   Isorropia::Operator::numElemsWithProperty() ";

%feature("docstring")  Isorropia::Partitioner::elemsInPart "virtual
void Isorropia::Partitioner::elemsInPart(int part, int *elementList,
int len) const =0

Fill user-allocated list (of length len) with the local element ids to
be located in the given part

Parameters:
-----------

part:  the part ID we consider

elementList:  array of elements that belongs to this part ID, must be
allocated by user with size at least len

len:  maximum number of elements we can put in the array. Usually, may
be the result of Isorropia::Partitioner::numElemsInPart(). .

See:   Isorropia::Operator::elemsWithProperty() ";

%feature("docstring")  Isorropia::Partitioner::extractPartsView "virtual int Isorropia::Partitioner::extractPartsView(int &size, const
int *&array) const

Give access of the part assignments array that is owned by the current
processor.

Parameters:
-----------

size:  Number of elements in the array.

array:  Pointer to the the part assignements array inside the object.

This pointer is only significant if the object still exists.
Otherwise, you must use

See:  Isorropia::Operator::extractPartsCopy()

Isorropia::Operator::extractPropertiesView() ";

%feature("docstring")  Isorropia::Partitioner::extractPartsCopy "virtual int Isorropia::Partitioner::extractPartsCopy(int len, int
&size, int *array) const

Copy a part of the part assignment array.

Parameters:
-----------

len:  of the array given by the user.

size:  Number of elements in the array.

array:  Array of part assignments. Allocated by the user with a size
of at least len elements.

Memory space which is not useful in the array is not initialized or
used in this method.

See:   Isorropia::Operator::extractPropertiesCopy() ";


// File: classIsorropia_1_1Partitioner2D.xml
%feature("docstring") Isorropia::Partitioner2D "

Interface (abstract base class) for computing a new 2D partitioning
and describing the layout of elements in the new partitions.

If the methods which describe the new partitioning (e.g.,
newPartitionNumber(), etc.) are called before compute_partitioning()
has been called, behavior is not well defined. Implementations will
either return empty/erroneous data, or throw an exception. In most
cases, implementations will probably call compute_partitioning()
internally in a constructor or factory method, so this won't usually
be an issue.

C++ includes: Isorropia_Partitioner2D.hpp ";

%feature("docstring")  Isorropia::Partitioner2D::~Partitioner2D "virtual Isorropia::Partitioner2D::~Partitioner2D()

Destructor ";

%feature("docstring")  Isorropia::Partitioner2D::partition "virtual
void Isorropia::Partitioner2D::partition(bool
force_repartitioning=false)=0

Method which does the work of computing a new partitioning.
Implementations of this interface will typically be constructed with
an object or information describing the existing ('old') partitioning.
This method computes a 'new' rebalanced partitioning for that input
data.

Parameters:
-----------

force_repartitioning:  Optional argument defaults to false. Depending
on the implementation, compute_partitioning() should only perform a
repartitioning the first time it is called, and subsequent repeated
calls are no-ops. If the user's intent is to re-compute the
partitioning (e.g., if parameters or other inputs have been changed),
then setting this flag to true will force a new partitioning to be
computed. ";

%feature("docstring")  Isorropia::Partitioner2D::numElemsInPart "virtual int Isorropia::Partitioner2D::numElemsInPart(int part) const
=0

Return the number of LOCAL elements in a given part.

Parameters:
-----------

part:  the part ID we want to know the number of local elements.

number of local elements that belongs to the given part.

See:   Isorropia::Operator::numElemsWithProperty() ";

%feature("docstring")  Isorropia::Partitioner2D::elemsInPart "virtual
void Isorropia::Partitioner2D::elemsInPart(int part, int *elementList,
int len) const =0

Fill user-allocated list (of length len) with the local element ids to
be located in the given part

Parameters:
-----------

part:  the part ID we consider

elementList:  array of elements that belongs to this part ID, must be
allocated by user with size at least len

len:  maximum number of elements we can put in the array. Usually, may
be the result of Isorropia::Partitioner::numElemsInPart(). .

See:   Isorropia::Operator::elemsWithProperty() ";


// File: classIsorropia_1_1Redistributor.xml
%feature("docstring") Isorropia::Redistributor "

Abstract base class for classes which are constructed with a
Partitioner instance, and define methods to redistribute their
objects.

C++ includes: Isorropia_Redistributor.hpp ";

%feature("docstring")  Isorropia::Redistributor::~Redistributor "virtual Isorropia::Redistributor::~Redistributor()

Destructor ";


// File: namespaceIsorropia.xml
%feature("docstring")  Isorropia::Utils::Isorropia_Version "std::string Isorropia::Isorropia_Version() ";


// File: namespaceIsorropia_1_1Utils.xml
%feature("docstring")  Isorropia::Utils::create_comm_plan "void
Isorropia::Utils::create_comm_plan(int myPID, const std::vector< int >
&all_proc_old_offsets, const std::vector< int > &all_proc_new_offsets,
std::vector< int > &send_info, std::vector< int > &recv_info)

Internal Isorropia implementation utility. Given a vector that
specifies all processors' old or current offsets into a global element
list, and another vector that specifies all processors' new or desired
offsets into a global element list, fill a send_info and recv_info
vector with data that can be unpacked/interpreted as follows.

while(i<send_info.size()) { send_info[i] == proc to send data to
send_info[i+1] == starting offset of local elements to be sent
send_info[i+2] == number of elements to send i += 3; }

while(i<recv_info.size()) { recv_info[i] == proc to recv from
recv_info[i+1] == offset at which incoming elements will be stored
recv_info[i+2] == number of incoming elements } ";

%feature("docstring")  Isorropia::Utils::cpu_time "double
Isorropia::Utils::cpu_time()

Return CPU time. To measure an elapsed time, take the difference
between two returned values. ";


// File: Isorropia__Colorer_8hpp.xml


// File: Isorropia__ConfigDefs_8hpp.xml


// File: Isorropia__CostDescriber_8hpp.xml


// File: Isorropia__Exception_8cpp.xml


// File: Isorropia__Exception_8hpp.xml


// File: Isorropia__LevelScheduler_8hpp.xml


// File: Isorropia__Operator_8hpp.xml


// File: Isorropia__Orderer_8hpp.xml


// File: Isorropia__Partitioner_8hpp.xml


// File: Isorropia__Partitioner2D_8hpp.xml


// File: Isorropia__Redistributor_8hpp.xml


// File: Isorropia__Utils_8cpp.xml


// File: Isorropia__Utils_8hpp.xml


// File: Isorropia__Version_8hpp.xml


// File: dir_c705b74b31d8ba556986ebaf028813d3.xml


// File: dir_400f949df9dcdcdf586638eea2faf104.xml

