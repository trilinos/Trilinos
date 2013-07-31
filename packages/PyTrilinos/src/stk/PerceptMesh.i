// -*- c++ -*-

%define %PerceptMesh_docstring 
"""
PyPercept is the Python interface for the Trilinos STK Percept package.

The following are examples of the actual PyPercept Python interface:

##################################################################################################################

Case 1: Norm of the difference of an Exact Solution and a Mesh Field

pMesh = PerceptMesh() #Create a PerceptMesh database
#Specify the geometric hex mesh to be a 3x3x3 box between 0,0,0 and 2,2,2
pMesh.new_mesh(GMeshSpec('3x3x3|bbox:0,0,0,2,2,2')) 
field = pMesh.add_field('coordinates', 1) #Add a coordinates field
pMesh.commit() #Commits the mesh to the mesh database

input_array = numpy.array([1.0, 0.5, 0.5]) #Create numpy arrays to be evaluated 
input_array_2 = numpy.array([1.0, 1.5, 1.5])
'
#Create a field function with a name, field base, and the mesh and the domain and co-domain dimensions
ff = FieldFunction('ff', field, pMesh, 3, 3) 
ff.addAlias('myalias') #Define another name for the field function
ff_output = ff.evaluate(input_array) #Evaluate the field function with the given array

f2 = FieldFunction('f2', field, pMesh, 3, 3)
f2_output = f2.evaluate(input_array_2)

#Create a function from a string 'x+y+z', with a name, the domain and co-domain dimensions
sf = StringFunction('x+y+z', 'myname', 3, 1) 
sf_output = sf.evaluate(input_array) #Evaluate the function with the given array

sf_diff = StringFunction('ff-f2', 'myname') #Create a function of the difference of the field functions

norm = L1Norm(pMesh.getBulkData()) #Create a norm object to evaluate the norm of the mesh field
value = norm.evaluate(ff) #Evaluate the norm
diffnorm = norm.evaluate(sf_diff) #Evaluate the norm of the difference of the exact solution and the mesh field

#Now use a helper function to evaluate the norm: eval_norm(bulkData, Function, int) 
value1 = eval_norm(pMesh.getBulkData(), ff, 1) 
diffnorm1 = eval_norm(pMesh.getBulkData(), sf_diff, 1)

#################################################################################################################

Case 2: Refine a mesh uniformly

pMesh = PerceptMesh()
pMesh.open('tet-mesh-brick-no-block.e') #Read an existing exodus file

uniform_refiner = Refiner(pMesh, TET4_TET4_8) #Create a refiner on the specified mesh of the specified pattern.
pMesh.commit()

i = 0
while i < 3:
  uniform_refiner.doBreak() #Call the doBreak() method to refine the mesh
  i = i + 1

pMesh.saveAs('tet-mesh-refined-3-times.e') #Save out the refined mesh

#################################################################################################################

Case 3: Norm of the difference of results from two separate simulations

subprocess.call('sierra aria -i tet_mesh.e -o result_0.e') #Run an external program from Python

pMesh = PerceptMesh()
pMesh.open('tet-mesh.e')
uniform_refiner = Refiner(pMesh, TET4_TET4_8)
pMesh.commit()

uniform_refiner.doBreak()
pMesh.saveAs('tet-mesh_refined.e')

subprocess.call('sierra aria -i tet_mesh_refined.e -o result_1.e')

pMesh_0 = PerceptMesh()
pMesh_1 = PerceptMesh()
pMesh_0.openReadOnly('result_0.e')
pMesh_1.openReadOnly('result_1.e')

ff_0 = Field_Function(pMesh_0)
ff_1 = Field_Function(pMesh_1)
diff = StringFunction('ff_0 - ff_1')

diffnorm = eval_norm(pMesh.getBulkData, diff, 2)
print 'diffnorm = ', diffnorm

#################################################################################################################
"""
%enddef

%module(package   = "PyPercept",
        autodoc   = "1",
        docstring = %PerceptMesh_docstring) PerceptMesh

%{
  
//from Teuchos
#include "Teuchos_ENull.hpp"
#include "Teuchos_RCPDecl.hpp"	
#include "Teuchos_RCP.hpp"
  
//from stk_adapt
#include "stk_adapt/UniformRefinerPattern.hpp"
#include "stk_adapt/Refiner.hpp"
  
//from stk_percept	
  
#include "stk_percept/norm/Norm.hpp"

#include "stk_percept/fixtures/QuadFixture.hpp"
#include "stk_percept/fixtures/WedgeFixture.hpp"
#include "stk_percept/fixtures/BeamFixture.hpp"
#include "stk_percept/fixtures/HeterogeneousFixture.hpp"


#define SWIG_FILE_WITH_INIT

//namespaces used
using namespace stk;
using namespace stk::diag;
using namespace stk::percept;
using namespace stk::mesh;
using namespace stk::adapt;
using namespace Teuchos;
	
%}

%apply long long { uint64_t }

#pragma SWIG nowarn=509
#pragma SWIG nowarn=467
#pragma SWIG nowarn=314

#define STK_HAS_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%include numpy.i
%init
%{
  import_array();
%}

// handle exceptions thrown from library
%exception
{ 
  try
  { 
    $action
  } 
  catch (std::runtime_error e)
  { 
    PyErr_SetString(PyExc_RuntimeError, e.what()); 
    SWIG_fail; 
  } 
  catch(...)
  {
    PyErr_SetString(PyExc_RuntimeError, "Unexpected exception"); 
    SWIG_fail; 
  }
}

%include stl.i
namespace std 
{
  %template(vectorfieldbase) vector<stk::mesh::FieldBase * >;
  %template(vectorparts) vector<stk::mesh::Part * >;
  %template(vectors) vector<string>;
  %template(vectorvs) vector<vector<string> >;
  %template(vectori) vector< int >;
  %template(vectord) vector< double >;
}

%include cpointer.i
%include typemaps.i

//General ignore directives
%ignore *::operator();
%ignore *::operator=;
%ignore *::print;

%feature("autodoc", 1);

// Include Percept documentation
%include "Percept_dox.i"

//from Teuchos
%include Teuchos_config.h
%ignore Teuchos::RCP::count;
%include Teuchos_RCPDecl.hpp
%template(RCPFunction) Teuchos::RCP<stk::percept::Function>;

// from intrepid

//ignore overloaded operators and functions not supported by SWIG
%ignore Intrepid::FieldContainer::operator [] (const int address) const;
%ignore Intrepid::FieldContainer::operator [] (const int address);
%ignore Intrepid::FieldContainer::getData() const;
%include Intrepid_FieldContainer.hpp

//from stk_util
%include Parallel.hpp
%include stk_util/diag/Writer.hpp

//from stk_io
%ignore get_cell_topology(const mesh::Part &part);
%ignore set_cell_topology(mesh::Part &part, const CellTopologyData * const cell_topology);
%include stk_io/IossBridge.hpp

//from stk_io/util

//ignore overloaded functions not supported by SWIG
%ignore stk::io::util::Gmesh_STKmesh_Fixture::getFEMMetaData() const;
%ignore stk::io::util::Gmesh_STKmesh_Fixture::getBulkData() const;
%include Gmesh_STKmesh_Fixture.hpp

//from stk_mesh/base

//ignore overloaded operators not supported by SWIG
%ignore stk::mesh::operator << ( std::ostream & , const FieldBase & );
%include stk_mesh/base/FieldBase.hpp

%include stk_mesh/base/Part.hpp

%include stk_mesh/base/Types.hpp

  //ignore overloaded operators not supported by SWIG
%ignore stk::mesh::Selector::operator ! () const;
%ignore stk::mesh::Selector::operator << (std::ostream & out,
                                          const stk::mesh::Selector & selector);  //FIXME
%ignore stk::mesh::operator!(const Part & A);
%include stk_mesh/base/Selector.hpp

%include stk_mesh/base/BulkData.hpp

//ignore overloaded operators not supported by SWIG
%ignore stk::mesh::BucketIterator::operator++();
%ignore stk::mesh::BucketIterator::operator--();
%ignore stk::mesh::BucketIterator::operator++(int);
%ignore stk::mesh::BucketIterator::operator--(int);
%ignore stk::mesh::BucketIterator::operator[]( const intType & n ) const;
%ignore stk::mesh::Bucket::operator[] ( size_t i ) const;
%ignore stk::mesh::operator << ( std::ostream & , const Bucket & );
%include stk_mesh/base/Bucket.hpp

//from stk_mesh/baseImpl
%include stk_mesh/baseImpl/FieldRepository.hpp
%include stk_mesh/fem/CellTopology.hpp

//from stk_mesh/fem
%include stk_mesh/fem/FEMMetaData.hpp
%template(get_field) stk::mesh::fem::FEMMetaData::get_field<stk::mesh::FieldBase>;

//from stk_percept/function
%include MDArray.hpp
%template(MDArray) Intrepid::FieldContainer<double>;

%typemap(in) stk::percept::MDArray & (stk::percept::MDArray mdarray)
{
  if (!PyArray_Check($input))
  {
    PyErr_SetString(PyExc_TypeError, "ERROR: input must be numpy array");
    SWIG_fail;
  }
  PyArrayObject* array = (PyArrayObject*) $input;
  if (PyArray_TYPE(array) != NPY_DOUBLE)
  {
    PyErr_SetString(PyExc_TypeError, "ERROR: input must be numpy array of type double");
    SWIG_fail;
  }
  int nd = PyArray_NDIM(array);
  if (nd > 5)
  {
    PyErr_SetString(PyExc_TypeError,
                    "ERROR: array dimensions must be less than or equal to 5");
    SWIG_fail;
  }
  if (nd == 1)
  {
    mdarray.resize(PyArray_DIMS(array)[0]);
    mdarray.setValues((double*) PyArray_DATA(array),
                      PyArray_DIMS(array)[0]*sizeof(double));
    $1 = &mdarray;
  }
  else if (nd == 2)
  {
    int dim1 = PyArray_DIMS(array)[0];
    int dim2 = PyArray_DIMS(array)[1];
    mdarray.resize(dim1, dim2);
    mdarray.setValues((double*) PyArray_DATA(array), dim1*dim2*sizeof(double));
    $1 = &mdarray;
  }
  else if (nd == 3)
  {
    int dim1 = PyArray_DIMS(array)[0];
    int dim2 = PyArray_DIMS(array)[1];
    int dim3 = PyArray_DIMS(array)[2];
    mdarray.resize(dim1, dim2, dim3);
    mdarray.setValues((double*) PyArray_DATA(array), dim1*dim2*dim3*sizeof(double));
    $1 = &mdarray;
  }
}

%typemap(typecheck) stk::percept::MDArray&
{
  void *vptr = 0;
  int res1 = SWIG_ConvertPtr($input, &vptr, $*descriptor, 0);
  int res2 = PyArray_Check($input);
  $1 = (SWIG_CheckState(res1) || res2);
}

%typemap(in, numinputs=0) stk::percept::MDArray& out (stk::percept::MDArray out)
{
  $1 = &out;
}

%typemap(out) stk::percept::MDArray
{
  int numdims = $1.rank();
  int dimensions[5];
  for (int i = 0; i < numdims; i++)
  {
    dimensions[i] = $1.dimension(i);  
  }  
  PyArrayObject* temp_array =
    (PyArrayObject*) PyArray_FromDims(numdims, dimensions, NPY_DOUBLE);
  double* buffer = (double*) PyArray_DATA(temp_array);
  Teuchos::ArrayRCP<double> rcp = $1.getData();
  double* md_temp = rcp.access_private_ptr(); 
  for (int i=0; i < $1.size(); i++) 
  {
    buffer[i] = md_temp[i];
  }
  $result = SWIG_Python_AppendOutput($result, (PyObject*) temp_array);
}

%module outarg

%typemap(argout) stk::percept::MDArray &MDOutVal
{
  int numdims = $1->rank();
  int dimensions[5];
  for (int i = 0; i < numdims; i++)
  {
    dimensions[i] = $1->dimension(i);  
  }  
  PyArrayObject* temp_array =
    (PyArrayObject*) PyArray_FromDims(numdims, dimensions, NPY_DOUBLE);
  double* buffer = (double*) PyArray_DATA(temp_array);
  Teuchos::ArrayRCP<double> rcp = $1->getData();
  double* md_temp = rcp.access_private_ptr(); 
  for (int i=0; i < $1->size(); i++) 
  {
    buffer[i] = md_temp[i];
  }
  $result = SWIG_Python_AppendOutput($result, (PyObject*)  temp_array);
  //$result =  (PyObject*) temp_array;
}

void stk::percept::Function::value(MDArray& domain,
                                   MDArray& MDOutVal,
                                   double time_value_optional=0.0);

%typemap(in) stk::percept::MDArrayString & (stk::percept::MDArrayString mdarray)
{
  if (!PyArray_Check($input))
  {
    PyErr_SetString(PyExc_TypeError, "ERROR: input must be numpy array");
    SWIG_fail;
  }
  PyArrayObject* array = (PyArrayObject*) $input;
  if (PyArray_TYPE(array) != NPY_STRING)
  {
    PyErr_SetString(PyExc_TypeError, "ERROR: input must be numpy array of type string");
    SWIG_fail;
  }
  int nd = PyArray_NDIM(array);
  bool debug = false;
  if(debug) std::cout << "nd= " << nd <<  std::endl;
  // FIXME
  if (nd > 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "ERROR: array dimensions must be less than or equal to 2");
    SWIG_fail;
  }
  if (nd == 1)
  {
    int dim1 = PyArray_DIMS(array)[0];
    if(debug) std::cout << "dim1= " << PyArray_DIMS(array)[0]
                        <<  std::endl;
    mdarray.resize(PyArray_DIMS(array)[0]);
    int itemsize = PyArray_ITEMSIZE(array);
    if(debug) std::cout << "dim1= " << PyArray_DIMS(array)[0]
                        <<  std::endl;
    //		mdarray.setValues((std::string*) PyArray_DATA(array) );
    MDArrayString tmp(dim1);
    char buf[itemsize+1];
    for (int i1=0; i1 < dim1; i1++)
    {
      // void* PyArray_GETPTR2(PyObject* obj, <npy_intp> i, <npy_intp> j)
      char *str = (char *)PyArray_GETPTR1(array, i1);
      for (int ii=0; ii < itemsize; ii++) buf[ii] = str[ii];
      buf[itemsize] = '\0';
      if(debug) std::cout << "str= " << buf
                          <<     std::endl;
      mdarray(i1) = buf;
    }
    if(debug) std::cout << "dim1= " << PyArray_DIMS(array)[0]
                        <<  std::endl;
    $1 = &mdarray;
  }
  else if (nd == 2)
  {
    int dim1 = PyArray_DIMS(array)[0];
    int dim2 = PyArray_DIMS(array)[1];
    if(debug) std::cout << "dim1,2= " << dim1 << " " << dim2
                        <<  std::endl;
    mdarray.resize(dim1, dim2);
    int itemsize = PyArray_ITEMSIZE(array);
    if(debug) std::cout << "dim1,2= " << dim1 << " " << dim2
                        << " itemsize = " << itemsize
                        <<  std::endl;
    
    MDArrayString tmp(dim1,dim2);
    char buf[itemsize+1];
    for (int i1=0; i1 < dim1; i1++)
      for (int i2=0; i2 < dim2; i2++)
      {
        // void* PyArray_GETPTR2(PyObject* obj, <npy_intp> i, <npy_intp> j)
        char *str = (char *)PyArray_GETPTR2(array, i1, i2);
        for (int ii=0; ii < itemsize; ii++) buf[ii] = str[ii];
        buf[itemsize] = '\0';
        if(debug) std::cout << "str= " << buf
                            <<     std::endl;
        mdarray(i1,i2) = buf;
      }
    //		mdarray.setValues((std::string*) PyArray_DATA(array) );
    if(debug) std::cout << "dim1,2= " << dim1 << " " << dim2
                        <<  std::endl;
    $1 = &mdarray;
  }
/*
  else if (nd == 3)
  {
    int dim1 = PyArray_DIMS(array)[0];
    int dim2 = PyArray_DIMS(array)[1];
    int dim3 = PyArray_DIMS(array)[2];
    mdarray.resize(dim1, dim2, dim3);
    mdarray.setValues((std::string*) PyArray_DATA(array) );
    $1 = &mdarray;
  }
*/
}

%typemap(typecheck) stk::percept::MDArrayString&
{
  void *vptr = 0;
  int res1 = SWIG_ConvertPtr($input, &vptr, $*descriptor, 0);
  int res2 = PyArray_Check($input);
  $1 = (SWIG_CheckState(res1) || res2);
}

// from stk_percept

%include stk_percept/Name.hpp
  //Allow user to pass in a string instead of a Name object
%typemap(in) stk::percept::Name (std::string name)
{
  if (!PyString_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "ERROR: name must be a string");
    SWIG_fail;
  }
  name = PyString_AsString($input);
  stk::percept::Name new_name(name);
  $1 = new_name;
}

%typemap(typecheck) stk::percept::Name
{
  int res = PyString_Check($input);
  $1 = res;
}

%feature("docstring") PerceptMesh "Create a PerceptMesh.";
%feature("docstring") GMeshSpec "Specify the mesh to be created by passing in a string, e.g. '3x3x3|bbox:0,0,0,1,1,1'.";
%feature("docstring") new_mesh "Create a new mesh object by passing in a GMeshSpec.";


%include PerceptMesh.hpp

//from stk_percept/function/internal
%include Dimensions.hpp

//from stk_percept/fixture
%feature("docstring") generate_mesh "This method generates the mesh.";
%include QuadFixture.hpp

%feature("immutable","1") stk::percept::QuadFixture<double>::bulk_data;
%feature("immutable","1") stk::percept::QuadFixture<double>::coord_field;
%feature("immutable","1") stk::percept::QuadFixture<double>::coord_gather_field;
%feature("immutable","1") stk::percept::QuadFixture<double>::meta_data;
%feature("immutable","1") stk::percept::QuadFixture<double>::quad_part;
%feature("docstring") QuadFixture<double, shards::Quadrilateral<4> > "This class creates a quad mesh.";
%template(QuadFixture_4) stk::percept::QuadFixture<double, shards::Quadrilateral<4> >;

%feature("immutable","1") stk::percept::QuadFixture<double, shards::Triangle<3> >::bulk_data;
%feature("immutable","1") stk::percept::QuadFixture<double, shards::Triangle<3> >::coord_field;
%feature("immutable","1") stk::percept::QuadFixture<double, shards::Triangle<3> >::coord_gather_field;
%feature("immutable","1") stk::percept::QuadFixture<double, shards::Triangle<3> >::meta_data;
%feature("immutable","1") stk::percept::QuadFixture<double, shards::Triangle<3> >::quad_part;
%feature("docstring") QuadFixture<double, shards::Triangle<3> > "This class creates a tri mesh.";
%template(QuadFixture_3) stk::percept::QuadFixture<double, shards::Triangle<3> >;

%include WedgeFixture.hpp

%feature("immutable","1") stk::percept::BeamFixture::m_bulkData;
%feature("immutable","1") stk::percept::BeamFixture::m_coordinates_field;
%feature("immutable","1") stk::percept::BeamFixture::m_metaData;
%feature("immutable","1") stk::percept::BeamFixture::m_block_beam;
%feature("immutable","1") stk::percept::BeamFixture::m_centroid_field;
%feature("immutable","1") stk::percept::BeamFixture::m_temperature_field;
%feature("immutable","1") stk::percept::BeamFixture::m_volume_field;
%feature("immutable","1") stk::percept::BeamFixture::m_element_node_coordinates_field;
%include BeamFixture.hpp

%feature("immutable","1") stk::percept::HeterogeneousFixture::m_bulkData;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_metaData;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_block_hex;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_block_wedge;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_block_tet;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_block_pyramid;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_coordinates_field;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_centroid_field;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_temperature_field;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_volume_field;
%feature("immutable","1") stk::percept::HeterogeneousFixture::m_element_node_coordinates_field;
%include HeterogeneousFixture.hpp

// from stk_percept/function  (depends on stk_percept)

//ignore overloaded function not supported by SWIG
%ignore stk::percept::GenericFunction::getNewCodomain() const;
%include GenericFunction.hpp

//eval is a Python command so we need to rename it to eval_func
%ignore stk::percept::operator<<(std::ostream& out,  Function& func);

%rename(eval_func) *::eval(double x, double y, double z, double t, Function& func);
%rename(eval_func) *::eval(double x, double y, double z, double t, Teuchos::RCP<Function>& func);
%rename(eval_func2) *::eval2(double x, double y, double t, Function& func);
%rename(eval_func2) *::eval2(double x, double y, double t, Teuchos::RCP<Function>& func);
%feature("docstring") eval_func "Evaluate the function by passing in values for x, y, z, t, and the Function. Returns a double.";
%include Function.hpp


%feature("docstring") FieldFunction "Create a field function by passing in a name, a field base, a PerceptMesh and input and output dimensions.";
%feature("docstring") evaluate "Evaluate the function with the given array.  Returns an array."
%include FieldFunction.hpp
 
%feature("notabstract") ConstantFunction;
%feature("notabstract") ConstantFunctionVec;
%include ConstantFunction.hpp

%feature("docstring") StringFunction "Create a function by passing in a function string, a name, and input and output dimensions";
%feature("docstring") evaluate "Evaluate the function with the given array.  Returns an array."
%include StringFunction.hpp
%extend stk::percept::StringFunction
{
  StringFunction __add__(StringFunction *rhs)
  {
    std::string newFuncString = "(" + $self->getFunctionString() + ") + (" +
      rhs->getFunctionString() + ")";
    return StringFunction(newFuncString.c_str());
  }
  StringFunction __div__(StringFunction *rhs)
  {
    std::string newFuncString = "(" + $self->getFunctionString() + ") / (" +
      rhs->getFunctionString() + ")";
    return StringFunction(newFuncString.c_str());
  }
  StringFunction __mul__(StringFunction *rhs)
  {
    std::string newFuncString = "(" + $self->getFunctionString() + ") * (" +
      rhs->getFunctionString() + ")";
    return StringFunction(newFuncString.c_str());
  }
  StringFunction __sub__(StringFunction *rhs)
  {
    std::string newFuncString = "(" + $self->getFunctionString() + ") - (" +
      rhs->getFunctionString() + ")";
    return StringFunction(newFuncString.c_str());
  }
};

//from stk_percept/norm
//%feature("docstring") eval_norm "Evaluate the norm by passing in a bulkData, a Function, and a 1 or 2 for L1 or L2 norm. Returns a double.";
%include Norm.hpp

%feature("docstring") Norm "Create a Norm object by passing in a bulkData.";
%feature("docstring") evaluate "Evaluate the norm by passing in a Function. Returns a double.";

%feature("docstring") Norm<int Power=2> "This class will evaluate the L2 norm of a function.";
%template(L2Norm) stk::percept::Norm<int Power=2>;

%feature("docstring") Norm<int Power=1> "This class will evaluate the L1 norm of a function.";
%template(L1Norm) stk::percept::Norm<int Power=1>;

%feature("docstring") Norm<int Power=-1> "This class will evaluate the L-inf norm of a function.";
%template(LInfNorm) stk::percept::Norm<int Power=-1>;

//from stk_adapt
%include stk_adapt/UniformRefinerPattern.hpp

%feature("docstring") Refiner "This class allows for refinement of meshes.

Create a Refiner by passing in a PerceptMesh and a Pattern.";
%feature("docstring") doBreak "Refine the mesh.";
%feature("docstring") __basePatterns__ "
The following is a list of patterns that may be used to create a Refiner object:
BEAM2_BEAM2_2    
BEAM2_BEAM3_1    
BEAM3_BEAM3_2
HEX20_HEX20_8   
HEX27_HEX27_8   
HEX8_HEX20_1
HEX8_HEX27_1   
HEX8_HEX8_8      
HEX8_TET4_24
HEX8_TET4_6_12    
LINE2_LINE2_2     
LINE3_LINE3_2
QUAD4_QUAD4_4   
QUAD4_QUAD4_4_OLD   
QUAD4_QUAD4_4_SIERRA
QUAD4_QUAD8_1   
QUAD4_QUAD9_1      
QUAD4_TRI3_2
QUAD4_TRI3_4    
QUAD4_TRI3_6        
QUAD8_QUAD8_4
QUAD9_QUAD9_4
SHELLLINE2_SHELLLINE2_2            
SHELLLINE3_SHELLLINE3_2
SHELLLINE3_SHELLLINE3_2            
SHELLQUAD4_SHELLQUAD8_1
SHELLQUAD8_SHELLQUAD8_4            
SHELLTRI3_SHELLTRI3_4
SHELLTRI6_SHELLTRI6_4
TET10_TET10_8    
TET4_TET10_1    
TET4_TET4_8
TRI3_TRI3_4      
TRI3_TRI6_1     
TRI6_TRI6_4
WEDGE15_WEDGE15_8    
WEDGE18_WEDGE18_8
WEDGE6_WEDGE15_1     
WEDGE6_WEDGE18_1
WEDGE6_WEDGE6_8";
%include stk_adapt/Refiner.hpp
%extend stk::adapt::Refiner
{
  void __basePatterns__()
  {  //here we add a dummy function for documentation purposes
  }
}
