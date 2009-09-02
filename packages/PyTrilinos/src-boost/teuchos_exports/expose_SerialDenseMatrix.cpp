
#include <boost/python.hpp>
using namespace boost::python;

#include <Teuchos_SerialDenseMatrix.hpp>
#include <numpy/arrayobject.h>


template<typename OrdinalType, typename ScalarType,char typekind>
class expose_sdmatrix {
public:
    typedef Teuchos::SerialDenseMatrix< OrdinalType, ScalarType > sdm;
    static void wrap(char * name)
    {
        class_<sdm>(name,init<OrdinalType,OrdinalType>())
            .def("numRows",&sdm::numRows)
            .def("numCols",&sdm::numCols)
		    .add_property("__array_struct__", array_struct__)
		    .def("__setitem__",&setitem)
		    .def("__getitem__",&getitem)
		;
		
		
	}
	
    static ScalarType getitem(sdm& self,object tup)
    {
        int rowIndex = extract<int>(tup[0]);
        int colIndex = extract<int>(tup[1]);
        
        return self(rowIndex,colIndex);
        
    }

    static void setitem(sdm& self,object tup,ScalarType val)
    {
        int rowIndex = extract<int>(tup[0]);
        int colIndex = extract<int>(tup[1]);
        
        self(rowIndex,colIndex) = val;
        
    }

    
    static void cleanup(void* obj)

    {

        PyArrayInterface* ai = (PyArrayInterface*) obj;

        delete [] ai->shape;

        delete [] ai->strides;

        delete ai;

    }

    static PyObject* array_struct__(sdm const& self)

    {

        // http://numpy.scipy.org/array_interface.shtml

        PyArrayInterface* ai = new PyArrayInterface;



        // contains the integer 2 as a sanity check

        ai->two = 2;

        // number of dimensions: two for a matrix

        ai->nd = 2;

        // kind in array --- character code of typestr

        ai->typekind = typekind;

        // size of each element

        ai->itemsize = sizeof(ScalarType);



        // how should be data interpreted

        ai->flags = NPY_CONTIGUOUS | NPY_ALIGNED | NPY_NOTSWAPPED | NPY_WRITEABLE;



        // A length-nd array of shape information

        ai->shape = new npy_intp[2];

        ai->shape[1] = self.numRows();

        ai->shape[0] = self.numCols();



        // A length-nd array of stride information

        ai->strides = new npy_intp[2];

        ai->strides[0] = ai->shape[1] * ai->itemsize;

        ai->strides[1] = ai->itemsize;



        // A pointer to the first element of the array

        ai->data = (void*) self.values();



        ai->descr = NULL;



        return PyCObject_FromVoidPtr(ai, cleanup);

    }
    
};


