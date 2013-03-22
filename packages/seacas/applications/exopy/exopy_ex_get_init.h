#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"

static PyObject * exopy_ex_get_init(PyObject *self, PyObject *args) {
  /* (title,num_dim,num_nodes,num_elem,num_elem_blk,num_node_sets,num_side_sets) = ex_get_init(exoid); */

  int exoid;
  char title[MAX_LINE_LENGTH+1];

  if ( !PyArg_ParseTuple(args, "i", &exoid) ) {
    return NULL;
  }

  int num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  int error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets);

  if ( error < 0 ) {
    PyErr_SetString(PyExc_RuntimeError, "error in exopy_ex_get_init()");
    return NULL;
  }

  return Py_BuildValue("(siiiiii)", title,num_dim,num_nodes,num_elem,num_elem_blk,num_node_sets,num_side_sets);
}
