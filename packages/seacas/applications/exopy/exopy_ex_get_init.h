#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"
#include "exopy_ref.h"

static PyObject * exopy_ex_get_init(PyObject *self, PyObject *args) {
  /* error = ex_get_init (exoid, title, &num_dim, &num_nodes,
                          &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets); */
  int exoid;
  char title[MAX_LINE_LENGTH+1];

  ref *title_ref, *num_dim_ref, *num_nodes_ref, *num_elem_ref, *num_elem_blk_ref, *num_node_sets_ref, *num_side_sets_ref;

  if ( !PyArg_ParseTuple(args, "iOOOOOOO", &exoid, &title_ref,
                         &num_dim_ref, &num_nodes_ref, &num_elem_ref,
                         &num_elem_blk_ref, &num_node_sets_ref, &num_side_sets_ref) ) {
    return NULL;
  }

  int num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  int error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets);

  title_ref->value = Py_BuildValue("s",title);
  num_dim_ref->value = PyInt_FromLong(num_dim);
  num_nodes_ref->value   = PyInt_FromLong(num_nodes);
  num_elem_ref->value   = PyInt_FromLong(num_elem);
  num_elem_blk_ref->value   = PyInt_FromLong(num_elem_blk);
  num_node_sets_ref->value   = PyInt_FromLong(num_node_sets);
  num_side_sets_ref->value   = PyInt_FromLong(num_side_sets);

  // Do this so error is raised pointing to this function
  if ( PyErr_Occurred() ) { return NULL; }

  return Py_BuildValue("i", error);
}
