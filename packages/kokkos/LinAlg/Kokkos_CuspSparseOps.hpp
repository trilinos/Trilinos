/*

notes:
* multiply requires encapsulation of multivector using make_array2d_view
* crs construction requires encapsulation of CSR matrix data via make_csr_matrix_view
  this object can be used for multiplication or conversion to another format.
* cusp::transpose() and cusp::convert() will come into play
* we may need to add numCols to our graph interface

*/
