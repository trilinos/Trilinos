The design directory contains some ideas and tests we should keep for
informational purposes.

The ViewOfView files show three different designs. The first two were
too restrictive due to added template parameters for the inner view to
force it to be unmanaged. If codes wanted to grab an inner view and
pass through functions, then they needed to have an exact match on the
template arguments or template their functions on the view type. This
was difficult for some codes to support. The final version constructs
unmanaged inner views at runtime by passing a pointer. This disables
padding so might not be a performance hit, but we have not really seen
this.

The Darma serialization file contains the code to serialize a
view-of-views using the Darma checkpoint library (maestro). Until
Darma is added to Trilinos as a TPL, we can't add this code to
trilinos to test against.
