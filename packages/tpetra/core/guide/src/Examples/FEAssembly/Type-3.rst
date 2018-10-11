.. _fem_assembly_type3:

Type 3 - Total Element Loop / Element Loop with Ghost Elements
##############################################################

This method constructs the graph without having the need for communication
because it already has a copy of its ghost elements.  This means that every
element knows enough about its neighbors so they don't have to communicate to
get the neighbors' contributions to the owned nodes.


