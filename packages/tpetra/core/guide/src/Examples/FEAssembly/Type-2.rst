.. _fem_assembly_type2:

Type 2 - Export Owned Graph to Shared Graph
###########################################

An optimization of type-1, here we separate the 'owned' nodes from the 'shared' 
elements and and communicate just the needed information across processors directly.
Owned nodes go into the `owned graph` and Shared nodes go into the `shared graph`



