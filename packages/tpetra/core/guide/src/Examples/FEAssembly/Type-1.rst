.. _fem_assembly_type1:

Type 1 - Owned Element Loop / FECrs{Graph,Matrix} Construction 
##############################################################

This method constructs the graph by looping over the OWNED elements 
and inserting the clique made of up the nodes associated with each 
of these elements into the graph using their global ids without regard 
to ownership of the nodes. FillComplete will communicate and sort out
ownership.

