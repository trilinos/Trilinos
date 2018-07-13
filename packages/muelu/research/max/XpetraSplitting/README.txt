######################################################################################################
################################## CONTACT INFO ######################################################
######################################################################################################

Authors: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
				 Raymond S. Tuminaro (rstumin@sandia.gov)

######################################################################################################
################################### COMPILING ########################################################
######################################################################################################

In order to compile this directory, insert the following line in your configure script to build Trilinos:

  -D MueLu_ENABLE_REGION_SPLITTING:BOOL=ON \

#######################################################################################################
########################## INPUT REQUIREMENTS #########################################################
#######################################################################################################	

Code basically reads two files and produces region matrices that it can print.
One of the files contains information on how global node ids are distributed
across regions. This file format is  

          %Global Number of Nodes             %Global Number of Regions
          6                                   2
          Node                                Region
          1                                   1
          2                                   1
          3                                   1
          4                                   1
          3                                   2
          4                                   2
          5                                   2
          6                                   2

In this case, this corresponds to 

          1---------3----------5
          |         |          |
          |  region |  region  |
          |    1    |    2     |
          |         |          |
          2---------4----------6


The second file is a matrix market file describing the global matrix.

#######################################################################################################
#################################### TOY PROBLEMS CREATION ############################################
#######################################################################################################

Matlab routines are provided in MueLu/research/max/XpetraSplitting to construct toy problems.
These can be used to test the performance of the code.
The rouintes essentially allow to build multiregional domains and construct a discretized Laplace operator on it.
The files are:
			
						create_problem.m 
						create_matrix.m
						create_regions.m

MueLu/research/max/MatrixMarket_routines: contains support routines to read/write matrices in MatrixMarket format.
The main file create_problem.m generates two files:

					A_multiregional.mm
					node_multiregional.mm

More documentation is provided in the .m files themselves.

#######################################################################################################
########################################### CODE EXECUTION ############################################
#######################################################################################################

These files A_multiregional.mm and node_multiregional.mm must be moved to the corresponding build version of MueLu/research/max/XpetraSplitting.
To run the code, do something like to following: 

mpirun -np number_processes ./MueLu_Test_xpetra.exe file.xml A_multiregional.mm node_multiregional.txt

#######################################################################################################
################################## CLASSES DESCRIPTION ################################################
#######################################################################################################


The code does this by creating an 'Xpetra_MatrixSplitting' which extends
an Xpetra::Matrix. Buried inside this class, the main data members are

      regionHandler_, compositeMatrixData_, regionMatrixData_

compositeMatrixData_ is pretty much the input matrix. regionMatrixData_ is essentially the desired output, i.e. the region matrices. More on regionHandler_ below.

Most of the Xpetra::Matrix methods operated on compositeMatrixData_ as opposed
to regionMatrixData_. We would need to look through each of these individual
methods to see if we need a version that works on regionMatrixData_ or not
(e.g., apply()). It might be necessary (or useful) to have 2 versions of 
some methods: one appropriate for compositeMatrixData_ and the other appropriate
for regionMatrixData_. 

The Xpetra_MatrixSplitting constructor basically does the following:
     1) creates a regionHandler_
     2) pulls things out of regionHandler_
     3) makes a map suitable for the composite matrix using the RowMap buried 
        in regionHandler_ (actually an array). could just make map inside of regionHandler_?
     4) reads MatrixMarket file shoving composite matrix into compositeMatrixData_
     5) builds the region matrices using the region maps (actually arrays)
        from regionHandler_. These region matrices are really tpetra that has to be 
        wrapped to xpetra and then put into the array regionMatrixData_.
        For each region, it seems like we need to get the list of 
        region global indices, find the corresponding map entry for the 
        composite matrix, map this to a local id, grab this matrix row from the 
        overlapped composite matrix, map these local col ids to global col ids,
        and finally map these to region global ids so that we can do an
        insertGlobalValues().

regionHandler_ has the following data members:

    nodes_                array of tuples storing (global composite id, region id)
    procs_per_region_     array of tuples storing (region id, list of proc ids
                                                   assigned to this region)
    nodesToRegion_        array of tuples storing (composite id, list of regions
                                                   containing this composite id)
    interfaceNodes_       array of tuples storing (composite id, list of regions
                                                   contianing this comp. id). 
                          This is a subset of nodesToRegion_
    regionToAll_          array of array of tuples. Each region has a list of
                          nodes where each node is a tuple: (region, comp id.).
    composite_map_        An array of gids assigned to me for composite matrix
    region_maps_          An array of arrays. For each region, this gives a 
                          list of my gids assigned to region.

##############################################################################################################################################

The class Xpetra::RegionAMG inherits from Xpetra::Operator and aims to operate as a preconditioner for a Belos-type linear solver.
RegionAMG essentially owns an Xpetra::MatrixSplitting object as a private variable. The functionality of RegionAMG is to 
use the stored Xpetra::MatrixSplitting to implement a new version of the multigrid V-cycle. The new version of V-cycle aims 
to confine as much computation as possible upon separate regions. This would allow one to minimize the 
communication across processes associated with separate reigons. 

Xpetra::RegionAMG stores Xpetra::Level objects. Xpetra::Level is currently not used in the code, but its goal is to 
store quantities needed at each level of the multigrid hierarchy without necessarily relying on MueLu::Hierarchy objects.
Quantities stored in Xpetra::Level should be needed when Xpetra::RegionAMG::apply() reaches its final stage. In fact, this method should run
independently of MueLu::Hierarchy objects.

The main method of this class is the public method apply(). 
N.B.: The implementation of this method is started but NOT finished yet.

The final goal is to have the apply() method do the same region-wise V-cycle as it is currently implemented in Matlab
Peter's code. 
As for now, the apply method executes a slightly different task. What is already implemented is the 
splitting of the composite INPUT multivector into region views, as well as the reassimbling of the composite OUTPUT multivector
from the region views of it.

The sequence of instructions run in the current Xpetra::RegionAMG::apply() is the following:

1) splitting of the composite INPUT multivector into region views
2) calculation od the region multigrid V-cycle (no interactions between different regions is accounted now, but it is needed to make the code do what we want)
3) reassimbling of the composite OUTPUT multivector from the region views of it
	 (entries associated with shared mesh nodes are ALREADY weighed)

Step 1) is accomplished via the private method computeRegionX().
Step 2 is executed via the MueLu::Hierarchy::Iterate() method for each region.
Step 3) is executed via the private method computeCompositeY() and the rescaling of shared entries is done with rescaleInterfaceEntries().

#################################################################################################################################
################################## STATE-OF-THE-ART AND FUTURE DEVELOPMENTS #####################################################
#################################################################################################################################


              fine level input							fine level output
							(already taken care of)				(already taken care of)

												\										/

													o								o    (here regions do not currently communicate, but they should on a similar fashion as "sharedSum" in Peter's code)

														\						/

														...				...   (here regions do not currently communicate, but they should on a similar fashion as "sharedSum" in Peter's code)    

																\		/

																	o   (here regions do not currently communicate, but they should on a similar fashion as "sharedSum" in Peter's code)

				
ALREADY DONE: communication between regions is already in place at the fine level. This is done with Import/Export routines
							that transfer info from composite-overlapping maps to composite-nonoverlapping maps or the other way around


THINGS TO DO: the same kind of communications across regions must be carried out at each coarsening level, so that 
							values attained by region multivectors at shared mesh nodes are properly averaged.
							To this goal, intermediate quantities needed are alreay calculated/stored but not at the desired final stage.
							For instance, the structure regionToAll is already injected at each level and stored in a Xpetra::Level object. 
							Thsi should be used to identify entries of region multivectors associated with corse mesh nodes shared by multiplve regions.
							

           





