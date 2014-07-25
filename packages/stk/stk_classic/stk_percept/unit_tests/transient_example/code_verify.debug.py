#!/usr/netpub/python-2.7/bin/python

from PerceptMesh import *
from math import *

import argparse

"""Entry point for Percept script."""
parser = argparse.ArgumentParser()

parser.add_argument("-nf", "--nodal_field",       nargs=1,                  help="name of nodal field")
parser.add_argument("-sf", "--exact_soln_string", nargs=1,                  help="name of nodal field")
parser.add_argument("results_files",              nargs=argparse.REMAINDER, help="list of Exodus results files")

args = parser.parse_args()

# A) user inputs 

# INPUT: list of exodus output files
mesh_files=args.results_files

#mesh_files=["output1.gold.e"]
num_meshes=len(mesh_files)

nodal_fields=args.nodal_field

#eps = 0.05
#sf_eps = StringFunction(str(eps), "eps", Dimensions(2), Dimensions(1));

#p0 = 35651.28116
#sf_p0 = StringFunction(str(p0), "p0", Dimensions(2), Dimensions(1));

#exact_soln_string="p0*(1.0 - eps*sin(PI*x)*cos(PI*y))"
exact_soln_string="35651.28116*(1.0 - 0.05*sin(PI*x)*cos(PI*y))"

exact_soln_name = [nodal_fields[0]+"_ex"]
sf_ex = [StringFunction(exact_soln_string, exact_soln_name[0], Dimensions(2), Dimensions(1))]

num_norms=2
errors=[0]*num_meshes*num_norms

expected_rates=[2,2]

tolerances=[1e-2]*num_norms

dofs_are_elems=True;
dofs=[0]*num_meshes

rates=[0]*(num_meshes-1)*num_norms
print "rates= ", rates

#####################################################################################

# B) loop over meshes and compute errors

for i in range(0,num_meshes):
    pMesh = PerceptMesh(2)
    pMesh.open(mesh_files[i])
    pMesh.commit()
    print "mesh_files[i]= " , i, mesh_files[i]

    spatial_dim=pMesh.get_spatial_dim()
    # TODO: check that this matches the StringFunction exact solution

    metaData = pMesh.get_fem_meta_data()
    bulkData = pMesh.get_bulk_data()
    
    nodal_field = metaData.get_field(nodal_fields[0])
    ff_Tnd = FieldFunction(nodal_fields[0], nodal_field, bulkData, Dimensions(spatial_dim), Dimensions(1))

    error_string = [exact_soln_name[0]+" - "+nodal_fields[0]];
    error_name = [nodal_fields[0]+"_err"]
    print "error_string= ", error_string
    sf_Terr = StringFunction(error_string[0], error_name[0], Dimensions(spatial_dim), Dimensions(1))

    numSteps = pMesh.get_database_time_step_count()
    print "numSteps= ", numSteps, " nodal_fields[0]= ", nodal_fields[0]
    pMesh.read_database_at_step(numSteps)


    p_field = pMesh.get_field("P")
    node = pMesh.get_node(1)
    p1 = pMesh.get_field_data(p_field, node)
    print "p1= " , p1

    p21 = eval_func2(2,1,0,ff_Tnd)
    p21_ex = eval_func2(2,1,0,sf_ex[0])
    err =  eval_func2(2,1,0,sf_Terr)
    print "p21 = ", p21, p21_ex, err

    # DEBUG
    print numSteps, pMesh.get_current_database_step(), pMesh.get_current_database_time()

    cubDegree = 2
    lInfNorm = LInfNorm(bulkData)
    lInfNorm.setCubDegree(cubDegree)

    l2Norm = L2Norm(bulkData)
    l2Norm.setCubDegree(cubDegree)

    errors[i*2]=l2Norm.evaluate(sf_Terr)
    errors[i*2+1]=lInfNorm.evaluate(sf_Terr)

    if dofs_are_elems:
        dofs[i] = float(pMesh.get_number_elements())
    else:
        dofs[i] = float(pMesh.get_number_nodes())
        
    print "done: ", numSteps, dofs[i], errors[i*2], errors[i*2+1]

for i in range(0,num_meshes-1):
    mesh_ratio=pow(dofs[i]/dofs[i+1],-1.0/spatial_dim)
    for n in range(0,2):
        rates[i*2+n]=log(errors[i*2+n]/errors[(i+1)*2+n])/log(mesh_ratio)
        print rates[i*2+n], rates[i*2+n]-expected_rates[n], tolerances[n]

# OUTPUT:
#   pass: if rates are within tolerance of expected rates
#   fail: otherwise

