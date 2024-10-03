#! /usr/bin/env python3

"""
Script for analyzing Panzer kernel performance on next-generation
architectures.  Runs hierarchic parallelism and generates plots from
data.
"""

__version__ = "1.0"
__author__  = "Roger Pawlowski"
__date__    = "Dec 2018"

# Import python modules for command-line options, the operating system, regular
# expressions, and system functions
import argparse
import os
import re
import sys
import datetime

#############################################################################

def main():

    """Script for analyzing Phalanx performance on next-generation architectures."""

    # Initialization
    print '****************************************'
    print '* Starting Panzer Analysis'
    print '****************************************'

    parser = argparse.ArgumentParser(description='Panzer hierarchic parallelism analysis script')
    parser.add_argument('-r', '--run', action='store_true', help='Run the executable to generate data and output to files.')
    parser.add_argument('-a', '--analyze', action='store_true', help='Analyze the data from files generated with --run.')
    parser.add_argument('-p', '--prefix', help='Add a prefix string to all output filenames.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print more data to screen.')
    parser.add_argument('-o', '--basis-order', type=int, required=True, help='FE basis order.')
    parser.add_argument('-ts', '--team-size', type=int, required=True, help='Team size for hierarchic parallelism.')
    parser.add_argument('-vs', '--vector-size', type=int, required=True, help='Vector size for hierarchic parallelism.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-b', '--base-exe', action='store_true', default=False, help="Use the Base executable")
    group.add_argument('-m', '--mixed-exe', action='store_true', default=False, help="Use the Mixed Field Type executable")
    group.add_argument('-d', '--device-dag-exe', action='store_true', default=True, help="Use the Device DAG executable")
    args = parser.parse_args()

    nx = 20
    ny = 20
    nz = 20
    order = args.basis_order
    ts = args.team_size
    vs = args.vector_size
    print "basis order = %d, team size = %d, vector size = %d\n" % (order, ts, vs)

    executable = "./PanzerAdaptersSTK_MixedPoissonExample.exe"

    print "Starting Workset Analysis"

    ws_step_size = 100
    workset_range = range(100,2000+ws_step_size,ws_step_size)
    print "workset range = "+str(workset_range)

    timings = {}
    if args.analyze:
        import numpy as np
        timings["panzer::AssemblyEngine::evaluate_volume(panzer::Traits::Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] GatherSolution (Tpetra): GRADPHI_FIELD (panzer::Traits::Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] DOFDiv: DIV_GRADPHI_FIELD (panzer::Traits::Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] Integrator_DivBasisTimesScalar (EVALUATES):  RESIDUAL_GRADPHI_FIELD"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] Sine Source"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] Integrator_DivBasisTimesScalar (CONTRIBUTES):  RESIDUAL_GRADPHI_FIELD"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] SCATTER_GRADPHI_FIELD Scatter Residual (Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] DOF: GRADPHI_FIELD accel_jac  (panzer::Traits::Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] Integrator_GradBasisDotVector (EVALUATES):  RESIDUAL_PHI_MASS_OP"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] GatherSolution (Tpetra): PHI (panzer::Traits::Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] DOFGradient: GRAD_PHI (panzer::Traits::Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] Integrator_GradBasisDotVector (EVALUATES):  RESIDUAL_PHI_DIFFUSION_OP"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] SumStatic Rank 2 Evaluator"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[panzer::Traits::Jacobian] SCATTER_PHI Scatter Residual (Jacobian)"] = np.zeros(len(workset_range),dtype=np.float64)
        

    #print dir(np)
        
    for i in range(len(workset_range)):

        ws = workset_range[i]
        
        filename = "fea_nx_%i_ny_%i_nz_%i_order_%i_ws_%i_ts_%i_vs_%i.dat" % (nx, ny, nz , order, ws, ts, vs)
        if args.prefix:
            filename = args.prefix+filename
        command = executable+" --x-elements=%i --y-elements=%i --z-elements=%i --hgrad-basis-order=%i --hdiv-basis-order=%i --workset-size=%i --use-shared-mem-for-ad --no-check-order" % (nx, ny, nz, order, order, ws)  +" >& "+filename

        if args.run:
            #print 'generating data...'
            if args.verbose:
                print "  Running \""+command+"\" ...",
                sys.stdout.flush()
            os.system(command);
            if args.verbose:
                print "completed!"
                sys.stdout.flush()

        if args.analyze:
            f = open(filename, mode='r')                
            lines = f.readlines()
            for line in lines:
                if args.verbose:
                    print line,
                for key,value in timings.iteritems():
                    if key in line:
                        split_line = line.split()
                        timings[key][i] += float(split_line[-4])
                        if args.verbose:
                            print "  found key: "+key+" = "+str(split_line[-4])
                        break
            f.close()

    do_jac = True
            
    if args.analyze:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.semilogy()
        # maroon = #990033, light blue = #00ffff
        #plt.plot(workset_range,timings["Jacobian Evaluation Time <<Host DAG>>"],label="Jac Total Time (Host DAG)",marker="o",color="#990033",markersize=8)
        #plt.plot(workset_range,timings["Jacobian Evaluation Time <<Device DAG>>"],label="Jac Total Time (Device DAG)",marker="s",color="r",markersize=8)
        #plt.plot(workset_range,timings["Residual Evaluation Time <<Host DAG>>"],label="Res Total Time (Host DAG)",marker="o",color="b",markersize=8)
        plt.plot(workset_range,timings["panzer::AssemblyEngine::evaluate_volume(panzer::Traits::Jacobian)"],label="Jacobian Volume Assembly Total Time",marker="s",color="#00ffff",markersize=8)
        plt.xlabel("Workset Size",fontsize=16)
        plt.ylabel("Time (s)",fontsize=16)
        plt.tick_params(labelsize=16)
        title = "nel=%i,order=%i" % (nx*ny*nz,order)
        plt.title(title)
        #plt.legend(bbox_to_anchor=(1,1))
        plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.0),ncol=2,fancybox=True,shadow=True, prop={'size': 12})
        plt.grid()
        dag_timings_filename = "total_time_nx_%i_ny_%i_nz_%i_order_%i_ts_%i_vs_%i.png" % (nx, ny, nz ,order, ts, vs)
        fig.savefig(dag_timings_filename)
        #plt.show()
            
        fig = plt.figure(2)
        #plt.clf()
        plt.semilogy()
        plt.plot(workset_range,timings["[panzer::Traits::Jacobian] Integrator_DivBasisTimesScalar (EVALUATES):  RESIDUAL_GRADPHI_FIELD"],label="Integrator DivBasisTimesScalar (eval)",marker='s')
        plt.plot(workset_range,timings["[panzer::Traits::Jacobian] Integrator_DivBasisTimesScalar (CONTRIBUTES):  RESIDUAL_GRADPHI_FIELD"],label="Integrator DivBasisTimesScalar (contrib)",marker='^')
        plt.plot(workset_range,timings["[panzer::Traits::Jacobian] Integrator_GradBasisDotVector (EVALUATES):  RESIDUAL_PHI_MASS_OP"],label="Integrator GradBasisDotVector (mass op)",marker='*')
        plt.plot(workset_range,timings["[panzer::Traits::Jacobian] Integrator_GradBasisDotVector (EVALUATES):  RESIDUAL_PHI_DIFFUSION_OP"],label="Integrator GradBasisDotVector (diff op)",marker='D')
        plt.plot(workset_range,timings["[panzer::Traits::Jacobian] DOF: GRADPHI_FIELD accel_jac  (panzer::Traits::Jacobian)"],label="DOF (GradPHI)",marker='+')
        plt.plot(workset_range,timings["[panzer::Traits::Jacobian] DOFGradient: GRAD_PHI (panzer::Traits::Jacobian)"],label="DOFGradient (GradPhi)",marker='x')
        #plt.plot(workset_range,timings["[panzer::Traits::Jacobian] DOFDiv: DIV_GRADPHI_FIELD (panzer::Traits::Jacobian)"],label="DOF Div (GradPhi)",marker='o')
        #plt.plot(workset_range,timings[""],label="Res Scatter",marker='.',color="#ff6600")
        plt.xlabel("Workset Size",fontsize=16)
        plt.ylabel("Time (s)",fontsize=16)
        plt.tick_params(labelsize=16)
        plt.ylim(1.0e-4,1.0e1)
        title = "nel=%i,order=%i" % (nx*ny*nz,order)
        plt.title(title)
        #plt.legend(bbox_to_anchor=(1,1))
        plt.legend(loc='upper center', bbox_to_anchor=(0.5,0.25),ncol=2,fancybox=True,shadow=True, prop={'size': 10})
        #plt.axis([0,2000,1.0e-4,0.1])
        plt.grid()
        res_evaluator_timings_filename = "kernel_timings_nx_%i_ny_%i_nz_%i_order_%i_ts_%i_vs_%i.png" % (nx, ny, nz, order, ts, vs)
        fig.savefig(res_evaluator_timings_filename)

        #print dir(plt)

        # Plot to assess savings
        count = 0;
        for key,value in timings.iteritems():
            filename_f = "raw_data_output_timer_%i_nx_%i_ny_%i_nz_%i_order_%i_ws_%i_ts_%i_vs_%i.csv" % (count, nx, ny, nz, order, ws, ts, vs)
            write_file = open(filename_f,'w')
            count += 1;
            write_file.write(str(key)+"\n")
            for i in range(len(workset_range)):
                write_file.write(str(workset_range[i])+", "+str(timings[key][i])+"\n")
        
    print "Finished Workset Analysis"

    if args.verbose:
        print timings
    

        # f = open(filename, mode='r')
        # lines = f.readlines()
        # for line in lines:
        #     print line,
        #     split_line = line.split(" ")
        #     print split_line[1]
        # f.close()
        

    
    #os.chdir('/Users')

    # Write timestamp for backup
    #os.chdir('/Users/rppawlo')
    #timestamp_file = open('BACKUP_DATE', 'w')
    #today = datetime.datetime.today()
    #date = today.strftime("YYYY.MM.DD: %Y.%m.%d at HH.MM.SS: %H.%M.%S")
    #timestamp_file.write(date)
    #timestamp_file.close()
    
    print '****************************************'
    print '* Finished Panzer Analysis!'
    print '****************************************'
    

#############################################################################
# If called from the command line, call main()
#############################################################################

if __name__ == "__main__":
    main()
