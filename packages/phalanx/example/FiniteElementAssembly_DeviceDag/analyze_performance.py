#! /usr/bin/env python

"""
Script for analyzing Phalanx performance on next-generation
architectures.  Runs hierarchic parallelism and Host/Device DAG
analysis. Generates plots from data.
"""

__version__ = "1.0"
__author__  = "Roger Pawlowski"
__date__    = "Dec 1 2017"

# Import python modules for command-line options, the operating system, regular
# expressions, and system functions
import commands
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
    print '* Starting Phalanx Analysis'
    print '****************************************'

    parser = argparse.ArgumentParser(description='Phalanx hierarchic parallelism and Host/Device DAG analysis script')
    parser.add_argument('-r', '--run', action='store_true', help='Run the executable to generate data and output to files.')
    parser.add_argument('-a', '--analyze', action='store_true', help='Analyze the data from files generated with --run.')
    parser.add_argument('-p', '--prefix', help='Add a prefix string to all output filenames.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print more data to screen.')
    parser.add_argument('-ne', '--num-eqns', type=int, default=16, help='Number of equations to solve for.')
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
    ne = args.num_eqns
    ts = args.team_size
    vs = args.vector_size
    print "number of equations = %d, team size = %d, vector size = %d\n" % (ne, vs, vs)

    executable = "./Phalanx_example_finite_element_assembly.exe"
    if args.mixed_exe:
        executable = "./Phalanx_example_finite_element_assembly_mixed_field_types.exe"
    elif args.device_dag_exe:
        executable = "./Phalanx_example_finite_element_assembly_device_dag.exe"

    print "Starting Workset Analysis"

    ws_step_size = 100
    workset_range = range(100,2000+ws_step_size,ws_step_size)
    print "workset range = "+str(workset_range)

    timings = {}
    if args.analyze:
        import numpy as np
        timings["Jacobian Evaluation Time <<Host DAG>>"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["Residual Evaluation Time <<Host DAG>>"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["Jacobian Evaluation Time <<Device DAG>>"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["Residual Evaluation Time <<Device DAG>>"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] Gather Solution"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] Gather Solution"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] ZeroContributedField"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] ZeroContributedField"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] ProjectValueToQP"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] ProjectValueToQP"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] ProjectGradientToQP"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] ProjectGradientToQP"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] IntegrateDiffusionTerm"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] IntegrateDiffusionTerm"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] IntegrateSourceTerm"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] IntegrateSourceTerm"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Jacobian] Scatter Jacobian"] = np.zeros(len(workset_range),dtype=np.float64)
        timings["[PHX::MyTraits::Residual] Scatter Residual"] = np.zeros(len(workset_range),dtype=np.float64)

    #print dir(np)
        
    for i in range(len(workset_range)):

        ws = workset_range[i]
        
        filename = "fea_nx_%i_ny_%i_nz_%i_ne_%i_ws_%i_ts_%i_vs_%i.dat" % (nx, ny, nz ,ne, ws, ts, vs)
        if args.prefix:
            filename = args.prefix+filename
        command = executable+" -nx %i -ny %i -nz %i -ne %i -ws %i -ts %i -vs %i" % (nx, ny, nz ,ne, ws, ts, vs)  +" >& "+filename

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
                        timings[key][i] += float(split_line[-2])
                        if args.verbose:
                            print "  found key: "+key+" = "+str(split_line[-2])
                        break
            f.close()

    do_jac = True
            
    if args.analyze:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.semilogy()
        if do_jac:
            # maroon = #990033, light blue = #00ffff
            plt.plot(workset_range,timings["Jacobian Evaluation Time <<Host DAG>>"],label="Jac Total Time (Host DAG)",marker="o",color="#990033",markersize=8)
            plt.plot(workset_range,timings["Jacobian Evaluation Time <<Device DAG>>"],label="Jac Total Time (Device DAG)",marker="s",color="r",markersize=8)
        plt.plot(workset_range,timings["Residual Evaluation Time <<Host DAG>>"],label="Res Total Time (Host DAG)",marker="o",color="b",markersize=8)
        plt.plot(workset_range,timings["Residual Evaluation Time <<Device DAG>>"],label="Res Total Time (Device DAG)",marker="s",color="#00ffff",markersize=8)
        plt.xlabel("Workset Size",fontsize=16)
        plt.ylabel("Time (s)",fontsize=16)
        plt.tick_params(labelsize=16)
        title = "nel=%i,neq=%i,nderiv=%i,ts=%i,vs=%i" % (nx*ny*nz,ne,8*ne,ts,vs)
        plt.title(title)
        #plt.legend(bbox_to_anchor=(1,1))
        plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.0),ncol=2,fancybox=True,shadow=True, prop={'size': 12})
        plt.grid()
        dag_timings_filename = "dag_timings_nx_%i_ny_%i_nz_%i_ne_%i_ts_%i_vs_%i.png" % (nx, ny, nz ,ne, ts, vs)
        fig.savefig(dag_timings_filename)
        #plt.show()

        if do_jac:
            fig = plt.figure(2)
            #plt.clf()
            plt.semilogy()
            plt.plot(workset_range,timings["Jacobian Evaluation Time <<Host DAG>>"],label="Jac Total Time (Host DAG)",marker='o')
            plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] Gather Solution"],label="Jac Gather",marker='s')
            plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] ZeroContributedField"],label="Jac Zero",marker='^')
            #plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] ProjectValueToQP"],label="Jac ProjToQP",marker='+')
            plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] ProjectGradientToQP"],label="Jac ProjGradToQP",marker='*')
            plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] IntegrateDiffusionTerm"],label="Jac Int Diff Term",marker='x')
            plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] IntegrateSourceTerm"],label="Jac Int Source Term",marker='D')
            plt.plot(workset_range,timings["[PHX::MyTraits::Jacobian] Scatter Jacobian"],label="Jac Scatter",marker='.',color="#ff6600")
            plt.xlabel("Workset Size",fontsize=16)
            plt.ylabel("Time (s)",fontsize=16)
            plt.tick_params(labelsize=16)
            title = "nel=%i,neq=%i,nderiv=%i,ts=%i,vs=%i" % (nx*ny*nz,ne,8*ne,ts,vs)
            plt.title(title)
            #plt.legend(bbox_to_anchor=(1,1))
            #plt.legend(loc=1, prop={'size': 10})
            plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.0),ncol=3,fancybox=True,shadow=True, prop={'size': 10})
            plt.grid()
            jac_evaluator_timings_filename = "jac_evalautor_timings_nx_%i_ny_%i_nz_%i_ne_%i_ts_%i_vs_%i.png" % (nx, ny, nz ,ne, ts, vs)
            fig.savefig(jac_evaluator_timings_filename)
            
        fig = plt.figure(3)
        #plt.clf()
        plt.semilogy()
        plt.plot(workset_range,timings["Residual Evaluation Time <<Host DAG>>"],label="Res Total Time (Host DAG)",marker='o')
        plt.plot(workset_range,timings["[PHX::MyTraits::Residual] Gather Solution"],label="Res Gather",marker='s')
        plt.plot(workset_range,timings["[PHX::MyTraits::Residual] ZeroContributedField"],label="Res Zero",marker='^')
        #plt.plot(workset_range,timings["[PHX::MyTraits::Residual] ProjectValueToQP"],label="Res ProjToQP",marker='+')
        plt.plot(workset_range,timings["[PHX::MyTraits::Residual] ProjectGradientToQP"],label="Res ProjGradToQP",marker='*')
        plt.plot(workset_range,timings["[PHX::MyTraits::Residual] IntegrateDiffusionTerm"],label="Res Int Diff Term",marker='x')
        plt.plot(workset_range,timings["[PHX::MyTraits::Residual] IntegrateSourceTerm"],label="Res Int Source Term",marker='D')
        plt.plot(workset_range,timings["[PHX::MyTraits::Residual] Scatter Residual"],label="Res Scatter",marker='.',color="#ff6600")
        plt.xlabel("Workset Size",fontsize=16)
        plt.ylabel("Time (s)",fontsize=16)
        plt.tick_params(labelsize=16)
        title = "nel=%i,neq=%i,nderiv=%i,ts=%i,vs=%i" % (nx*ny*nz,ne,8*ne,ts,vs)
        plt.title(title)
        #plt.legend(bbox_to_anchor=(1,1))
        plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.0),ncol=3,fancybox=True,shadow=True, prop={'size': 10})
        #plt.axis([0,2000,1.0e-4,0.1])
        plt.grid()
        res_evaluator_timings_filename = "res_evalautor_timings_nx_%i_ny_%i_nz_%i_ne_%i_ts_%i_vs_%i.png" % (nx, ny, nz ,ne, ts, vs)
        fig.savefig(res_evaluator_timings_filename)

        #print dir(plt)


        # Plot to assess savings
        filename_f = "raw_data_output_f_nx_%i_ny_%i_nz_%i_ne_%i_ws_%i_ts_%i_vs_%i.csv" % (nx, ny, nz ,ne, ws, ts, vs)
        write_file = open(filename_f,'w')
        for i in range(len(workset_range)):
            write_file.write(str(workset_range[i])+", "+str(timings["Residual Evaluation Time <<Host DAG>>"][i])+"\n")

        filename_J = "raw_data_output_J_nx_%i_ny_%i_nz_%i_ne_%i_ws_%i_ts_%i_vs_%i.csv" % (nx, ny, nz ,ne, ws, ts, vs)
        write_file = open(filename_J,'w')
        for i in range(len(workset_range)):
            write_file.write(str(workset_range[i])+", "+str(timings["Jacobian Evaluation Time <<Host DAG>>"][i])+"\n")

        
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
    print '* Finished Phalanx Analysis!'
    print '****************************************'
    

#############################################################################
# If called from the command line, call main()
#############################################################################

if __name__ == "__main__":
    main()
