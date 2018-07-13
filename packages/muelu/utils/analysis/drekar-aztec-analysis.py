#!/usr/bin/env python
import matplotlib.cm     as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
import optparse
import re
import yaml


def nonlinear_history_iterations(yaml_data, mode, ax = None):
    """Show number of linear iterations per nonlinear step across the whole simulation"""
    ticks = []

    itit   = []
    maxits = 0

    offset = 0
    mx     = 0
    i      = 0
    for step in sorted(yaml_data['Steps']):
        c_step = yaml_data['Steps'][step]

        # check whether nl_its is contained in list
        # if YES, it is a dynamic simulation, i.e. Steps contains time step (transient steps)
        # if NO, it is a static simulation with nonlinear iterations only
	bDynamicSim = True
	try:
	  c_step['nl_its']
	except KeyError:
	  bDynamicSim = False

	if bDynamicSim == False:
	  raise RuntimeError("Static simulations not supported, yet")
	
        its = []
        for nlstep in sorted(c_step):
            if nlstep == 'nl_its':
                continue
            ticks.append(str(len(its)))
            its.append(c_step[nlstep]['its'])
        if maxits < len(its):
	    maxits = len(its)

	itit.append(its)

        i += 1

    print itit
    
    ikk = []
    for k in range(maxits):
	ik = []
	for l in range(len(itit)):
	    its = itit[l]
	    if len(its) > k:
	        ik.append(its[k])
	    else:
	        ik.append(0)
	ikk.append(ik)

    ind = np.arange(len(itit))
    
    bottomvals = []
    for l in range(len(itit)):
        bottomvals.append(0)
    for k in range(maxits):
	p1 = plt.bar(ind, ikk[k], 1.0, color='r',bottom = bottomvals)
	for l in range(len(itit)):
	    bottomvals[l] = bottomvals[l] + ikk[k][l]

    ax.set_ylabel('Number of accumulated linear iterations')
    ax.set_xlabel('Nonlinear iteration index')
    
    plt.show()

if __name__ == '__main__':
    p = optparse.OptionParser()

    # action arguments
    p.add_option('-i', '--input-file',  dest='input_file')
    p.add_option('-o', '--output-file', dest='output_file')
    p.add_option('-d', '--display',     dest='display',     default='print')
    p.add_option('-a', '--analysis',    dest='analysis',    default='mode1')

    # parse
    options, arguments = p.parse_args()

    # validate options
    if options.input_file == None:
        raise RuntimeError("Please specify the input file")

    with open(options.input_file) as data_file:
        yaml_data = yaml.safe_load(data_file)

    f, ax = plt.subplots(1)

    analysis = options.analysis
    display  = options.display
    
    nonlinear_history_iterations(yaml_data, mode=display, ax=ax)

    if options.output_file != None:
        plt.savefig(options.output_file, bbox_inches='tight')
    elif display == 'display':
        plt.show()
