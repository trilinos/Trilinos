#!/usr/bin/env python
import optparse
import re
import yaml
import logging
from   transitions      import Machine, logger

# YAML scheme structure
#
# scheme:                   <string>            # Particular YAML scheme ['muelu' | 'albany' | 'drekar']
# Steps:
#     c_step_1:
#         nl_its:           <int>               # number of nonlinear iterations
#         nl_step_1:
#             its:          <int>               # number of linear iterations
#             res_hist:     <double list>       # residual history
#             setup_time:   <double>            # setup time in s
#             solve_time:   <double>            # solve time in s
#         nl_step_2:
#             its:          <int>
#             res_hist:     <double list>
#             setup_time:   <double>            # setup time in s
#             solve_time:   <double>            # solve time in s
# Timer names:              <string list>
# Total times:
#    <timer_name>:          {MaxOverProcs: <double>, MeanOverCallCounts: <double>, MeanOverProcs: <double>, MinOverProcs: <double>}
# Call counts:
#    <timer_name>:          {MaxOverProcs: <double>, MeanOverCallCounts: <double>, MeanOverProcs: <double>, MinOverProcs: <double>}
# Number of processes:      1
# Statistics collected:     [MinOverProcs, MeanOverProcs, MaxOverProcs, MeanOverCallCounts]
# Time unit:                's'

class Drekar(object):
    pass

def log2yaml(filename, presvec):
    # construct the YAML string
    # 2nd version, using FSM (Finite State Machine)
    states = ['none', 'transient', 'nonlinear', 'linear', 'l_setup', 'l_solutiontime', 'l_iteration']
    transitions = [
        { 'trigger': 'n2t',              'source': 'none',                        'dest': 'transient'  },
        { 'trigger': 't2n',              'source': 'transient',                   'dest': 'none'       },
        { 'trigger': 't2nl',             'source': 'transient',                   'dest': 'nonlinear'  },
        { 'trigger': 'nl2t',             'source': 'nonlinear',                   'dest': 'transient'  },
        { 'trigger': 'nl2l',             'source': 'nonlinear',                   'dest': 'linear'     },
        { 'trigger': 'l2l_setup',        'source': 'linear',                      'dest': 'l_setup'    },
        { 'trigger': 'l_setup2l',        'source': 'l_setup',                     'dest': 'linear'     },
        { 'trigger': 'l2l_solutiontime', 'source': 'linear',                      'dest': 'l_solutiontime'    },
        { 'trigger': 'l_solutiontime2l_iteration', 'source': 'l_solutiontime',    'dest': 'l_iteration'},        
        { 'trigger': 'l_iteration2t',    'source': 'l_iteration',                 'dest': 'transient'  },
        { 'trigger': 'n2nl',             'source': 'none',                        'dest': 'nonlinear'  },
        { 'trigger': 'n2l_setup',        'source': 'none',                        'dest': 'l_setup'  }
    ]

    logger.setLevel(logging.INFO)

    machine = Drekar()
    M = Machine(model=machine, states=states, transitions=transitions, initial='none')

    mode = 'drekar'
    start_stepper                   = 'Entering Rythmos::.*::advanceStepperToTime'
    end_stepper                     = 'Leaving Rythmos::.*::advanceStepperToTime'
    start_c_step                    = 'Entering Rythmos::.*::takeStep'
    end_c_step                      = 'Leaving Rythmos::.*::takeStep'
    start_nl_step                   = '(?<=Nonlinear Solver Step )\d*'
    end_nl_step_good                = '\(Converged!\)'
    end_nl_step_bad                 = '\(Failed!\)'
    start_l_step                    = 'CALCULATING FORCING TERM'
    start_residual_aztecoo          = 'Problem: Thyra::DefaultBlockedLinearOp'
    start_setuptime                 = '(?<=TimeMonitor results over )\d*'
    end_setuptime                   = '(?<=MueLu setup time )\d*'
    end_simulation                  = '(?<=Drekar run completed.)'
    start_solutiontime              = '(?<=Solution time: ).*(?= \(sec.\))'
    end_l_step                      = '(?<=total iterations: )\d*'
    close_file                      = 'CLOSEFILE'
    
    #end_residual_belos_good         = '(?<=returned a solve status of "SOLVE_STATUS_CONVERGED" in )\d*'
    #end_residual_belos_bad          = '(?<=returned a solve status of "SOLVE_STATUS_UNCONVERGED" in )\d*'
    
    #mid_timers                      = '.* \(.*\)\s*$'
    #end_timers                      = '==========='

    timer_names  = []
    timer_times  = {}
    timer_calls  = {}
    timer_serial = None

    nl_step = 0
    t_step  = 0
    lineno  = 0
    
    l_setuptimes = []

    cur_setup_time = 0.0

    yaml_string = '{"scheme":"' + mode + '","Steps":{'
    with open(filename) as f:
        try:
            for line in f:
                lineno += 1
                if   re.search(start_stepper, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    t_step  = 0
                    
                elif re.search(end_stepper, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)

                elif re.search(start_c_step, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.n2t()

                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"c_step_' + str(t_step) + '":{'

                elif re.search(end_c_step, line) != None:
                    if machine.state != 'transient':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.t2n()

                    yaml_string += ', "nl_its":' + str(nl_step-1) + '}'

                    nl_step  = 0
                    t_step  += 1
                
                # start NONLINEAR state
                elif re.search(start_nl_step, line) != None:
		    print 'start nonlinear state...'
                    if machine.state != 'transient' and machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    
                    if machine.state == 'none':
		      machine.n2nl()
                    elif machine.state == 'transient':
		      machine.t2nl()
		    

                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"nl_step_' + str(nl_step) + '":{'

                    nl_step += 1
                    
                elif re.search(end_nl_step_good, line) != None or re.search(end_nl_step_bad, line) != None:
		    print 'end nonlinear state...'
                    if machine.state != 'nonlinear':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.nl2t()

		    print yaml_string

                    # Get rid of ",nl_step_?:{}
                    i = 1
                    while (yaml_string[-i] != ','):
                         i += 1
                    yaml_string = yaml_string[:-i]
              
                    print yaml_string
                # start LINEAR state
                elif re.search(start_l_step, line) != None:
		    print 'start linear state...'
                    if machine.state != 'nonlinear':
                        raise RuntimeError('Wrong state: ' + machine.state)
		    machine.nl2l()

		    l_setuptimes[:] = [] # empty list of setup times
		    	    
                # Collect setup time data for MueLu             
                elif re.search(start_setuptime, line) != None:
		    print 'collect setup time...'
		    if machine.state == 'none' or machine.state == 'transient':
			continue
                    if machine.state != 'linear' and machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
		      
		    machine.l2l_setup()

		elif re.search(end_setuptime, line) != None:
		    print 'collect setup time done...'
                    if machine.state != 'l_setup':
                        raise RuntimeError('Wrong state: ' + machine.state)
		                   
                    m = re.search(end_setuptime, line)
                    if m != None:
                        tt = m.group()
                    else:
                        m = re.search(end_timers_aztecoo, line)
                        tt = m.group()
                    cur_setup_time = line.split(" ")[6]
                    
                    l_setuptimes.append(float(cur_setup_time))
                    
                    machine.l_setup2l()

		# extract linear solution time
                elif re.search(start_solutiontime, line) != None:
		    print 'extract solution time...'
                    if machine.state != 'linear':
                        raise RuntimeError('Wrong state: ' + machine.state)

		    machine.l2l_solutiontime()
		    
                    m = re.search(start_solutiontime, line)
                    tt = m.group()
                    yaml_string += '  "solve_time":' + str(tt)   
                   
                    machine.l_solutiontime2l_iteration()

		# extract linear iterations
		# end LINEAR state
                elif re.search(end_l_step, line) != None:
		    print 'extract iterations...'
                    if machine.state != 'l_iteration':
                        raise RuntimeError('Wrong state: ' + machine.state)

                    m = re.search(end_l_step, line)
                    tt = m.group()
                    yaml_string += ', "its":' + str(tt)   
                   
		    yaml_string += ', "setup_time":' + str(sum(l_setuptimes)) + '}'
		    l_setuptimes[:] = []
                   
                    machine.l_iteration2t()                             
                    print 'end linear state...'
                elif re.search(end_simulation, line) != None:
		  
		    print 'end simulation'
		    
                elif re.search(close_file, line) != None:
		    if machine.state == 'linear':
			yaml_string += ' "its": -1' + ' }'
                        machine.l2t()
                       
		    if machine.state == 'nonlinear':
                        machine.nl2t()

			# Get rid of ",nl_step_?:{}
			i = 1
			while (yaml_string[-i] != ','):
			    i += 1
			yaml_string = yaml_string[:-i]
                    
                    if machine.state == 'transient':
		      	i = 1
			while (yaml_string[-i] != ','):
			    i += 1
			yaml_string = yaml_string[:-i]
			yaml_string += ', "nl_its":' + str(nl_step) + '}'
			nl_step  = 0
			machine.t2n()

		  
        except RuntimeError as e:
            raise RuntimeError("Caught an error while parsing on line " + str(lineno) + ": " + e.args[0])

    if timer_serial == None:
        # We did not encounter any timers
        yaml_string += '}}'

    print yaml_string
    print yaml.load(yaml_string)

    try:
        yaml_data = yaml.load(yaml_string)
    except yaml.parser.ParserError:
        raise RuntimeError('Could not parse YAML out. Did you select the right mode?')

    
    return yaml_data
    




if __name__ == '__main__':
    p = optparse.OptionParser()

    # action arguments
    p.add_option('-i', '--input-file',    dest='input_file')
    p.add_option('-o', '--output-file',   dest='output_file')
    p.add_option('-m', '--mode',          dest='mode',       default='muelu')
    p.add_option('-r', '--residual',      dest='resvec', default='0')

    # parse
    options, arguments = p.parse_args()

    # validate options
    if options.input_file == None:
        raise RuntimeError("Please specify an input file")
    filename = options.input_file

    mode = options.mode
    assert(mode == 'drekar')

    resvec = options.resvec

    yaml_data = log2yaml(filename, presvec=resvec)

    # dump the data
    output_file = options.output_file
    if output_file == None:
      output_file = filename + '.yaml'

    f = open(output_file, 'w')
    yaml.dump(yaml_data, stream=f)
