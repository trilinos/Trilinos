#!/usr/bin/env python3
import yaml
from log import Log
import logging
import re
from transitions import Machine, logger

class DrekarBelosLog(Log):
    """Log parser of transient Drekar simulations with Belos."""

    states = ['none', 'transient', 'nonlinear', 'linear', 'timers']

    transitions = [
        { 'trigger': 'n2m',     'source': 'none',       'dest': 'timers'     },
        { 'trigger': 'n2t',     'source': 'none',       'dest': 'transient'  },
        { 'trigger': 't2n',     'source': 'transient',  'dest': 'none'       },
        { 'trigger': 't2nl',    'source': 'transient',  'dest': 'nonlinear'  },
        { 'trigger': 'nl2t',    'source': 'nonlinear',  'dest': 'transient'  },
        { 'trigger': 'nl2l',    'source': 'nonlinear',  'dest': 'linear'     },
        { 'trigger': 'l2nl',    'source': 'linear',     'dest': 'nonlinear'  },
        { 'trigger': 'l2t',     'source': 'linear',     'dest': 'transient'  },
        { 'trigger': 'm2n',     'source': 'timers',     'dest': 'none'       }
    ]

    def run(self, filename):
        logger.setLevel(logging.INFO)

        machine = self
        M = Machine(model=machine, states=self.states, transitions=self.transitions, initial='none')

        start_stepper                   = 'Entering Rythmos::.*::advanceStepperToTime'
        end_stepper                     = 'Leaving Rythmos::.*::advanceStepperToTime'
        start_c_step                    = 'Entering Rythmos::.*::takeStep'
        end_c_step                      = 'Leaving Rythmos::.*::takeStep'
        start_nl_step                   = '(?<=Nonlinear Solver Step )\d*'
        end_nl_step_good                = '\(Converged!\)'
        end_nl_step_bad                 = '\(Failed!\)'
        start_residual_belos            = '\*\*\*\*\* Belos Iterative Solver: '
        mid_residual                    = 'Iter.*, \[.*\] :\s*.*'
        end_residual_belos_good         = '(?<=returned a solve status of "SOLVE_STATUS_CONVERGED" in )\d*'
        end_residual_belos_bad          = '(?<=returned a solve status of "SOLVE_STATUS_UNCONVERGED" in )\d*'
        start_timers                    = '(?<=TimeMonitor results over )\d*'
        mid_timers                      = '.* \(.*\)\s*$'
        end_timers                      = '==========='

        timer_names  = []
        timer_times  = {}
        timer_calls  = {}
        timer_serial = None

        nl_step = 0
        t_step  = 0
        lineno  = 0

        yaml_string = '{"scheme":"drekar","Steps":{'
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

                        yaml_string += ', "nl_its":' + str(nl_step) + '}'

                        nl_step  = 0
                        t_step  += 1

                    elif re.search(start_nl_step, line) != None:
                        if machine.state != 'transient':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.t2nl()

                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"nl_step_' + str(nl_step) + '":{'

                        nl_step += 1

                    elif re.search(end_nl_step_good, line) != None or re.search(end_nl_step_bad, line) != None:
                        if machine.state != 'nonlinear':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.nl2t()

                        # Get rid of ",nl_step_?:{}
                        i = 1
                        while (yaml_string[-i] != ','):
                             i += 1
                        yaml_string = yaml_string[:-i]

                    elif re.search(start_residual_belos, line) != None:
                        if machine.state != 'nonlinear':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.nl2l()

                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"res_hist":['

                    elif re.search(end_residual_belos_good, line) != None or re.search(end_residual_belos_bad, line) != None:
                        if machine.state != 'linear':
                            raise RuntimeError('Wrong state: ' + machine.state)

                        machine.l2t()

                        m = re.search(end_residual_belos_good, line)
                        if m != None:
                            its = m.group()
                        else:
                            m = re.search(end_residual_belos_bad, line)
                            its = m.group()
                        yaml_string += '], "its":' + its

                        m = re.search('(?<=with total CPU time of ).*(?=\ sec)', line)
                        belos_time = m.group()
                        yaml_string += ', "solve_time":' + belos_time + '}'

                    elif re.search(mid_residual, line) != None:
                        if machine.state != 'linear':
                            raise RuntimeError('Wrong state: ' + machine.state)

                        m = re.search('[^\s]*$', line)
                        res = m.group()
                        if yaml_string[-1] != '[':
                            yaml_string += ','
                        yaml_string += res

                    elif re.search(start_timers, line) != None:
                        if machine.state != 'none':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.n2m()

                        m = re.search(start_timers, line)
                        nprocs = m.group()

                        if nprocs == "1":
                            timer_serial = True
                        else:
                            timer_serial = False

                        # Finalize stepping
                        yaml_string += '}'

                        yaml_string += ',"Number of processes":' + nprocs
                        yaml_string += ',"Time unit":s'
                        yaml_string += ',"Statistics collected":["MinOverProcs","MeanOverProcs","MaxOverProcs","MeanOverCallCounts"]'

                    elif re.search(end_timers, line) != None:
                        if machine.state != 'timers':
                            # there could be other ======== lines
                            continue
                        machine.m2n()

                        # Timer names
                        yaml_string += ',"Timer names":['
                        for name in timer_names:
                            if yaml_string[-1] != '[':
                                yaml_string += ','
                            yaml_string += '"' + name + '"'
                        yaml_string += ']'

                        # Total times
                        yaml_string += ',"Total times":{'
                        for name in timer_names:
                            if yaml_string[-1] != '{':
                                yaml_string += ','
                            yaml_string += '"' + name + '":{'
                            yaml_string += '"MinOverProcs":'        + timer_times[name]['MinOverProcs']
                            yaml_string += ',"MeanOverProcs":'      + timer_times[name]['MeanOverProcs']
                            yaml_string += ',"MaxOverProcs":'       + timer_times[name]['MaxOverProcs']
                            yaml_string += ',"MeanOverCallCounts":' + timer_times[name]['MeanOverCallCounts']
                            yaml_string += '}'
                        yaml_string += '}'

                        # Call counts
                        yaml_string += ',"Call counts":{'
                        for name in timer_calls:
                            if yaml_string[-1] != '{':
                                yaml_string += ','
                            yaml_string += '"' + name + '":{'
                            yaml_string += '"MinOverProcs":'        + timer_times[name]['MinOverProcs']
                            yaml_string += ',"MeanOverProcs":'      + timer_times[name]['MeanOverProcs']
                            yaml_string += ',"MaxOverProcs":'       + timer_times[name]['MaxOverProcs']
                            yaml_string += ',"MeanOverCallCounts":' + timer_times[name]['MeanOverCallCounts']
                            yaml_string += '}'
                        yaml_string += '}'

                        yaml_string += '}'

                    elif re.search(mid_timers, line) != None:
                        if machine.state != 'timers':
                            # there could be other matching lines
                            continue

                        if re.search('Timer Name', line) != None:
                            # Skip header (in serial it matches the pattern)
                            continue

                        splits = line.split()
                        if timer_serial == True:
                            name = ' '.join(splits[0:-2])
                            name = name.replace('"', '\\"')
                            timer_names.append(name)
                            timer_times[name] = {}
                            timer_times[name]['MinOverProcs']       = splits[-2];
                            timer_times[name]['MeanOverProcs']      = splits[-2];
                            timer_times[name]['MaxOverProcs']       = splits[-2];
                            timer_times[name]['MeanOverCallCounts'] = splits[-2];
                            timer_calls[name] = {}
                            timer_calls[name]['MinOverProcs']       = splits[-1][1:-1];
                            timer_calls[name]['MeanOverProcs']      = splits[-1][1:-1];
                            timer_calls[name]['MaxOverProcs']       = splits[-1][1:-1];
                            timer_calls[name]['MeanOverCallCounts'] = splits[-1][1:-1];
                        else:
                            name = ' '.join(splits[0:-8])
                            name = name.replace('"', '\\"')
                            timer_names.append(name)
                            timer_times[name] = {}
                            timer_times[name]['MinOverProcs']       = splits[-8];
                            timer_times[name]['MeanOverProcs']      = splits[-6];
                            timer_times[name]['MaxOverProcs']       = splits[-4];
                            timer_times[name]['MeanOverCallCounts'] = splits[-2];
                            timer_calls[name] = {}
                            timer_calls[name]['MinOverProcs']       = splits[-7][1:-1];
                            timer_calls[name]['MeanOverProcs']      = splits[-5][1:-1];
                            timer_calls[name]['MaxOverProcs']       = splits[-3][1:-1];
                            timer_calls[name]['MeanOverCallCounts'] = splits[-1][1:-1];
            except RuntimeError as e:
                raise RuntimeError("Caught an error while parsing on line " + str(lineno) + ": " + e.args[0])

        if timer_serial == None:
            # We did not encounter any timers
            yaml_string += '}}'

        try:
            yaml_data = yaml.load(yaml_string)
        except yaml.parser.ParserError:
            raise RuntimeError('Could not parse YAML out. Did you select the right mode?')

        return yaml_data
