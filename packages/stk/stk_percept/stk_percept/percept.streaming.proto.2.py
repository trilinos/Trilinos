#!/usr/netpub/python-2.7/bin/python

# TODO: how to ensure python 2.7 is installed everywhere?
#       then use !/usr/bin/env python
# TODO: make a python .egg and have users install it?

"""Main function for the Percept command line tool.

This program provides functionality for a number of Percept tasks from
the command line.

Example usage:
# Uniformly refine a mesh containing different element types.\"
percept refine -n \"Summer Vacation 2009\" -t Vermont ~/photos/vacation2009/*

# blah blah
percept refine

# Convert a mesh of hexahedral elements to a mesh of tetrahedra
percept convert

# Convert a mesh of linear elements, e.g., HEX8,  TET4, QUAD4, TRI3,
# into a mesh of quadratic elements, e.g., HEX27, TET8, QUAD9, TRI6
percept enrich

Some terminology in use:
  task: What the client wants done (e.g. refine, enrich, convert).

"""

import os
import sys
import argparse
import textwrap
import subprocess

from multiprocessing import Pool


def setup_parser():
    """Set up the parser.

    Returns:
      argparse.ArgumentParser with all options configured.
    """

    epilog_txt = textwrap.dedent('''\
        Some commonly used %(prog)s commands are:
           refine    Divide blocks of elements uniformly
           convert   Convert blocks of elements to another specified type of element
           enrich    Add more nodes to blocks of elements to increase polynomial order

        See "%(prog)s COMMAND --help" for more information on a specific command.

        ''')

    parser = argparse.ArgumentParser(prog='percept', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Percept is a tool for verifying modeling and simulation codes and solutions.',
                                     epilog=epilog_txt)

    # version number of this script
    # TODO: get the version from percept
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    subparsers = parser.add_subparsers()

    # the block argument can be used with multiple tasks, predefine string here
    blocks_arg_list = ('-b', '--blocks')
    blocks_arg_dict = {'metavar':'BLOCK_LIST', 'help':'include or exclude block numbers in a quoted, comma delimited string.'}
        #(1) empty string or option not specified: convert all blocks in the input mesh file
        #(2) file:my_filename.my_ext (e.g. file:filelist.dat) which will read input block names from the given file
        #(3) [+]block_name_1,[+]block_name_2, etc ,block_name_n to include only these blocks, plus sign is optional
        #(4) a single input block name (e.g. block_3) to be converted
        #(5) -block_3,-block_5 to exclude blocks from those included (all blocks or include-only blocks), minus sign is mandatory
        #(6) block_1..block_10 include the range of blocks #1 to #10
        #(7) any combination of [+] and - options and range (..) option can be specified
        #Note: wherever you specify block_number this can be replaced with just the #, e.g. "1,2,4,5"

    blocks_file_arg_list = ('-B','--blocks-from')
    blocks_file_arg_dict = {'metavar':'BLOCK_FILE', 'help':'include or exclude block numbers listed in a file.'}

    # hide the --remove-original-elements argument, for now
    #original_arg_list = ('--keep-original-elements',)
    #original_arg_dict = {'action':'store_true', 'default':False, 'help':'Keep the original mesh embedded in the output file.'}

    # Tasks:
    # each task has a subparser and possibly task specific options

    # refine
    parser_refine = subparsers.add_parser('refine')
    parser_refine.add_argument('-n', '--number-refines', type=int, default=1, help='number of refinements (default: %(default)s)')
    parser_refine.add_argument('-t', '--refine-type', default='DEFAULT', metavar='REFINE_TYPE', help='type of elements to refine into another type (one of: Quad4_Quad4_4, Tri3_Tri3_4, Tet4_Tet4_8, Hex8_Hex8_8, Wedge6_Wedge6_8, Tri6_Tri6_4, Quad9_Quad9_4, Hex27_Hex27_8, Tet10_Tet10_8, Wedge18_Wedge18_8, ShellTri3_ShellTri3_4, ShellQuad4_ShellQuad4_4)')
    parser_refine.add_argument(*blocks_arg_list, **blocks_arg_dict)
    parser_refine.add_argument(*blocks_file_arg_list, **blocks_file_arg_dict)
    #parser_refine.add_argument(*original_arg_list, **original_arg_dict)

    # convert
    parser_convert = subparsers.add_parser('convert')
    parser_convert.add_argument('-t', '--convert-type', default="", metavar='CONVERT_TYPE', help='type of elements to convert to another type (one of: Hex8_Tet4_24, Hex8_Tet4_6, Quad4_Tri3_2, Quad4_Tri3_6, Quad4_Tri3_4)')
    parser_convert.add_argument(*blocks_arg_list, **blocks_arg_dict)
    parser_convert.add_argument(*blocks_file_arg_list, **blocks_file_arg_dict)
    #parser_convert.add_argument(*original_arg_list, **original_arg_dict)

    # enrich
    parser_enrich = subparsers.add_parser('enrich')
    parser_enrich.add_argument('-t', '--enrich-type', default='DEFAULT', metavar='ENRICH_TYPE', help='type of elements to convert to another type (one of: Quad4_Quad8_1, Quad4_Quad9_1, Tri3_Tri6_1, Tet4_Tet10_1, Hex8_Hex20_1, Hex8_Hex27_1, Wedge6_Wedge15_1, Wedge6_Wedge18_1)')
    parser_enrich.add_argument(*blocks_arg_list, **blocks_arg_dict)
    parser_enrich.add_argument(*blocks_file_arg_list, **blocks_file_arg_dict)
    #parser_enrich.add_argument(*original_arg_list, **original_arg_dict)

    # infile and outfile
    parser.add_argument('infile', nargs=1, help='input mesh file name')
    parser.add_argument('outfile', nargs='?', help='output mesh file name')

    # Other options
    parser.add_argument('-C','--directory', metavar='DIR', help='change working directory')
    
    # Parallel operation
    parallel = parser.add_argument_group('parallel control')
    # NOTE: --load-balance is taken care of by number of processes argument
    parallel.add_argument('-j', '--num-proc', type=int, default=1, metavar='NP', help='number of processors to use (default: %(default)s)')
    parallel.add_argument('--proc-rank-field', action='store_true', default=False, help='store the processor-rank in an element field on mesh output')

    # Streaming operation
    streaming = parser.add_argument_group('streaming control')
    streaming.add_argument('-M', '--streaming-size', type=int, default=0, metavar='NM', help='number of virtual processors in streaming mode = M in file.e.M.iM decomposition (default: %(default)s)')
    streaming.add_argument('-W', '--streaming-W', type=int, default=0, metavar='NW', help='number of actual parallel worker processes in streaming mode - must have (M mod W) == 0 (default: %(default)s)')

    # Informative output
    informative = parser.add_argument_group('informative output')
    # NOTE: --print-info N is now verbose
    informative.add_argument('-v', '--verbose', type=int, default=0, help='print more verbose information int >= 0  (default: %(default)s)')
    informative.add_argument('--echo-command-line', default=False, action='store_true', help='output original command line and continue as normal')
    informative.add_argument('-o','--output-log', default=('%s.log' % os.path.basename(sys.argv[0])), help='path where output log file will be written (default: %(prog)s.log)')
    informative.add_argument('--show-defaults', action='store_true', default=False, help='show %(prog)s defaults')

    # diagnostic and debugging options
    diagnostic = parser.add_argument_group('diagnostic and debugging')
    diagnostic.add_argument('--pause-for-debugging', default=False, action='store_true', help='wait for your input to allow attaching debugger')
    diagnostic.add_argument('--pout', default='-', help='path where per processor log files will be written (default: %(default)s)')
    diagnostic.add_argument('--dout', default='cout', help='path, or stream (one of cout, cerr), where diagnostic log file will be written (default: %(default)s)')
    diagnostic.add_argument('-D', '--debug', default=False, action='store_true', help='for internal debugging (default: %(default)s)')
    # NOTE: -dw is now --writer
    diagnostic.add_argument('--writer', metavar='WRITER', help='name of the diagnostic writer, for example "all"')
    diagnostic.add_argument('--timer', metavar="TIMER", help='name of the diagnostic timer, for example "mesh"')
    # NOTE: for now hide more esoteric diagnostics
    #diagnostic.add_argument('--test-memory-elements', metavar='NUM-ELEMS', type=int, default=0, help='specify number of elements (default: %(default)s)')
    #diagnostic.add_argument('--test-memory-nodes', metavar='NUM-NODES', type=int, default=0, help='specify number of nodes (default: %(default)s)')

    return parser

def convert_to_stk_adapt_option(args, arg_name):
    arg_value = getattr(args, arg_name)
    if (arg_name == "block"):
        if (arg_value is not None):
            return "--block-name=%s" % arg_value
    elif (arg_name == "blocks_from"):
        if (arg_value is not None):
            return "--block-name=file:%s" % arg_value
    elif (arg_name == "directory"):
        if (arg_value is not None):
            return "--directory=%s" % arg_value
    #elif (arg_name == "dout"):
    #    if (arg_value is not None):
    #        return "--dout=%s" % arg_value
    elif (arg_name == "echo_command_line"):
        if (arg_value):
            return "--echo-command-line"
    elif (arg_name == "infile"):
        assert type(arg_value) == list
        assert len(arg_value) == 1
        return "--input_mesh=%s" % arg_value[0]
    elif (arg_name == "number_refines"):
        return "--number_refines=%d" % arg_value
    elif (arg_name == "outfile"):
        if (arg_value is None):
            arg_value = getattr(args, "infile")[0].split(".")[0] + "_out.e"
        return "--output_mesh=%s" % arg_value
    #elif (arg_name == "output_log"):
    #    return "--output-log=%s" % arg_value
    elif (arg_name == "pause_for_debugging"):
        if (arg_value):
            return "--pause-for-debugging"
    #elif (arg_name == "pout"):
    #    if (arg_value is not None):
    #        return "--pout=%s" % arg_value
    elif (arg_name == "proc_rank_field"):
        if (arg_value):
            return "--proc_rank_field=%d" % 1
    elif (arg_name == "refine_type"):
        return "--refine=%s" % arg_value
    elif (arg_name == "convert_type"):
        return "--convert=%s" % arg_value
    elif (arg_name == "enrich_type"):
        return "--enrich=%s" % arg_value
    elif (arg_name == "show_defaults"):
        if (arg_value):
            return "" # TODO
    elif (arg_name == "timer"):
        if (arg_value is not None):
            return "--timer=%s" % arg_value
    elif (arg_name == "verbose"):
        if (arg_value > 0):
            return "--print_info=%d" % arg_value
    elif (arg_name == "writer"):
        if (arg_value is not None):
            return "--writer=%s" % arg_value

    return ""

def get_cmd_str(args, close_string=True):
    streaming_size = args.streaming_size
    streaming_W = args.streaming_W

    # Construct the sierra command from args
    cmd_str = "sierra"

    # TODO: use release version by omitting this later
    #       for now we default to head
    #cmd_str += " -V head"

    # Find what args/options need to go to the sierra command itself (not the app)
    skip_list = ["num_proc"]
    num_proc = args.num_proc

    # for now, only allow serial streaming (later we can do M --> N by P workers, P==num_proc)
    if streaming_size != 0:
        msg = ""
        if num_proc > 1:
            msg = " [Note: num_proc reduced from %d to 1 for serial streaming mode - parallel streaming not yet supported] " % num_proc
        num_proc = 1
        print "streaming_size= " , streaming_size, " num_proc= " , num_proc, msg

    cmd_str += " -j %d" % num_proc

    # Find what args/options need to go to stk_adapt
    cmd_str += ' stk_adapt_exe -O "'
    for arg in dir(args):
        if (not arg.startswith("_") and arg not in skip_list):
            cmd_str += " %s" % convert_to_stk_adapt_option(args, arg)

    if (num_proc > 1 and streaming_size == 0):
        cmd_str += " --load_balance=1"
    else:
        cmd_str += " --load_balance=0"

    if streaming_size > 0:
        cmd_str += " --streaming_size=%d" % streaming_size
    if streaming_W > 0:
        cmd_str += " --streaming_W=%d" % streaming_W

    if close_string:    
        cmd_str += '"'
    return cmd_str

def run_cmd(cmd_str):
    print "run_cmd: ", cmd_str

    proc = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, errput = proc.communicate()
    stat = proc.wait()

    if (stat != 0):
        raise SystemExit("%s failed:\n\nError: %s" % (os.path.basename(sys.argv[0]), errput))
    return output

def run_if(cmd_str, debug=False):
    if debug:
        print "debug print cmd_str= " , cmd_str
    else:
        output = run_cmd(cmd_str)
    return output

def create_scripts(cmd_str, debug=False):
    if debug:
        print "debug print cmd_str= " , cmd_str
    else:
        output = run_cmd(cmd_str)
    return output

def add_args(i_pass, streaming_W, iW):
    add = " --streaming_pass_start=%d --streaming_pass_end=%d --streaming_W=%d --streaming_iW=%d" % (i_pass, i_pass, streaming_W, iW)
    add += '"'
    return add

def run_if_par_0(cmd_str, i_pass, streaming_W, iW, debug, timeout=60):
    cmd_str1 = cmd_str + add_args(i_pass, streaming_W, iW)
    run_if(cmd_str1, debug)

def run_if_par(cmd_str, i_pass, streaming_W, debug, timeout=600):

    pool = Pool(processes=streaming_W)

    # the first time we create the sierra.W.iW.sh script, then we execute it below
    add_opt = " --script-only; mv sierra.sh sierra.W.%d.%d.sh"
    argsl = []
    for iW in range(streaming_W):
        cmd_str1 = cmd_str + add_args(i_pass, streaming_W, iW)
        cmd_str1 = cmd_str1 + (add_opt % (streaming_W,iW))
        argsl.append(cmd_str1)

    # this must be done in serial
    result = map(create_scripts, argsl)

    if True:
        for iW in range(streaming_W):
            print "create_scripts result[", iW, "]= ", result[iW]

    # now create the final sierra command
    add_opt = " --exec-script=sierra.W.%d.%d.sh"
    argsl = []
    for iW in range(streaming_W):
        cmd_str1 = cmd_str + add_args(i_pass, streaming_W, iW)
        cmd_str1 = cmd_str1 + (add_opt % (streaming_W,iW))
        argsl.append(cmd_str1)

    if False:
        result = map(run_if, argsl)
    else:
        # parallel run
        result = pool.map(run_if, argsl)

    if True:
        for iW in range(streaming_W):
            print "result[", iW, "]= ", result[iW]

def run(args, debug):
    streaming_W = args.streaming_W
 
    if streaming_W == 0:
        cmd_str = get_cmd_str(args)
        run_if(cmd_str, debug)
    else:
        cmd_str = get_cmd_str(args, False)

        for i_pass in range(-1, 4):

            # initial sync point
            iW = -1
            cmd_str1 = cmd_str + add_args(i_pass, streaming_W, iW)
            run_if(cmd_str1, debug)

            # parallel operation
            if False:
                for iW in range(0, streaming_W):
                    cmd_str1 = cmd_str + add_args(i_pass, streaming_W, iW)
                    run_if(cmd_str1, debug)
            else:
                print "startparallel ..."
                run_if_par(cmd_str, i_pass, streaming_W, debug)
                print "end parallel"

            # final sync point
            iW = streaming_W
            cmd_str1 = cmd_str + add_args(i_pass, streaming_W, iW)
            run_if(cmd_str1, debug)


def main():
    """Entry point for Percept script."""
    parser = setup_parser()
    #(options, args) = parser.parse_args()
    args = parser.parse_args()

    debug = args.debug
    #debug = True

    # TODO: setup_logger(options)
    if not args:
        # TODO: run_interactive(parser)
        print 'Interactive mode not yet working.'
    else:

        run(args, debug)

def exit_from_int(*args):
    """Handler for SIGINT signal."""
    print ''
    exit(0)


if __name__ == '__main__':
    import signal
    signal.signal(signal.SIGINT, exit_from_int)
    main()
