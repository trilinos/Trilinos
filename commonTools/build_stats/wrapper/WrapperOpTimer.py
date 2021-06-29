import subprocess
import csv
import os
from WrapperCommandLineParser import WrapperCommandLineParser

def get_full_header(fields_list,full_header_map):
  return ','.join([ full_header_map[f] for f in fields_list ])


# the values are
usr_bin_time_csv_map = {
  "E":
    "elapsed_real_time_fmt",
  "e":
    "elapsed_real_time_sec",
  "S":
    "cpu_sec_kernel_mode",
  "U":
    "cpu_sec_user_mode",
  "P":
    "perc_cpu_used",
  "M":
    "max_resident_size_Kb",
  "t":
    "avg_resident_size_Kb",
  "K":
    "avg_total_memory_used_Kb",
  "D":
    "avg_size_unshared_data_area_Kb",
  "p":
    "avg_size_unshared_stack_area_Kb",
  "X":
    "avg_size_unshared_text_area_Kb",
  "Z":
    'page_size_bytes',
  "F":
    "num_major_page_faults",
  "R":
    "num_minor_page_faults",
  "W":
    "num_swapped",
  "c":
    "num_involuntary_context_switch",
  "w":
    "num_waits",
  "I":
    "num_filesystem_inputs",
  "O":
    "num_filesystem_outputs",
  "r":
    "num_socket_msg_recv",
  "s":
    "num_socket_msg_sent",
  "k":
    "num_signals",
  "x":
    "exit_status",
}

usr_bin_time_desc_map = {
  "E":
    "Elapsed real time ([h:]m:s)",
  "e":
    "Elapsed real time (s)",
  "S":
    "Total number of CPU-seconds that the process spent in kernel mode",
  "U":
    "Total number of CPU-seconds that the process spent in user mode",
  "P":
    "Percentage of the CPU that this job got",
  "M":
    "Maximum resident set size of the process during its lifetime (Kb)",
  "t":
    "(Not in tcsh.) Average resident set size of the process (Kb)",
  "K":
    "Average total (data+stack+text) memory use of the process (Kb)",
  "D":
    "Average size of unshared data area (Kb)",
  "p":
    "Average size of unshared stack space (Kb)",
  "X":
    "Average size of shared text space (Kb)",
  "Z":
    "System page size (bytes)",
  "F":
    "Number of major page faults",
  "R":
    "Number of minor or recoverable page faults",
  "W":
    "Number of times the process was swapped out of main memory",
  "c":
    "Number of times the process was context-switched involuntarily",
  "w":
    "Number of waits",
  "I":
    "Number of file system inputs by the process",
  "O":
    "Number of file system outputs by the process",
  "r":
    "Number of socket messages received by the process",
  "s":
    "Number of socket messages sent by the process",
  "k":
    "Number of signals delivered to the process",
  "x":
    "(Not in tcsh.) Exit status of the command",
}

default_fields = [
  "e",
  "M",
  "K",
  "D",
  "X",
  "F",
  "R",
  "W",
  "w",
  "c",
  "S",
  "U",
  "P",
  "I",
  "O",
  "r",
  "s",
  "k",
  "x",
  ]

field_header_full = get_full_header(default_fields, usr_bin_time_csv_map) #','.join([ WrapperOpTimer.usr_bin_time_csv_map[f] for f in default_fields ])
field_header_short = ','.join(default_fields)
field_arg = '--format=' + field_header_full + '\n' + ','.join([ '%{}'.format(f) for f in default_fields] )

class WrapperOpTimer:

  @staticmethod
  def run_cmd(cmd):
    p = subprocess.Popen(cmd)
    p.communicate()
    returncode = p.returncode
    return returncode

  @staticmethod
  def time_op(wcp):
    """
      evaluate 'op' with 'op_args', and gather stats into output_stats_file
    """
    # if os.path.exists(output_stats_file) and os.path.getsize(output_stats_file) > 0:
    #   print("WARNING: File '"+output_stats_file+"' exists and will be overwritten")
    #   print("op='"+op+"'")
    #   print("op_args='"+str(op_args)+"'")
    #   print("op_output_file='"+op_output_file+"'")

    # initializing the titles and rows list
    fields = []
    csv_row = {}

    cmdcount = 0
    returncode = 0
    for cmd in wcp.commands:
      if cmdcount == 0:
        cmd = [ wcp.time_cmd,
                # '--append',
                '--output=' + wcp.output_stats_file,
                field_arg,
               ] + cmd
      cmdcount += 1
      returncode |= WrapperOpTimer.run_cmd(cmd)

    # reading csv file
    with open(wcp.output_stats_file, 'r') as csvfile:
      # creating a csv reader object
      csvreader = csv.reader(csvfile)

      # extracting field names through first row
      fields = next(csvreader)

      # extracting each data row one by one
      # we effectively retain only the last row.
      # it isn't clear if we should expect multiple rows per file
      #
      # In the bash version of this I was able to handle multiple rows per file
      # We could do that here, but it would require returning a list of csv maps
      # On the system side of things, it is very murky.  We would need to ensure
      # file integrity (concurrent reads/writes).  For now, it's
      # best to enforce 1 file per operation performed. (which should happen if we
      # name things correctly) - That is invalid is there is a cycle in the Build graph,
      # but that is a larger problem.
      for row in csvreader:
        csv_row = dict(zip(fields, row))

    # FileSize
    csv_row['FileSize'] = WrapperOpTimer.get_file_size(wcp.op_output_file)

    # add a field with the short op
    csv_row['op'] = os.path.basename(wcp.op)

    # FileName
    if wcp.base_build_dir:
      abs_base_build_dir = os.path.abspath(wcp.base_build_dir)
      current_working_dir = os.path.abspath(os.getcwd())
      rel_path_to_base_build_dir = os.path.relpath(
        current_working_dir, start=abs_base_build_dir)
      rel_op_output_file = os.path.join(rel_path_to_base_build_dir, wcp.op_output_file)
    else:
      rel_op_output_file = wcp.op_output_file
    csv_row['FileName'] = rel_op_output_file

    # Remove the build stats output file if the build failed
    if returncode != 0 and os.path.exists(wcp.output_stats_file):
      os.remove(wcp.output_stats_file)

    return (csv_row, returncode)


  # returns the file size in bytes
  @staticmethod
  def get_file_size(filename):
    sz = -1
    try:
      sz = os.stat(filename).st_size
    except:
      pass
    return sz


