#!/usr/bin/perl -w
#
# This perl script create an input file for the meshreader(...) function
# for only a single processor given the direct output from the Triangle
# program run as 'triangle -pe meshname.poly'.
#
use strict;
use Getopt::Long;
#
my $g_use_msg =
  "Use: make-serial-mesh.pl --input-base=base_name\n";
if( scalar(@ARGV) < 1 ) {
  print STDERR $g_use_msg;
  exit(1);
}
my $input_mesh_base_fn = "";
my $debug = 0;
GetOptions(
  "input-base=s"     => \$input_mesh_base_fn,
  "debug=i"          => \$debug
  );
$input_mesh_base_fn ne "" || die;
#
my $line;
my $i;
open MESH_OUT, ">${input_mesh_base_fn}.000";
#
print "\nRead in the ${input_mesh_base_fn}.node node file and output nodes ...\n\n";
#
open MESH_IN_NODE, "<${input_mesh_base_fn}.node" || die;
$line = <MESH_IN_NODE>;
#print "line = $line" if $debug;
my ($num_nodes, $dim, $num_attr, $num_bdry_markers) = split(" ",$line);
print "num_nodes = $num_nodes\n";
print "dim = $dim\n";
print "num_attr = $num_attr\n";
print "num_bdry_markers = $num_bdry_markers\n";
print MESH_OUT "${num_nodes} 0\n";
for( $i = 0; $i < $num_nodes; ++$i ) {
  $line = <MESH_IN_NODE>;
  my ($node_idx, $node_x, $node_y) = split(" ",$line);
  print "(node_idx,node_x,node_y) = ($node_idx,$node_x,$node_y)\n" if $debug;
  printf MESH_OUT "%d %e %e %d\n", $node_idx, $node_x, $node_y, 0;
}
close MESH_IN_NODE || die;
#
print "\nWrite the bogus node lines that metis would have written ...\n";
#
print MESH_OUT "0 0\n";
for( $i = 0; $i < $num_nodes; ++$i ) {
  printf MESH_OUT "%d %e %e %d\n", 0, 0, 0, 0;
}
#
print "\nRead the ${input_mesh_base_fn}.ele element file and write the elements ...\n\n";
#
open MESH_IN_ELE, "<${input_mesh_base_fn}.ele" || die;
$line = <MESH_IN_ELE>;
#print "line = $line" if $debug;
my ($num_ele, $num_nodes_per_ele, $num_attr2) = split(" ",$line);
print "num_ele = $num_ele\n";
print "num_nodes_per_ele = $num_nodes_per_ele\n";
print "num_attr2 = $num_attr2\n";
print MESH_OUT "${num_ele}\n";
for( $i = 0; $i < $num_ele; ++$i ) {
  $line = <MESH_IN_ELE>;
  my ($ele_idx,$node_1,$node_2,$node_3) = split(" ",$line);
  print "(ele_idx,node_1,node_2,node_3) = ($ele_idx,$node_1,$node_2,$node_3)\n" if $debug;
  printf MESH_OUT "%d %d %d\n", $node_1, $node_2, $node_3;
}
close MESH_IN_ELE || die;
#
print "\nRead the ${input_mesh_base_fn}.poly file which contains boundry edges and boundry markers ...\n\n";
#
open MESH_IN_POLY, "<${input_mesh_base_fn}.poly" || die;
$line = <MESH_IN_POLY>;
my ($num_nodes2, $dim2, $num_attr3, $num_bdry_markers2) = split(" ",$line);
print "num_nodes2 = $num_nodes2\n";
print "dim2 = $dim2\n";
print "num_attr3 = $num_attr3\n";
print "num_bdry_markers2 = $num_bdry_markers2\n";
$line = <MESH_IN_POLY>;
my ($num_bdry_edges, $num_bdry_markers3) = split(" ",$line);
print "num_bdry_edges = $num_bdry_edges\n";
print "num_bdry_markers3 = $num_bdry_markers3\n";
print MESH_OUT "${num_bdry_edges}\n";
for( $i = 0; $i < $num_bdry_edges; ++$i ) {
  $line = <MESH_IN_POLY>;
  my ($edge_idx,$node_1,$node_2,$bdry_marker) = split(" ",$line);
  print "(edge_idx,node_1,node_2,bdry_marker) = ($edge_idx,$node_1,$node_2,$bdry_marker)\n" if $debug;
  printf MESH_OUT "%d %d %d\n", $node_1, $node_2, $bdry_marker;
}
close MESH_IN_POLY || die;
#
close MESH_OUT || die;
print "\nWrote file ${input_mesh_base_fn}.000\n\n";
