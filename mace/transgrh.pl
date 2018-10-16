#!/usr/bin/perl

# transgrh, transform a graph file of edge list type to adjacet matrix type
# i.e., convert a file
# (0,1) <== edge between vertex 0 and vertex 1
# (2,1) <== edge between vertex 2 and vertex 1
# (2,4) <== edge between vertex 3 and vertex 4
# ...
# to the file
# 1 2 5 9    <=  vertex 0 is adjacent to 1, 2, 5, 9
# 2        <=  vertex 1 is adjacent to 2
# 3 4 6 7  <=  vertex 2 is adjacent to 3, 4, 6, 7
# ...
# actually, vertex 1 and 2 are adjacent to vertex 0, but they are omitted
# since in line 0, there are 1 and 2. Similary, for an edge (x, y) 
# we put y on xth line if x<y, and vice versa.

# the graph file is input from standard input, and output to standard output
# if B or b is specified in the first parameter, the input graph is considered
# as a bipartite graph
# in B option, edge (0,1) is converted so that the 0th line has vertex 1
# in b option, edge (0,1) is converted so that the 1st line has vertex 0
# In both cases, the number to be written in a line will be added m, 
# which is the maximum number among all first numbers in all lines (all 
# second numbers, in B option
 
# if D or d is specified in the first parameter, the input graph is considered
# as a directed graph
# in D option, edge (0,1) is verted so that the 0th line has vertex 1
# in d option, edge (0,1) is converted so that the 1st line has vertex 0


$ARGC = @ARGV;
if ( $ARGC < 0 ){
  printf ("transgrp.pl: [BbdD] [separator] < input-file > output-file\n");
  exit (1);
}
@lines = <STDIN>
$count = 0;
%numbers = ();


$m = 0;
$m1 = 0;
$m2 = 0;
$sep = " ";
$c = 0;
$b = 0;
$d = 0;
if ( $ARGC > $c ){
  if ( index ( $ARGV[0], "b") >= 0 ){ $b = 1; $c = 1;}
  if ( index ( $ARGV[0], "B") >= 0 ){ $b = 2; $c = 1;}
  if ( index ( $ARGV[0], "d") >= 0 ){ $d = 1; $c = 1;}
  if ( index ( $ARGV[0], "D") >= 0 ){ $d = 2; $c = 1;}
  if ( $ARGC > $c ){ $sep = $ARGV[$c]; } 
}

if ( $b >0 ){
  foreach $trans( @lines ) {
    chomp ($trans);
    @eles = split($sep, $trans);
    if ( $eles[0] > $m1 ){ $m1 = $eles[0]; }
    if ( $eles[1] > $m2 ){ $m2 = $eles[1]; }
  }
}

if ( $b==1 ){ $m = $m1; }
if ( $b==2 ){ $m = $m2; }

foreach $trans( @lines ){
  chomp ($trans);
  $_ =~ s/$sep$sep/$sep/g;
  @eles = split($sep, $trans);
  if ( $b == 1 ){ push ( @{$t[$eles[0]]}, $eles[1]+$m1+1 ); }
  elsif ( $b == 2 ){ push ( @{$t[$eles[1]]}, $eles[0]+$m2+1 ); }
  else {
    if ( $d == 1 ){ push ( @{$t[$eles[0]]}, $eles[1] ); }
    elsif ( $d == 2 ){ push ( @{$t[$eles[1]]}, $eles[0] ); }
    elsif ( $eles[0] > $eles[1] ){ push ( @{$t[$eles[1]]}, $eles[0] ); }
    else { push ( @{$t[$eles[0]]}, $eles[1] ); }
    if ( $eles[0] > $m ){ $m = $eles[0]; }
    if ( $eles[1] > $m ){ $m = $eles[1]; }
  }
}

for ( $i=0 ; $i<=$m ; $i++ ){ print "@{$t[$i]}\n";}


