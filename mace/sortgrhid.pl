#!/usr/bin/perl

#sortgrhid.pl sorts the lines of graph file according to the id.
#In the input graph file, a line corresponds to the list of a vertex, and the
# ID of the vertex is the first number in the line.
#sortgrhid.pl sorts makes the file so that the i-th line of the output
# graph file is the adjacent vertex list of vertex i, and remove the frst column, 
# which is the id of the line. the id is removed in the output file.

# One can specify the separator of the input file by the first parameter.

$ARGC = @ARGV;
if ( $ARGC < 0 ){
  printf ("sortgrhid.pl: [separator] < input-file > output-file\n");
  exit (1);
}
$count = 0;
%numbers = ();

$m=0;
$c=0;
$sep = " ";
if ( $ARGC > $c ){ $sep = $ARGV[$c]; } 

while (<STDIN>){
  chomp;
  @eles = split($sep, $_);
  $t[$eles[0]] = $_;
  if ( $m<$eles[0] ){ $m = $eles[0]; }
}

for ( $i=0 ; $i<=$m ; $i++ ){
  @eles = split($sep, $t[$i]);
  shift @eles;
  print "@eles";
}


