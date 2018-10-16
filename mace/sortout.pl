#!/usr/bin/perl

# input: a file each whose line is a sequence of numbers
# output: the file obtained by (1) sorting each line in increasing order, 
#   (2) sort the lines in the increasing order (as strings )
# the input file is read from standard input, and the result is output to 
# the standard output
# We can specify the separator of the number, by the first parameter

$ARGC = @ARGV;
if ( $ARGC < 0 ){
  printf ("sortout.pl: [separator] < input-file > output-file\nif separator is '-', do not sort the lines");
  exit (1);
}
$count = 0;
%numbers = ();
$m=0;
$linesort = 1;

$sep = " ";
if ( $ARGC > 0 ){
  if ( $ARGV[0] eq "-" ){ $linesort = 0; }
  else { $sep = $ARGV[0]; }
} 
while (<STDIN>){
  chomp;
  @eles = split($sep, $_);
  @eles = sort { $a <=> $b } (@eles);
  foreach $cell(@eles){
    if ( $cell == 0 ){
      if (index ( $cell, "0") >= 0 ){ $lines2[$m] .= $cell." "; }
    } else { $lines2[$m] .= $cell." "; }
  }
  $m++;
}
if ( $linesort == 1 ){ @lines2 = sort (@lines2); }
for ( $mm=0 ; $mm<$m ; $mm++ ){
  print "$lines2[$mm]\n";
}


