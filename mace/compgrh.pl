#!/usr/bin/perl

# output the complement graph input from the standard input, to the
# standard output
# if b command is given, the graph is treated as a bipartite graph

$ARGC = @ARGV;
if ( $ARGC < 0 ){
  printf ("compgrh [b] [separator] < input-file > output-file\n");
  exit (1);
}
$count = 0;
$m = 0;

$sep = " ";
$bipartite = 0;
$c = 0;

if ( $ARGC > 0 ){
  if ( $ARGV[0] eq "b" ){ $bipartite = 1; $c++; }
  if ( $ARGC > $c ){ $sep = $ARGV[$c]; } 
}

while (<STDIN>){
  chomp;
  @eles = split($sep, $_);
  foreach $item( @eles ){
#    if ( $item == 0 ){
#      if (index ( $cell, "0") < 0 ){ next; }
#    }
    if ( $bipartite == 0 ){ push ( @{$t[$item]}, $count ); }
    push ( @{$t[$count]}, $item );
    if ( $item>$m ){ $m = $item;}
  }
  $count++;
}
$m++;
if ( $bipartite == 0 ){ $m = $count; }

for ( $i=0 ; $i<$count ; $i++ ){
  $jj=0;
  $flag =0;
  @all = sort { $a <=> $b }(@{$t[$i]});
  for ( $j=0 ; $j<$m ; $j++ ){
    if ( $all[$jj] == $j ){ while ( $all[$jj] == $j ){$jj++; }}
    else {
      if ( $bipartite==1 || $j > $i ){ 
        if ( $flag == 0 ){$flag = 1;}
        else { print " ";}
        print "$j";
      }
    }
  }
  print "\n";
}


