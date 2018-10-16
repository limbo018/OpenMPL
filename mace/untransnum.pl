#!/usr/bin/perl

# transform the numbers to the strings according to the table file
# basically, this script is for untransform the numbers transformed by 
# the transnum.pl

use Text::ParseWords;

$ARGC = @ARGV;
if ( $ARGC < 1 ){
      # error routine
  printf ("untransnum.pl: output-table-file [separator] < input-file > output-file\n");
  exit (1);
}
    # initialization
open (TABLEFILE, "<$ARGV[0]" );
%numbers = ();
$c=0;
$sep=" "; $sep2=" ";
if ( $ARGC >1 ){ $sep = $ARGV[1]; }  # separator
if ( $ARGC >2 ){ $sep2 = $ARGV[2]; }  # separator

    # read transform-table
while (<TABLEFILE>){
  chomp;
  $e0 = $_; $e1 = $_;
  $e0 =~ s/ .*$//;
  $e1 =~ s/^[^ ]* //;
  $numbers{$e0} = $e1;
  $c++;
}

    # read transform-file
while (<STDIN>){
  $_ =~ s/[\r\n]//g;
  chomp;
  $_ =~ s/$sep$sep/$sep/g;
  @eles = split($sep, $_);
  $all = @eles;
  $c = 0;
  foreach $item( @eles ){
	if ( $item < 0 ){
      print "*";
    } elsif (!exists $numbers{$item}) { 
      print "$item";
    } else { print "$numbers{$item}"; }
    $c++;
    if ($c<$all){ print $sep2;}
  }
  print "\n"
}
