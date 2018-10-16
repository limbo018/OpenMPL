#!/usr/bin/perl

# transnum assigns unique number beginning from 1 to each string appearing in
# the inpu-file, and convert all the strings to the numbers so that the same
# strings are given the same number
# the table of correspondence between numbers and strings is written
# to the table-file
# if -1, -2,..., or -9 is/are specified, ignore the specified columns, i.e., 
# the specified columns are not converted to numbers.

# if -n is given, all the numbers (words including no-alphabet) will be unified into one number

$ARGC = @ARGV;
if ( $ARGC < 1 ){
  printf ("transnum.pl: output-table-file [ignore-columns(-0,-1,-2,...)] [-n] [separator] < input-file > output-file\n");
  exit (1);
}
open ( TABLEFILE, ">$ARGV[0]" );
#@lines = <STDIN>;
$count = 1;
%numbers = ();
$c=1;
$numflag=0;

$sep = " ";
$ignore[1] = $ignore[2] = $ignore[3] = $ignore[4] = $ignore[5] = $ignore[6] = $ignore[7] = $ignore[8] = $ignore[9] = 0;

while ( $ARGC > $c ){
  if ( $ARGV[$c] eq "-1" ){ $ignore[1] = 1; }
  elsif ( $ARGV[$c] eq "-2" ){ $ignore[2] = 1; }
  elsif ( $ARGV[$c] eq "-3" ){ $ignore[3] = 1; }
  elsif ( $ARGV[$c] eq "-4" ){ $ignore[4] = 1; }
  elsif ( $ARGV[$c] eq "-5" ){ $ignore[5] = 1; }
  elsif ( $ARGV[$c] eq "-6" ){ $ignore[6] = 1; }
  elsif ( $ARGV[$c] eq "-7" ){ $ignore[7] = 1; }
  elsif ( $ARGV[$c] eq "-8" ){ $ignore[8] = 1; }
  elsif ( $ARGV[$c] eq "-9" ){ $ignore[9] = 1; }
  elsif ( $ARGV[$c] eq "-n" ){
    $numflag = 1;
    $count = 2;
    print TABLEFILE "1 __numbers__\n";
  }
  else { $sep = $ARGV[$c]; }
  $c++;
} 

while (<STDIN>){
  $_ =~ s/[\r\n]//g;
  chomp;
  $_ =~ s/$sep$sep/$sep/g;
  $_ =~ s/\t/$sep/g;
  @eles = split($sep, $_);
  $all = @eles;
  $c = 1;
  foreach $item( @eles ) {
    if ( $c < 10 && $ignore[$c] == 1 ){
      if ( $item ne "" ){
        print $item;
        if ($c<$all){ print " ";}
        $c++;
      }
    } else {
      if ( $item ne "" ){
        if ( $numflag==1 ){
          $item2 = $item;
          $item2 =~ s/[0-9]//g;
          if ( $item2 eq "" ){
            print "1";
            if ($c<$all){ print " ";}
            $c++;
            next;
          }
        }
        if (!exists $numbers{$item}) { 
          $numbers{$item} = $count;
          print TABLEFILE "$count $item\n";
          $count++;
        }
        print "$numbers{$item}";
        if ($c<$all){ print " ";}
        $c++;
      }
    }
  }
  print "\n"
}
