#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#
# checkTE.pl is a program that checks coordinates and translates formats
# of TEnest LTR output files. The use of this program is twofold, first is
# to translate *.LTR to *.gff (version 3) format to use TEnest repeat coordinates
# on  genome annotation display programs.  The second is to translate *.gff
# to *.LTR format to use third party or user identified repeat annotations
# with the TEnest svg_ltr display program, or to add gene annotations to the
# svg_ltr display.  
#
# Use of the display program svg_ltr requires exactly placed TE and gene
# basepair coordinates.  This script will check each coordinate to make sure
# there are no problems and prompt the user to make decisions for any problems.
#
# Any alteration of the *.LTR file requires a check of the coordinates to 
# make sure there are no basepair overlaps.  This program should be used
# anytime a *.LTR file is altered, even if two *.LTR files are combined
# use this file to translate from *.LTR to *.LTR to check all coordinates.
#
# If you input a gff format it is expected to be in the following format.
# There are 9 columns: project name, annotation program, type, start, end, .,
# direction, ., group. Project name and annotation program can be anything
# and will be changed to your output file name and TEnest when printing a gff
# output.  Start and end are number values, direction is a + or -.  Direction
# can also be shown by flipping the start and end and changing the + or -, ie
# 100 to 200 direction + is the same as 200 to 100 direction -.  Type is the 
# gff feature type, group contains the labels for each entry as well as the
# relational values for type.  Because it displays segmented (spliced) 
# annotations, it expects at least two lines for each entry.  Type 'TE' is 
# the whole TE annotation, type 'part' are each segmented annotation, ie if
# the parts of an annotation goes from 300-400 and 600-700, TE is 300-700.
# Group contains type relationships and notes. Here you will need to give
# TE annotation class (pair, solo, frag, nltr), TE annotation number, 
# TE type, and TE based coordinates.
# To differentiate between LTRs and internal retrotransposon regions PAIR
# needs to have one of three designations associated with it in the note section: for L for
# left LTR, R for right LTR and M for the middle region.
# The age of the LTR retrotransposon is put before the TE coordinates in BSR format : diguus 0.05-1-8000

#example  TEnest part  2      3109   .  +  .   TE FRAG f8; Note "huck 1-3000"
#example  TEnest TE    2      3109   .  +  .   TE FRAG f8; Note "huck 1-3000"
#example  TEnest part  3111   5111   .  -  .   TE PAIR p7; Note "diguus L-1-2000"
#example  TEnest part  12063  16100  .  -  .   TE PAIR p7; Note "diguus M-2000-6000"
#example  TEnest part  16101  18101  .  -  .   TE PAIR p7; Note "diguus R-6000-8000"
#example  TEnest TE    3111   18101  .  -  .   TE PAIR p7; Note "diguus 0.05-1-8000"
#
############################# OPTIONS ########################################

my $input = '';   # --input gff or ltr
my $output = '';  # --output gff or ltr

GetOptions('input=s' => \$input, 'output=s' => \$output);

##############################################################################

if($#ARGV < 1 || $input eq '' || $output eq '')
  {
  print "usage: LTR_checker.pl INPUT OUTPUT --input (gff or ltr) --output (gff or ltr)\n";
  exit;
  }
my $file = $ARGV[0];
my $out = $ARGV[1];
my $largest;
open(INPUT,$file), or die "$file not found\n";
my @ltr = ();
my @input = ();
while(<INPUT>)
  {
  my $line = $_;
  chomp($line);
  $line =~ s/\s+/ /g;
  $line =~ s/\;//g;
  $line =~ s/\"//g;
  if($line !~ m/^\#/)
    {push @input, [split(/ /,$line)];}
  }
close(INPUT);
if($input eq 'ltr')
  {@ltr = @input;}
if($input eq 'gff')
  {
  #TRANSLATE TO LTR FORMAT
  for(my $x=0;$x<@input;$x++)
    {
    if($input[$x][4] < $input[$x][3])
      {
      my $hold = $input[$x][4];
      $input[$x][4] = $input[$x][3];
      $input[$x][3] = $hold;
      if($input[$x][6] eq '-')
        {$input[$x][6] = '+';}
       else
        {$input[$x][6] = '-';}
      }
    if($input[$x][6] eq '-')
      {$input[$x][6] = 1;}
     else
      {$input[$x][6] = 0;}
    if($input[$x][9] eq 'PAIR')
      {($input[$x][13],$input[$x][14],$input[$x][15]) = split(/-/,$input[$x][13]);}
     else
      {($input[$x][13],$input[$x][14]) = split(/-/,$input[$x][13]);}
    }
  @input = sort {$a->[2] cmp $b->[2]} @input;
  @input = sort {$a->[10] cmp $b->[10]} @input;
  $largest = 0;
  for(my $x=0;$x<@input;$x++)
    {
    if($input[$x][3] > $largest)
      {$largest = $input[$x][3];}
    if($input[$x][4] > $largest)
      {$largest = $input[$x][4];}
    }
  $largest++;
  my $ltr_line = 0;
  my $end = @input;
  $ltr[$ltr_line][0] = $input[0][9];
  $ltr[$ltr_line][1] = $input[0][10];
  $ltr[$ltr_line][2] = $input[0][12];
  $ltr[$ltr_line][3] = $input[0][6];
  $ltr[$ltr_line][4] = 0;
  $ltr[$ltr_line][5] = 0;
  $ltr[$ltr_line][6] = 0;
  $ltr_line++;
  $ltr[$ltr_line][0] = $input[1][10];
  for(my $x=1;$x<@input;$x++)
    {
    my $set = 0;
    do{
      if($input[$x][9] eq 'FRAG' || $input[$x][9] eq 'NLTR' || $input[$x][9] eq 'SOLO')
        {
        $ltr[$ltr_line][($set*4)+1] = $input[$x][3];
        $ltr[$ltr_line][($set*4)+2] = $input[$x][4];
        $ltr[$ltr_line][($set*4)+3] = $input[$x][13];
        $ltr[$ltr_line][($set*4)+4] = $input[$x][14];
        $set++;
        }
      $x++;
      }until($x == $end || $input[$x][2] eq 'TE');
    if($x < $end)
      {
      if($input[$x][9] eq 'FRAG' || $input[$x][9] eq 'NLTR' || $input[$x][9] eq 'SOLO')
        {
        $ltr_line++;
        $ltr[$ltr_line][0] = $input[$x][9];
        $ltr[$ltr_line][1] = $input[$x][10];
        $ltr[$ltr_line][2] = $input[$x][12];
        $ltr[$ltr_line][3] = $input[$x][6];
        $ltr[$ltr_line][4] = 0;
        $ltr[$ltr_line][5] = 0;
        $ltr[$ltr_line][6] = 0;
        $ltr_line++;
        $ltr[$ltr_line][0] = $input[$x+1][10];
        }
      }
    }

#PAIRS
  my @pair = ();
  for(my $x=0;$x<@input;$x++)
    {
    if($input[$x][9] eq 'PAIR')
      {
      my $move = "$input[$x][0] $input[$x][1] $input[$x][2] $input[$x][3] $input[$x][4] $input[$x][5] $input[$x][6] $input[$x][7] $input[$x][8] $input[$x][9] $input[$x][10] $input[$x][11] $input[$x][12] $input[$x][13] $input[$x][14] $input[$x][15]";
      push @pair, [split(/ /,$move)];
      }
    }
  my $pair_count = @pair;
  my $last_pair = $pair[$pair_count-1][10];
  $last_pair =~ s/p//;
  for(my $y=0;$y<$last_pair+1;$y++)
    {
    my @pair_L = ();
    my @pair_R = ();
    my @pair_M = ();
    my $p = 'p'.$y;
    for(my $x=0;$x<@pair;$x++)
      {
      if($pair[$x][10] eq $p && $pair[$x][2] eq 'TE')
        {
        my $pair_line = "PAIR $p $pair[$x][12] $pair[$x][6] $pair[$x][13] 0 0 0";
        push @ltr, [split(/ /,$pair_line)];
        }
       elsif($pair[$x][10] eq $p)
        {
        if($pair[$x][13] eq 'L')
          {
          my $move = "$pair[$x][0] $pair[$x][1] $pair[$x][2] $pair[$x][3] $pair[$x][4] $pair[$x][5] $pair[$x][6] $pair[$x][7] $pair[$x][8] $pair[$x][9] $pair[$x][10] $pair[$x][11] $pair[$x][12] $pair[$x][13] $pair[$x][14] $pair[$x][15]";
          push @pair_L, [split(/ /,$move)];
          }
        if($pair[$x][13] eq 'R')
          {
          my $move = "$pair[$x][0] $pair[$x][1] $pair[$x][2] $pair[$x][3] $pair[$x][4] $pair[$x][5] $pair[$x][6] $pair[$x][7] $pair[$x][8] $pair[$x][9] $pair[$x][10] $pair[$x][11] $pair[$x][12] $pair[$x][13] $pair[$x][14] $pair[$x][15]";
          push @pair_R, [split(/ /,$move)];
          }
        if($pair[$x][13] eq 'M')
          {
          my $move = "$pair[$x][0] $pair[$x][1] $pair[$x][2] $pair[$x][3] $pair[$x][4] $pair[$x][5] $pair[$x][6] $pair[$x][7] $pair[$x][8] $pair[$x][9] $pair[$x][10] $pair[$x][11] $pair[$x][12] $pair[$x][13] $pair[$x][14] $pair[$x][15]";
          push @pair_M, [split(/ /,$move)];
          }
        }
      }
    my $pair_line = "$p L";
    for(my $x=0;$x<@pair_L;$x++)
      {$pair_line = $pair_line . " $pair_L[$x][3] $pair_L[$x][4] $pair_L[$x][14] $pair_L[$x][15]";}
    push @ltr, [split(/ /,$pair_line)];
    $pair_line = "$p R";
    for(my $x=0;$x<@pair_R;$x++)
      {$pair_line = $pair_line . " $pair_R[$x][3] $pair_R[$x][4] $pair_R[$x][14] $pair_R[$x][15]";}
    push @ltr, [split(/ /,$pair_line)];
    $pair_line = "$p M";
    for(my $x=0;$x<@pair_M;$x++)
      {$pair_line = $pair_line . " $pair_M[$x][3] $pair_M[$x][4] $pair_M[$x][14] $pair_M[$x][15]";}
    push @ltr, [split(/ /,$pair_line)];
    } 
  }

#Either:
#move gff @input to @ltr format
#move @ltr to hash format
#check coords while in hash

#or

#move gff @input to @ltr format
#check coords


#CHECK @ltr COORDS
 # decide whether to move to hash and check, or check in array

#PRINT OUTPUT
if($output eq 'ltr')
  {
  my $outfile = $out . ".LTR";
  open(OUT,">$outfile");
  print OUT "$largest\n";
  print OUT "$out\n";
  for(my $x=0;$x<@ltr;$x++)
    {
    print OUT "$ltr[$x][0]";
    for(my $y=1;$y<scalar @{$ltr[$x]};$y++)
      {print OUT " $ltr[$x][$y]";}
    print OUT "\n";
    }
  close(OUT);
  }
if($output eq 'gff')
  {
  my $outfile = $out . ".gff";
  open(OUT,">$outfile");
  print OUT "# file: $out\n";
  my $type;
  my $name;
  my $startSEQ;
  my $endSEQ;
  my $startTE;
  my $endTE;
  my $largeSEQ;
  my $smallSEQ;
  my $largeTE;
  my $smallTE;
  my $class;
  my $dir;
  for(my $x=0;$x<@ltr;$x++)
    {
    if($ltr[$x][0] eq 'SOLO' || $ltr[$x][0] eq 'FRAG' || $ltr[$x][0] eq 'NLTR')
      {
      $name = $ltr[$x][1];
      $type = $ltr[$x][2];
      $class = $ltr[$x][0];
      if($ltr[$x][3] == 0)
        {$dir = '+';}
       elsif($ltr[$x][3] == 1)
        {$dir = '-';}
      $x++;
      $smallSEQ = $ltr[$x][1];
      $smallTE = $ltr[$x][3];
      for(my $y=1;$y<scalar @{$ltr[$x]};$y=$y+4)
        {
        $startSEQ = $ltr[$x][$y];
        $endSEQ = $ltr[$x][$y+1];
        $startTE = $ltr[$x][$y+2];
        $endTE = $ltr[$x][$y+3];
        print OUT "$out\tTEnest\tpart\t$startSEQ\t$endSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type $startTE-$endTE\"\n";
        }
      $largeTE = $endTE;
      $largeSEQ = $endSEQ;
      print OUT "$out\tTEnest\tTE\t$smallSEQ\t$largeSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type $smallTE-$largeTE\" \n";
      }
    if($ltr[$x][0] eq 'PAIR')
      {
      my $age = $ltr[$x][4];
      $type = $ltr[$x][2];
      $name = $ltr[$x][1];
      $class = $ltr[$x][0];
      if($ltr[$x][3] == 0)
        {$dir = '+';}
       elsif($ltr[$x][3] == 1)
        {$dir = '-';}
      $x++;
      $smallSEQ = $ltr[$x][2];
      $smallTE = $ltr[$x][4];
      for(my $y=2;$y<scalar @{$ltr[$x]};$y=$y+4)
        {
        $startSEQ = $ltr[$x][$y];
        $endSEQ = $ltr[$x][$y+1];
        $startTE = $ltr[$x][$y+2];
        $endTE = $ltr[$x][$y+3];
        print OUT "$out\tTEnest\tpart\t$startSEQ\t$endSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type L-$startTE-$endTE\"\n";
        }
      $x++;
      for(my $y=2;$y<scalar @{$ltr[$x]};$y=$y+4)
        {
        $startSEQ = $ltr[$x][$y];
        $endSEQ = $ltr[$x][$y+1];
        $startTE = $ltr[$x][$y+2];
        $endTE = $ltr[$x][$y+3];
        print OUT "$out\tTEnest\tpart\t$startSEQ\t$endSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type R-$startTE-$endTE\"\n";
        }
      $largeTE = $endTE;
      $largeSEQ = $endSEQ;
      $x++;
      for(my $y=2;$y<scalar @{$ltr[$x]};$y=$y+4)
        {
        $startSEQ = $ltr[$x][$y];
        $endSEQ = $ltr[$x][$y+1];
        $startTE = $ltr[$x][$y+2];
        $endTE = $ltr[$x][$y+3];
        print OUT "$out\tTEnest\tpart\t$startSEQ\t$endSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type M-$startTE-$endTE\"\n";
        }
      print OUT "$out\tTEnest\tTE\t$smallSEQ\t$largeSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type $age-$smallTE-$largeTE\" \n";
      }
    if($ltr[$x][0] eq 'GENE' || $ltr[$x][0] eq 'PSDO')
      {
      $name = $ltr[$x][1];
      $type = $ltr[$x][2];
      $class = $ltr[$x][0];
      if($ltr[$x][3] == 0)
        {$dir = '+';}
       elsif($ltr[$x][3] == 1)
        {$dir = '-';}
      $x++;
      $smallSEQ = $ltr[$x][1];
      $smallTE = $ltr[$x][3];
      for(my $y=1;$y<scalar @{$ltr[$x]};$y=$y+4)
        {
        $startSEQ = $ltr[$x][$y];
        $endSEQ = $ltr[$x][$y+1];
        $startTE = $ltr[$x][$y+2];
        $endTE = $ltr[$x][$y+3];
        print OUT "$out\tTEnest\tpart\t$startSEQ\t$endSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type $startTE-$endTE\"\n";
        }
      $largeTE = $endTE;
      $largeSEQ = $endSEQ;
      print OUT "$out\tTEnest\tTE\t$smallSEQ\t$largeSEQ\t.\t$dir\t.\tTE $class $name; Note \"$type $smallTE-$largeTE\" \n";
      }
    }
  close(OUT);
  }

