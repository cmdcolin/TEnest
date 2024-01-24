#!/usr/bin/perl -w

#############################################################################
#                                LICENSE                                    #
#############################################################################

# The TE nest software package, including TE_nest.pl, svg_ltr.pl and the
# associated repeat databases were written and are maintained by
# Brent Kronmiller (bak@iastate.edu).  Please email me with any questions,
# suggestions, or problems.

# TE nest is distributed under the GNU general public license.
# Please see http://www.gnu.org/licenses/gpl.txt for more information.

############################ UPDATES ########################################

# VERSION TE_nest_1.0
#5/17/07 - Added GNU license and version number

#5/4/07 added option to display just portions of *.LTR file
# also fixed thickness of DNA line, was too big in long sequences

#3/19/07 fixed naming for new LTR (actual TE name in *.LTR), and now correctly sorts
# legend, and correctly orders colors

############################# OPTIONS ########################################
use strict;
use POSIX qw(ceil floor);
use Getopt::Long;

#MAP FULL TEs WITH LTRs - Default, 0, MAPs each group.  Change to 1 or use option to turn off 
my $map_pair = 'T';  # --map_pair
#MAP SOLO LTRs
my $map_solo = 'T';  # --map_solo
#MAP FULL NON-LTR TEs
my $map_nltr = 'T';  # --map_nltr
#MAP FRAGMENT HITS
my $map_frag = 'T';  # --map_frag
#MAP GENES
my $map_gene = 'T';  # --map_gene
#MAP PSEUDO-GENES
my $map_psdo = 'T';  # --map_psdo
#PRINT COORDINATES ON TE INSERTION (HELPS TO SEE IF WHOLE TE IS PRESENT)
my $print_coords = 'F'; # Default is OFF, options are TE based (TE) or input sequence based (SEQ)
                      # while print_coords is OFF it over-rides individual options for group types
                      # --print_coords (TE or SEQ)
 #FOR PAIRs
 my $pair_coords = 'T'; # --pair_coords
 #FOR SOLOs
 my $solo_coords = 'T'; # --solo_coords
 #FOR FRAGs
 my $frag_coords = 'T'; # --frag_coords
 #FOR NLTRs
 my $nltr_coords = 'T'; # --nltr_coords
 #FOR GENEs
 my $gene_coords = 'T'; # --gene_coords
 #FOR PSUDO-GENESs
 my $psdo_coords = 'T'; # --psdo_coords
#SHOW UNIDENTIFIED AREAS IN TE HITS
my $white_out = 'T';    # --white_out
#SHOW BSR OR MYA - Default is mya (option --mya) select option --bsr to change 
my $mya = 'T';       # --mya
my $bsr = 'F';       # --bsr
#ONLY MAKE DISPLAY FOR A CERTAIN AREA
my $start_coord; # --start
my $end_coord;   # --end
my $coord_split = 'T';  # --split  F means if your entered display coords cut a TE, it won't be displayed
#LENGTH OF AREA WITHIN GENE ANNOTATION TO SHOW DIRECTION ARROW
my $gene_arrow_length = 250;
#DISTANCE BETWEEN SLASHES FOR PSUDEOGENES
my $slash = 150;
#LEGEND WIDTH
my $legend_width = 5;  # --width
#SCALE LENGTH
my $scale_size = 10000; # --scale
#TIC AMOUNTS
my $tic = 10;          # --tic
GetOptions('map_pair=s' => \$map_pair, 'map_solo=s' => \$map_solo, 'map_nltr=s' => \$map_nltr,
           'map_frag=s' => \$map_frag, 'white_out=s' => \$white_out, 'width=i' => \$legend_width,
           'print_coords=s' => \$print_coords, 'pair_coords=s' => \$pair_coords,
           'solo_coords=s' => \$solo_coords, 'nltr_coords=s' => \$nltr_coords, 'tic=i' => \$tic,
           'frag_coords=s' => \$frag_coords, 'scale=i' => \$scale_size, 'mya=s' => \$mya, 'bsr=s' => \$bsr,
           'start=i' => \$start_coord, 'end=i' => \$end_coord, 'split=s' => \$coord_split,
           'gene_coords=s' => \$gene_coords, 'psdo_coords=s' => \$psdo_coords, 'map_gene=s' => \$map_gene,
           'map_psdo=s' => \$map_psdo);
if($bsr eq 'T')  # If both were selected in command line, goes with BSR
  {$mya='F';}
##############################################################################

#LOAD HASH FILE
if ($#ARGV == -1)
  {
  print "enter *.LTR to graph\n";
  exit;
  }
my $hash_file = $ARGV[0];
my %retro = ();
my $total_inserts = 0;
if($start_coord || $end_coord)
  {($hash_file,$start_coord,$end_coord) = COORD_PRE($hash_file,$start_coord,$end_coord,$coord_split);}
open(HASH, $hash_file) or die "$hash_file not found\n";
my $line = <HASH>;
chomp($line);
my $total_seq = $line;
$total_seq++;
$line = <HASH>;
chomp($line);
my $fasta = $line;
if(!$start_coord)
  {$start_coord = 1;}
if(!$end_coord)
  {$end_coord = $total_seq;}
while(<HASH>)
  {
  my $line = $_;
  chomp($line);
  my @line = split(/ /,$line);
  if($line[0] eq 'PAIR')
    {
    if($map_pair eq 'T')
      {
      $total_inserts++;
      $retro{pair}{$line[1]} = {type => $line[2], dir => $line[3], bsr => $line[4]};
      for(my $y=0;$y<3;$y++)
        {
        $line = <HASH>;
        chomp($line);
        @line = split(/ /,$line);
        my $pair_count = @line;
        $pair_count = ($pair_count-2)/4;
        my $section_size = 0;
        for(my $x=0;$x<$pair_count;$x++)
          {
          $retro{pair}{$line[0]}{$line[1]}{$x} = {SEQ_start => $line[($x*4)+2], SEQ_end => $line[($x*4)+3], TE_start => $line[($x*4)+4], TE_end => $line[($x*4)+5]};
          $section_size = $section_size + $retro{pair}{$line[0]}{$line[1]}{$x}{SEQ_end} - $retro{pair}{$line[0]}{$line[1]}{$x}{SEQ_start};
          }
        if($section_size != 0)
          {
          $retro{pair}{$line[0]}{size}{$line[1]} = $section_size;
          }
        if($y==0)
          {
          $retro{pair}{$line[0]}{start} = $retro{pair}{$line[0]}{$line[1]}{0}{SEQ_start};
          }
         elsif($y==1)
          {
          $retro{pair}{$line[0]}{end} = $retro{pair}{$line[0]}{$line[1]}{$pair_count-1}{SEQ_end};
          }
        }
      }
    }
  }
close(HASH);
my $pair_amt = scalar keys %{$retro{pair}};
my @seq_spots = ();
my $hash_pair = scalar keys %{$retro{pair}};
my $x = my $found = 0;
do{
  my $p = 'p'.$x;
  if($retro{pair}{$p})
    {
    my $amt = 0;
    my $hit_amt = scalar keys %{$retro{pair}{$p}{L}};
    for(my $z=0;$z<$hit_amt;$z++)
      {
      $retro{pair}{$p}{coords}{$amt}{SEQ_start} = $retro{pair}{$p}{L}{$z}{SEQ_start};
      $retro{pair}{$p}{coords}{$amt}{SEQ_end} = $retro{pair}{$p}{L}{$z}{SEQ_end};
      $retro{pair}{$p}{coords}{$amt}{TE_start} = $retro{pair}{$p}{L}{$z}{TE_start};
      $retro{pair}{$p}{coords}{$amt}{TE_end} = $retro{pair}{$p}{L}{$z}{TE_end};
      $amt++;
      }
    $hit_amt = scalar keys %{$retro{pair}{$p}{M}};
    for(my $z=0;$z<$hit_amt;$z++)
      {
      $retro{pair}{$p}{coords}{$amt}{SEQ_start} = $retro{pair}{$p}{M}{$z}{SEQ_start};
      $retro{pair}{$p}{coords}{$amt}{SEQ_end} = $retro{pair}{$p}{M}{$z}{SEQ_end};
      $retro{pair}{$p}{coords}{$amt}{TE_start} = $retro{pair}{$p}{M}{$z}{TE_start};
      $retro{pair}{$p}{coords}{$amt}{TE_end} = $retro{pair}{$p}{M}{$z}{TE_end};
      $amt++;
      }
    $hit_amt = scalar keys %{$retro{pair}{$p}{R}};
    for(my $z=0;$z<$hit_amt;$z++)
      {
      $retro{pair}{$p}{coords}{$amt}{SEQ_start} = $retro{pair}{$p}{R}{$z}{SEQ_start};
      $retro{pair}{$p}{coords}{$amt}{SEQ_end} = $retro{pair}{$p}{R}{$z}{SEQ_end};
      $retro{pair}{$p}{coords}{$amt}{TE_start} = $retro{pair}{$p}{R}{$z}{TE_start};
      $retro{pair}{$p}{coords}{$amt}{TE_end} = $retro{pair}{$p}{R}{$z}{TE_end};
      $amt++;
      }
    for(my $z=0;$z<$amt;$z++)
      {
      my $seq_spots = "$retro{pair}{$p}{coords}{$z}{SEQ_start} $retro{pair}{$p}{coords}{$z}{SEQ_end}";
      push @seq_spots, [split(/ /,$seq_spots)];
      }
    $found++;
    }
  $x++;
  }until($found == $hash_pair);
my $seq_spots_count = @seq_spots;
open(HASH, $hash_file);
while(<HASH>)
  {
  my $line = $_;
  chomp($line);
  my @line = split(/ /,$line);
  if($line[0] eq 'SOLO')
    {
    if($map_solo eq 'T')
      {
      my $unique = 0;
      my $number = $line[1];
      my $type = $line[2];
      my $dir = $line[3];
      my $line = <HASH>;
      chomp($line);
      my @line = split(/ /,$line);
      my $solo_count = @line;
      $solo_count = ($solo_count-1)/4; 
      my $count = 0;
      for(my $x=0;$x<$solo_count;$x++)
        {
        my $dont_use = 0;
        for(my $y=0;$y<$seq_spots_count;$y++)
          {
          if($dont_use == 0)
            {
            if($line[($x*4)+1] < $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1] && $line[($x*4)+2] > $seq_spots[$y][0])
              {
              if($dir == 0)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+2] - $seq_spots[$y][0]);} 
               elsif($dir == 1)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+2] - $seq_spots[$y][0]);}
              $line[($x*4)+2] = $seq_spots[$y][0];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] > $seq_spots[$y][1] && $line[($x*4)+1] < $seq_spots[$y][1])
              {
              if($dir == 0)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+1] - $seq_spots[$y][1]);}
               elsif($dir == 1)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+1] - $seq_spots[$y][1]);}
              $line[($x*4)+1] = $seq_spots[$y][1];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1])
              {
              $dont_use = 1;
              }
            }
          }
        if($dont_use == 0 && ($line[($x*4)+2] - $line[($x*4)+1]) > 50)
          {
          $retro{solo}{$number}{coords}{$count} = {SEQ_start => $line[($x*4)+1], SEQ_end => $line[($x*4)+2], TE_start => $line[($x*4)+3], TE_end => $line[($x*4)+4]};
          $unique = 1;
          $count++;
          my $seq_spots = "$line[($x*4)+1] $line[($x*4)+2]";
          push @seq_spots, [split(/ /, $seq_spots)];          
          }
        }
      if($unique == 1)
        {
        $total_inserts++;
        $retro{solo}{$number}{type} = $type;
        $retro{solo}{$number}{dir} = $dir;
        $retro{solo}{$number}{start} = $retro{solo}{$number}{coords}{0}{SEQ_start};
        $retro{solo}{$number}{end} = $retro{solo}{$number}{coords}{$count-1}{SEQ_end};
        }
      }
    }
   elsif($line[0] eq 'FRAG')
    {
    if($map_frag eq 'T')
      {
      my $seq_spots_count = @seq_spots;
      my $unique = 0;
      my $number = $line[1];
      my $type = $line[2];
      my $dir = $line[3];
      my $line = <HASH>;
      chomp($line);
      my @line = split(/ /,$line);
      my $frag_count = @line;
      $frag_count = ($frag_count-1)/4;
      my $count = 0;
      for(my $x=0;$x<$frag_count;$x++)
        {
        my $dont_use = 0;
        for(my $y=0;$y<$seq_spots_count;$y++)
          {
          if($dont_use == 0)
            {
            if($line[($x*4)+1] < $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1] && $line[($x*4)+2] > $seq_spots[$y][0])
              {
              if($dir == 0)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+2] - $seq_spots[$y][0]);}
               elsif($dir == 1)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+2] - $seq_spots[$y][0]);}
              $line[($x*4)+2] = $seq_spots[$y][0];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] > $seq_spots[$y][1] && $line[($x*4)+1] < $seq_spots[$y][1])
              {
              if($dir == 0)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+1] - $seq_spots[$y][1]);}
               elsif($dir == 1)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+1] - $seq_spots[$y][1]);}
              $line[($x*4)+1] = $seq_spots[$y][1];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1])
              {
              $dont_use = 1;
              }
            }
          }
        if($dont_use == 0 && ($line[($x*4)+2] - $line[($x*4)+1]) > 50)
          {
          $retro{frag}{$number}{coords}{$count} = {SEQ_start => $line[($x*4)+1], SEQ_end => $line[($x*4)+2], TE_start => $line[($x*4)+3], TE_end => $line[($x*4)+4]};
          $unique = 1;
          $count++;
          }
        }
      if($unique == 1)
        {
        $total_inserts++;
        $retro{frag}{$number}{type} = $type;
        $retro{frag}{$number}{dir} = $dir;
        $retro{frag}{$number}{start} = $retro{frag}{$number}{coords}{0}{SEQ_start};
        $retro{frag}{$number}{end} = $retro{frag}{$number}{coords}{$count-1}{SEQ_end};
        }
      }
    }
   elsif($line[0] eq 'NLTR')
    {
    if($map_nltr eq 'T')
      {
      my $seq_spots_count = @seq_spots;
      my $unique = 0;
      my $number = $line[1];
      my $type = $line[2];
      my $dir = $line[3];
      my $line = <HASH>;
      chomp($line);
      my @line = split(/ /,$line);
      my $nltr_count = @line;
      $nltr_count = ($nltr_count-1)/4;
      my $count = 0;
      for(my $x=0;$x<$nltr_count;$x++)
        {
        my $dont_use = 0;
        for(my $y=0;$y<$seq_spots_count;$y++)
          {
          if($dont_use == 0)
            {
            if($line[($x*4)+1] < $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1] && $line[($x*4)+2] > $seq_spots[$y][0])
              {
              if($dir == 0)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+2] - $seq_spots[$y][0]);}
               elsif($dir == 1)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+2] - $seq_spots[$y][0]);}
              $line[($x*4)+2] = $seq_spots[$y][0];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] > $seq_spots[$y][1] && $line[($x*4)+1] < $seq_spots[$y][1])
              {
              if($dir == 0)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+1] - $seq_spots[$y][1]);}
               elsif($dir == 1)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+1] - $seq_spots[$y][1]);}
              $line[($x*4)+1] = $seq_spots[$y][1];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1])
              {
              $dont_use = 1;
              }
            }
          }
        if($dont_use == 0 && ($line[($x*4)+2] - $line[($x*4)+1]) > 50)
          {
          $retro{nltr}{$number}{coords}{$count} = {SEQ_start => $line[($x*4)+1], SEQ_end => $line[($x*4)+2], TE_start => $line[($x*4)+3], TE_end => $line[($x*4)+4]};
          $unique = 1;
          $count++;
          }
        }
      if($unique == 1)
        {
        $total_inserts++;
        $retro{nltr}{$number}{type} = $type;
        $retro{nltr}{$number}{dir} = $dir;
        $retro{nltr}{$number}{start} = $retro{nltr}{$number}{coords}{0}{SEQ_start};
        $retro{nltr}{$number}{end} = $retro{nltr}{$number}{coords}{$count-1}{SEQ_end};
        }
      }
    }
  }
close(HASH);
my %gene = ();
open(HASH, $hash_file) or die "$hash_file not found\n";
my $total_gene = 0;
while(<HASH>)
  {
  my $line = $_;
  chomp($line);
  my @line = split(/ /,$line);
  if($line[0] eq 'GENE')
    {
    if($map_gene eq 'T')
      {
      my $seq_spots_count = @seq_spots;
      my $unique = 0;
      my $number = $line[1];
      my $dir = $line[2];
      my $line_count = @line;
      my $type = '';
      for(my $x=3;$x<$line_count;$x++)
        {$type = $type . " " . $line[$x];}
      $type =~ s/^ //g;
      my $line = <HASH>;
      chomp($line);
      my @line = split(/ /,$line);
      my $gene_count = @line;
      $gene_count = ($gene_count-1)/4;
      my $count = 0;
      for(my $x=0;$x<$gene_count;$x++)
        {
        my $dont_use = 0;
        for(my $y=0;$y<$seq_spots_count;$y++)
          {
          if($dont_use == 0)
            {
            if($line[($x*4)+1] < $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1] && $line[($x*4)+2] > $seq_spots[$y][0])
              {
              if($dir == 0)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+2] - $seq_spots[$y][0]);}
               elsif($dir == 1)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+2] - $seq_spots[$y][0]);}
              $line[($x*4)+2] = $seq_spots[$y][0];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] > $seq_spots[$y][1] && $line[($x*4)+1] < $seq_spots[$y][1])
              {
              if($dir == 0)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+1] - $seq_spots[$y][1]);}
               elsif($dir == 1)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+1] - $seq_spots[$y][1]);}
              $line[($x*4)+1] = $seq_spots[$y][1];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1])
              {
              $dont_use = 1;
              }
            }
          }
        if($dont_use == 0 && ($line[($x*4)+2] - $line[($x*4)+1]) > 50)
          {
          $gene{gene}{$number}{coords}{$count} = {SEQ_start => $line[($x*4)+1], SEQ_end => $line[($x*4)+2], TE_start => $line[($x*4)+3], TE_end => $line[($x*4)+4]};
          $unique = 1;
          $count++;
          my $seq_spots = "$line[($x*4)+1] $line[($x*4)+2]";
          push @seq_spots, [split(/ /, $seq_spots)];
          }
        }
      if($unique == 1)
        {
        $total_gene++;
        $gene{gene}{$number}{type} = $type;
        $gene{gene}{$number}{dir} = $dir;
        $gene{gene}{$number}{start} = $gene{gene}{$number}{coords}{0}{SEQ_start};
        $gene{gene}{$number}{end} = $gene{gene}{$number}{coords}{$count-1}{SEQ_end};
        }
      }
    }
   elsif($line[0] eq 'PSDO')
    {
    if($map_gene eq 'T')
      {
      my $seq_spots_count = @seq_spots;
      my $unique = 0;
      my $number = $line[1];
      my $dir = $line[2];
      my $line_count = @line;
      my $type = '';
      for(my $x=3;$x<$line_count;$x++)
        {$type = $type . " " . $line[$x];}
      $type =~ s/^ //g;
      my $line = <HASH>;
      chomp($line);
      my @line = split(/ /,$line);
      my $psdo_count = @line;
      $psdo_count = ($psdo_count-1)/4;
      my $count = 0;
      for(my $x=0;$x<$psdo_count;$x++)
        {
        my $dont_use = 0;
        for(my $y=0;$y<$seq_spots_count;$y++)
          {
          if($dont_use == 0)
            {
            if($line[($x*4)+1] < $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1] && $line[($x*4)+2] > $seq_spots[$y][0])
              {
              if($dir == 0)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+2] - $seq_spots[$y][0]);}
               elsif($dir == 1)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+2] - $seq_spots[$y][0]);}
              $line[($x*4)+2] = $seq_spots[$y][0];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] > $seq_spots[$y][1] && $line[($x*4)+1] < $seq_spots[$y][1])
              {
              if($dir == 0)
                {$line[($x*4)+3] = $line[($x*4)+3] + ($line[($x*4)+1] - $seq_spots[$y][1]);}
               elsif($dir == 1)
                {$line[($x*4)+4] = $line[($x*4)+4] - ($line[($x*4)+1] - $seq_spots[$y][1]);}
              $line[($x*4)+1] = $seq_spots[$y][1];
              }
             elsif($line[($x*4)+1] > $seq_spots[$y][0] && $line[($x*4)+2] < $seq_spots[$y][1])
              {
              $dont_use = 1;
              }
            }
          }
        if($dont_use == 0 && ($line[($x*4)+2] - $line[($x*4)+1]) > 50)
          {
          $gene{psdo}{$number}{coords}{$count} = {SEQ_start => $line[($x*4)+1], SEQ_end => $line[($x*4)+2], TE_start => $line[($x*4)+3], TE_end => $line[($x*4)+4]};
          $unique = 1;
          $count++;
          }
        }
      if($unique == 1)
        {
        $total_gene++;
        $gene{psdo}{$number}{type} = $type;
        $gene{psdo}{$number}{dir} = $dir;
        $gene{psdo}{$number}{start} = $gene{psdo}{$number}{coords}{0}{SEQ_start};
        $gene{psdo}{$number}{end} = $gene{psdo}{$number}{coords}{$count-1}{SEQ_end};
        }
      }
    }
  }
close(HASH);

#DETERMINE GROUP, ORDER, LEVEL AGAIN FOR ALL DATA TYPES
#SOLO
my $hash_solo = scalar keys %{$retro{solo}};
my @gol = ();
$x = $found = 0;
do{
  my $s = 's'.$x;
  if($retro{solo}{$s})
    {
    my $solo_count = scalar keys %{$retro{solo}{$s}{coords}};
    my $gol = "solo $s $retro{solo}{$s}{coords}{0}{SEQ_start} $retro{solo}{$s}{coords}{$solo_count-1}{SEQ_end}";
    push @gol, [split(/ /,$gol)];
    $found++;
    }
  $x++;
  }until($found == $hash_solo);
#PAIR
$hash_pair = scalar keys %{$retro{pair}};
$x = $found = 0;
do{
  my $p = 'p'.$x;
  if($retro{pair}{$p})
    {
    my $pair_count = scalar keys %{$retro{pair}{$p}{R}};
    my $gol = "pair $p $retro{pair}{$p}{L}{0}{SEQ_start} $retro{pair}{$p}{R}{$pair_count-1}{SEQ_end}";
    push @gol, [split(/ /,$gol)];
    $found++;
    }
  $x++;
  }until($found == $hash_pair);
#FRAG
my $hash_frag = scalar keys %{$retro{frag}};
$x = $found = 0;
do{
  my $f = 'f'.$x;
  if($retro{frag}{$f})
    {
    my $frag_count = scalar keys %{$retro{frag}{$f}{coords}};
    my $gol = "frag $f $retro{frag}{$f}{coords}{0}{SEQ_start} $retro{frag}{$f}{coords}{$frag_count-1}{SEQ_end}";
    push @gol, [split(/ /,$gol)];
    $found++;
    }
  $x++;
  }until($found == $hash_frag);
#NLTR
my $hash_nltr = scalar keys %{$retro{nltr}};
$x = $found = 0;
do{
  my $n = 'n'.$x;
  if($retro{nltr}{$n})
    {
    my $nltr_count = scalar keys %{$retro{nltr}{$n}{coords}};
    my $gol = "nltr $n $retro{nltr}{$n}{coords}{0}{SEQ_start} $retro{nltr}{$n}{coords}{$nltr_count-1}{SEQ_end}";
    push @gol, [split(/ /,$gol)];
    $found++;
    }
  $x++;
  }until($found == $hash_nltr);
@gol = sort{$a->[2] <=> $b->[2]} @gol;
my $gol_count = @gol;
$gol[0][4] = $gol[0][5] = $gol[0][6] = my $zero = 0; my $largest_level = my $largest_order = my $largest_group= 0;
for(my $x=1;$x<$gol_count;$x++)
  {
  my $do = 0;
  my $y = $x - 1;
  while($gol[$y][0] && $do==0 && $y>=0)
    {
    if($gol[$x][3] < $gol[$y][3])
      {
      $gol[$x][4] = $gol[$y][4];
      $gol[$x][6] = $gol[$y][6] + 1;
      $do=1;
      }
     else
      {$y--;}
    }
  if($do==0)
    {
    $gol[$x][4] = $gol[$x-1][4] + 1;
    $gol[$x][6] = 0;
    }
  }
for(my $x=1;$x<$gol_count;$x++)
  {
  if($gol[$x][6] > $largest_level)
    {
    $largest_level = $gol[$x][6];
    }
  if($gol[$x][4] > $largest_group)
    {
    $largest_group = $gol[$x][4];
    }
  }
for(my $x=0;$x<$largest_level+1;$x++)
  {
  my $order = 0;
  for(my $y=0;$y<$gol_count;$y++)
    {
    if($gol[$y][6] eq $x)
      {
      $gol[$y][5] = $order;
      $order++;
      }
    }
  }
for(my $x=1;$x<$gol_count;$x++)
  {
  if($gol[$x][5] > $largest_order)
    {
    $largest_order = $gol[$x][5];
    }
  }
for(my $x=0;$x<$gol_count;$x++)
  {
  $retro{$gol[$x][0]}{$gol[$x][1]}{group} = $gol[$x][4];
  $retro{$gol[$x][0]}{$gol[$x][1]}{order} = $gol[$x][5];
  $retro{$gol[$x][0]}{$gol[$x][1]}{level} = $gol[$x][6];
  }

#WRITE HASH TO ARRAY
my @retro = ();
my $retro;
$retro{dna}{d0} = {group => -1, level => -1, order => 0, ipoint => ($total_seq-($start_coord-1))/2, start => 0, end => $total_seq};
my $for_count = scalar keys %{$retro{dna}};
for(my $p=0;$p<$for_count;$p++)
  {
  my $x = 'd'.$p;
  $retro = "dna $x . . $retro{dna}{$x}{group} $retro{dna}{$x}{order} $retro{dna}{$x}{level} $retro{dna}{$x}{start} $retro{dna}{$x}{end}";
  push @retro, [split(/ /, $retro)];
  $retro[0][9] = ($total_seq-($start_coord-1))/2;
  }
$for_count = scalar keys %{$retro{solo}};
my $done = my $t = 0;
if($for_count > 0)
  {
  do{
    my $x = 's'.$t;
    if($retro{solo}{$x})
      {
      $retro = "solo $x $retro{solo}{$x}{type} $retro{solo}{$x}{dir} $retro{solo}{$x}{group} $retro{solo}{$x}{order} $retro{solo}{$x}{level} $retro{solo}{$x}{start} $retro{solo}{$x}{end}";
      push @retro, [split(/ /, $retro)];
      $done++;
      }
    $t++;
    }until($done == $for_count);
  }
$for_count = scalar keys %{$retro{frag}};
$done = $t = 0;
if($for_count > 0)
  {
  do{
    my $x = 'f'.$t;
    if($retro{frag}{$x})
      {
      $retro = "frag $x $retro{frag}{$x}{type} $retro{frag}{$x}{dir} $retro{frag}{$x}{group} $retro{frag}{$x}{order} $retro{frag}{$x}{level} $retro{frag}{$x}{start} $retro{frag}{$x}{end}";
      push @retro, [split(/ /, $retro)];
      $done++;
      }
    $t++;
    }until($done == $for_count);
  }
$for_count = scalar keys %{$retro{nltr}};
$done = $t = 0;
if($for_count > 0)
  {
  do{
    my $x = 'n'.$t;
    if($retro{nltr}{$x})
      {
      $retro = "nltr $x $retro{nltr}{$x}{type} $retro{nltr}{$x}{dir} $retro{nltr}{$x}{group} $retro{nltr}{$x}{order} $retro{nltr}{$x}{level} $retro{nltr}{$x}{start} $retro{nltr}{$x}{end}";
      push @retro, [split(/ /, $retro)];
      $done++;
      }
    $t++;
    }until($done == $for_count);
  }
$for_count = scalar keys %{$retro{pair}};
$done = $t = 0;
if($for_count > 0)
  {
  do{
    my $x = 'p'.$t;
    if($retro{pair}{$x})
      {
      $retro = "pair $x $retro{pair}{$x}{type} $retro{pair}{$x}{dir} $retro{pair}{$x}{group} $retro{pair}{$x}{order} $retro{pair}{$x}{level} $retro{pair}{$x}{start} $retro{pair}{$x}{end}";
      push @retro, [split(/ /, $retro)];
      $done++;
      }
    $t++;
    }until($done == $for_count);
  }
my @gene = ();
$for_count = scalar keys %{$gene{gene}};
$done = $t = 0;
if($for_count > 0)
  {
  do{
    my $x = 'g'.$t;
    if($gene{gene}{$x})
      {
      my $amts = scalar keys %{$gene{gene}{$x}{coords}};
      for(my $y=0;$y<$amts;$y++)
        {
        $retro = "gene $x $gene{gene}{$x}{dir} $gene{gene}{$x}{coords}{$y}{SEQ_start} $gene{gene}{$x}{coords}{$y}{SEQ_end}";
        push @gene, [split(/ /, $retro)];
        }
      $done++;
      }
    $t++;
    }until($done == $for_count);
  }
$for_count = scalar keys %{$gene{psdo}};
$done = $t = 0;
if($for_count > 0)
  {
  do{
    my $x = 'u'.$t;
    if($gene{psdo}{$x})
      {
      my $amts = scalar keys %{$gene{psdo}{$x}{coords}};
      for(my $y=0;$y<$amts;$y++)
        {
        $retro = "psdo $x $gene{psdo}{$x}{dir} $gene{psdo}{$x}{coords}{$y}{SEQ_start} $gene{psdo}{$x}{coords}{$y}{SEQ_end}";
        push @gene, [split(/ /, $retro)];
        }
      $done++;
      }
    $t++;
    }until($done == $for_count);
  }

#DETERMINE INSERTION POINT IN PREVIOUS TE
@retro = sort{$a->[5] <=> $b->[5]} @retro;
@retro = sort{$a->[6] <=> $b->[6]} @retro;
for(my $x=-1;$x<=$largest_level;$x++)
  {
  my $zero_p = 0;
  for(my $y=0;$y<=$largest_group;$y++)
    {
    my $previous = 0;
    for(my $z=0;$z<=$largest_order;$z++)
      {
      for(my $a=0;$a<$total_inserts+1;$a++)
        {
        my $ipoint = 0;
        if($retro[$a][6] == $x && $retro[$a][4] == $y && $retro[$a][5] == $z)
          {
          if($retro[$a][6] == 0)
            {
            $ipoint = $retro[$a][7] - $zero_p;
            $retro[$a][9] = $ipoint;
            $retro{$retro[$a][0]}{$retro[$a][1]}{ipoint} = $ipoint;
            }
          for(my $b=0;$b<$total_inserts+1;$b++)
            {
            if($retro[$b][4] == $retro[$a][4] && $retro[$b][6]+1 == $retro[$a][6] &&
               $retro[$a][7] > $retro[$b][7] && $retro[$a][8] < $retro[$b][8])
              {
              for(my $c=0;$c<$total_inserts+1;$c++)
                {
                if($retro[$c][4] == $retro[$a-1][4] && $retro[$c][6]+1 == $retro[$a-1][6] &&
                   $retro[$a-1][7] > $retro[$c][7] && $retro[$a-1][8] < $retro[$c][8])
                  {
                  if($retro[$b][1] ne $retro[$c][1])
                    {$previous = 0;}
                  }
                }
              $ipoint = $retro[$a][7] - $retro[$b][7] - $previous;
              $retro[$a][9] = $ipoint;
              $retro{$retro[$a][0]}{$retro[$a][1]}{ipoint} =  $ipoint;
              }
            }
          $previous = $previous + $retro[$a][8] - $retro[$a][7];
          $zero_p = $zero_p + $retro[$a][8] - $retro[$a][7];
          }
        }
      }
    }
  }   

=cut
#PRINT OUT HASH TO MAKE SURE IT LOADED CORRECTLY
$for_count = scalar keys %{$retro{dna}};
for(my $p=0;$p<$for_count;$p++)
  {
  my $x = 'd'.$p;
  print "dna $x . . $retro{dna}{$x}{group} $retro{dna}{$x}{order} $retro{dna}{$x}{level} $retro{dna}{$x}{start} $retro{dna}{$x}{end} $retro{dna}{$x}{ipoint}\n";
  }
$done = $t = 0;
$for_count = scalar keys %{$retro{solo}};
do{
  my $x = 's'.$t;
  if($retro{solo}{$x})
    {
    print "solo $x $retro{solo}{$x}{type} $retro{solo}{$x}{dir} $retro{solo}{$x}{group} $retro{solo}{$x}{order} $retro{solo}{$x}{level} $retro{solo}{$x}{start} $retro{solo}{$x}{end} $retro{solo}{$x}{ipoint}\n";
    $done++;
    my $coord_count = scalar keys %{$retro{solo}{$x}{coords}};
    for(my $y=0;$y<$coord_count;$y++)
      {
      print "     $retro{solo}{$x}{coords}{$y}{SEQ_start} $retro{solo}{$x}{coords}{$y}{SEQ_end} $retro{solo}{$x}{coords}{$y}{TE_start} $retro{solo}{$x}{coords}{$y}{TE_end}\n";
      }
    }
  $t++;
  }until($done == $for_count);
$done = $t = 0;
$for_count = scalar keys %{$retro{frag}};
do{
  my $x = 'f'.$t;
  if($retro{frag}{$x})
    {
    print "frag $x $retro{frag}{$x}{type} $retro{frag}{$x}{dir} $retro{frag}{$x}{group} $retro{frag}{$x}{order} $retro{frag}{$x}{level} $retro{frag}{$x}{start} $retro{frag}{$x}{end} $retro{frag}{$x}{ipoint}\n";
    $done++;
    my $coord_count = scalar keys %{$retro{frag}{$x}{coords}};
    for(my $y=0;$y<$coord_count;$y++)
      {
      print "     $retro{frag}{$x}{coords}{$y}{SEQ_start} $retro{frag}{$x}{coords}{$y}{SEQ_end} $retro{frag}{$x}{coords}{$y}{TE_start} $retro{frag}{$x}{coords}{$y}{TE_end}\n";
      }
    }
  $t++; 
  }until($done == $for_count);
$done = $t = 0;
$for_count = scalar keys %{$retro{nltr}};
do{
  my $x = 'n'.$t;
  if($retro{nltr}{$x})
    {
    print "nltr $x $retro{nltr}{$x}{type} $retro{nltr}{$x}{dir} $retro{nltr}{$x}{group} $retro{nltr}{$x}{order} $retro{nltr}{$x}{level} $retro{nltr}{$x}{start} $retro{nltr}{$x}{end} $retro{nltr}{$x}{ipoint}\n";
    $done++;
    my $coord_count = scalar keys %{$retro{nltr}{$x}{coords}};
    for(my $y=0;$y<$coord_count;$y++)
      {
      print "     $retro{nltr}{$x}{coords}{$y}{SEQ_start} $retro{nltr}{$x}{coords}{$y}{SEQ_end} $retro{nltr}{$x}{coords}{$y}{TE_start} $retro{nltr}{$x}{coords}{$y}{TE_end}\n";
      }
    }
  $t++;
  }until($done == $for_count);
$done = $t = 0;
$for_count = scalar keys %{$retro{pair}};
do{
  my $x = 'p'.$t;
  if($retro{pair}{$x})
    {
    print "pair $x $retro{pair}{$x}{type} $retro{pair}{$x}{dir} $retro{pair}{$x}{bsr} $retro{pair}{$x}{group} $retro{pair}{$x}{order} $retro{pair}{$x}{level} $retro{pair}{$x}{start} $retro{pair}{$x}{end} $retro{pair}{$x}{ipoint}\n";
    $done++;
    my $coord_count = scalar keys %{$retro{pair}{$x}{coords}};
    for(my $y=0;$y<$coord_count;$y++)
      {
      print "     $retro{pair}{$x}{coords}{$y}{SEQ_start} $retro{pair}{$x}{coords}{$y}{SEQ_end} $retro{pair}{$x}{coords}{$y}{TE_start} $retro{pair}{$x}{coords}{$y}{TE_end}\n";
      }
    }
  $t++;
  }until($done == $for_count);
=cut

#DETERMINE REAL SIZE OF TE AND WHITE COORDs
for(my $x=0;$x<$total_inserts+1;$x++)
  {
  my $amt_hits = scalar keys %{$retro{$retro[$x][0]}{$retro[$x][1]}{coords}};
  for(my $z=0;$z<$amt_hits;$z++)
    {
    $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{WHITE_start} = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{SEQ_start};
    $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{WHITE_end} = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{SEQ_end};
    }
  my $size = $retro[$x][8] - $retro[$x][7];
  for(my $z=0;$z<$amt_hits;$z++)
    {
    for(my $a=0;$a<2;$a++)
      {
      my $subtract = 0; 
      my $hit;
      if($a==0)
        {$hit = 'WHITE_start';}
       elsif($a==1)
        {$hit = 'WHITE_end';}
      for(my $y=0;$y<$total_inserts+1;$y++)
        {
        if($retro[$x][6] == $retro[$y][6]-1 && $retro[$y][7] > $retro[$x][7] && $retro[$y][8] < $retro[$x][8])
          {
          if($retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{$hit} > (($retro[$y][8]-$retro[$y][7])/2)+$retro[$y][7])
            {
            $subtract = $subtract + ($retro[$y][8] - $retro[$y][7]);
            }
          }
        }
      $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{$hit} = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$z}{$hit} - $subtract;
      }
    }
  for(my $y=0;$y<$total_inserts+1;$y++)
    {
    if($retro[$x][6] == $retro[$y][6]-1 && $retro[$y][7] > $retro[$x][7] && $retro[$y][8] < $retro[$x][8])
      {
      $size = $size - ($retro[$y][8] - $retro[$y][7]);
      }
    }
  $retro[$x][10] = $size;
  }        
#GET COORDs and FIX OVERLAPPING TRIANGLES
@retro = sort{$a->[5] <=> $b->[5]} @retro; #sort by order 
@retro = sort{$a->[4] <=> $b->[4]} @retro; #sort by group
@retro = sort{$a->[6] <=> $b->[6]} @retro; #sort by level
$retro[0][13] = 0;
$retro[0][11] = $retro[0][13] - $retro[0][10]/2;
$retro[0][12] = $retro[0][13] + $retro[0][10]/2;
for(my $x=0;$x<$largest_level+1;$x++)
  {
  #GET COORDs, TRIANGLE START, END, MID
  for(my $y=0;$y<$total_inserts+1;$y++)
    {
    if($retro[$y][6] == $x)
      {
      for(my $z=$y;$z>=0;$z--)
        {
        my $z_begin = $retro[$z][7];
        my $z_end = $retro[$z][8];
        if($retro[$z][6] == -1 && $total_inserts == 1)
          {
          $z_begin = $retro[$z][7] + $start_coord;
          $z_end = $retro[$z][8] + $start_coord;
          }
#        print "y $y z $z $retro[$y][6]-1 == $retro[$z][6] && $retro[$y][7] > $z_begin && $retro[$y][8] < $z_end\n";
        if($retro[$y][6]-1 == $retro[$z][6] && $retro[$y][7] > $z_begin && $retro[$y][8] < $z_end)
          {
          $retro[$y][13] = $retro[$z][11] + $retro[$y][9];
          $retro[$y][11] = $retro[$y][13] - $retro[$y][10]/2; 
          $retro[$y][12] = $retro[$y][13] + $retro[$y][10]/2;
          }
        }
      }
    }
  #FIX OVERLAPPING REGIONS
  my $any_overlap;
  my $any_overlap2;
  do{
    $any_overlap = 0;
    for(my $y=0;$y<$total_inserts+1;$y++)
      {
      if($retro[$y][6] == $x)
        {
        do{
          $any_overlap2 = 0;
          for(my $z=0;$z<$total_inserts+1;$z++)
            {
            if($retro[$z][6] == $x)
              {
#print "$retro[$y][12] + 500 > $retro[$z][11] && $retro[$y][13] < $retro[$z][13] && $retro[$z][1] ne $retro[$y][1]\n";
              if($retro[$y][12] + 500 > $retro[$z][11] && $retro[$y][13] < $retro[$z][13] && $retro[$z][1] ne $retro[$y][1])
                {
                $any_overlap = $any_overlap2 = 1;
                my $overlap = $retro[$y][12] - $retro[$z][11];
                $overlap = 50;
                $retro[$y][11] = $retro[$y][11] - $overlap;
                $retro[$y][12] = $retro[$y][12] - $overlap;
                $retro[$z][11] = $retro[$z][11] + $overlap;
                $retro[$z][12] = $retro[$z][12] + $overlap;
                }
              }
            }
          }until($any_overlap2 == 0);
        }
      }
    }until($any_overlap == 0);
  }

#FIND BIGEST LEVEL, USE IT AS REFERENCE FOR SIZE OF PICTURE
@retro = sort{$a->[11] <=> $b->[11]} @retro;
my $window_start = $retro[0][11];
@retro = sort{$b->[12] <=> $a->[12]} @retro;
my $window_end = $retro[0][12];
my $length = $window_end - $window_start;
my $offset = -$window_start;
#UPDATE COORDS WITH OFFSET VALUE
for(my $x=0;$x<$total_inserts+1;$x++)
  {
  $retro[$x][11] = $retro[$x][11] + $offset;
  $retro[$x][12] = $retro[$x][12] + $offset;
  $retro[$x][13] = $retro[$x][13] + $offset;
  my $amt_hits = scalar keys %{$retro{$retro[$x][0]}{$retro[$x][1]}{coords}};
  for(my $y=0;$y<$amt_hits;$y++)
    {
    $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_start} = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_start} - $retro[$x][7] + $retro[$x][11];
    $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_end} = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_end} - $retro[$x][7] + $retro[$x][11];
    }
  }
@retro = sort{$a->[6] <=> $b->[6]} @retro;
for(my $x=0;$x<$total_inserts+1;$x++)
  {
  $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} = $retro[$x][11];
  $retro{$retro[$x][0]}{$retro[$x][1]}{Tend} = $retro[$x][12];
  $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid} = $retro[$x][13];
  }
my $height = $largest_level + 1; 
$height = $height * 5000;

my @gene_locs = ();
for(my $x=1;$x<$total_inserts+1;$x++)
  {
  my $amt_hits = scalar keys %{$retro{$retro[$x][0]}{$retro[$x][1]}{coords}};
  for(my $y=0;$y<$amt_hits;$y++)
    {
    my $ws = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_start};
    my $we = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_end};
    my $end = 0;
    if($y == 0)
      {$ws = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart};}
    if($y == $amt_hits-1)
      {
      $we = $retro{$retro[$x][0]}{$retro[$x][1]}{Tend};
      $end = 1;
      }
    my $gene_loc =  "$retro[$x][1] $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{SEQ_start} $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{SEQ_end} $retro{$retro[$x][0]}{$retro[$x][1]}{level} $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid} $ws $we $end";
    push @gene_locs, [split(/ /,$gene_loc)];
    }
  }
push @gene_locs, [split(/ /,"$retro[0][1] $retro[0][7] $retro[0][7] $retro[0][6] $retro[0][11] $retro[0][11] $retro[0][11] 0")];
push @gene_locs, [split(/ /,"$retro[0][1] $retro[0][8] $retro[0][8] $retro[0][6] $retro[0][11] $retro[0][12] $retro[0][12] 0")];
@gene_locs = sort{$a->[1] <=> $b->[1]} @gene_locs;
#print "$gene_locs[0][0] $gene_locs[0][1] $gene_locs[0][2] $gene_locs[0][3] $gene_locs[0][4] $gene_locs[0][5] $gene_locs[0][6]\n";
for(my $x=1;$x<@gene_locs;$x++)
  {
  for(my $y=0;$y<@gene;$y++)
    {
    if($gene[$y][3] > $gene_locs[$x-1][1] && $gene[$y][4] < $gene_locs[$x][2])
      {
      $gene[$y][5] = $gene_locs[$x-1][3];
      $gene[$y][6] = ($gene[$y][3] - $gene_locs[$x-1][2]) + $gene_locs[$x-1][6];
      $gene[$y][7] = ($gene[$y][4] - $gene_locs[$x-1][2]) + $gene_locs[$x-1][6];
      if($gene_locs[$x-1][7] == 1)
        {
        $gene[$y][5]--;
        $gene[$y][6] = ($gene[$y][3] - $gene_locs[$x-1][2]) + $gene_locs[$x-1][4];
        $gene[$y][7] = ($gene[$y][4] - $gene_locs[$x-1][2]) + $gene_locs[$x-1][4];
        }
#      print "  $gene[$y][1] $gene[$y][2] $gene[$y][3] $gene[$y][4] $gene[$y][5] $gene[$y][6] $gene[$y][7]\n";
      }
    }
#  print "$gene_locs[$x][0] $gene_locs[$x][1] $gene_locs[$x][2] $gene_locs[$x][3] $gene_locs[$x][4] $gene_locs[$x][5] $gene_locs[$x][6]\n";
  }

@gene = sort{$a->[3] <=> $b->[3]} @gene;
@gene = sort{$a->[1] cmp $b->[1]} @gene;
my $gene = '';
my @g_set = ();
my @g_mod = ();
for(my $y=0;$y<@gene;$y++)
  {
  if($gene[$y][1] eq $gene)
    {
    my $g_set = "$gene[$y][1] $gene[$y][2] $gene[$y][3] $gene[$y][4] $gene[$y][5] $gene[$y][6] $gene[$y][7]";
    push @g_set, [split(/ /,$g_set)];
    }
   else
    {
    if(@g_set > 0)
      {
      if($g_set[0][1] == 0)
        {
        my $g_set_count = @g_set;
        my $last = $g_set_count-1;
        my $g_mod = "modify $last $g_set[$last][0] $g_set[$last][1] $g_set[$last][2] $g_set[$last][3] $g_set[$last][4] $g_set[$last][5] $g_set[$last][6]";
        push @g_mod, [split(/ /,$g_mod)];
        }
      elsif($g_set[0][1] == 1)
        {
        my $last = 0;
        my $g_mod = "modify $last $g_set[$last][0] $g_set[$last][1] $g_set[$last][2] $g_set[$last][3] $g_set[$last][4] $g_set[$last][5] $g_set[$last][6]";
        push @g_mod, [split(/ /,$g_mod)];
        }
      }
    @g_set = ();
    my $g_set = "$gene[$y][1] $gene[$y][2] $gene[$y][3] $gene[$y][4] $gene[$y][5] $gene[$y][6] $gene[$y][7]";
    push @g_set, [split(/ /,$g_set)];
    }
  $gene = $gene[$y][1];
  }
if(@g_set > 0)
  {
  if($g_set[0][1] == 0)
    {
    my $g_set_count = @g_set;
    my $last = $g_set_count-1;
    my $g_mod = "modify $last $g_set[$last][0] $g_set[$last][1] $g_set[$last][2] $g_set[$last][3] $g_set[$last][4] $g_set[$last][5] $g_set[$last][6]";
    push @g_mod, [split(/ /,$g_mod)];
    }
  elsif($g_set[0][1] == 1)
    {
    my $last = 0;
    my $g_mod = "modify $last $g_set[$last][0] $g_set[$last][1] $g_set[$last][2] $g_set[$last][3] $g_set[$last][4] $g_set[$last][5] $g_set[$last][6]";
    push @g_mod, [split(/ /,$g_mod)];
    }
  }
my $g_set_count = @g_set;
my $g_mod_count = @g_mod;
#Change each in @g_mod to a trapizoid - rectangle with arrow - or triangle if less than arrow length
for(my $x=0;$x<@g_mod;$x++)
  {
  #Remove original exon from gene list
  my $remove = '1000000000000000000000000111';
  for(my $y=0;$y<@gene;$y++)
    {
    if($g_mod[$x][2] eq $gene[$y][1] && $g_mod[$x][4] == $gene[$y][3])
      {$remove = $y;}
    }
  my @gene_tmp = ();
  for(my $y=0;$y<@gene;$y++)
    { 
    if($y != $remove)
      {
      my $len = scalar @{$gene[$y]};
      my $gene = '';
      for(my $z=0;$z<$len;$z++)
        {$gene .= " $gene[$y][$z]";}
      $gene =~ s/^\s//g;
      push @gene_tmp, [split(/ /,$gene)];
      }
    }
  @gene = @gene_tmp;
  #Alter new exon to have arrow - add back to gene list
  my $type;
  if($g_mod[$x][2] =~ m/g/)
    {$type = 'gene'}
   elsif($g_mod[$x][2] =~ m/u/)
    {$type = 'psdo'}
  if($g_mod[$x][3] == 0)
    {
    if($g_mod[$x][5] - $g_mod[$x][4] > $gene_arrow_length)
      {
      my $g_sub = $g_mod[$x][8] - $gene_arrow_length;
      my $g_mod = "$type $g_mod[$x][2] $g_mod[$x][3] $g_mod[$x][4] $g_mod[$x][5] $g_mod[$x][6] $g_mod[$x][7] $g_mod[$x][8] $g_mod[$x][7] $g_sub $g_mod[$x][8]";
      push @gene, [split(/ /,$g_mod)];
      }
     else
      {
      my $g_mod = "$type $g_mod[$x][2] $g_mod[$x][3] $g_mod[$x][4] $g_mod[$x][5] $g_mod[$x][6] $g_mod[$x][7] $g_mod[$x][8] $g_mod[$x][7] $g_mod[$x][7] $g_mod[$x][8]";
      push @gene, [split(/ /,$g_mod)];
      }
    }
   elsif($g_mod[$x][3] == 1)
    {
    if($g_mod[$x][5] - $g_mod[$x][4] > $gene_arrow_length)
      {
      my $g_sub = $g_mod[$x][7] + $gene_arrow_length; 
      my $g_mod = "$type $g_mod[$x][2] $g_mod[$x][3] $g_mod[$x][4] $g_mod[$x][5] $g_mod[$x][6] $g_mod[$x][7] $g_mod[$x][8] $g_mod[$x][8] $g_sub $g_mod[$x][7]";
      push @gene, [split(/ /,$g_mod)];
      }
     else
      {
      my $g_mod = "$type $g_mod[$x][2] $g_mod[$x][3] $g_mod[$x][4] $g_mod[$x][5] $g_mod[$x][6] $g_mod[$x][7] $g_mod[$x][8] $g_mod[$x][8] $g_mod[$x][8] $g_mod[$x][7]";
      push @gene, [split(/ /,$g_mod)];
      }
    }
  }
@gene = sort{$a->[3] <=> $b->[3]} @gene;
@gene = sort{$a->[1] cmp $b->[1]} @gene;
for(my $x=0;$x<@gene;$x++)
  {$gene[$x][11] = $gene{$gene[$x][0]}{$gene[$x][1]}{type};}

#for(my $x=0;$x<$total_inserts+1;$x++)
#  {
#  print "$retro[$x][0] $retro[$x][1] $retro[$x][2] $retro[$x][3] $retro[$x][4] $retro[$x][5] $retro[$x][6] $retro[$x][7] $retro[$x][8] $retro[$x][9] $retro[$x][10] $retro[$x][11] $retro[$x][12] $retro[$x][13]\n";
#  }
#COLOR ARRAY
my @color_list = ('green', 'blue', 'red', 'orange', 'purple', 'firebrick', 'darkgrey', 'olive', 'slateblue', 'peru', 'teal', 'khaki', 'forestgreen', 'magenta', 'navy', 'lime', 'aqua', 'blueviolet', 'chartreuse', 'darkgoldenrod', 'deeppink', 'gold', 'lightgrey', 'yellow', 'indianred', 'darkorange', 'coral', 'crimson', 'cyan', 'green', 'blue', 'red', 'orange', 'purple', 'firebrick', 'darkgrey', 'olive', 'slateblue', 'peru', 'teal', 'khaki', 'forestgreen', 'magenta', 'navy', 'lime', 'aqua', 'blueviolet', 'chartreuse', 'darkgoldenrod', 'deeppink', 'gold', 'lightgrey', 'indianred', 'darkorange', 'coral', 'crimson', 'cyan', 'green', 'blue', 'red', 'orange', 'purple', 'firebrick', 'darkgrey', 'olive', 'slateblue', 'peru', 'teal', 'khaki', 'forestgreen', 'magenta', 'navy', 'lime', 'aqua', 'blueviolet', 'chartreuse', 'darkgoldenrod', 'deeppink', 'gold', 'lightgrey', 'yellow', 'indianred', 'darkorange', 'coral', 'crimson', 'cyan', 'green', 'blue', 'red', 'orange', 'purple', 'firebrick', 'darkgrey', 'olive', 'slateblue', 'peru', 'teal', 'khaki', 'forestgreen', 'magenta', 'navy', 'lime', 'aqua', 'blueviolet', 'chartreuse', 'darkgoldenrod', 'deeppink', 'gold', 'lightgrey', 'indianred', 'darkorange', 'coral', 'crimson', 'cyan', 'green', 'blue', 'red', 'orange', 'purple', 'firebrick', 'darkgrey', 'olive', 'slateblue', 'peru', 'teal', 'khaki', 'forestgreen', 'magenta', 'navy', 'lime', 'aqua', 'blueviolet', 'chartreuse', 'darkgoldenrod', 'deeppink', 'gold', 'lightgrey', 'yellow', 'indianred', 'darkorange', 'coral', 'crimson', 'cyan', 'green', 'blue', 'red', 'orange', 'purple', 'firebrick', 'darkgrey', 'olive', 'slateblue', 'peru', 'teal', 'khaki', 'forestgreen', 'magenta', 'navy', 'lime', 'aqua', 'blueviolet', 'chartreuse', 'darkgoldenrod', 'deeppink', 'gold', 'lightgrey', 'indianred', 'darkorange', 'coral', 'crimson', 'cyan');
#COLORS for TEs
my $color_found = 0;
my @te_color = ();
if($total_inserts > 0)
  {
  my @color_setA = ();
  my @color_set1 = ();
  for(my $x=1;$x<$total_inserts+1;$x++)
    {
    my $color_set = "$retro[$x][1] $retro[$x][2]";
    if($retro[$x][2] =~ m/^[0-9]/)
      {push @color_set1, [split(/ /, $color_set)];}
     elsif(lc($retro[$x][2]) =~ m/^[a-z]/)
      {push @color_setA, [split(/ /, $color_set)];}
    }
  @color_setA = sort{$a->[1] cmp $b->[1]} @color_setA;
  @color_set1 = sort{$a->[1] <=> $b->[1]} @color_set1;
  my @color_set = ();
  for(my $x=0;$x<@color_setA;$x++)
    {
    my $color_set = "$color_setA[$x][0] $color_setA[$x][1]";
    push @color_set, [split(/ /,$color_set)];
    }
  for(my $x=0;$x<@color_set1;$x++)
    {
    my $color_set = "$color_set1[$x][0] $color_set1[$x][1]";
    push @color_set, [split(/ /,$color_set)];
    }
  $te_color[0][0] = $color_set[0][1];
  $te_color[0][1] = $color_list[0];
  $color_found = 1;
  for(my $x=1;$x<$total_inserts;$x++)
    {
    if($color_set[$x-1][1] ne $color_set[$x][1])
      {
      $te_color[$color_found][0] = $color_set[$x][1];
      $te_color[$color_found][1] = $color_list[$color_found];
      $color_found++;
      }
    }
  }
#COLORS for genes
my $gene_color_found = 0;
my @gene_color = ();
if(@gene > 0)
  {
  my @color_set = ();
  my $color_set_count = 0;
  for(my $x=0;$x<@gene;$x++)
    {
    $color_set[$color_set_count][0] = $gene[$x][1];
    $color_set[$color_set_count][1] = $gene[$x][11];
    $color_set_count++;
    }
  @color_set = sort{$a->[1] cmp $b->[1]} @color_set;
  $gene_color[0][0] = $color_set[0][1];
  $gene_color[0][1] = $color_list[0];
  $gene_color_found = 1;
  for(my $x=0;$x<@gene;$x++)
    {
    if($color_set[$x-1][1] ne $color_set[$x][1])
      {
      $gene_color[$gene_color_found][0] = $color_set[$x][1];
      $gene_color[$gene_color_found][1] = $color_list[$gene_color_found];
      $gene_color_found++;
      }
    }
  }

#MAKE PICTURE LEGEND
my $legend_depth = sprintf(ceil($color_found / $legend_width));
my $add = $color_found * 1500 / $legend_width;
my $legend = $height + 4000 + $add;
#MAKE SVG FILE
my $svg_file = "$fasta.svg";
open(SVG, ">$svg_file");
print SVG "<?xml version=\"1.0\"  standalone=\"yes\"?>\n";
print SVG "<svg width=\"100%\" height=\"100%\" version=\"1.1\" viewBox = \"0,-1500 $length,$legend\" xmlns=\"http://www.w3.org/2000/svg\">\n";
print SVG "<rect x=\"0\" y=\"-500\" width=\"$length\" height=\"$legend\" fill=\"white\" stroke=\"black\" />\n";
my $x_loc_box = 1000;
my $x_loc_txt = 2100;
for(my $x=0;$x<$legend_width;$x++)
  {
  my $y_loc_box = $height + 1500;
  my $y_loc_txt = $y_loc_box + 700;
  for(my $y=0;$y<$legend_depth;$y++)
    {
    if(($x*$legend_depth)+$y < $color_found)
      {
      print SVG "<rect fill=\"$te_color[($x*$legend_depth)+$y][1]\" x=\"$x_loc_box\" y=\"$y_loc_box\" width=\"1000\" height=\"1000\"/>\n";
      print SVG "<text x=\"$x_loc_txt\" y=\"$y_loc_txt\" fill=\"black\" font-size=\"1000\"> $te_color[($x*$legend_depth)+$y][0]</text>\n";
      }
    $y_loc_box = $y_loc_box + 1500;
    $y_loc_txt = $y_loc_txt + 1500;
    }
  $x_loc_box = $x_loc_box + 7500;
  $x_loc_txt = $x_loc_txt + 7500;
  }
my $scale_start = $x_loc_box + 2000;
my $scale_end = $scale_start + $scale_size;
my $scale_height = $height + 1500;
my $bar_up = $scale_height + 500;
my $bar_down = $scale_height - 500;
print SVG "<line stroke=\"black\" stroke-width=\"200\" x1=\"$scale_start\" y1=\"$scale_height\" x2=\"$scale_end\" y2=\"$scale_height\"/>\n";
print SVG "<line stroke=\"black\" stroke-width=\"200\" x1=\"$scale_start\" y1=\"$bar_up\" x2=\"$scale_start\" y2=\"$bar_down\"/>\n";
print SVG "<line stroke=\"black\" stroke-width=\"200\" x1=\"$scale_end\" y1=\"$bar_up\" x2=\"$scale_end\" y2=\"$bar_down\"/>\n";
my $tic_dist = $scale_size / $tic;
for(my $x=0;$x<$tic;$x++)
  {
  my $tic_loc = $scale_start + ($x*$tic_dist);
  print SVG "<line stroke=\"black\" stroke-width=\".1%\" x1=\"$tic_loc\" y1=\"$bar_up\" x2=\"$tic_loc\" y2=\"$bar_down\"/>\n";
  }
$bar_down = $bar_down + 2000;
print SVG "<text x=\"$scale_start\" y=\"$bar_down\" fill=\"black\" font-size=\"1000\">$scale_size bp</text>\n";
print SVG "<line stroke=\"black\" stroke-width=\"250\" x1=\"$retro{dna}{d0}{Tstart}\" y1=\"$height\" x2=\"$retro{dna}{d0}{Tend}\" y2=\"$height\"/>\n";
for(my $x=1;$x<$total_inserts+1;$x++)
  {
  my $top = ($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000;
  my $bottom = ($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}+1) * 5000;
  print SVG "<polygon stroke=\"black\" stroke-width=\"150\" points=\"$retro{$retro[$x][0]}{$retro[$x][1]}{Tstart},$top $retro{$retro[$x][0]}{$retro[$x][1]}{Tend},$top $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid},$bottom\"/>\n";
  for(my $y=0;$y<$color_found;$y++)
    {
    if($te_color[$y][0] eq $retro[$x][2])
      {
      print SVG "<polygon fill=\"$te_color[$y][1]\" points=\"$retro{$retro[$x][0]}{$retro[$x][1]}{Tstart},$top $retro{$retro[$x][0]}{$retro[$x][1]}{Tend},$top $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid},$bottom\"/>\n";
      }
    }
  }
#WHITE OUT NON-HITS
if($white_out eq 'T')
  {
  for(my $x=1;$x<$total_inserts+1;$x++)
    {
    my $amt_hits = scalar keys %{$retro{$retro[$x][0]}{$retro[$x][1]}{coords}};
    my $top = ($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000;
    my $bottom = ($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}+1) * 5000;
    for(my $y=1;$y<$amt_hits;$y++)
      {
      print SVG "<polygon fill=\"white\" points=\"$retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y-1}{WHITE_end},$top $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_start},$top $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid},$bottom\"/>\n";
      }    
    }
  }

#DRAW GENES AND PSUEDOGENES
for(my $x=0;$x<@gene;$x++)
  {
  for(my $y=0;$y<$gene_color_found;$y++)
    {
    if($gene_color[$y][0] eq $gene[$x][11])
      {
      my $top = (($largest_level - $gene[$x][5]) * 5000) - 300;
      my $size = $gene[$x][7] - $gene[$x][6];
      my $mid = $top + 300;
      my $bottom = $top + 600;
      if($gene[$x][8])
        {
        print SVG "<polygon stroke=\"black\" stroke-width=\"150\" points=\"$gene[$x][8],$top $gene[$x][9],$top $gene[$x][10],$mid $gene[$x][9],$bottom $gene[$x][8],$bottom\"/>\n";
        print SVG "<polygon fill=\"$gene_color[$y][1]\" points=\"$gene[$x][8],$top $gene[$x][9],$top $gene[$x][10],$mid $gene[$x][9],$bottom $gene[$x][8],$bottom\"/>\n";
        if($gene[$x][0] eq 'psdo')
          {
          my $slash_count = $size/$slash;
          my $for_start = 0;
          my $for_end = $slash_count;
          if($gene[$x][2] == 0)
            {$for_end--;}
           elsif($gene[$x][2] == 1)
            {$for_start++;}
          for(my $z=$for_start;$z<$for_end;$z++)
            {
            my $start = $gene[$x][6] + ($z * $slash);
            my $end = $gene[$x][6] + $slash + ($z * $slash);
            print SVG "<line stroke=\"white\" stroke-width=\"35\" x1=\"$start\" y1=\"$top\" x2=\"$end\" y2=\"$bottom\"/>\n";
            }
          }
        }
       else
        {
        print SVG "<rect stroke=\"black\" stroke-width=\"150\" x=\"$gene[$x][6]\" y=\"$top\" width=\"$size\" height=\"600\"/>\n";
        print SVG "<rect fill=\"$gene_color[$y][1]\" x=\"$gene[$x][6]\" y=\"$top\" width=\"$size\" height=\"600\"/>\n";
        if($gene[$x][0] eq 'psdo')
          {
          my $slash_count = $size/$slash;
          for(my $z=0;$z<$slash_count;$z++)
            {
            my $start = $gene[$x][6] + ($z * $slash);
            my $end = $gene[$x][6] + $slash + ($z * $slash);
            print SVG "<line stroke=\"white\" stroke-width=\"35\" x1=\"$start\" y1=\"$top\" x2=\"$end\" y2=\"$bottom\"/>\n";
            }
          }
        }
      }
    }
  }
 #CONNECT GENE EXONS
@gene = sort{$a->[6] <=> $b->[6]} @gene;
@gene = sort{$a->[1] cmp $b->[1]} @gene;
$gene = $gene[0][1];
my $start = $gene[0][7];  
for(my $x=1;$x<@gene;$x++)
  {
  if($gene[$x][1] eq $gene)
    {
    my $bottom = (($largest_level - $gene[$x-1][5]) * 5000) - 300;
    my $top = (($largest_level - $gene[$x-1][5]) * 5000) - 800;
    my $end = $gene[$x][6];
    my $half = ($end - $start)/2 + $start;
    print SVG "<polyline stroke=\"black\" stroke-width=\"150\"  points=\"$start,$bottom $half,$top $end,$bottom\" style=\"fill:none\"/>\n";
    }
  $gene = $gene[$x][1];
  $start = $gene[$x][7];  
  }

 #PRINT GENE COORDS
#if($print_coords eq 'TE' || $print_coords eq 'SEQ')
#  {
#  my $print_coords_start = $print_coords . '_start';
#  my $print_coords_end = $print_coords . '_end';
#  my @to_print = ('','','','');
#  if($pair_coords eq 'T')
#    {$to_print[0] = 'pair';}
#  if($nltr_coords eq 'T')
#    {$to_print[1] = 'nltr';}
#  if($solo_coords eq 'T')
#    {$to_print[2] = 'solo';}
#  if($frag_coords eq 'T')
#    {$to_print[3] = 'frag';}
#  for(my $x=0;$x<$total_inserts+1;$x++)
#    {
#    for(my $y=0;$y<4;$y++)
#      {
#      if($retro[$x][0] eq $to_print[$y])
#        {
#        my $textys = (($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000) - 500;
#        my $textye = (($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000) - 1000;
#        my $amt_hits = scalar keys %{$retro{$retro[$x][0]}{$retro[$x][1]}{coords}};
#        for(my $y=0;$y<$amt_hits;$y++)
#          {
#          my $textxs = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_start} - 500;
#          my $textxe = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_end} - 500;
#          print SVG "<text x=\"$textxs\" y=\"$textys\" fill=\"black\" font-size=\"500\"> $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{$print_coords_start}</text>\n";
#          print SVG "<text x=\"$textxe\" y=\"$textye\" fill=\"black\" font-size=\"500\"> $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{$print_coords_end}</text>\n";
#          }
#        }
#      }
#    }
#  }

#PRINT COORDs
if($print_coords eq 'TE' || $print_coords eq 'SEQ')
  {
  my $print_coords_start = $print_coords . '_start';
  my $print_coords_end = $print_coords . '_end';
  my @to_print = ('','','','');
  if($pair_coords eq 'T')
    {$to_print[0] = 'pair';}
  if($nltr_coords eq 'T')
    {$to_print[1] = 'nltr';}
  if($solo_coords eq 'T')
    {$to_print[2] = 'solo';}
  if($frag_coords eq 'T')
    {$to_print[3] = 'frag';}
  for(my $x=0;$x<$total_inserts+1;$x++)
    {
    for(my $y=0;$y<4;$y++)
      {
      if($retro[$x][0] eq $to_print[$y])
        {
        my $textys = (($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000) - 500;
        my $textye = (($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000) - 1000;
        my $amt_hits = scalar keys %{$retro{$retro[$x][0]}{$retro[$x][1]}{coords}};
        for(my $y=0;$y<$amt_hits;$y++)
          {
          my $textxs = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_start} - 500;
          my $textxe = $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{WHITE_end} - 500;
          print SVG "<text x=\"$textxs\" y=\"$textys\" fill=\"black\" font-size=\"500\"> $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{$print_coords_start}</text>\n";
          print SVG "<text x=\"$textxe\" y=\"$textye\" fill=\"black\" font-size=\"500\"> $retro{$retro[$x][0]}{$retro[$x][1]}{coords}{$y}{$print_coords_end}</text>\n";
          }
        }
      }
    }
  }

#PUT IN TEXT AND BOXES FOR BSR (MYA)
for(my $x=1;$x<$total_inserts+1;$x++)
  {
  if($retro[$x][0] eq 'pair')
    {
    my $top = ($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000;
    my $textx = 0;
    if(($retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart})/2 + $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} < $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid})
      {
      $textx = (($retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart}) /2) + $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart};
      }
     elsif(($retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart})/2 + $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} > $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid})
      {
      $textx = (($retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart}) /2) + $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} - 2000;
      }
     elsif(($retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart})/2 + $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} == $retro{$retro[$x][0]}{$retro[$x][1]}{Tmid})
      {
      $textx = (($retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart}) /2) + $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} - 1500;
      }
    my $texty = $top + 250;
    print SVG "<rect stroke-width=\"150\" x=\"$textx\" y=\"$texty\" width=\"3000\" height=\"1000\" fill=\"black\" stroke=\"black\" />\n";
    for(my $y=0;$y<$color_found;$y++)
      {
      if($te_color[$y][0] eq $retro[$x][2])
        {
        print SVG "<rect x=\"$textx\" y=\"$texty\" width=\"3000\" height=\"1000\" fill=\"$te_color[$y][1]\" stroke=\"black\" />\n";
        }
      }
    $textx = $textx + 250;
    $texty = $texty + 875;
    my $bsr_mya = $retro{$retro[$x][0]}{$retro[$x][1]}{bsr};
    if($mya eq 'T')
      {
      $bsr_mya = $bsr_mya/.026;
      $bsr_mya = sprintf("%.3f", $bsr_mya);
      }
    print SVG "<text x=\"$textx\" y=\"$texty\" fill=\"white\" font-size=\"1000\"> $bsr_mya </text>\n";
    #LTR ARROWS
    my $arrow_top = $top + 500;
    my $arrow_bot = $top - 500;
    my $arrow_start;
    my $arrow_end;
    my $Lltr_size = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro{$retro[$x][0]}{$retro[$x][1]}{size}{L};
    my $Rltr_size = $retro{pair}{$retro[$x][1]}{Tend} - $retro{pair}{$retro[$x][1]}{size}{R};
    #LEFT ARROW
    print SVG "<line stroke=\"black\" stroke-width=\"200\" x1=\"$retro{$retro[$x][0]}{$retro[$x][1]}{Tstart}\" y1=\"$top\" x2=\"$Lltr_size\" y2=\"$top\"/>\n";
    if($retro[$x][3] == 0)
      {
      $arrow_start = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro{$retro[$x][0]}{$retro[$x][1]}{size}{L} - 400;
      $arrow_end = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro{$retro[$x][0]}{$retro[$x][1]}{size}{L};
      print SVG "<polygon fill=\"black\" points=\"$arrow_start,$arrow_top $arrow_start,$arrow_bot $arrow_end,$top\"/>\n";
      }
     elsif($retro[$x][3] == 1)
      {
      $arrow_start = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart};
      $arrow_end = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + 400;
      print SVG "<polygon fill=\"black\" points=\"$arrow_start,$top $arrow_end,$arrow_bot $arrow_end,$arrow_top\"/>\n";
      }
    #RIGHT ARROW
    print SVG "<line stroke=\"black\" stroke-width=\"200\" x1=\"$Rltr_size\" y1=\"$top\" x2=\"$retro{$retro[$x][0]}{$retro[$x][1]}{Tend}\" y2=\"$top\"/>\n";
    if($retro[$x][3] == 0)
      {
      $arrow_start = $retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - 400;
      $arrow_end = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro{$retro[$x][0]}{$retro[$x][1]}{size}{L};
      print SVG "<polygon fill=\"black\" points=\"$arrow_start,$arrow_top $arrow_start,$arrow_bot $retro{$retro[$x][0]}{$retro[$x][1]}{Tend},$top\"/>\n";
      }
     elsif($retro[$x][3] == 1)
      {
      $arrow_start = $retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{size}{R};
      $arrow_end = $retro{$retro[$x][0]}{$retro[$x][1]}{Tend} - $retro{$retro[$x][0]}{$retro[$x][1]}{size}{R} + 400;
      print SVG "<polygon fill=\"black\" points=\"$arrow_start,$top $arrow_end,$arrow_bot $arrow_end,$arrow_top\"/>\n";
      } 
    }
  }
for(my $x=1;$x<$total_inserts+1;$x++)
  {
  if($retro[$x][0] eq 'solo')
    {
    #SOLO LTR ARROWS
    my $top = ($largest_level - $retro{$retro[$x][0]}{$retro[$x][1]}{level}) * 5000;
    my $arrow_top = $top + 500;
    my $arrow_bot = $top - 500;
    my $arrow_start;
    my $arrow_end;
    my $ltr_size = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro[$x][10];
    print SVG "<line stroke=\"black\" stroke-width=\"200\" x1=\"$retro{$retro[$x][0]}{$retro[$x][1]}{Tstart}\" y1=\"$top\" x2=\"$ltr_size\" y2=\"$top\"/>\n";
    if($retro[$x][3] == 0)
      {
      $arrow_start = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro[$x][10] - 400;
      $arrow_end = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + $retro[$x][10];
      print SVG "<polygon fill=\"black\" points=\"$arrow_start,$arrow_top $arrow_start,$arrow_bot $arrow_end,$top\"/>\n";
      }
     elsif($retro[$x][3] == 1)
      {
      $arrow_start = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart};
      $arrow_end = $retro{$retro[$x][0]}{$retro[$x][1]}{Tstart} + 400;
      print SVG "<polygon fill=\"black\" points=\"$arrow_start,$top $arrow_end,$arrow_bot $arrow_end,$arrow_top\"/>\n";
      }
    }
  }
print SVG "</svg>\n";
close(SVG);
####################################################################################################

sub COORD_PRE
{
my $hash_file = $_[0];
my $start_coord = $_[1];
my $end_coord = $_[2];
my $coord_split = $_[3];
my $new_hash = "$hash_file.new";
open(OUT, ">$new_hash");
open(HASH, $hash_file) or die "$hash_file not found\n";
my $line = <HASH>;
chomp($line);
my $total_seq = $line;
my $new_frag = $total_seq;
if(!$start_coord)
  {$start_coord = 1;}
if(!$end_coord)
  {$end_coord = $total_seq;}
$total_seq = $end_coord - ($start_coord-1);
print OUT "$total_seq\n";
$line = <HASH>;
chomp($line);
print OUT "$line\n";
while(<HASH>)
  {
  $line = $_;
  chomp($line);
  if($line =~ m/^SOLO/ || $line =~ m/^FRAG/ || $line =~ m/^NLTR/)
    {
    my $coords = <HASH>;
    chomp($coords);
    my @coords = split(/ /,$coords);
    my $coord_count = @coords;
    my $coord_count_s = ($coord_count-1)/4;
    my $new_coords = "$coords[0]";
    for(my $x=0;$x<$coord_count_s;$x++)
      {
      if($coords[($x*4)+1] > $start_coord && $coords[($x*4)+2] < $end_coord)
        {
        $coords[($x*4)+1] = $coords[($x*4)+1] - ($start_coord - 1);
        $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
        $new_coords = $new_coords . " $coords[($x*4)+1] $coords[($x*4)+2] $coords[($x*4)+3] $coords[($x*4)+4]";
        }
       elsif($coords[($x*4)+1] < $start_coord && $coords[($x*4)+2] > $start_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
          $new_coords = $new_coords . " 1 $coords[($x*4)+2] $coords[($x*4)+3] $coords[($x*4)+4]";
          }
        }
       elsif($coords[($x*4)+1] < $end_coord && $coords[($x*4)+2] > $end_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+1] = $coords[($x*4)+1] - ($start_coord - 1);
          $new_coords = $new_coords . " $coords[($x*4)+1] $end_coord $coords[($x*4)+3] $coords[($x*4)+4]";
          }
        }
      }
    my @new_coords = split(/ /,$new_coords);
    my $new_count = @new_coords;
    if($new_count > 3)
      {
      if($coord_split eq 'T')
        {
        print OUT "$line\n";
        print OUT "$new_coords\n";
        }
       elsif($coord_split eq 'F' && $coord_count == $new_count)
        {
        print OUT "$line\n";
        print OUT "$new_coords\n";
        }
      }
    }
  if($line =~ m/^PAIR/)
    {
    my $coords = <HASH>;
    chomp($coords);
    my @coords = split(/ /,$coords);
    my $coord_count_L = @coords;
    my $coord_count_s = ($coord_count_L-2)/4;
    my $new_coords_L = "$coords[0] $coords[1]";
    for(my $x=0;$x<$coord_count_s;$x++)
      {
      if($coords[($x*4)+2] > $start_coord && $coords[($x*4)+3] < $end_coord)
        {
        $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
        $coords[($x*4)+3] = $coords[($x*4)+3] - ($start_coord - 1);
        $new_coords_L = $new_coords_L . " $coords[($x*4)+2] $coords[($x*4)+3] $coords[($x*4)+4] $coords[($x*4)+5]";
        }
       elsif($coords[($x*4)+2] < $start_coord && $coords[($x*4)+3] > $start_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+3] = $coords[($x*4)+3] - ($start_coord - 1);
          $new_coords_L = $new_coords_L . " 1 $coords[($x*4)+3] $coords[($x*4)+4] $coords[($x*4)+5]";
          }
        }
       elsif($coords[($x*4)+2] < $end_coord && $coords[($x*4)+3] > $end_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
          $new_coords_L = $new_coords_L . " $coords[($x*4)+2] $end_coord $coords[($x*4)+4] $coords[($x*4)+5]";
          }
        }
      }
    $coords = <HASH>;
    chomp($coords);
    @coords = split(/ /,$coords);
    my $coord_count_R = @coords;
    $coord_count_s = ($coord_count_R-2)/4;
    my $new_coords_R = "$coords[0] $coords[1]";
    for(my $x=0;$x<$coord_count_s;$x++)
      {
      if($coords[($x*4)+2] > $start_coord && $coords[($x*4)+3] < $end_coord)
        {
        $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
        $coords[($x*4)+3] = $coords[($x*4)+3] - ($start_coord - 1);
        $new_coords_R = $new_coords_R . " $coords[($x*4)+2] $coords[($x*4)+3] $coords[($x*4)+4] $coords[($x*4)+5]";
        }
       elsif($coords[($x*4)+2] < $start_coord && $coords[($x*4)+3] > $start_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+3] = $coords[($x*4)+3] - ($start_coord - 1);
          $new_coords_R = $new_coords_R . " 1 $coords[($x*4)+3] $coords[($x*4)+4] $coords[($x*4)+5]";
          }
        }
       elsif($coords[($x*4)+2] < $end_coord && $coords[($x*4)+3] > $end_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
          $new_coords_R = $new_coords_R . " $coords[($x*4)+2] $end_coord $coords[($x*4)+4] $coords[($x*4)+5]";
          }
        }
      }
    $coords = <HASH>;
    chomp($coords);
    @coords = split(/ /,$coords);
    my $coord_count_M = @coords;
    $coord_count_s = ($coord_count_M-2)/4;
    my $new_coords_M = "$coords[0] $coords[1]";
    for(my $x=0;$x<$coord_count_s;$x++)
      {
      if($coords[($x*4)+2] > $start_coord && $coords[($x*4)+3] < $end_coord)
        {
        $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
        $coords[($x*4)+3] = $coords[($x*4)+3] - ($start_coord - 1);
        $new_coords_M = $new_coords_M . " $coords[($x*4)+2] $coords[($x*4)+3] $coords[($x*4)+4] $coords[($x*4)+5]";
        }
       elsif($coords[($x*4)+2] < $start_coord && $coords[($x*4)+3] > $start_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+3] = $coords[($x*4)+3] - ($start_coord - 1);
          $new_coords_M = $new_coords_M . " 1 $coords[($x*4)+3] $coords[($x*4)+4] $coords[($x*4)+5]";
          }
        }
       elsif($coords[($x*4)+2] < $end_coord && $coords[($x*4)+3] > $end_coord)
        {
        if($coord_split eq 'T')
          {
          $coords[($x*4)+2] = $coords[($x*4)+2] - ($start_coord - 1);
          $new_coords_M = $new_coords_M . " $coords[($x*4)+2] $end_coord $coords[($x*4)+4] $coords[($x*4)+5]";
          }
        }
      }
    my @new_coords_L = split(/ /,$new_coords_L);
    my $new_count_L = @new_coords_L;
    my @new_coords_R = split(/ /,$new_coords_R);
    my $new_count_R = @new_coords_R;
    my @new_coords_M = split(/ /,$new_coords_M);
    my $new_count_M = @new_coords_M;
    if($coord_split eq 'T')
      {
      if($new_count_L < 3 || $new_count_R < 3)
        {
        if($new_count_L > 3 || $new_count_M > 3 || $new_count_R > 3)
          {
          #MOVE TO FRAG
          my @new_frag = split(/ /,$line);
          print OUT "FRAG f$new_frag $new_frag[2] $new_frag[3] $new_frag[4] $new_frag[5] $new_frag[6] $new_frag[7]\n";
          print OUT "f$new_frag";
          for(my $x=0;$x<(@new_coords_L-2)/4;$x++)
            {print OUT " $new_coords_L[($x*4)+2] $new_coords_L[($x*4)+3] $new_coords_L[($x*4)+4] $new_coords_L[($x*4)+5]";}
          for(my $x=0;$x<(@new_coords_M-2)/4;$x++)
            {print OUT " $new_coords_M[($x*4)+2] $new_coords_M[($x*4)+3] $new_coords_M[($x*4)+4] $new_coords_M[($x*4)+5]";}
          for(my $x=0;$x<(@new_coords_R-2)/4;$x++)
            {print OUT " $new_coords_R[($x*4)+2] $new_coords_R[($x*4)+3] $new_coords_R[($x*4)+4] $new_coords_R[($x*4)+5]";}
          print OUT "\n";
          $new_frag++;
          }
        }
       else
        {
        print OUT "$line\n";
        print OUT "$new_coords_L\n";
        print OUT "$new_coords_R\n";
        print OUT "$new_coords_M\n";
        }
      }
     elsif($coord_split eq 'F' && $coord_count_L == $new_count_L && $coord_count_R == $new_count_R && $coord_count_M == $new_count_M)
      {
      print OUT "$line\n";
      print OUT "$new_coords_L\n";
      print OUT "$new_coords_R\n";
      print OUT "$new_coords_M\n";
      }
    }
  }
close(HASH);
close(OUT);
return($new_hash,$start_coord,$end_coord);
}

