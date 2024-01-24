#!/usr/bin/perl -w

###########################################################################
#                                LICENSE                                  #
###########################################################################
#
# The TE nest software package, including TE_nest.pl, svg_ltr.pl and the
# associated repeat databases were written and are maintained by
# Brent Kronmiller (bak@iastate.edu).  Please email me with any questions,
# suggestions, or problems.
#
# TE nest is distributed under the GNU general public license.
# Please see http://www.gnu.org/licenses/gpl.txt for more information.
#
###########################################################################
#                        PROGRAM DESCRIPTION                              #
###########################################################################
#
# VERSION TE_nest_1.0
#5/17/07 - Added GNU license and version number
#
#5/6/07  - Combined GDCB version and local version to keep it all the same
#
#4/27/07 - fixed custom DB directory locations
#        - updated with new FRAG/NLTR system
#
#2/7/07  - fixed discrepency and overlap process to pick best hit in entire region, not first
#
#12/9/06 - made command line options
#        - added ltr_list to .LTR so not needed for svg_ltr.pl
#        - added organism switch, changed directory structure (ltr_list found in org dir)
#
###########################################################################
#                             PARAMETERS                                  #

use strict;
use Getopt::Long;
#
my $processors = 2;         # --proc
my $db_directory = '/home/bak/TE_nest/'; # --current 
my $blast_directory = '/home/bak/bin/';                       # --blast
my $lal_directory = '/home/bak/bin/FASTA/';                         # --lal
my $output_directory = '/home/bak/TE_nest/';
my $input_ltr = 'ltr_list'; # --te_list
my $organism = 'maize';                                       # --org
my $pair_only = 'F';   # use T to only return PAIR information, S to return PAIR and SOLO
my $frag_only = 'F';   # use T to only return FRAG and NLTR information
our $running;
#
### LTR PARAMETERS ###
#
# LTR alignment lalign parameters
 # gap open
our $ltr_lal_f = 30;  # --ltr_gap_open
 # gap extension
our $ltr_lal_g = 15;  # --ltr_gap_ext
 # alignments to report
our $ltr_lal_k = 7;   # --ltr_gap_rep
#
# LALIGN parse score
our $ltr_lal_score = 20; # --ltr_score
#
# Amount of overlap ignored when combining overlapping LTR alignments.  Value
# of 0 does not allow any overlap, value of 5 allows slight overlap between
# sections.
our $ltr_unique_overlap = 25; # --ltr_overlap
#
# Amount of bases to allow for overlapping LTR alignments. Value of 0 requires
# exact LTR alignment coordinates (no overlapping), however this may result
# in some LTR hits not reported (especially fragmented or very small LTRs). 
# Value of 5 allows some wiggle room, coords will be corrected later.
# Parameter is tied to $ltr_pwr_offset, should be similar value.
our $ltr_ucoord_fudge = 30; # --ltr_ucoord
#
# Amount of bases for offset of LTR powerset joining.  Value of 0 allows no overlap
# between combined sections.  Value of 10 allows slight overlap.
# Parameter is tied to $ltr_ucoord_fudge, should be similar value.
our $ltr_pwr_offset = 30; # --ltr_pwr_offset
#
# Smallest size of LTR section allowed after powerset join and trimming.
our $ltr_small_size = 25; # --ltr_small
#
# LTR powerset gap distance division value
our $ltr_pwr_gap = 100; # --ltr_pwr_div
#
# Max length between LTR sections reconstituted with LTR PowerSet.
our $ltr_pwr_max = 100000; # --ltr_pwr_max
#
#BSR alignment lalign parameters
#
 #gap open
our $bsr_lal_f = 12; # --bsr_gap_open
 #gap extension
our $bsr_lal_g = 4; # --bsr_gap_ext
#
# Amount to drop BSR compare each round (starts at 0.6, value of 0.1 will comapre
# 0.6, 0.5, 0.4, 0.3, 0.2, 0.1.
our $bsr_compare_drop = 0.1; # --bsr_drop
#
### MID PARAMETERS ###
#
# MID alignment lalign parameters
 # gap open
our $mid_lal_f = 75;  # --mid_gap_open
 # gap extension
our $mid_lal_g = 75;  # --mid_gap_ext
 # amount of alignments to report
our $mid_lal_k = 7;   # --mid_gap_rep
#
# LALIGN parse score
our $mid_lal_score = 20; #--mid_score
#
### FRAG/NLTR
#
# F/N alignment lalign parameters
 # gap open
our $fn_lal_f = 75;  # --fn_gap_open
 # gap extension
our $fn_lal_g = 75;  # --fn_gap_ext
 # amount of alignments to report
our $fn_lal_k = 7;   # --fn_gap_rep
#
# LALIGN parse score
our $fn_lal_score = 20; #--fn_score
#
GetOptions('ltr_gap_open=i' => \$ltr_lal_f, 'ltr_gap_ext=i' => \$ltr_lal_g, 'ltr_gap_rep=i' => \$ltr_lal_k,
           'ltr_score=i' => \$ltr_lal_score,
           'mid_gap_open=i' => \$mid_lal_f, 'mid_gap_ext=i' => \$mid_lal_g, 'mid_gap_rep=i' => \$mid_lal_k,
           'mid_score=i' => \$mid_lal_score,
           'fn_gap_open=i' => \$fn_lal_f, 'fn_gap_ext=i' => \$fn_lal_g, 'fn_gap_rep=i' => \$fn_lal_k,
           'fn_score=i' => \$fn_lal_score,
           'bsr_gap_open=i' => \$bsr_lal_f, 'bsr_gap_ext=i' => \$bsr_lal_g,
           'ltr_overlap=i' => \$ltr_unique_overlap,
           '--ltr_ucoord=i' => \$ltr_ucoord_fudge,
           '--ltr_pwr_offset=i' => \$ltr_pwr_offset,
           '--ltr_small=i' => \$ltr_small_size,
           '--ltr_pwr_div=i' => \$ltr_pwr_gap,
           '--ltr_pwr_max=i' => \$ltr_pwr_max,
           '--bsr_drop=i' => \$bsr_compare_drop,
           '--proc=i' => \$processors, 
           '--te_list=s' => \$input_ltr,
           '--pair_only=s' => \$pair_only,
           '--frag_only=s' => \$frag_only,
           '--org=s' => \$organism,
           '--output=s' => \$output_directory);
my $org_directory = $organism;
if($organism !~ m/^\//)
  {
  $organism =~ tr/[a-z]/[A-Z]/;
  $org_directory = $db_directory . $organism;
  }
#
#
#                                                                         #
###########################################################################
#                                MAIN                                     #
###########################################################################

#GET DATE FOR DIRECTORY AND FILE NAMES
my $date = DATE();
#GET FASTA FILE TO ANALIZE
if ($#ARGV == -1)
  {
  print "gimme some commands\n";
  exit;
  }
my $input_fasta = $ARGV[0];
my @directory = split(/\//,$input_fasta);
my $directory;
if(@directory > 1)
  {
  for(my $x=0;$x<@directory-1;$x++)
    {$directory = $directory . $directory[$x] . "/";}
  }
 else
  {$directory = './';}
$input_fasta = $directory[@directory-1];
if(!(-e "$output_directory/$input_fasta"))
  {system("cp $directory/$input_fasta $output_directory");}
system("$blast_directory/pressdb $output_directory/$input_fasta");
$date = 'TE_DIR_'. $date . '_' . $input_fasta;   #FOR LOCAL VERSION
#$date = 'TE_DIR_' . $input_fasta; #FOR WEB-CGI VERSION
print "Date is $date process will occur in this directory, files will be named with date\n";
system("mkdir $output_directory/$date");
my $bug_file = "$output_directory/$date/debug.out";
open(BUG, ">$bug_file");
print BUG "--org $organism\n";
#OPEN ltr_list AND WRITE TO @LTR_LIST
if($input_ltr eq 'ltr_list')
  {open(LTR_LIST, "$org_directory/$input_ltr") or die "$input_ltr not found\n";}
 else
  {open(LTR_LIST, "$input_ltr") or die "$input_ltr not found\n";}
our @ltr_list = ();
our @te_list = ();
while(<LTR_LIST>)
  {
  my $ltr_list = $_;
  chomp($ltr_list);
  $ltr_list =~ s/\s+/ /g;
  my @list_split = split(/ /, $ltr_list);
  push(@te_list, $list_split[0]);
  if($list_split[1] eq '0')
    {
    push(@ltr_list, $list_split[0]);
    }
  }
close(LTR_LIST);
#OPEN FASTA LOAD SEQUENCE INTO SEQ
my $seq_ref = FASTA($input_fasta);
our @seq = @$seq_ref;
our $seq_size = @seq;
#FOR EACH LTR IN LIST, MULTI PROCESS BLAST, ALIGN EACH BLAST HIT, WRITE RESULTS TO FILE
my $total = @ltr_list;
my $total_te = @te_list;
print BUG "$total LTRs SENT TO LTR ALIGNMENT ANALYSIS\n";
my @ltr = ();
my $lal_file = "$output_directory/$date/lal_total";
open(OUT_LAL, ">$lal_file");
my $ltr_size = "$output_directory/$date/ltr_size";
open(SIZE, ">$ltr_size");
if($frag_only eq 'F')
  {MULTI($total, '0', '0');}
close(OUT_LAL);
close(SIZE);
#OPEN LTR ALIGNMENTS FILE, WRITE TO ARRAY @lal_total
open(LAL_TOTAL, $lal_file) or die "$lal_file not found\n";
my @lal_total = ();
while(<LAL_TOTAL>)
  {
  my $line = $_;
  chomp($line);
  push @lal_total, [split(/ /, $line)];
  }
system("rm $output_directory/$date/lal_total");
#OPEN AND GET LTR SIZES
open(SIZE, $ltr_size) or die "$ltr_size not found\n";
my @ltr_size = ();
while(<SIZE>)
  {
  my $line = $_;
  chomp($line);
  push @ltr_size, [split(/ /, $line)];
  }
close(SIZE);
#system("rm $output_directory/$date/ltr_size");
#UNIQUE @lal_total
my $ltr_u_ref = UNIQUE(@lal_total);
my @ltr_u = @$ltr_u_ref;
my $ltr_u_count = @ltr_u;
#PRINT UNIQUE LTRs
print BUG "print UNIQUE ltrs\n";
for(my $x=0;$x<$ltr_u_count;$x++)
  {
  print BUG "$ltr_u[$x][0] $ltr_u[$x][1] $ltr_u[$x][2] $ltr_u[$x][3] $ltr_u[$x][4] $ltr_u[$x][5] $ltr_u[$x][6]\n";
  }
#PUT EACH FULL LTR IN POWERSET FILE, SEND PARTIALS THROUGH POWERSET
my $pwr_file = "$output_directory/$date/pwr_total";
open(OUT_PWR, ">$pwr_file");
our @ltr_unique = ();
for(my $x=0;$x<$ltr_u_count;$x++)
  {
  if(($ltr_u[$x][3] - $ltr_u[$x][2]) / $ltr_size[$ltr_u[$x][0]][0] > .99)
    {
    print OUT_PWR "$ltr_u[$x][0] $ltr_u[$x][1] $ltr_u[$x][2] $ltr_u[$x][3] $ltr_u[$x][4] $ltr_u[$x][5] $ltr_u[$x][6]\n";
    }
   else
    {
    my $ltr_unique = "$ltr_u[$x][0] $ltr_u[$x][1] $ltr_u[$x][2] $ltr_u[$x][3] $ltr_u[$x][4] $ltr_u[$x][5] $ltr_u[$x][6]";
    push @ltr_unique, [split(/ /, $ltr_unique)];
    } 
  }   
my $ltr_unique_count = @ltr_unique;
print BUG "NOW SENDING TO POWERSET $ltr_unique_count from an original $ltr_u_count\n";
#FOR EACH partial LTR SEND @ltr_unique TO POWERSET, JOIN PARTIAL LTRs
MULTI($total, '0', '1');
close(OUT_PWR);
#OPEN POWERSET FILE, WRITE TO ARRAY @pwr_total
open(PWR_TOTAL, $pwr_file) or die "$pwr_file not found\n";
our @pwr_total = ();
my $x = 0;
while(<PWR_TOTAL>)
  {
  my $line = $_;
  chomp($line);
  $line = "s$x " . $line;
  push @pwr_total, [split(/ /, $line)];
  $x++;
  }
system("rm $output_directory/$date/pwr_total");
my $pwr_total_count = @pwr_total;
#TRIM @pwr_total SECTIONS FOR ANY SLIGHT OVERLAP
for(my $x=0;$x<$pwr_total_count;$x++)
  {
  my $power_lenX = scalar @{$pwr_total[$x]};
  for(my $y=0;$y<($power_lenX-1)/7;$y++)
    {
    for(my $i=0;$i<$pwr_total_count;$i++)
      {
      my $power_lenI = scalar @{$pwr_total[$i]};
      for(my $j=0;$j<($power_lenI-1)/7;$j++)
        {
        my $x_check = 0;
        my $i_check = 0;
        if($x == $i)
          {$x_check = 1;}
        if($y != $j)
          {$i_check = 1;}
        if($i_check == 0 || $x_check == 0)
          {
          if($pwr_total[$x][$y*7+5] < $pwr_total[$i][$j*7+5])
            {
            if($pwr_total[$x][$y*7+6] < $pwr_total[$i][$j*7+6])
              {
              if($pwr_total[$x][$y*7+6] >= $pwr_total[$i][$j*7+5])
                {
                if($pwr_total[$x][$y*7+7] >= $pwr_total[$i][$j*7+7])
                  {
                  $pwr_total[$i][$j*7+5] = $pwr_total[$x][$y*7+6] + 1;
                  }
                 else
                  {
                  $pwr_total[$x][$y*7+6] = $pwr_total[$i][$j*7+5] - 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
my @ltr_temp = ();
my $z=0;
for(my $x=0;$x<$pwr_total_count;$x++)
  {
  my $power_lenX = scalar @{$pwr_total[$x]};
  for(my $y=0;$y<($power_lenX-1)/7;$y++)
    {
    if(($pwr_total[$x][($y*7)+6]-$pwr_total[$x][($y*7)+5]) < $ltr_small_size)
      {
      splice(@{$pwr_total[$x]},($y*7)+1,7);
      $power_lenX = scalar @{$pwr_total[$x]};
      $y--;
      }
    }
  $power_lenX = scalar @{$pwr_total[$x]};
  if($power_lenX > 2)
    {
    my $s = "s".$z;
    $z++;
    my $ltr_temp = "$s";
    for(my $y=0;$y<($power_lenX-1)/7;$y++)
      {
      $ltr_temp = $ltr_temp . " " . "$pwr_total[$x][($y*7)+1] $pwr_total[$x][($y*7)+2] $pwr_total[$x][($y*7)+3] $pwr_total[$x][($y*7)+4] $pwr_total[$x][($y*7)+5] $pwr_total[$x][($y*7)+6] $pwr_total[$x][($y*7)+7]"; 
      }
    push @ltr_temp, [split(/ /,$ltr_temp)];
    }
  }
$pwr_total_count = @ltr_temp;
@pwr_total = sort{$a->[5] <=> $b->[5]} @ltr_temp;
@pwr_total = sort{$a->[1] <=> $b->[1]} @pwr_total;
#WRITE POWERSET LTRs TO %retro
my %retro = ();
$pwr_total_count = @pwr_total;
for(my $x=0;$x<$pwr_total_count;$x++)
  {
  $retro{solo}{$pwr_total[$x][0]} = {type => $pwr_total[$x][1], dir => $pwr_total[$x][2], bsr => 100, order => 100, level => 100, group => 100};
  my $power_len = scalar @{$pwr_total[$x]};
  $power_len = ($power_len-1)/7;
  for(my $y=0;$y<$power_len;$y++)
    {
    $retro{solo}{$pwr_total[$x][0]}{coords}{$y} = {SEQ_start => $pwr_total[$x][($y*7)+5],
                                    SEQ_end => $pwr_total[$x][($y*7)+6],
                                    TE_start => $pwr_total[$x][($y*7)+3],
                                    TE_end => $pwr_total[$x][($y*7)+4]};
    }
  }
#PRINT POWERSET UNIQUE
print BUG "print POWERSET ltrs\n";
for(my $x=0;$x<$pwr_total_count;$x++)
  {
  print BUG "$pwr_total[$x][0]";
  my $power_lenX = scalar @{$pwr_total[$x]};
  for(my $y=1;$y<$power_lenX;$y++)
    {
    print BUG " $pwr_total[$x][$y]";
    }
  print BUG "\n";
  }

#FOR EACH LTR ALIGN CUT LTR TO CUT LTR, DETERMINE BSR AND PAIR
my $bsr_file = "$output_directory/$date/bsr_total";
open(OUT_BSR, ">$bsr_file");
MULTI($total, '0', '2');
close(OUT_BSR);
#OPEN BSR FILE, WRITE TO ARRAY @bsr_total
open(BSR_TOTAL, $bsr_file) or die "$bsr_file not found\n";
our @bsr_total = ();
$x = 0;
while(<BSR_TOTAL>)
  {
  my $line = $_;
  chomp($line);
  #ADD UNIQUE IDENTIFIER TO EACH LTR PAIR
  $line = "p$x " . $line;
  push @bsr_total, [split(/ /, $line)];
  $x++;
  }
system("rm $output_directory/$date/bsr_total");
@bsr_total = sort{$a->[4] <=> $b->[4]} @bsr_total;
my $bsr_total_count = @bsr_total;
#DETERMINE NESTING GROUP, LEVEL and ORDER
$bsr_total[0][11] = $bsr_total[0][12] = 0;
for(my $x=1;$x<$bsr_total_count;$x++)
  {
  my $do = 0;
  my $y = $x-1;
  while($bsr_total[$y][0] && $do==0 && $y>=0)
    {
    if($bsr_total[$x][9] < $bsr_total[$y][9])
      {
      $bsr_total[$x][11] = $bsr_total[$y][11];
      $bsr_total[$x][12] = $bsr_total[$y][12] + 1;
      $do = 1;
      }
     else
      {$y--;}
    }
  if($do==0)
    {
    $bsr_total[$x][11] = $bsr_total[$x-1][11] + 1;
    $bsr_total[$x][12] = 0;
    }
  }
my $nest_group = $bsr_total[$bsr_total_count-1][11];
my $zero = $bsr_total[0][13] = 0;
for(my $x=1;$x<$bsr_total_count;$x++)
  {
  my $y = $x-1;
  if($bsr_total[$x][12] == 0)
    {
    $zero++;
    $bsr_total[$x][13] = $zero;
    }
   elsif($bsr_total[$y][12] < $bsr_total[$x][12])
    {$bsr_total[$x][13] = 0;}
   else
    {
    my $do = 0;
    while($y>0 && $bsr_total[$y][11] == $bsr_total[$x][11] && $do==0)
      {
      if($bsr_total[$y][12] == $bsr_total[$x][12])
        {
        $bsr_total[$x][13] = $bsr_total[$y][13] + 1;
        $do=1;
        }
      $y--;
      }
    }
  }
#WRITE PAIRS TO HASH, REMOVE SOLO LTRs THAT ARE IN PAIR SETS
my $solo_amt = scalar keys %{$retro{solo}};
for(my $x=0;$x<$bsr_total_count;$x++)
  {
  $retro{pair}{$bsr_total[$x][0]} = {type => $bsr_total[$x][1], dir => $retro{solo}{$bsr_total[$x][3]}{dir},
                                     bsr => $bsr_total[$x][10], group => $bsr_total[$x][11],
                                     level => $bsr_total[$x][12], order => $bsr_total[$x][13]};
  my $ltr_split = scalar keys %{$retro{solo}{$bsr_total[$x][3]}{coords}};
  for(my $y=0;$y<$ltr_split;$y++)
    {
    $retro{pair}{$bsr_total[$x][0]}{coords}{L}{$y} = 
                                 {SEQ_start => $retro{solo}{$bsr_total[$x][3]}{coords}{$y}{SEQ_start},
                                  SEQ_end => $retro{solo}{$bsr_total[$x][3]}{coords}{$y}{SEQ_end},
                                  TE_start => $retro{solo}{$bsr_total[$x][3]}{coords}{$y}{TE_start},
                                  TE_end => $retro{solo}{$bsr_total[$x][3]}{coords}{$y}{TE_end}};
    }
  $ltr_split = scalar keys %{$retro{solo}{$bsr_total[$x][7]}{coords}};
  for(my $y=0;$y<$ltr_split;$y++)
    {
    $retro{pair}{$bsr_total[$x][0]}{coords}{R}{$y} = 
                                 {SEQ_start => $retro{solo}{$bsr_total[$x][7]}{coords}{$y}{SEQ_start},
                                  SEQ_end => $retro{solo}{$bsr_total[$x][7]}{coords}{$y}{SEQ_end},
                                  TE_start => $retro{solo}{$bsr_total[$x][7]}{coords}{$y}{TE_start},
                                  TE_end => $retro{solo}{$bsr_total[$x][7]}{coords}{$y}{TE_end}};
    }
  }
$pwr_total_count--;
my @solo_tmp = ();
for(my $x=0;$x<scalar keys %{$retro{solo}};$x++)
  {
  my $bad = 0;
  my $p = "s".$x;
  for(my $z=0;$z<$bsr_total_count;$z++)
    {
    if($p eq $bsr_total[$z][3] || $p eq $bsr_total[$z][7])
      {$bad = 1;}
    }
  if($bad == 0)
    {
    my $solo_tmp = "$retro{solo}{$p}{type} $retro{solo}{$p}{dir} $retro{solo}{$p}{bsr} $retro{solo}{$p}{group} $retro{solo}{$p}{level} $retro{solo}{$p}{order}";
    for(my $y=0;$y<scalar keys %{$retro{solo}{$p}{coords}};$y++)
      {
      $solo_tmp = $solo_tmp . " $retro{solo}{$p}{coords}{$y}{SEQ_start} $retro{solo}{$p}{coords}{$y}{SEQ_end} $retro{solo}{$p}{coords}{$y}{TE_start} $retro{solo}{$p}{coords}{$y}{TE_end}";
      }
    push @solo_tmp, [split(/ /,$solo_tmp)];
    }
  }
delete($retro{solo});
my $count = @solo_tmp;
for(my $x=0;$x<$count;$x++)
  {
  my $p = "s".$x;
  $retro{solo}{$p} = {type => $solo_tmp[$x][0], dir => $solo_tmp[$x][1], bsr => $solo_tmp[$x][2],
                      group => $solo_tmp[$x][3], level => $solo_tmp[$x][4], order => $solo_tmp[$x][5]};
  my $power_len = ((scalar @{$solo_tmp[$x]})-6)/4;
  for(my $y=0;$y<$power_len;$y++)
    {
    $retro{solo}{$p}{coords}{$y} = {SEQ_start => $solo_tmp[$x][($y*4)+6],
                                    SEQ_end => $solo_tmp[$x][($y*4)+7],
                                    TE_start => $solo_tmp[$x][($y*4)+8],
                                    TE_end => $solo_tmp[$x][($y*4)+9]};
    }
  }
print BUG "BSR PAIRS\n";
for(my $x=0;$x<$bsr_total_count;$x++)
  {
  print BUG "$bsr_total[$x][0] $bsr_total[$x][1] $bsr_total[$x][2] $bsr_total[$x][3] $bsr_total[$x][4] $bsr_total[$x][5] $bsr_total[$x][6] $bsr_total[$x][7] $bsr_total[$x][8] $bsr_total[$x][9] $bsr_total[$x][10] $bsr_total[$x][11] $bsr_total[$x][12] $bsr_total[$x][13]\n";
  }

#FOR EACH NEST GROUP, FIND MID LOCATIONS
my $mid_file = "$output_directory/$date/mid_total";
open(OUT_MID, ">$mid_file");
MULTI($nest_group+1, '0', '3');
close(OUT_MID);
#OPEN MID FILE, WRITE TO ARRAY @mid_total
open(MID_TOTAL, $mid_file) or die "$mid_file not found\n";
our @mid_total = ();
while(<MID_TOTAL>)
  {
  my $line = $_;
  chomp($line);
  push @mid_total, [split(/ /, $line)];
  }
close(MID_TOTAL);
system("rm $output_directory/$date/mid_total");
my $mid_total_count = @mid_total;
#WRITE POWERSET MID HITS TO HASH
for(my $x=0;$x<$mid_total_count;$x++)
  {
  my $mid_split = scalar @{$mid_total[$x]};
  $mid_split = ($mid_split/6);
  for(my $y=0;$y<$mid_split;$y++)
    {
    $retro{pair}{$mid_total[$x][0]}{coords}{M}{$y} = {SEQ_start => $mid_total[$x][($y*6)+3],
                                                      SEQ_end => $mid_total[$x][($y*6)+4], 
                                                      TE_start => $mid_total[$x][($y*6)+1],
                                                      TE_end => $mid_total[$x][($y*6)+2]};
    }
  }
#PRINT POWERSET UNIQUE
print BUG "print POWERSET mids\n";
for(my $x=0;$x<$mid_total_count;$x++)
  {
  print BUG "$mid_total[$x][0]";
  my $power_len = scalar @{$mid_total[$x]};
  for(my $y=1;$y<$power_len;$y++)
    {
    print BUG " $mid_total[$x][$y]";
    }
  print BUG "\n";
  }

my @first_dis = ();
my @second_dis = ();
my @first_head = ();
my @second_head = ();
our @blank = ();
if($pair_only eq 'F')
{
#MAKE SURE THERE ARE NO COORDINATE NESTING DISCREPANCIES BETWEEN PAIRS
DISCREPANCY('PAIR','PAIR');
#MAKE SURE THERE ARE NO COORDINATE NESTING DISCREPANCIES BETWEEN SOLOS
DISCREPANCY('SOLO','SOLO');
#MAKE SURE THERE ARE NO COORDINATE NESTING DISCREPANCIES BETWEEN PAIRS AND SOLOS
DISCREPANCY('SOLO','PAIR');
#MAKE SURE THERE ARE NO SOLOs IN MID LOCATIONS, IF SO REMOVE THE SOLO
OVERLAP_S_M('SOLO','PAIR');
#MAKE SURE THERE ARE NO PAIRS IN MID LOCATIONS, IF SO REMOVE THE PAIR
OVERLAP_P_M('PAIR','PAIR');

#GET FRAG LOCATIONS
#READ HASH GET PAIR AND SOLO COORDS
my $pair_amt = scalar keys %{$retro{pair}};
my @blank_all= ();
print BUG "PAIR AMT $pair_amt\n";
for(my $x=0;$x<$pair_amt;$x++)
  {
  for(my $y=0;$y<3;$y++)
    {
    my $section = '';
    if($y==0)
      {$section = 'L';}
     elsif($y==1)
      {$section = 'R';}
     elsif($y==2)
      {$section = 'M';}
    my $p = 'p'.$x;
    my $coord_amt = scalar keys %{$retro{pair}{$p}{coords}{$section}};
    for(my $z=0;$z<$coord_amt;$z++)
      {
      my $blank_all = '';
      if($retro{pair}{$p}{coords}{$section}{$z}{SEQ_start} < $retro{pair}{$p}{coords}{$section}{$z}{SEQ_end})
        {
        $blank_all = "$retro{pair}{$p}{coords}{$section}{$z}{SEQ_start} $retro{pair}{$p}{coords}{$section}{$z}{SEQ_end}";
        }
       else
        {
        $blank_all = "$retro{pair}{$p}{coords}{$section}{$z}{SEQ_end} $retro{pair}{$p}{coords}{$section}{$z}{SEQ_start}";
        }
#print "blank_all $blank_all-$retro{pair}{$p}{coords}{$section}{$z}{SEQ_start} < $retro{pair}{$p}{coords}{$section}{$z}{SEQ_end}\n";
      push @blank_all, [split(/ /, $blank_all)];
      }
    }
  }
#UNIQUE SEQUENCE LOCATIONS
my $blank_all_count = @blank_all;
my @blank_unique = (); 
if($blank_all_count > 0)
  {
  @blank_all = sort {$a->[0] <=> $b->[0]} @blank_all;
  my $blank_unique = "$blank_all[0][0] $blank_all[0][1]";
  push @blank_unique, [split(/ /, $blank_unique)];
  my $unique_raise=1;
  for(my $x=1;$x<$blank_all_count;$x++)
    {
    if(($blank_all[$x-1][1]+20) >= $blank_all[$x][0])
      {
      $blank_unique[$unique_raise-1][1] = $blank_all[$x][1];
      }
     else
      {
      $blank_unique = "$blank_all[$x][0] $blank_all[$x][1]";
      push @blank_unique, [split(/ /, $blank_unique)];
      $unique_raise++;
      }
    }
  }
#GET BLANK LOCATIONS
my $blank_uni_count = @blank_unique;
my $blank = '';
my @blank_u = ();
if($blank_uni_count > 0)
  {
  if ($blank_unique[0][0] > 20)
    {
    $blank = "1 $blank_unique[0][0]";
    push @blank_u, [split(/ /, $blank)];
    }
  for($x=1;$x<$blank_uni_count;$x++)
    {
    $blank = "$blank_unique[$x-1][1] $blank_unique[$x][0]";
    push @blank_u, [split(/ /, $blank)];
    }
  if($blank_unique[$blank_uni_count-1][1]+20 < $seq_size)
    {
    $blank = "$blank_unique[$blank_uni_count-1][1] $seq_size";
    push @blank_u, [split(/ /, $blank)];
    }
  }
 else
  {
  $blank_u[0][0] = 1;
  $blank_u[0][1] = $seq_size;
  }
my $blank_count = @blank_u;
#SPLIT BLANK LOCATIONS
for(my $x=0;$x<$blank_count;$x++)
  {
  my $blank_size = $blank_u[$x][1] - $blank_u[$x][0];
  if($blank_size <= 10000)
    {
    $blank = "$blank_u[$x][0] $blank_u[$x][1]";
    push @blank, [split(/ /,$blank)];
    }
   else
    {
    my $divide_amt = int($blank_size / 10000) + 1;
    my $add_amt = int(($blank_size / $divide_amt) + .5);
    my $start = $blank_u[$x][0];
    my $end;
    for(my $y=1;$y<$divide_amt;$y++)
      {
      $end = $start + $add_amt; 
      $blank = "$start $end";
      push @blank, [split(/ /,$blank)];
      $start = $end
      }
    $blank = "$end $blank_u[$x][1]";
    push @blank, [split(/ /,$blank)];
    }
  }
$blank_count = @blank;
for(my $x=0;$x<$blank_count;$x++)
  {
  print BUG "BLANK $x $blank[$x][0] $blank[$x][1]\n";
  }
my $frag_file = "$output_directory/$date/frag_total";
open(OUT_FRAG, ">$frag_file");
print BUG "total frag locs is $blank_count\n";
MULTI($blank_count, '0', '4');
close(OUT_FRAG);
system("rm $output_directory/$date/rev*");
#OPEN FRAG FILE, WRITE TO ARRAY @frag_total
open(FRAG_TOTAL, $frag_file) or die "$frag_file not found\n";
our @frag_set = ();
my $frag_set_count = 0;
while(<FRAG_TOTAL>)
  {
  my $line = $_;
  chomp($line);
  push @frag_set, [split(/ /, $line)];
  for(my $x=0;$x<$total_te;$x++)
    {
    $frag_set[$frag_set_count][5] =~ s/\.fasta//;
    $frag_set[$frag_set_count][5] =~ s/\.con//;
    $frag_set[$frag_set_count][5] =~ s/_con//;
    $frag_set[$frag_set_count][5] =~ s/_R//;
    if($frag_set[$frag_set_count][5] eq $te_list[$x])
      {
      $frag_set[$frag_set_count][5] = $x;
      }
    }
  $frag_set_count++;
  }
close(FRAG_TOTAL);
system("rm $output_directory/$date/frag_total");
#FIND INSERTION GROUPS FOR EACH FRAG
if($frag_set_count > 1)
  {
  @frag_set = sort {$a->[0] <=> $b->[0]} @frag_set;
  }
for(my $w=0;$w<$frag_set_count;$w++)
  {
print BUG "$w FRAG_SET $frag_set[$w][0] $frag_set[$w][1] $frag_set[$w][2] $frag_set[$w][3] $frag_set[$w][4] $frag_set[$w][5] $frag_set[$w][6] $frag_set[$w][7]\n";
  $frag_set[$w][1] = $frag_set[$w][1]-1;
  my $coordA;
  my $coordB;
  my $coordLoc;
  my $small_size = 100000000;
  my $made_coord = 0;
  my $coord_amt;
  for(my $x=0;$x<$pair_amt;$x++)
    {
    $coord_amt = scalar keys %{$retro{pair}{"p$x"}{coords}{L}};
    for(my $y=0;$y<$coord_amt;$y++)
      {
      if($frag_set[$w][0]+5 >= $retro{pair}{"p$x"}{coords}{L}{$y}{SEQ_end})
        {
        if($retro{pair}{"p$x"}{coords}{L}{$y+1})
          {
          if($frag_set[$w][1]-5 <= $retro{pair}{"p$x"}{coords}{L}{$y+1}{SEQ_start})
            {
            $coordA = $retro{pair}{"p$x"}{coords}{L}{$y}{SEQ_end};
            $coordB = $retro{pair}{"p$x"}{coords}{L}{$y+1}{SEQ_start};
            $coordLoc = "p$x";
            $made_coord = 1;
            }
          }
         elsif($retro{pair}{"p$x"}{coords}{M}{0})
          {
          if($frag_set[$w][1]-5 <= $retro{pair}{"p$x"}{coords}{M}{0}{SEQ_start})
            {
            $coordA = $retro{pair}{"p$x"}{coords}{L}{$y}{SEQ_end};
            $coordB = $retro{pair}{"p$x"}{coords}{M}{0}{SEQ_start};
            $coordLoc = "p$x";
            $made_coord = 1;
            }
          }
         else
          {
          if($frag_set[$w][1]-5 <= $retro{pair}{"p$x"}{coords}{R}{0}{SEQ_start})
            {
            $coordA = $retro{pair}{"p$x"}{coords}{L}{$y}{SEQ_end};
            $coordB = $retro{pair}{"p$x"}{coords}{R}{0}{SEQ_start};
            $coordLoc = "p$x";
            $made_coord = 1;
            }
          }
        }
      }
    $coord_amt = scalar keys %{$retro{pair}{"p$x"}{coords}{M}};
    for(my $y=0;$y<$coord_amt;$y++)
      {
      if($frag_set[$w][0]+5 >= $retro{pair}{"p$x"}{coords}{M}{$y}{SEQ_end})
        {
        if($retro{pair}{"p$x"}{coords}{M}{$y+1})
          {
          if($frag_set[$w][1]-5 <= $retro{pair}{"p$x"}{coords}{M}{$y+1}{SEQ_start})
            {
            $coordA = $retro{pair}{"p$x"}{coords}{M}{$y}{SEQ_end};
            $coordB = $retro{pair}{"p$x"}{coords}{M}{$y+1}{SEQ_start};
            $coordLoc = "p$x";
            $made_coord = 1;
            }
          }
         else
          {
          if($frag_set[$w][1]-5 <= $retro{pair}{"p$x"}{coords}{R}{0}{SEQ_start})
            {
            $coordA = $retro{pair}{"p$x"}{coords}{M}{$y}{SEQ_end};
            $coordB = $retro{pair}{"p$x"}{coords}{R}{0}{SEQ_start};
            $coordLoc = "p$x";
            $made_coord = 1;
            }
          }
        }
      }
    $coord_amt = scalar keys %{$retro{pair}{"p$x"}{coords}{R}};
    for(my $y=0;$y<$coord_amt;$y++)
      {
      if($frag_set[$w][0]+5 >= $retro{pair}{"p$x"}{coords}{R}{$y}{SEQ_end})
        {
        if($retro{pair}{"p$x"}{coords}{R}{$y+1})
          {
          if($frag_set[$w][1]-5 <= $retro{pair}{"p$x"}{coords}{R}{$y+1}{SEQ_start})
            {
            $coordA = $retro{pair}{"p$x"}{coords}{R}{$y}{SEQ_end};
            $coordB = $retro{pair}{"p$x"}{coords}{R}{$y+1}{SEQ_start};
            $coordLoc = "p$x";
            $made_coord = 1;
            }
          }
        }
      }
    if($made_coord == 1 && $coordB - $coordA < $small_size)
      {
      $small_size = $coordB - $coordA;
      $frag_set[$w][7] = $coordA;
      $frag_set[$w][8] = $coordB;
      $frag_set[$w][9] = $coordLoc;
      }
    }
  if($made_coord == 0)
    {
    $frag_set[$w][7] = 0;
    $frag_set[$w][8] = $seq_size;
    $frag_set[$w][9] = 'dna';
    }
  }
if($frag_set_count > 1)
  {
  @frag_set = sort {$a->[5] cmp $b->[5]} @frag_set;
  @frag_set = sort {$a->[8] <=> $b->[8]} @frag_set;
  @frag_set = sort {$a->[9] cmp $b->[9]} @frag_set;
  @frag_set = sort {$a->[7] <=> $b->[7]} @frag_set;
  }
my @frag_group = ();
$x = 0;
my $frag = my $nltr = 0;
my $push_last = 0;
for(my $x=0;$x<$frag_set_count;$x++)
  {
  my $frag_group = "$frag_set[$x][0] $frag_set[$x][1] $frag_set[$x][2] $frag_set[$x][3] $frag_set[$x][4] $frag_set[$x][5] $frag_set[$x][6] $frag_set[$x][7] $frag_set[$x][8] $frag_set[$x][9]";
  }
my @frag_left = ();
for(my $x=0;$x<$frag_set_count;$x++)
  {
  my $frag_left = "$frag_set[$x][0] $frag_set[$x][1] $frag_set[$x][2] $frag_set[$x][3] $frag_set[$x][4] $frag_set[$x][5] $frag_set[$x][6] $frag_set[$x][7] $frag_set[$x][8] $frag_set[$x][9]";
  push @frag_left, [split(/ /, $frag_left)];
  }
my $frag_left_count = @frag_left;
my @frag_hold = ();
if($frag_left_count > 0)
  {
do{
  my $frag_group = "$frag_left[0][0] $frag_left[0][1] $frag_left[0][2] $frag_left[0][3] $frag_left[0][4] $frag_left[0][5] $frag_left[0][6] $frag_left[0][7] $frag_left[0][8] $frag_left[0][9]";
  push @frag_group, [split(/ /, $frag_group)];
  for(my $z=1;$z<$frag_left_count;$z++)
    {
    my $frag_group_count = @frag_group;
    if($frag_left[$z][7] == $frag_group[0][7] && $frag_left[$z][8] == $frag_group[0][8] && $frag_left[$z][5] eq $frag_group[0][5] && $frag_left[$z][9] eq $frag_group[0][9] && $frag_group[$frag_group_count-1][3]-100 < $frag_left[$z][2])
      {
      my $frag_group = "$frag_left[$z][0] $frag_left[$z][1] $frag_left[$z][2] $frag_left[$z][3] $frag_left[$z][4] $frag_left[$z][5] $frag_left[$z][6] $frag_left[$z][7] $frag_left[$z][8] $frag_left[$z][9]";
      push @frag_group, [split(/ /, $frag_group)];
      }
     else
      {
      my $frag_hold = "$frag_left[$z][0] $frag_left[$z][1] $frag_left[$z][2] $frag_left[$z][3] $frag_left[$z][4] $frag_left[$z][5] $frag_left[$z][6] $frag_left[$z][7] $frag_left[$z][8] $frag_left[$z][9]";
      push @frag_hold, [split(/ /, $frag_hold)];
      }
    }
  my $frag_group_count = @frag_group;
  my $frag_pct = 0;
  for(my $y=0;$y<$frag_group_count;$y++)
    {
    $frag_pct = $frag_pct + ($frag_group[$y][3] - $frag_group[$y][2]);
    }
  if($frag_group_count > 0)
    {
    if($frag_pct/$frag_group[0][6] > .8)
      {
      my $n = 'n'.$nltr;
      my $fn_dir;
      if($frag_group[0][2] < $frag_group[0][3])
        {$fn_dir = 0;}
       elsif($frag_group[0][2] > $frag_group[0][3])
        {$fn_dir = 1;}
      $retro{nltr}{$n} = {type => $frag_group[0][5], dir => $fn_dir, bsr => 100, order => 100, level => 100, group => 100};
      for(my $y=0;$y<$frag_group_count;$y++)
        {
        $retro{nltr}{$n}{coords}{$y} = {SEQ_start => $frag_group[$y][0],
                                        SEQ_end => $frag_group[$y][1],
                                        TE_start => $frag_group[$y][2],
                                        TE_end => $frag_group[$y][3]};
        }
      $nltr++;
      }
     else
      {
      my $f = 'f'.$frag;
      my $fn_dir;
      if($frag_group[0][2] < $frag_group[0][3])
        {$fn_dir = 0;}
       elsif($frag_group[0][2] > $frag_group[0][3])
        {$fn_dir = 1;}
      $retro{frag}{$f} = {type => $frag_group[0][5], dir => $fn_dir, bsr => 100, order => 100, level => 100, group => 100};
      for(my $y=0;$y<$frag_group_count;$y++)
        {
        $retro{frag}{$f}{coords}{$y} = {SEQ_start => $frag_group[$y][0],
                                        SEQ_end => $frag_group[$y][1],
                                        TE_start => $frag_group[$y][2],
                                        TE_end => $frag_group[$y][3]};
        }
      $frag++;
      }
    }
  @frag_left = @frag_hold;
  $frag_left_count = @frag_left;
  @frag_group = ();
  @frag_hold = ();
  }until($frag_left_count == 0);
  }

#MAKE SURE THERE ARE NO COORDINATE NESTING DISCREPANCIES BETWEEN FRAGs and NLTRs
DISCREPANCY('FRAG-NLTR','FRAG-NLTR');
#MAKE SURE THERE ARE NO COORDINATE NESTING DISCREPANCIES BETWEEN FRAGs NLTRs AND SOLOs
DISCREPANCY('SOLO','FRAG-NLTR');
#FIX OVERLAPPING FRAG-NLTR AND SOLOS
OVERLAP_S_FN('SOLO','FRAG-NLTR');
}

#DETERMINE GROUP, ORDER, LEVEL AGAIN FOR ALL DATA TYPES
#SOLO
my $hash_solo = scalar keys %{$retro{solo}};
my @gol = ();
for(my $x=0;$x<$hash_solo;$x++)
  {
  my $s = 's'.$x;
  my $solo_count = scalar keys %{$retro{solo}{$s}{coords}};
  my $gol = "solo $s $retro{solo}{$s}{coords}{0}{SEQ_start} $retro{solo}{$s}{coords}{$solo_count-1}{SEQ_end}";
  push @gol, [split(/ /,$gol)];
  }
#PAIR
my $hash_pair = scalar keys %{$retro{pair}};
for(my $x=0;$x<$hash_pair;$x++)
  {
  my $p = 'p'.$x;
  my $pair_count = scalar keys %{$retro{pair}{$p}{coords}{R}};
  my $gol = "pair $p $retro{pair}{$p}{coords}{L}{0}{SEQ_start} $retro{pair}{$p}{coords}{R}{$pair_count-1}{SEQ_end}";
  push @gol, [split(/ /,$gol)];
  }
#FRAG
my $hash_frag = scalar keys %{$retro{frag}};
for(my $x=0;$x<$hash_frag;$x++)
  {
  my $f = 'f'.$x;
  my $frag_count = scalar keys %{$retro{frag}{$f}{coords}};
  my $gol = "frag $f $retro{frag}{$f}{coords}{0}{SEQ_start} $retro{frag}{$f}{coords}{$frag_count-1}{SEQ_end}";
  push @gol, [split(/ /,$gol)];
  }
#NLTR
my $hash_nltr = scalar keys %{$retro{nltr}};
for(my $x=0;$x<$hash_nltr;$x++)
  {
  my $n = 'n'.$x;
  my $nltr_count = scalar keys %{$retro{nltr}{$n}{coords}};
  my $gol = "nltr $n $retro{nltr}{$n}{coords}{0}{SEQ_start} $retro{nltr}{$n}{coords}{$nltr_count-1}{SEQ_end}";
  push @gol, [split(/ /,$gol)];
  }
@gol = sort{$a->[2] <=> $b->[2]} @gol;
my $gol_count = @gol;
$gol[0][4] = $gol[0][5] = $gol[0][6] = $zero = my $largest_level = 0;
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
  }
$nest_group = $gol[$gol_count-1][4];
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
for(my $x=0;$x<$gol_count;$x++)
  {
  $retro{$gol[$x][0]}{$gol[$x][1]}{group} = $gol[$x][4];
  $retro{$gol[$x][0]}{$gol[$x][1]}{order} = $gol[$x][5];
  $retro{$gol[$x][0]}{$gol[$x][1]}{level} = $gol[$x][6];
  }

#PRINT HASH
print BUG "PRINT HASH\n";
#SOLO
$hash_solo = scalar keys %{$retro{solo}};
print BUG "SOLOs $hash_solo\n";
for(my $x=0;$x<$hash_solo;$x++)
  {
  my $s = 's'.$x;
  print BUG "solo-$x type-$retro{solo}{$s}{type} dir-$retro{solo}{$s}{dir}";
  my $solo_count = scalar keys %{$retro{solo}{$s}{coords}};
  print BUG " parts-$solo_count\n";
  }
#PAIR
$hash_pair = scalar keys %{$retro{pair}};
print BUG "PAIRs\n";
for(my $x=0;$x<$hash_pair;$x++)
  {
  my $p = 'p' . $x;
  print BUG "pair-$x type-$retro{pair}{$p}{type} dir-$retro{pair}{$p}{dir} rate-$retro{pair}{$p}{bsr}";
  print BUG " group-$retro{pair}{$p}{group} level-$retro{pair}{$p}{level} order-$retro{pair}{$p}{order}\n";
  for(my $y=0;$y<3;$y++)
    {
    my $part = '';
    if($y==0)
      {$part = 'L';}
     elsif($y==1)
      {$part = 'R';}
     elsif($y==2)
      {$part = 'M';}
    my $pair_count = scalar keys %{$retro{pair}{$p}{coords}{$part}};
    print BUG " $part-parts-$pair_count";
    for(my $z=0;$z<$pair_count;$z++)
      {
      print BUG "  $retro{pair}{$p}{coords}{$part}{$z}{SEQ_start} $retro{pair}{$p}{coords}{$part}{$z}{SEQ_end}";
      }
    print BUG "\n";
    }
  }
#FRAG
$hash_frag = scalar keys %{$retro{frag}};
print BUG "FRAGs $hash_frag\n";
for(my $x=0;$x<$hash_frag;$x++)
  {
  my $f = 'f' . $x;
  print BUG "frag-$f type-$retro{frag}{$f}{type} dir-$retro{frag}{$f}{dir}";
  my $frag_count = scalar keys %{$retro{frag}{$f}{coords}};
  print BUG " parts-$frag_count\n";
  }
#NLTR
$hash_nltr = scalar keys %{$retro{nltr}};
print BUG "NLTRs $hash_nltr\n";
for(my $x=0;$x<$hash_nltr;$x++)
  {
  my $n = 'n'.$x;
  print BUG "nltr-$x type-$retro{nltr}{$n}{type} dir-$retro{nltr}{$n}{dir}";
  my $nltr_count = scalar keys %{$retro{nltr}{$n}{coords}};
  print BUG " parts-$nltr_count\n";
  }

# PRINT .LTR FILE
my $LTR = $input_fasta;
$LTR =~ s/.fasta//g;
$LTR = $LTR.'.LTR';
my $hash_file = "$output_directory/$date/$LTR";
open(OUT_HASH, ">$hash_file");
#PRINT HASH
print OUT_HASH "$seq_size\n";
print OUT_HASH "$input_fasta\n";
#SOLO
if($pair_only ne 'T')
  {
  $hash_solo = scalar keys %{$retro{solo}};
  for(my $x=0;$x<$hash_solo;$x++)
    {
    my $s = 's'.$x;
    print OUT_HASH "SOLO $s $te_list[$retro{solo}{$s}{type}] $retro{solo}{$s}{dir} $retro{solo}{$s}{level} $retro{solo}{$s}{order} $retro{solo}{$s}{level}\n";
    my $solo_count = scalar keys %{$retro{solo}{$s}{coords}};
    print OUT_HASH "$s";
    for(my $z=0;$z<$solo_count;$z++)
      {
      print OUT_HASH " $retro{solo}{$s}{coords}{$z}{SEQ_start} $retro{solo}{$s}{coords}{$z}{SEQ_end} $retro{solo}{$s}{coords}{$z}{TE_start} $retro{solo}{$s}{coords}{$z}{TE_end}";
      for(my $y=$retro{solo}{$s}{coords}{$z}{SEQ_start};$y<$retro{solo}{$s}{coords}{$z}{SEQ_end}+1;$y++)
        {
        $seq[$y] = 'N';
        }
      }
    print OUT_HASH "\n";
    }
  }
#PAIR
$hash_pair = scalar keys %{$retro{pair}};
for(my $x=0;$x<$hash_pair;$x++)
  {
  my $p = 'p' . $x;
  print OUT_HASH "PAIR $p $te_list[$retro{pair}{$p}{type}] $retro{pair}{$p}{dir} ";
  printf OUT_HASH ("%.3f", $retro{pair}{$p}{bsr});
  print OUT_HASH " $retro{pair}{$p}{group} $retro{pair}{$p}{order} $retro{pair}{$p}{level}\n";
  for(my $y=0;$y<3;$y++)
    {
    my $part = '';
    if($y==0)
      {$part = 'L';}
     elsif($y==1)
      {$part = 'R';}
     elsif($y==2)
      {$part = 'M';}
    my $pair_count = scalar keys %{$retro{pair}{$p}{coords}{$part}};
    print OUT_HASH "$p $part";
    for(my $z=0;$z<$pair_count;$z++)
      {
      print OUT_HASH " $retro{pair}{$p}{coords}{$part}{$z}{SEQ_start} $retro{pair}{$p}{coords}{$part}{$z}{SEQ_end} $retro{pair}{$p}{coords}{$part}{$z}{TE_start} $retro{pair}{$p}{coords}{$part}{$z}{TE_end}";
      for(my $y=$retro{pair}{$p}{coords}{$part}{$z}{SEQ_start};$y<$retro{pair}{$p}{coords}{$part}{$z}{SEQ_end}+1;$y++)
        {
        $seq[$y] = 'N';
        }
      }
    print OUT_HASH "\n";
    }
  }
#FRAG
$hash_frag = scalar keys %{$retro{frag}};
for(my $x=0;$x<$hash_frag;$x++)
  {
  my $f = 'f' . $x;
  print OUT_HASH "FRAG $f $te_list[$retro{frag}{$f}{type}] $retro{frag}{$f}{dir} $retro{frag}{$f}{group} $retro{frag}{$f}{order} $retro{frag}{$f}{level}\n";
  my $frag_count = scalar keys %{$retro{frag}{$f}{coords}};
  print OUT_HASH "$f";
  for(my $z=0;$z<$frag_count;$z++)
    {
    print OUT_HASH " $retro{frag}{$f}{coords}{$z}{SEQ_start} $retro{frag}{$f}{coords}{$z}{SEQ_end} $retro{frag}{$f}{coords}{$z}{TE_start} $retro{frag}{$f}{coords}{$z}{TE_end}";
    for(my $y=$retro{frag}{$f}{coords}{$z}{SEQ_start};$y<$retro{frag}{$f}{coords}{$z}{SEQ_end}+1;$y++)
      {
      $seq[$y] = 'N';
      }
    }
  print OUT_HASH "\n";
  }
#NLTR
$hash_nltr = scalar keys %{$retro{nltr}};
for(my $x=0;$x<$hash_nltr;$x++)
  {
  my $n = 'n'.$x;
  print OUT_HASH "NLTR $n $te_list[$retro{nltr}{$n}{type}] $retro{nltr}{$n}{dir} $retro{nltr}{$n}{group} $retro{nltr}{$n}{order} $retro{nltr}{$n}{level}\n";
  my $nltr_count = scalar keys %{$retro{nltr}{$n}{coords}};
  print OUT_HASH "$n";
  for(my $z=0;$z<$nltr_count;$z++)
    {
    print OUT_HASH " $retro{nltr}{$n}{coords}{$z}{SEQ_start} $retro{nltr}{$n}{coords}{$z}{SEQ_end} $retro{nltr}{$n}{coords}{$z}{TE_start} $retro{nltr}{$n}{coords}{$z}{TE_end}";
    for(my $y=$retro{nltr}{$n}{coords}{$z}{SEQ_start};$y<$retro{nltr}{$n}{coords}{$z}{SEQ_end}+1;$y++)
      {
      $seq[$y] = 'N';
      }
    }
  print OUT_HASH "\n";
  }
close(OUT_HASH);
close(BUG);
#PRINT MASKED SEQUENCE FILE
my $seq_out = "$output_directory/$date/$input_fasta.mask";
open(OUT_SEQ, ">$seq_out");
for($x=0;$x<($seq_size/70)-1;$x++)
  {
  for(my $y=0;$y<70;$y++)
    {
    print OUT_SEQ "$seq[($x*70)+$y]";
    }
  print OUT_SEQ "\n";
  }
for(my $y=$x*70;$y<$seq_size;$y++)
  {
  print OUT_SEQ "$seq[$y]";
  }
print OUT_SEQ "\n";
close(OUT_SEQ);

###########################################################################
#                             SUBROUTINES                                 #
###########################################################################

sub UNIQUE
{
my @lal_total = @_;
my $lal_total_count = @lal_total;
#SORT HITS BY DIRECTION AND LOCATION
if ($lal_total_count > 1)
  {
  @lal_total = sort {$a->[4] <=> $b->[4]} @lal_total;
  }
#UNIQUE LALIGN HITS, IF OVERLAPPING TAKE ONE WITH BIGGER SIZE AND BETTER SCORE
my @ltr_unique = ();
my $ltr_unique_count = 0;
### ltr type / dir / te start / te end / seq start / seq end / score ###
if ($lal_total_count > 1)
  {
  my $align = "$lal_total[0][0] $lal_total[0][1] $lal_total[0][2] $lal_total[0][3] $lal_total[0][4] $lal_total[0][5] $lal_total[0][6]";
  push @ltr_unique, [split(/ /,$align)];
  $ltr_unique_count++;
  for(my $y=0;$y<$lal_total_count;$y++)
    {
    if ($ltr_unique[$ltr_unique_count-1][5] - $ltr_ucoord_fudge >= $lal_total[$y][4])
      {
      my $u_size = $ltr_unique[$ltr_unique_count-1][5] - $ltr_unique[$ltr_unique_count-1][4];
      my $s_size = $lal_total[$y][5] - $lal_total[$y][4];
      my $u_score = $ltr_unique[$ltr_unique_count-1][6];
      my $s_score = $lal_total[$y][6];
      if(($s_size + 2 > $u_size && $s_score > $u_score) ||(($s_size + ($s_size * .2)) > $u_size && $s_score > $u_score))
        {
        $ltr_unique[$ltr_unique_count-1][0] = $lal_total[$y][0];
        $ltr_unique[$ltr_unique_count-1][1] = $lal_total[$y][1];
        $ltr_unique[$ltr_unique_count-1][2] = $lal_total[$y][2];
        $ltr_unique[$ltr_unique_count-1][3] = $lal_total[$y][3];
        $ltr_unique[$ltr_unique_count-1][4] = $lal_total[$y][4];
        $ltr_unique[$ltr_unique_count-1][5] = $lal_total[$y][5];
        $ltr_unique[$ltr_unique_count-1][6] = $lal_total[$y][6];
        }
      }
     else
      {
      $align = "$lal_total[$y][0] $lal_total[$y][1] $lal_total[$y][2] $lal_total[$y][3] $lal_total[$y][4] $lal_total[$y][5] $lal_total[$y][6]";
      push @ltr_unique, [split(/ /,$align)];
      $ltr_unique_count++;
      }
    }
  }
#IF ONLY ONE YOU CAN'T UNIQUE IT SILLY
if ($lal_total_count == 1)
  {
  my $align = "$lal_total[0][0] $lal_total[0][1] $lal_total[0][2] $lal_total[0][3] $lal_total[0][4] $lal_total[0][5] $lal_total[0][6]";
  push @ltr_unique, [split(/ /,$align)];
  $ltr_unique_count++;
  }
#SORT HITS BY LTR
if ($ltr_unique_count > 1)
  {
  @ltr_unique = sort {$a->[0] <=> $b->[0]} @ltr_unique;
  }
return(\@ltr_unique);
}
###########################################################################

sub LTR
#OPEN AND LOAD LTR SEQUENCE INTO @LTR
{
my $ltr_loc = $_[0];
my $amt_run = $_[1];
open(LTR_SEQ, "$org_directory/$ltr_loc") or die "$ltr_loc not found YO!\n";
my $ltr_seq = '';
while(<LTR_SEQ>)
  {
  my $line = $_;
  chomp($line);
  if ($line !~ /\>/)
    {$ltr_seq .= $line;}
  }
$ltr[$amt_run] = [split(//, $ltr_seq)];
my $ltr_size = scalar @{$ltr[$amt_run]};
print SIZE "$ltr_size\n";
return($ltr_size);
}
###########################################################################

sub FASTA
#OPEN AND LOAD FASTA FILE INTO @SEQ
{
my $fasta = $_[0];
my $fasta_string = '';
open(INFILE, "$output_directory/$fasta") or die "$fasta not found\n";
while(<INFILE>)
  {
  my $line = $_;
  chomp($line);
  if ($line !~ /\>/)
    {
    $fasta_string .= $line;
    }
  }
my @seq = split(//, $fasta_string);
return (\@seq);
}
###########################################################################

sub MULTI
#MULTI PROCESS SECTION
{
$SIG{CHLD} = \&REAPER;
$total = $_[0];
$running = $_[1];
my $amt_run = 0;
my $group = $_[2];
while ($running >= 0 && $amt_run < $total)
  {
  while ($running < $processors)
    {
    if (fork() == 0)
      {
      if($group == 0)
        {
        print BUG "Sending $ltr_list[$amt_run] to LTR_ANAL\n";
        LTR_ANAL($amt_run);
        }
       elsif($group == 1)
        {
        LTR_PWR($amt_run);
        }
       elsif($group == 2)
        {
        LTR_BSR($amt_run);
        }
       elsif($group == 3)
        {
        MID($amt_run);
        }
       elsif($group == 4)
        {
        FRAGS($amt_run);
        }
      exit 0;
      }
     else
      {
      $running++;
      $amt_run++;
      }
    }
  sleep;
  }
wait;
}
###########################################################################

sub BLAST_LTR
#BLAST THE @LTR TO @SEQ, PARSE RESULTS
{
my $ltr_loc = $_[0];
my $ltr_size = $_[1];
my $amt_run = $_[2];
system("$blast_directory/blastn $output_directory/$input_fasta $org_directory/$ltr_loc -noseqs > $output_directory/$date/ltr$amt_run.blast");
my $blast = "$output_directory/$date/ltr$amt_run.blast";
open(BLAST, $blast) or die "$blast not found\n";
my @b_hit = ();
#PARSE BLAST RESULTS SEND TO @b_hit(query start, query end, subj start, subj end)
while(<BLAST>)
  {
  my $b_line = $_;
  chomp($b_line);
  my @b_line = split(/\=/, $b_line);
  if ($b_line[3])
    {
    my $eval = $b_line[3];
    if ($eval eq ' 0.')
      {
      $eval = '0-1000';
      }
    my @eval = split(/-/, $eval);
    if ($eval[1])
      {
      if ($eval[1] >= 10)
        {
        $b_line = <BLAST>;
        $b_line = <BLAST>;
        $b_line = <BLAST>;
        chomp($b_line);
        my @query = split(/ /, $b_line);
        $b_line = <BLAST>;
        chomp($b_line);
        my @subj = split(/ /, $b_line);
        my $b_hit = "$query[1] $query[3] $subj[1] $subj[3]";
        push @b_hit, [split(/ /, $b_hit)];
#print BUG "Blast $b_hit\n";
        }
      }
    }
  }
system("rm $output_directory/$date/ltr$amt_run.blast");
#COMBINE OVERLAPPING BLAST HITS
@b_hit = sort {$a->[2] <=> $b->[2]} @b_hit;
my @b_combine = ();
my $b_hit_count = @b_hit;
my $b_combine_count = 0;
if ($b_hit_count > 0)
  {
  my $b_hit = "$b_hit[0][0] $b_hit[0][1] $b_hit[0][2] $b_hit[0][3]";
  push @b_combine, [split(/ /, $b_hit)];
  $b_combine_count = 1;
  for(my $x=1;$x<$b_hit_count;$x++)
    {
    if ($b_combine[$b_combine_count - 1][3] < $b_hit[$x][3])
      {
      if ($b_combine[$b_combine_count - 1][3] >= $b_hit[$x][2])
        {
        $b_combine[$b_combine_count - 1][3] = $b_hit[$x][3];
        }
       else
        {
        $b_hit = "$b_hit[$x][0] $b_hit[$x][1] $b_hit[$x][2] $b_hit[$x][3]";
        push @b_combine, [split(/ /, $b_hit)];
        $b_combine_count++;
        }
      }
    }
  }
#INCREASE BLAST HITS BY SIZE OF LTR
for(my $x=0;$x<$b_combine_count;$x++)
  {
  my $center = int(($b_combine[$x][3] - $b_combine[$x][2])/2) + $b_combine[$x][2];
  $b_combine[$x][2] = $center - $ltr_size;
  $b_combine[$x][3] = $center + $ltr_size;
  if ($b_combine[$x][2] < 1)
    {
    $b_combine[$x][2] = 1;
    }
  if ($b_combine[$x][3] > $seq_size)
    {
    $b_combine[$x][3] = $seq_size;
    }
  }
my @b_unique = ();
my $b_unique = '';
my $b_unique_count = 0;
if($b_combine_count>0)
  {
  $b_unique = "$b_combine[0][0] $b_combine[0][1] $b_combine[0][2] $b_combine[0][3]";
  push @b_unique, [split(/ /, $b_unique)];
  $b_unique_count++;
  }
for(my $x=1;$x<$b_combine_count;$x++)
  {
  if($b_combine[$x][2] <= $b_unique[$b_unique_count-1][3])
    {
    $b_unique[$b_unique_count-1][3] = $b_combine[$x][3];
    }
   else
    {
    $b_unique = "$b_combine[$x][0] $b_combine[$x][1] $b_combine[$x][2] $b_combine[$x][3]";
    push @b_unique, [split(/ /, $b_unique)];
    $b_unique_count++;
    }
  }
return (\@b_unique);
}
###########################################################################

sub ALIGN
{
my $amt_run = $_[0];
my $cut_start = $_[1];
my $cut_end = $_[2];
my $ltr_size = $_[3];
my $seq_align = "$output_directory/$date/seq_align$amt_run";
open(OUT_SEQ, ">$seq_align");
print OUT_SEQ ">SEQ\n";
#MAKE SEQUENCE FROM INPUT FASTA
for(my $y=$cut_start;$y<$cut_end;$y++)
  {
  print OUT_SEQ $seq[$y];
  }
close(OUT_SEQ);
my @lal_score = ();
my $lal_score_count = 0;
for(my $y=0;$y<2;$y++)
  {
  my $ltr_align = "$output_directory/$date/ltr_align$amt_run";
  open(OUT_LTR, ">$ltr_align");
  print OUT_LTR ">LTR\n";
##MAKE FORWARD LTR
  if($y==0)
    {
    for(my $z=0;$z<$ltr_size;$z++)
      {
      print OUT_LTR "$ltr[$amt_run][$z]";
      }
    }
##MAKE REVERSE LTR
  if ($y==1)
    {
    for(my $z=$ltr_size-1;$z>=0;$z--)
      {
      if ($ltr[$amt_run][$z] eq 'A' || $ltr[$amt_run][$z] eq 'a')
        {
        print OUT_LTR 'T';
        }
      if ($ltr[$amt_run][$z] eq 'T' || $ltr[$amt_run][$z] eq 't')
        {
        print OUT_LTR 'A';
        }
      if ($ltr[$amt_run][$z] eq 'C' || $ltr[$amt_run][$z] eq 'c')
        {
        print OUT_LTR 'G';
        }
      if ($ltr[$amt_run][$z] eq 'G' || $ltr[$amt_run][$z] eq 'g')
        {
        print OUT_LTR 'C';
        }
      }
    }
  close(OUT_LTR);
  system("$lal_directory/lalign -f $ltr_lal_f -g $ltr_lal_g -K $ltr_lal_k -n $ltr_align $seq_align > $output_directory/$date/lal$amt_run.out");
  system("rm $output_directory/$date/ltr_align$amt_run");
  my $lal_out = "$output_directory/$date/lal$amt_run.out";
  open(LAL, $lal_out) or die "$lal_out not found\n";
##PARSE ALIGNMENT
  while(<LAL>)
    {
    my $line = $_;
    chomp($line);
    if ($line =~ /ident/)
      {
      $line =~ s/\s+/ /g;
      my @line = split(/ /, $line);
      my $coords = $line[7];
      my $eval = $line[11];
      $coords =~ s/\(//g;
      $coords =~ s/\;//g;
      $coords =~ s/\)//g;
      $coords =~ s/\-/:/g;
      my @coords = split(/:/, $coords);
      if ($eval !~ /\+/)
        {
        if ($eval =~ /e-/)
          {
          my @eval = split(/-/, $eval);
          $eval = $eval[1];
          }
        if ($eval == 0)
          {
          $eval = 1000;
          }
        if ($eval >= $ltr_lal_score)
          {
          if($coords[2] < $coords[3])
            {
            my $start = $coords[2];
            my $end = $coords[3];
            }
           else
            {
            my $start = $coords[3];
            my $end = $coords[2];
            }
          my $align = "$y $coords[0] $coords[1] $coords[2] $coords[3] $eval";
          push @lal_score, [split(/ /,$align)];
          $lal_score_count++;
          }
        }
      }
    }
  close(LAL);
  system("rm $output_directory/$date/lal$amt_run.out");
  }
system("rm $output_directory/$date/seq_align$amt_run");
#SORT HITS BY DIRECTION AND LOCATION
if ($lal_score_count > 1)
  {
  @lal_score = sort {$a->[3] <=> $b->[3]} @lal_score;
  }
for(my $y=0;$y<$lal_score_count;$y++)
  {
  $lal_score[$y][3] = $lal_score[$y][3] + $cut_start;
  $lal_score[$y][4] = $lal_score[$y][4] + $cut_start;
#  print BUG "LTR ALIGN $lal_score[$y][0] $lal_score[$y][1] $lal_score[$y][2] $lal_score[$y][3] $lal_score[$y][4] $lal_score[$y][5]\n";
  }
#UNIQUE LALIGN HITS, IF OVERLAPPING TAKE ONE WITH BIGGER SIZE AND BETTER SCORE
my @lal_fin = ();
my $lal_fin_count = 0;
if ($lal_score_count > 1)
  {
  my $align = "$lal_score[0][0] $lal_score[0][1] $lal_score[0][2] $lal_score[0][3] $lal_score[0][4] $lal_score[0][5]";
  push @lal_fin, [split(/ /,$align)];
  $lal_fin_count++;
  for(my $y=0;$y<$lal_score_count;$y++)
    {
    if ($lal_fin[$lal_fin_count-1][4] - $ltr_unique_overlap >= $lal_score[$y][3])
      {
      my $u_size = $lal_fin[$lal_fin_count-1][4] - $lal_fin[$lal_fin_count-1][3];
      my $s_size = $lal_score[$y][4] - $lal_score[$y][3];
      if($s_size > $u_size)
        {
        if ($lal_fin[$lal_fin_count-1][5] > $lal_score[$y][5])
          {print BUG "WHOA, SCORE WAS BETTER IN UNIQUING $lal_fin[$lal_fin_count-1][5] > $lal_score[$y][5]\n";}
        $lal_fin[$lal_fin_count-1][0] = $lal_score[$y][0];
        $lal_fin[$lal_fin_count-1][1] = $lal_score[$y][1];
        $lal_fin[$lal_fin_count-1][2] = $lal_score[$y][2];
        $lal_fin[$lal_fin_count-1][3] = $lal_score[$y][3];
        $lal_fin[$lal_fin_count-1][4] = $lal_score[$y][4];
        $lal_fin[$lal_fin_count-1][5] = $lal_score[$y][5];
        }
       else
        {if ($lal_fin[$lal_fin_count-1][5] < $lal_score[$y][5])
          {print BUG "WHOA, SCORE WAS BETTER IN UNIQUING $lal_fin[$lal_fin_count-1][5] < $lal_score[$y][5]\n";}
        }
      }
     else
      {
      $align = "$lal_score[$y][0] $lal_score[$y][1] $lal_score[$y][2] $lal_score[$y][3] $lal_score[$y][4] $lal_score[$y][5]";
      push @lal_fin, [split(/ /,$align)];
      $lal_fin_count++;
      }
    }
  }
#IF ONLY ONE YOU CAN'T UNIQUE IT SILLY
if ($lal_score_count == 1)
  {
  my $align = "$lal_score[0][0] $lal_score[0][1] $lal_score[0][2] $lal_score[0][3] $lal_score[0][4] $lal_score[0][5]";
  push @lal_fin, [split(/ /,$align)];
  $lal_fin_count++;
  }
return(\@lal_fin);
}
###########################################################################

sub LTR_ANAL 
{
#RUN SUB LTR TO LOAD INTO @LTR
my $amt_run = $_[0];
my $ltr = $ltr_list[$amt_run];
my $ltr_file = $ltr . "-0.fasta";
my $ltr_loc = "LTRs/$ltr_file";
my $ltr_size = LTR($ltr_loc, $amt_run);
#RUN SUB BLAST TO BLAST @LTR TO @SEQ
my $b_unique_ref = BLAST_LTR($ltr_loc, $ltr_size, $amt_run);
my @b_unique = @$b_unique_ref;
my $b_unique_count = @b_unique;
#RUN LOCAL ALIGNMENT FOR EACH UNIQUE BLAST HIT
my @all_lal = ();
for (my $x=0;$x<$b_unique_count;$x++)
  {
  my $lal_hits_ref = ALIGN($amt_run, $b_unique[$x][2], $b_unique[$x][3], $ltr_size);
  my @lal_hits = @$lal_hits_ref;
  my $lal_hits_count = @lal_hits;
  for(my $i=0;$i<$lal_hits_count;$i++)
    {
    my $lal_hits = "$amt_run $lal_hits[$i][0] $lal_hits[$i][1] $lal_hits[$i][2] $lal_hits[$i][3] $lal_hits[$i][4] $lal_hits[$i][5]";
    push @all_lal, [split(/ /, $lal_hits)];
    }
  }
#WRITE TO FILE TO SEND TO MAIN
my $all_lal_count = @all_lal;
for(my $x=0;$x<$all_lal_count;$x++)
  {
  print OUT_LAL "$all_lal[$x][0] $all_lal[$x][1] $all_lal[$x][2] $all_lal[$x][3] $all_lal[$x][4] $all_lal[$x][5] $all_lal[$x][6]\n";
  }
sleep 1;
}
###########################################################################

sub LTR_PWR
{
my $which_ltr = $_[0];
my $ltr_unique_count = @ltr_unique;
my @ltr_just = ();
print BUG "LTR_PWR $which_ltr\n";
for(my $x=0;$x<$ltr_unique_count;$x++)
  {
  if($ltr_unique[$x][0] eq $which_ltr)
    {
    my $just = "$ltr_unique[$x][0] $ltr_unique[$x][1] $ltr_unique[$x][2] $ltr_unique[$x][3] $ltr_unique[$x][4] $ltr_unique[$x][5] $ltr_unique[$x][6]";
    push @ltr_just, [split(/ /, $just)];
print BUG "$just\n";
    }
  }
#MAKE POWERSET
my $ltr_just_count = @ltr_just;
my @word = ();
for (my $x=0;$x<$ltr_just_count;$x++)
  {
  push(@word, $x);
  }
my $power_ref = POWERSET(@word);
my @power = @$power_ref;
my $power_count = @power;
my @ltr_power = ();
my $ltr_power_count = 0;
for(my $x=0;$x<$power_count-1;$x++)
  {
#print "$power[$x][0]";
  my $power_len = scalar @{$power[$x]};
  my $pwr_set_ok = 1;
  if($power_len > 0)
    {
    for(my $y=1;$y<$power_len;$y++)
      {
#print " $power[$x][$y]";
      if($ltr_just[$power[$x][$y-1]][3] - $ltr_pwr_offset > $ltr_just[$power[$x][$y]][2] || $ltr_just[$power[$x][$y]][5] - $ltr_just[$power[$x][$y-1]][4] > $ltr_pwr_max)
        {
        $pwr_set_ok = 0;
        }
      }
    }
#print "\n";
  if($pwr_set_ok == 1)
    {
#print "powerset ok\n";
    my $pwr_set = "$power[$x][0]";
    my $pwr_size = $ltr_just[$power[$x][0]][3] - $ltr_just[$power[$x][0]][2] + 1;
    if($power_len > 0)
      {
      for(my $y=1;$y<$power_len;$y++)
        {
        $pwr_set = $pwr_set . " $power[$x][$y]";
        $pwr_size = $pwr_size + $ltr_just[$power[$x][$y]][3] - $ltr_just[$power[$x][$y]][2];
        }
      }
    my $pwr_group = "$pwr_size " . $pwr_set;
#print "$pwr_group\n";
    push @ltr_power, [split(/ /,$pwr_group)];
    my $power_len = scalar @{$ltr_power[$ltr_power_count]};
    for(my $y=2;$y<$power_len;$y++)
      {
#print " $ltr_power[$ltr_power_count][0] = $ltr_power[$ltr_power_count][0] - ($ltr_just[$ltr_power[$ltr_power_count][$y]][4] - $ltr_just[$ltr_power[$ltr_power_count][$y-1]][5])\n";
      $ltr_power[$ltr_power_count][0] = $ltr_power[$ltr_power_count][0] - ($ltr_just[$ltr_power[$ltr_power_count][$y]][4] - $ltr_just[$ltr_power[$ltr_power_count][$y-1]][5]);
      }
#print BUG "$x $ltr_power[$ltr_power_count][0]\n";
    $ltr_power_count++;
    }
  }
$ltr_power_count = @ltr_power;
my @pwr_unique = ();
my @pwr_u_set = ();
my $pwr_u_set_count = 0;
@ltr_power = sort {$b->[0] <=> $a->[0]} @ltr_power;
for(my $w=0;$w<$ltr_power_count;$w++)
  {
  my $power_found = 0;
  my $power_lenW = scalar @{$ltr_power[$w]};
#print BUG "$ltr_power[$w][0] ";
  for(my $x=1;$x<$power_lenW;$x++)
    {
#print BUG " $ltr_power[$w][$x] ";
    for(my $y=0;$y<$pwr_u_set_count;$y++)
      {
      my $power_lenY = scalar @{$pwr_u_set[$y]};
      for(my $z=0;$z<$power_lenY;$z++)
        {
        if($ltr_power[$w][$x] == $pwr_u_set[$y][$z])
          {
#print BUG "(ltr_power[$w][$x] == pwr_u_set[$y][$z], $ltr_power[$w][$x] == $pwr_u_set[$y][$z])";
          $power_found = 1;
          }
        }
      }
    }
#print BUG "\n";
  if($power_found == 0)
    {
    my $pwr_unique = "$ltr_just[$ltr_power[$w][1]][0] $ltr_just[$ltr_power[$w][1]][1] $ltr_just[$ltr_power[$w][1]][2] $ltr_just[$ltr_power[$w][1]][3] $ltr_just[$ltr_power[$w][1]][4] $ltr_just[$ltr_power[$w][1]][5] $ltr_just[$ltr_power[$w][1]][6]";
    my $pwr_u_set = "$ltr_power[$w][1]";
    my $power_lenW = scalar @{$ltr_power[$w]};
    if($power_lenW-1 > 0)
      {
      for(my $x=2;$x<$power_lenW;$x++)
        {
        $pwr_unique = $pwr_unique . " $ltr_just[$ltr_power[$w][$x]][0] $ltr_just[$ltr_power[$w][$x]][1] $ltr_just[$ltr_power[$w][$x]][2] $ltr_just[$ltr_power[$w][$x]][3] $ltr_just[$ltr_power[$w][$x]][4] $ltr_just[$ltr_power[$w][$x]][5] $ltr_just[$ltr_power[$w][$x]][6]";
        $pwr_u_set = $pwr_u_set . " $ltr_power[$w][$x]";
        }
      }
#print BUG "writing $pwr_u_set\n";
    push @pwr_unique, [split(/ /, $pwr_unique)];
    push @pwr_u_set, [split(/ /, $pwr_u_set)];
    $pwr_u_set_count++;
    }
  }
my $pwr_unique_count = @pwr_unique;
#WRITE TO FILE TO SEND TO MAIN
for(my $x=0;$x<$pwr_unique_count;$x++)
  {
  print OUT_PWR "$pwr_unique[$x][0]";
  print BUG "$pwr_unique[$x][0]";
  my $power_lenX = scalar @{$pwr_unique[$x]};
  for(my $y=1;$y<$power_lenX;$y++)
    {
    print OUT_PWR " $pwr_unique[$x][$y]";
    print BUG " $pwr_unique[$x][$y]";
    }
  print OUT_PWR "\n";
  print BUG "\n";
  }
sleep 1;
}
###########################################################################

sub POWERSET
{
return [[]] if @_ == 0;
my $first = shift;
my $pow = &POWERSET;
[ map { [$first, @$_ ], [ @$_] } @$pow ];
}
###########################################################################

sub REAPER
  {
  wait;
  $running--;
  $SIG{CHLD} = \&REAPER;
  }
###########################################################################

sub LTR_BSR
{
my $amt_run = $_[0];
open(SIZE, $ltr_size) or die "$ltr_size not found\n";
my @ltr_size = ();
while(<SIZE>)
  {
  my $line = $_;
  chomp($line);
  push @ltr_size, [split(/ /, $line)];
  }
close(SIZE);
#MAKE ARRAY FOR THIS LTR OF EACH LTR MATCH FOR EACH DIRECTION
my $pwr_total_count = @pwr_total;
for(my $w=0;$w<2;$w++)
  {
  my @ltr_set = ();
  for(my $x=0;$x<$pwr_total_count;$x++)
    {
    if($pwr_total[$x][1] == $amt_run)
      {
      if($pwr_total[$x][2] == $w)
        {    
        my $hit_size = 0;
        my $power_len = scalar @{$pwr_total[$x]};
        my $ltr_set = $pwr_total[$x][0];
        for(my $y=1;$y<$power_len;$y++)
          {
          $ltr_set = $ltr_set . " $pwr_total[$x][$y]";
          }
        for(my $y=0;$y<($power_len-1)/7;$y++)
          {
          $hit_size = $hit_size + $pwr_total[$x][($y*7)+6] - $pwr_total[$x][($y*7)+5]; 
          }
        if($hit_size > ($ltr_size[$amt_run][0]/2) || $hit_size > 300)
          {
          push @ltr_set, [split(/ /,$ltr_set)];
          }
        }
      }
    }
#GET THE SEQ FOR THIS LTR DIRECTION, WRITE TO ARRAY
  my $ltr_set_count = @ltr_set;
  for(my $x=0;$x<$ltr_set_count;$x++)
    {
    my $ltr_len = scalar @{$ltr_set[$x]};
    my $ltr_parts = $ltr_len / 8;
    my $ltr_seq = '';
    for(my $y=0;$y<$ltr_parts;$y++)
      {
      my $ltr_start = $ltr_set[$x][5+($y*7)];
      my $ltr_end = $ltr_set[$x][6+($y*7)];
      if($ltr_set[$x][2] == 0)
        {
        for(my $z=$ltr_start;$z<$ltr_end;$z++)
          {
          $ltr_seq = $ltr_seq . "$seq[$z]";
          $seq[$z] = 'N';
          }
        }
       elsif($ltr_set[$x][2] == 1)
        {
        for(my $z=$ltr_end-1;$z>=$ltr_start;$z--)
          {
          if ($seq[$z] eq 'A' || $seq[$z] eq 'a')
            {
            $ltr_seq = $ltr_seq . 'T';
            }
          if ($seq[$z] eq 'T' || $seq[$z] eq 't')
            {
            $ltr_seq = $ltr_seq . 'A';
            }
          if ($seq[$z] eq 'C' || $seq[$z] eq 'c')
            {
            $ltr_seq = $ltr_seq . 'G';
            }
          if ($seq[$z] eq 'G' || $seq[$z] eq 'g')
            {
            $ltr_seq = $ltr_seq . 'C';
            }
          $seq[$z] = 'N';
          }
        }
      }
    $ltr_set[$x][$ltr_len] = $ltr_seq;
    }
#ALIGN THE TWO LTRs, EXAMINE RESULTS
  my @bsr = ();
  for(my $x=0;$x<$ltr_set_count;$x++)
    {
    my $ltr_lenX = scalar @{$ltr_set[$x]};
    my $first = "$output_directory/$date/first$amt_run";
    open(OUT_SEQ, ">$first");
    print OUT_SEQ ">SEQ\n";
    print OUT_SEQ $ltr_set[$x][$ltr_lenX-1];
    close(OUT_SET);
    for(my $y=0;$y<$ltr_set_count;$y++)
      {
      my $second = "$output_directory/$date/second$amt_run";
      open(OUT_SEQ, ">$second");
      print OUT_SEQ ">SEQ\n";
      my $ltr_lenY = scalar @{$ltr_set[$y]};
      print OUT_SEQ $ltr_set[$y][$ltr_lenY-1];
      close(OUT_SET);
      $bsr[$x][$y] = 1;
      my $dont_do = 0;
      if ($x != $y)
        {
#print BUG "$amt_run $w $x $y \n";
        my $ltr_len = scalar @{$ltr_set[$x]};
        my $ltr_parts = $ltr_len / 8;
        my $ltr_startX = $ltr_set[$x][5];
        my $ltr_endX = $ltr_set[$x][6+(($ltr_parts-1)*7)];
        $ltr_len = scalar @{$ltr_set[$y]};
        $ltr_parts = $ltr_len / 8;
        my $ltr_startY = $ltr_set[$y][5];
        my $ltr_endY = $ltr_set[$y][6+(($ltr_parts-1)*7)];
        if($ltr_startY > $ltr_startX && $ltr_startY < $ltr_endX)
          {
          $dont_do = 1;
          }
         elsif($ltr_startX > $ltr_startY && $ltr_startX < $ltr_endY)
          {
          $dont_do = 1;
          }
        if($dont_do == 0)
          {
#print BUG "dont_do ok\n";
          my $ltr_lal = "$output_directory/$date/ltr_lal-$amt_run-$ltr_set[$x][0]-$ltr_set[$y][0]";
          system("$lal_directory/lalign -f $bsr_lal_f -g $bsr_lal_g -K 1 -n $first $second > $ltr_lal");
          open(LAL, "$ltr_lal")  or die "$ltr_lal not found\n";
          my $all = '';
          my $line = <LAL>;
          $line = <LAL>;
          $line = <LAL>;
          $line = <LAL>;
          $line = <LAL>;
          $line = <LAL>;
          if($line)
            {
            my @match = split(/ /, $line);
            my $match = $match[5];
            $line = <LAL>;
            $line = <LAL>;
            $line = <LAL>;
            while(<LAL>)
              {
              $line = $_;
              chomp($line);
              $line =~ s/^       //;
              $all = $all.$line;
              $line = <LAL>;
              $line = <LAL>;
              $line = <LAL>;
              $line = <LAL>;
              $line = <LAL>;
              }
            $all =~ s/\s+/ /g;
            my @Asplit = split(/ /, $all);
            my $all_size = @Asplit;
            my $bsr = '';
            if($match eq 0)
              {
              $bsr = 0;
              }
             elsif($match eq 'nt')
              {
              $bsr = 0;
              }
             else
              {
              $bsr = $all_size / $match;
print BUG "$w $x $y match $match size $all_size bsr $bsr\n";
              }
            my $mya = $bsr / .013;
            $bsr[$x][$y] = $bsr;
            }
#           else
#            {
#            print BUG "NO LINE\n";
#            }
          close(LAL);
          system("rm $ltr_lal");
          }
        }
      system("rm $output_directory/$date/second$amt_run");
      }
    system("rm $output_directory/$date/first$amt_run");
    }
  if($ltr_set_count > 1)
    {
    print BUG "BSR compare on direction $w type $amt_run for $ltr_set_count\n";
    for(my $x=0;$x<$ltr_set_count;$x++)
      {
      print BUG "$ltr_set[$x][0] ";
      for(my $y=0;$y<$ltr_set_count;$y++)
        {
        print BUG "$bsr[$x][$y] ";
        }
      print BUG "\n";
      }
    }
  my @pair = ();
  for(my $x=0;$x<$ltr_set_count;$x++)
    {
    my $small = $bsr[$x][0];
    my $set = "$x 0";
    for(my $y=1;$y<$ltr_set_count;$y++)
      {
      if ($bsr[$x][$y] < $small)
        {
        $small = $bsr[$x][$y];
        $set = "$x $y";
        }
      }
    my $small_set = "$set $small";
    push @pair, [split(/ /, $small_set)];
    }
  my @bsr_pair = ();
  for(my $z=.6;$z>.001;$z=$z-$bsr_compare_drop)
    {
    for(my $x=0;$x<$ltr_set_count;$x++)
      {
      for(my $y=0;$y<$ltr_set_count;$y++)
        {
        my $roundx = sprintf("%${z}f", $pair[$x][2]);
        my $roundy = sprintf("%${z}f", $pair[$y][2]);
        if ($pair[$x][0] eq $pair[$y][1] && $pair[$x][1] eq $pair[$y][0] &&
            $roundx eq $roundy && $pair[$x][0] < $pair[$x][1] && $pair[$x][2] < 1)
          {
          $a = $pair[$x][0];
          $b = $pair[$x][1];
          my $ltr_lenA = scalar @{$ltr_set[$a]};
          my $ltr_lenB = scalar @{$ltr_set[$b]};
          my $bsr_pair =  "$amt_run $a $ltr_set[$a][0] $ltr_set[$a][5] $ltr_set[$a][$ltr_lenA-3] $b $ltr_set[$b][0] $ltr_set[$b][5] $ltr_set[$b][$ltr_lenB-3] $pair[$x][2]";
          push @bsr_pair, [split(/ /, $bsr_pair)];
print BUG "found a pair at $z, $pair[$x][0] $pair[$x][1]\n";
          $pair[$x][0] = 1;
          $pair[$x][1] = 1;
          }
        }
      }
    }
#WRITE TO FILE TO SEND TO MAIN
  my $bsr_pair_count = @bsr_pair;
  for(my $x=0;$x<$bsr_pair_count;$x++)
    {
    print OUT_BSR "$bsr_pair[$x][0]";
    my $bsr_len = scalar @{$bsr_pair[$x]};
    for(my $y=1;$y<$bsr_len;$y++)
      {
      print OUT_BSR " $bsr_pair[$x][$y]";
      }
    print OUT_BSR "\n";
    }
  }
sleep 1;
}
###########################################################################

sub MID
{
my $amt_run = $_[0];
#MAKE ARRAY OF LTR SETS WITHIN NEST GROUP, CHANGE COORDS FOR BEGINNING AND END OF MID SECTION
$bsr_total_count = @bsr_total;
my @mid_set = ();
my $mid_count = 0;
for(my $x=0;$x<$bsr_total_count;$x++)
  {
  if($bsr_total[$x][0])
    {
    if($bsr_total[$x][11] == $amt_run)
      {
      my $mid_set = "$bsr_total[$x][0] $bsr_total[$x][1] $bsr_total[$x][2] $bsr_total[$x][3] $bsr_total[$x][4] $bsr_total[$x][5] $bsr_total[$x][6] $bsr_total[$x][7] $bsr_total[$x][8] $bsr_total[$x][9] $bsr_total[$x][10]";
      push @mid_set, [split(/ /, $mid_set)];
      $mid_set[$mid_count][11] = $mid_set[$mid_count][9] - $mid_set[$mid_count][4];
      $mid_count++;
      }
    }
  }
my $mid_set_count = @mid_set;
#SORT BY SIZE - START WITH SMALLEST FIRST
@mid_set = sort {$a->[11] <=> $b->[11]} @mid_set;
#MAKE QUERY FILE
my $to_align = "$output_directory/$date/cut_te$amt_run";
my $align_out = "$output_directory/$date/consensus_lal$amt_run";
for(my $y=0;$y<$mid_set_count;$y++)
  {
  open(OUTFILE, ">$to_align");
  print OUTFILE ">TE $y\n";
#print BUG "SEQ W/N LTRs FOR $mid_set[$y][0] IS $mid_set[$y][5] to $mid_set[$y][8]\n";
  for(my $x=$mid_set[$y][5];$x<$mid_set[$y][8];$x++)
    {
    print OUTFILE $seq[$x];
    }
  close(OUTFILE);
#SCORE FORWARD AND REVERSE ALIGNMENT
  my $te = '';
  my @unique_sls_for = my @unique_sls_rev = ();
  my $align_sum_for = my $align_sum_rev = 0;
  for(my $x=0;$x<2;$x++)
    {
    if ($x==0)
      {
      $te = "$ltr_list[$mid_set[$y][1]]";
      }
     elsif($x==1)
      {
      $te = "$ltr_list[$mid_set[$y][1]]_R";
      }
    system("$lal_directory/lalign -f $mid_lal_f -g $mid_lal_g -K $mid_lal_k -n $org_directory/TEs/$te $to_align > $align_out");
    open(INFILE, "$align_out") or die "$align_out not found\n";
    my @lal_score = ();
    while(<INFILE>)
      {
      my $line = $_;
      chomp($line);
      if ($line =~ /ident/)
        {
        $line =~ s/\s+/ /g;
        my @line = split(/ /, $line);
        my $coords = $line[7];
        my $eval = $line[11];
        $coords =~ s/\(//g;
        $coords =~ s/\;//g;
        $coords =~ s/\)//g;
        $coords =~ s/\-/:/g;
        my @coords = split(/:/, $coords);
        if ($eval !~ /\+/)
          {
          if ($eval =~ /e-/)
            {
            my @eval = split(/-/, $eval);
            $eval = $eval[1];
            }
          if ($eval == 0)
            {
            $eval = 1000;
            }
          if ($eval >= $mid_lal_score)
            {
            my $align = "$coords[0] $coords[1] $coords[2] $coords[3] $eval";
            push @lal_score, [split(/ /,$align)];
#print BUG "eval $align\n";
            }
          }
        }
      }
    close(INFILE);
    my $align_sum = 0;
    my $unique_no = 0;
    my $array = @lal_score;
    my @unique_sls = ();
#print BUG "unique\n";
    if ($array > 0)
      {
      @lal_score = sort {$a->[2] <=> $b->[2]} @lal_score;
      my $unique = "$mid_set[$y][0] $lal_score[0][0] $lal_score[0][1] $lal_score[0][2] $lal_score[0][3] $lal_score[0][4]";
#print BUG "$unique\n";
      push @unique_sls, [split(/ /,$unique)];
      $unique_no++;
      for(my $z=1;$z<$array;$z++)
        {
        if($unique_sls[$unique_no-1][4] < $lal_score[$z][3])
          {
          if($unique_sls[$unique_no-1][4] >= $lal_score[$z][2] && $unique_sls[$unique_no-1][2] < $lal_score[$z][1] && $unique_sls[$unique_no-1][2] >= $lal_score[$z][0])
            {
            $unique_sls[$unique_no-1][4] = $lal_score[$z][3];
            $unique_sls[$unique_no-1][2] = $lal_score[$z][1];
            }
           else
            {
            $unique = "$mid_set[$y][0] $lal_score[$z][0] $lal_score[$z][1] $lal_score[$z][2] $lal_score[$z][3] $lal_score[$z][4]";
#print BUG "$unique\n";
            push @unique_sls, [split(/ /,$unique)];
            $unique_no++;
            }
          }
        }
      for(my $z=0;$z<$unique_no;$z++)
        {
        my $diff = $unique_sls[$z][4] - $unique_sls[$z][3];
        $align_sum = $align_sum + $diff;
        }
      }
    system("rm $output_directory/$date/consensus_lal$amt_run");
    if ($x == 0)
      {
      @unique_sls_for = @unique_sls;
      $align_sum_for = $align_sum;
      }
     elsif ($x == 1)
      {
      @unique_sls_rev = @unique_sls;
      $align_sum_rev = $align_sum;
      }
    }
#COMPARE FOR AND REV ALIGNEMENTS
  my @use_align = ();
  if ($align_sum_rev < $align_sum_for)
    {
    @use_align = @unique_sls_for;
    }
   else
    {
    @use_align = @unique_sls_rev;
    }
  system("rm $output_directory/$date/cut_te$amt_run");
  my $use_align_count = @use_align;
  my @mid_align = ();
#print BUG "mid_set $mid_set[$y][0] $mid_set[$y][4] $mid_set[$y][9]\n";
  for(my $x=0;$x<$use_align_count;$x++)
    {
    $use_align[$x][3] = $use_align[$x][3] + $mid_set[$y][5];
    $use_align[$x][4] = $use_align[$x][4] + $mid_set[$y][5];
    my $safe = 0;
#print BUG "use $use_align[$x][0] $use_align[$x][1] $use_align[$x][2] $use_align[$x][3] $use_align[$x][4] $use_align[$x][5]\n";
    for(my $z=0;$z<$y;$z++)
      {
#print BUG "compare $use_align[$x][3] > $mid_set[$z][4] && $use_align[$x][4] < $mid_set[$z][9]\n";
      if($use_align[$x][3] > $mid_set[$z][4] && $use_align[$x][4] < $mid_set[$z][9])
        {$safe = 1;}
      }
    if($safe == 0)
      {
      my $mid_align = "$use_align[$x][0] $use_align[$x][1] $use_align[$x][2] $use_align[$x][3] $use_align[$x][4] $use_align[$x][5]";
      push @mid_align, [split(/ /,$mid_align)];
#print BUG "to_pwr $mid_align\n";
      }
    }
  my $mid_pwr_ref = MID_PWR($y, \@mid_align);
  my @mid_pwr = @$mid_pwr_ref;
  my $mid_pwr_count = @mid_pwr;
#WRITE 'N' FOR MID HITS IN @SEQ
  my @to_n = ();
  my $to_n;
  for(my $x=0;$x<$mid_pwr_count;$x++)
    {
    my $mid_len = scalar @{$mid_pwr[$x]};
    $mid_len = ($mid_len)/6;
    for(my $z=0;$z<$mid_len;$z++)
      {
      $to_n = "$mid_pwr[$x][($z*6)+3] $mid_pwr[$x][($z*6)+4]";
      push @to_n, [split(/ /,$to_n)];
      }
    }
  $to_n = "$mid_set[$y][4] $mid_set[$y][5]";
  push @to_n, [split(/ /,$to_n)];
  $to_n = "$mid_set[$y][8] $mid_set[$y][9]";
  push @to_n, [split(/ /,$to_n)];
  my $to_n_count = @to_n;
  for(my $x=0;$x<$to_n_count;$x++)
    {
    for(my $z=$to_n[$x][0];$z<=$to_n[$x][1];$z++)
      {
      $seq[$z] = 'N';
      }
    }
  for(my $x=0;$x<$mid_pwr_count;$x++)
    {
    print OUT_MID "$mid_pwr[$x][0]";
    my $power_len = scalar @{$mid_pwr[$x]};
    for(my $z=1;$z<$power_len;$z++)
      {
      print OUT_MID " $mid_pwr[$x][$z]";
      }
    print OUT_MID "\n";
    }
  }
sleep 1;
}
###########################################################################

sub MID_PWR
{
my ($amt_run, $mid_group_ref) = @_;
my @mid_group = @$mid_group_ref;
my $mid_group_count = @mid_group;
#MAKE POWERSET
my @word = ();
for (my $x=0;$x<$mid_group_count;$x++)
  {
  push(@word, $x);
#print BUG "in_pwr $mid_group[$x][0] $mid_group[$x][1] $mid_group[$x][2] $mid_group[$x][3] $mid_group[$x][4] $mid_group[$x][5]\n";
  }
my $power_ref = POWERSET(@word);
my @power = @$power_ref;
my $power_count = @power;
my @mid_power = ();
for(my $x=0;$x<$power_count-1;$x++)
  {
  my $power_len = scalar @{$power[$x]};
  my $pwr_set_ok = 1;
  if($power_len > 0)
    {
    for(my $y=1;$y<$power_len;$y++)
      {
      my $mid_offset;
      if($mid_group[$power[$x][$y-1]][2] > $mid_group[$power[$x][$y]][1])
        {$mid_offset = $mid_group[$power[$x][$y-1]][2]/10;}
       elsif($mid_group[$power[$x][$y-1]][2] <= $mid_group[$power[$x][$y]][1])
        {$mid_offset = $mid_group[$power[$x][$y]][1]/10;}
      if($mid_group[$power[$x][$y-1]][2] > $mid_group[$power[$x][$y]][1]+$mid_offset)
        {
        $pwr_set_ok = 0;
        }
      }
    }
  if($pwr_set_ok == 1)
    {
    my $pwr_set = "$power[$x][0]";
    my $pwr_size = $mid_group[$power[$x][0]][2] - $mid_group[$power[$x][0]][1] + 1;
    if($power_len > 0)
      {
      for(my $y=1;$y<$power_len;$y++)
        {
        $pwr_set = $pwr_set . " $power[$x][$y]";
        $pwr_size = $pwr_size + $mid_group[$power[$x][$y]][2] - $mid_group[$power[$x][$y]][1];
        }
      }
    my $pwr_group = "$pwr_size " . $pwr_set;
#print BUG "pwr_group $pwr_group\n";
    push @mid_power, [split(/ /,$pwr_group)];
    }
  }
my $mid_power_count = @mid_power;
my @mid_unique = ();
my @mid_u_set = ();
my $mid_u_set_count = 0;
for(my $w=0;$w<$mid_power_count;$w++)
  {
  my $power_score = 0;
  my $power_lenW = scalar @{$mid_power[$w]};
  for(my $x=1;$x<$power_lenW;$x++)
    {
    for(my $y=0;$y<$mid_u_set_count;$y++)
      {
      my $power_lenY = scalar @{$mid_u_set[$y]};
      for(my $z=0;$z<$power_lenY;$z++)
        {
        if($mid_power[$w][$x] == $mid_u_set[$y][$z])
          {
          $power_score = 1;
          }
        }
      }
    }
  if($power_score == 0)
    {
#print BUG "power long unique $w\n";
    my $mid_unique = "$mid_group[$mid_power[$w][1]][0] $mid_group[$mid_power[$w][1]][1] $mid_group[$mid_power[$w][1]][2] $mid_group[$mid_power[$w][1]][3] $mid_group[$mid_power[$w][1]][4] $mid_group[$mid_power[$w][1]][5]";
    my $mid_u_set = "$mid_power[$w][1]";
    my $power_lenW = scalar @{$mid_power[$w]};
    if($power_lenW-1 > 0)
      {
      for(my $x=2;$x<$power_lenW;$x++)
        {
        $mid_unique = $mid_unique . " $mid_group[$mid_power[$w][$x]][0] $mid_group[$mid_power[$w][$x]][1] $mid_group[$mid_power[$w][$x]][2] $mid_group[$mid_power[$w][$x]][3] $mid_group[$mid_power[$w][$x]][4] $mid_group[$mid_power[$w][$x]][5]";
        $mid_u_set = $mid_u_set . " $mid_power[$w][$x]";
        }
      }
    push @mid_unique, [split(/ /, $mid_unique)];
    push @mid_u_set, [split(/ /, $mid_u_set)];
    $mid_u_set_count++;
    }
  }
### IF MORE THAN ONE MID PWRSET FOUND, PICK LARGEST ONLY
my $mid_unique_count = @mid_unique;
if($mid_unique_count > 1)
  {
  my @mid_new = ();
  my $mid_u_largest = 0;
  my $power_lenB = scalar @{$mid_unique[0]};
  for(my $y=0;$y<$power_lenB-1;$y++)
    {
    $mid_new[0][$y] = $mid_unique[0][$y];
    }
  for($x=1;$x<$mid_unique_count;$x++)
    {
    my $power_lenA = scalar @{$mid_unique[$x-1]};
    my $power_lenB = scalar @{$mid_unique[$x]};
    if($mid_unique[$x-1][$power_lenA-1] < $mid_unique[$x][$power_lenB-1])
      {
      $mid_u_largest = $x;
#print "largest now $mid_unique[$x][$power_lenB-1]\n";
      for(my $y=0;$y<$power_lenB-1;$y++)
        {
        $mid_new[0][$y] = $mid_unique[$x][$y];
        }
      }
    }
#print "mid_u_largest is $mid_u_largest\n";
  @mid_unique = @mid_new;
  }
### END
return(\@mid_unique);
}
###########################################################################

sub FRAGS
{
#FOR EACH BLANK, ADD SEQ TO THE ARRAY
my $amt_run = $_[0];
$blank[$amt_run][2] = '';
for(my $z=$blank[$amt_run][0];$z<$blank[$amt_run][1];$z++)
  {
  $blank[$amt_run][2] = $blank[$amt_run][2] . $seq[$z];
  }
#FOR EACH BLANK ENTRY, BLAST TO GET POSSIBLE MATCHES
my $blast_seq = "$output_directory/$date/blast_seq$amt_run";
open(BLAST_SEQ, ">$blast_seq");
print BLAST_SEQ ">SEQ\n";
print BLAST_SEQ $blank[$amt_run][2];
close(BLAST_SEQ);
system("$blast_directory/blastn $org_directory/unique_repeats $blast_seq -noseqs -V 10000 -B 10000 -S 500 > $output_directory/$date/frag.blast$amt_run");
system("rm $output_directory/$date/blast_seq$amt_run");
my $blast = "$output_directory/$date/frag.blast$amt_run";
open(BLAST, $blast) or die "$blast not found\n";
my @subj = ();
while(<BLAST>)
  {
  my $b_line = $_;
  chomp($b_line);
  if ($b_line =~ /^\>/)
    {
    $b_line =~ s/\>//;
    my @b_line = split(/ /, $b_line);
    push(@subj, $b_line[0]);
    }
  }
close(BLAST);
my $subj_count = @subj;
for(my $y=0;$y<$subj_count;$y++)
  {
  print BUG "BLANK $amt_run $subj[$y]\n";
  }
system("rm $output_directory/$date/frag.blast$amt_run");
#FOR EACH POSSIBLE MATCH FOR EACH BLANK, LALIGN TO FIND BEST MATCH
my @lal_hits = ();
for(my $y=0;$y<$subj_count;$y++)
  {
  my $te_size;
  my $second;
  my $ok_frag = 0;
  for(my $v=0;$v<@te_list;$v++)
    {
    if($te_list[$v] eq $subj[$y])
      {$ok_frag = 1;}
    }
  if($ok_frag == 1)
    {
    for(my $z=0;$z<2;$z++)
      {
      my $first = "$output_directory/$date/first$amt_run";
      open(OUT_SEQ, ">$first");
      print OUT_SEQ ">SEQ\n";
      print OUT_SEQ $blank[$amt_run][2];
      close(OUT_SEQ);
      if($z==0)
        {
        $second = "$org_directory/CON-1/".$subj[$y];
        open(SIZE, $second) or die "$second not found\n";
        my $te_seq = '';
        while(<SIZE>)
          {
          my $lsize = $_;
          chomp($lsize);
          if($lsize !~ />/)
            {
            $te_seq = $te_seq . $lsize;
            }
          }
        my @te_size = split(//, $te_seq);
        $te_size = @te_size;
        }
       elsif($z==1)
        {
        my $rev = "$output_directory/$date/rev$amt_run";
        open(REV, ">$rev");
        open(FOR, $second)  or die "$second not found\n";
        my $rev_seq;
        while(<FOR>)
          {
          my $line = $_;
          chomp($line);
          $rev_seq = $rev_seq . $line;
          }
        close(FOR);
        my @rev_seq = split(//, $rev_seq);
        my $rev_seq_count = @rev_seq;
        print REV ">REV\n";
        for(my $w=$rev_seq_count-1;$w>=0;$w--)
          {
          if ($rev_seq[$w] eq 'A' || $rev_seq[$w] eq 'a')
            {
            print REV 'T';
            }
          if ($rev_seq[$w] eq 'T' || $rev_seq[$w] eq 't')
            {
            print REV 'A';
            }
          if ($rev_seq[$w] eq 'C' || $rev_seq[$w] eq 'c')
            {
            print REV 'G';
            }
          if ($rev_seq[$w] eq 'G' || $rev_seq[$w] eq 'g')
            {
            print REV 'C';
            }
          }
        close(REV);
        $second = "$rev";
        }
      system("$lal_directory/lalign -f 30 -g 1 -n $first $second > $output_directory/$date/frag.lal$amt_run");
      system("rm $output_directory/$date/first$amt_run");
      my $lalign = "$output_directory/$date/frag.lal$amt_run";
      open(LAL, $lalign) or die "$lalign not found\n";
      my $lal_line;
      my @lal_line = ();
      my $eval;
      while(<LAL>)
        {
        $lal_line = $_;
        chomp($lal_line);
        if($lal_line =~ /ident/)
          {
          $lal_line =~ s/\s+/ /g;
          @lal_line = split(/ /, $lal_line);
          $eval = $lal_line[11];
          $lal_line[7] =~ s/\(//;
          $lal_line[7] =~ s/\)//;
          $lal_line[7] =~ s/\;//;
          my @lal_coord = split(/:/, $lal_line[7]);
          my @lal_subj = split(/-/,$lal_coord[0]);
          my @lal_query = split(/-/,$lal_coord[1]);
          if($eval !~ /\+/)
            {
            if ($eval =~ /e-/)
              {
              my @eval = split(/-/, $eval);
              $eval = $eval[1];
              }
            if ($eval == 0)
              {
              $eval = 1000;
              }
            if ($eval >= 30)
              {
              if($lal_subj[1] - $lal_subj[0] > 100)
                {
                my $hit_size = $lal_subj[1] - $lal_subj[0];
                my $lal_hits = "$lal_subj[0] $lal_subj[1] $lal_query[0] $lal_query[1] $eval $subj[$y] $te_size $hit_size";
                push @lal_hits, [split(/ /, $lal_hits)];
                }
              }
            }
          }
        }
      close(LAL);
      system("rm $output_directory/$date/frag.lal$amt_run");
      }
    }
  }

##REMOVE OVERLAPPING HITS IN @lal_hits
my $lal_hits_count = @lal_hits;
for(my $y=0;$y<$lal_hits_count;$y++)
  {
  print BUG "$amt_run $y $lal_hits[$y][0] $lal_hits[$y][1] $lal_hits[$y][2] $lal_hits[$y][3] $lal_hits[$y][4] $lal_hits[$y][5] $lal_hits[$y][6] $lal_hits[$y][7]\n";
  }
if($lal_hits_count > 1)
  {
  @lal_hits = sort {$a->[0] <=> $b->[0]} @lal_hits;
  }
for(my $y=0;$y<$lal_hits_count;$y++)
  {
  $lal_hits[$y][0] = $lal_hits[$y][0] + $blank[$amt_run][0];
  $lal_hits[$y][1] = $lal_hits[$y][1] + $blank[$amt_run][0]; 
  my $frag_size = $lal_hits[$y][1] - $lal_hits[$y][0];
#  print OUT_FRAG "$lal_hits[$y][0] $lal_hits[$y][1] $lal_hits[$y][2] $lal_hits[$y][3] $lal_hits[$y][4] $lal_hits[$y][5] $lal_hits[$y][6] $lal_hits[$y][7]\n";
  print BUG "$lal_hits[$y][0] $lal_hits[$y][1] $lal_hits[$y][2] $lal_hits[$y][3] $lal_hits[$y][4] $lal_hits[$y][5] $lal_hits[$y][6] $lal_hits[$y][7]\n";
  }
my @frag_set = ();
my @blank_sections = ();
my $blank_sections = "$blank[$amt_run][0] $blank[$amt_run][1]";
push @blank_sections, [split(/ /,$blank_sections)];
for(my $x=0;$x<@blank_sections;$x++)
  {
  my @within = ();
  for(my $y=0;$y<@lal_hits;$y++)
    {
    if(($lal_hits[$y][0] > $blank_sections[$x][0] && $lal_hits[$y][0] < $blank_sections[$x][1]) || ($lal_hits[$y][1] > $blank_sections[$x][0] && $lal_hits[$y][1] < $blank_sections[$x][1]) || ($lal_hits[$y][0] > $blank_sections[$x][0] && $lal_hits[$y][1] < $blank_sections[$x][1]))
      {
      my $within = "$lal_hits[$y][0] $lal_hits[$y][1] $lal_hits[$y][2] $lal_hits[$y][3] $lal_hits[$y][4] $lal_hits[$y][5] $lal_hits[$y][6] $lal_hits[$y][7]";
      push @within, [split(/ /,$within)];
      }
    }
  if(@within > 0)
    {
    @within = sort {$b->[7] <=> $a->[7]} @within;
    @within = sort {$b->[4] <=> $a->[4]} @within;
    if($within[0][0] < $blank_sections[$x][0])
      {
      $within[0][0] = $blank_sections[$x][0];
      }
    if($within[0][1] > $blank_sections[$x][1])
      {
      $within[0][1] = $blank_sections[$x][1];
      }
    my $frag_set = "$within[0][0] $within[0][1] $within[0][2] $within[0][3] $within[0][4] $within[0][5] $within[0][6] $within[0][7]";
    push @frag_set, [split(/ /,$frag_set)];
    if(($within[0][0] - $blank_sections[$x][0]) > 25)
      {
      my $end = $within[0][0]-1;
      my $blank_sections = "$blank_sections[$x][0] $end";
      push @blank_sections, [split(/ /,$blank_sections)];
      }
    if(($blank_sections[$x][1] - $within[0][1]) > 25)
      {
      my $start = $within[0][1]+1;
      my $blank_sections = "$start $blank_sections[$x][1]";
      push @blank_sections, [split(/ /,$blank_sections)];
      }
    }
  }
my $frag_set_count = @frag_set;
for(my $y=0;$y<$frag_set_count;$y++)
  {
  my $frag_size = $frag_set[$y][1] - $frag_set[$y][0];
  if($frag_size > 100 )
    {
    print OUT_FRAG "$frag_set[$y][0] $frag_set[$y][1] $frag_set[$y][2] $frag_set[$y][3] $frag_set[$y][4] $frag_set[$y][5] $frag_set[$y][6] $frag_size\n";
    }
  }
}
###########################################################################

sub DATE
{
my ($sec,$min,$hour,$mday,$mon,$year) = localtime(time);
$year = $year -100;
my $date = '';
if($year < 10)
  {$date .= 0;}
$date .= $year;
if($mon < 10)
  {$date .= 0;}
$date .= $mon+1;
if($mday < 10)
  {$date .= 0;}
$date .= $mday;
$date .= '-';
if($hour < 10)
  {$date .= 0;}
$date .= $hour;
if($min < 10)
  {$date .= 0;}
$date .= $min;
if($sec < 10)
  {$date .= 0;}
$date .= $sec;
return($date);
}
###########################################################################

sub LOAD_DISCREPANCY
{
my $first_dis = $_[0];
my $second_dis = $_[1];
@first_dis = ();
@second_dis = ();
@first_head = ();
@second_head = ();
if($first_dis eq 'SOLO')
  {
  for(my $x=0;$x<scalar keys %{$retro{solo}};$x++)
    {
    my $p = "s".$x;
    $first_head[$x][1] = $retro{solo}{$p}{type};
    $first_head[$x][2] = $retro{solo}{$p}{dir};
    $first_head[$x][3] = $retro{solo}{$p}{bsr};
    $first_head[$x][4] = $retro{solo}{$p}{group};
    $first_head[$x][5] = $retro{solo}{$p}{level};
    $first_head[$x][6] = $retro{solo}{$p}{order};
    $first_head[$x][0] = $p;
    for(my $y=0;$y<scalar keys %{$retro{solo}{$p}{coords}};$y++)
      {
      $first_dis[$x][$y*4] = $retro{solo}{$p}{coords}{$y}{SEQ_start};
      $first_dis[$x][($y*4)+1] = $retro{solo}{$p}{coords}{$y}{SEQ_end};
      $first_dis[$x][($y*4)+2] = $retro{solo}{$p}{coords}{$y}{TE_start};
      $first_dis[$x][($y*4)+3] = $retro{solo}{$p}{coords}{$y}{TE_end};
      }
    }
  }
 elsif($first_dis eq 'PAIR')
  {
  for(my $x=0;$x<scalar keys %{$retro{pair}};$x++)
    {
    my $p = "p".$x;
    $first_head[$x][1] = $retro{pair}{$p}{type};
    $first_head[$x][2] = $retro{pair}{$p}{dir};
    $first_head[$x][3] = $retro{pair}{$p}{bsr};
    $first_head[$x][4] = $retro{pair}{$p}{group};
    $first_head[$x][5] = $retro{pair}{$p}{level};
    $first_head[$x][6] = $retro{pair}{$p}{order};
    $first_head[$x][0] = $p;
    for(my $y=0;$y<scalar keys %{$retro{pair}{$p}{coords}{L}};$y++)
      {
      $first_dis[$x*3][$y*4] = $retro{pair}{$p}{coords}{L}{$y}{SEQ_start};
      $first_dis[$x*3][($y*4)+1] = $retro{pair}{$p}{coords}{L}{$y}{SEQ_end};
      $first_dis[$x*3][($y*4)+2] = $retro{pair}{$p}{coords}{L}{$y}{TE_start};
      $first_dis[$x*3][($y*4)+3] = $retro{pair}{$p}{coords}{L}{$y}{TE_end};
      }
    for(my $y=0;$y<scalar keys %{$retro{pair}{$p}{coords}{R}};$y++)
      {
      $first_dis[($x*3)+1][$y*4] = $retro{pair}{$p}{coords}{R}{$y}{SEQ_start};
      $first_dis[($x*3)+1][($y*4)+1] = $retro{pair}{$p}{coords}{R}{$y}{SEQ_end};
      $first_dis[($x*3)+1][($y*4)+2] = $retro{pair}{$p}{coords}{R}{$y}{TE_start};
      $first_dis[($x*3)+1][($y*4)+3] = $retro{pair}{$p}{coords}{R}{$y}{TE_end};
      }
    if(scalar keys %{$retro{pair}{$p}{coords}{M}} == 0)
      {
      $first_dis[($x*3)+2][0] = 'blank';
      }
    for(my $y=0;$y<scalar keys %{$retro{pair}{$p}{coords}{M}};$y++)
      {
      $first_dis[($x*3)+2][$y*4] = $retro{pair}{$p}{coords}{M}{$y}{SEQ_start};
      $first_dis[($x*3)+2][($y*4)+1] = $retro{pair}{$p}{coords}{M}{$y}{SEQ_end};
      $first_dis[($x*3)+2][($y*4)+2] = $retro{pair}{$p}{coords}{M}{$y}{TE_start};
      $first_dis[($x*3)+2][($y*4)+3] = $retro{pair}{$p}{coords}{M}{$y}{TE_end};
      }
    }
  }
 elsif($first_dis eq 'FRAG-NLTR')
  {
  my $z=0;
  for(my $x=0;$x<(scalar keys %{$retro{frag}})+(scalar keys %{$retro{nltr}});$x++)
    {
    my $type;
    my $p;
    if($x < (scalar keys %{$retro{frag}}))
      {
      $type = 'frag';
      $p = "f".$x;
      }
     elsif($x >= (scalar keys %{$retro{frag}}))
      {
      $type = 'nltr';
      $p = "n".$z;
      $z++;
      }
    $first_head[$x][1] = $retro{$type}{$p}{type};
    $first_head[$x][2] = $retro{$type}{$p}{dir};
    $first_head[$x][3] = $retro{$type}{$p}{bsr};
    $first_head[$x][4] = $retro{$type}{$p}{group};
    $first_head[$x][5] = $retro{$type}{$p}{level};
    $first_head[$x][6] = $retro{$type}{$p}{order};
    $first_head[$x][0] = $p;
    for(my $y=0;$y<scalar keys %{$retro{$type}{$p}{coords}};$y++)
      {
      $first_dis[$x][$y*4] = $retro{$type}{$p}{coords}{$y}{SEQ_start};
      $first_dis[$x][($y*4)+1] = $retro{$type}{$p}{coords}{$y}{SEQ_end};
      $first_dis[$x][($y*4)+2] = $retro{$type}{$p}{coords}{$y}{TE_start};
      $first_dis[$x][($y*4)+3] = $retro{$type}{$p}{coords}{$y}{TE_end};
      }
    }
  }
if($second_dis eq 'SOLO')
  {
  for(my $x=0;$x<scalar keys %{$retro{solo}};$x++)
    {
    my $p = "s".$x;
    $second_head[$x][1] = $retro{solo}{$p}{type};
    $second_head[$x][2] = $retro{solo}{$p}{dir};
    $second_head[$x][3] = $retro{solo}{$p}{bsr};
    $second_head[$x][4] = $retro{solo}{$p}{group};
    $second_head[$x][5] = $retro{solo}{$p}{level};
    $second_head[$x][6] = $retro{solo}{$p}{order};
    $second_head[$x][0] = $p;
    for(my $y=0;$y<scalar keys %{$retro{solo}{$p}{coords}};$y++)
      {
      $second_dis[$x][$y*4] = $retro{solo}{$p}{coords}{$y}{SEQ_start};
      $second_dis[$x][($y*4)+1] = $retro{solo}{$p}{coords}{$y}{SEQ_end};
      $second_dis[$x][($y*4)+2] = $retro{solo}{$p}{coords}{$y}{TE_start};
      $second_dis[$x][($y*4)+3] = $retro{solo}{$p}{coords}{$y}{TE_end};
      }
    }
  }
 elsif($second_dis eq 'PAIR')
  {
  for(my $x=0;$x<scalar keys %{$retro{pair}};$x++)
    {
    my $p = "p".$x;
    $second_head[$x][1] = $retro{pair}{$p}{type};
    $second_head[$x][2] = $retro{pair}{$p}{dir};
    $second_head[$x][3] = $retro{pair}{$p}{bsr};
    $second_head[$x][4] = $retro{pair}{$p}{group};
    $second_head[$x][5] = $retro{pair}{$p}{level};
    $second_head[$x][6] = $retro{pair}{$p}{order};
    $second_head[$x][0] = $p;
    for(my $y=0;$y<scalar keys %{$retro{pair}{$p}{coords}{L}};$y++)
      {
      $second_dis[$x*3][$y*4] = $retro{pair}{$p}{coords}{L}{$y}{SEQ_start};
      $second_dis[$x*3][($y*4)+1] = $retro{pair}{$p}{coords}{L}{$y}{SEQ_end};
      $second_dis[$x*3][($y*4)+2] = $retro{pair}{$p}{coords}{L}{$y}{TE_start};
      $second_dis[$x*3][($y*4)+3] = $retro{pair}{$p}{coords}{L}{$y}{TE_end};
      }
    for(my $y=0;$y<scalar keys %{$retro{pair}{$p}{coords}{R}};$y++)
      {
      $second_dis[($x*3)+1][$y*4] = $retro{pair}{$p}{coords}{R}{$y}{SEQ_start};
      $second_dis[($x*3)+1][($y*4)+1] = $retro{pair}{$p}{coords}{R}{$y}{SEQ_end};
      $second_dis[($x*3)+1][($y*4)+2] = $retro{pair}{$p}{coords}{R}{$y}{TE_start};
      $second_dis[($x*3)+1][($y*4)+3] = $retro{pair}{$p}{coords}{R}{$y}{TE_end};
      }
    if(scalar keys %{$retro{pair}{$p}{coords}{M}} == 0)
      {
      $second_dis[($x*3)+2][0] = 'blank';
      }
    for(my $y=0;$y<scalar keys %{$retro{pair}{$p}{coords}{M}};$y++)
      {
      $second_dis[($x*3)+2][$y*4] = $retro{pair}{$p}{coords}{M}{$y}{SEQ_start};
      $second_dis[($x*3)+2][($y*4)+1] = $retro{pair}{$p}{coords}{M}{$y}{SEQ_end};
      $second_dis[($x*3)+2][($y*4)+2] = $retro{pair}{$p}{coords}{M}{$y}{TE_start};
      $second_dis[($x*3)+2][($y*4)+3] = $retro{pair}{$p}{coords}{M}{$y}{TE_end};
      }
    }
  }
 elsif($second_dis eq 'FRAG-NLTR')
  {
  my $fn_size = (scalar keys %{$retro{frag}})+(scalar keys %{$retro{nltr}});
  my $z=0;
  for(my $x=0;$x<$fn_size;$x++)
    {
    my $type;
    my $p;
    if($x < (scalar keys %{$retro{frag}}))
      {
      $type = 'frag';
      $p = "f".$x;
      }
     elsif($x >= (scalar keys %{$retro{frag}}))
      {
      $type = 'nltr';
      $p = "n".$z;
      $z++;
      }
    $second_head[$x][1] = $retro{$type}{$p}{type};
    $second_head[$x][2] = $retro{$type}{$p}{dir};
    $second_head[$x][3] = $retro{$type}{$p}{bsr};
    $second_head[$x][4] = $retro{$type}{$p}{group};
    $second_head[$x][5] = $retro{$type}{$p}{level};
    $second_head[$x][6] = $retro{$type}{$p}{order};
    $second_head[$x][0] = $p;
    for(my $y=0;$y<scalar keys %{$retro{$type}{$p}{coords}};$y++)
      {
      $second_dis[$x][$y*4] = $retro{$type}{$p}{coords}{$y}{SEQ_start};
      $second_dis[$x][($y*4)+1] = $retro{$type}{$p}{coords}{$y}{SEQ_end};
      $second_dis[$x][($y*4)+2] = $retro{$type}{$p}{coords}{$y}{TE_start};
      $second_dis[$x][($y*4)+3] = $retro{$type}{$p}{coords}{$y}{TE_end};
      }
    }
  }
if($first_dis eq 'SOLO')
  {delete($retro{solo});}
 elsif($first_dis eq 'PAIR')
  {delete($retro{pair});}
 elsif($first_dis eq 'FRAG-NLTR')
  {delete($retro{frag});
   delete($retro{nltr});}
if($second_dis eq 'SOLO')
  {delete($retro{solo});}
 elsif($second_dis eq 'PAIR')
  {delete($retro{pair});}
 elsif($second_dis eq 'FRAG-NLTR')
  {delete($retro{frag});
   delete($retro{nltr});}
}
###########################################################################

sub UNLOAD_DISCREPANCY
{
my $first_dis = $_[0];
if($first_dis eq 'SOLO')
  {
  my $count = @first_head;
  for(my $x=0;$x<$count;$x++)
    {
    my $p = "s".$x;
    $retro{solo}{$p} = {type => $first_head[$x][1], dir => $first_head[$x][2], bsr => $first_head[$x][3],
                        group => $first_head[$x][4], level => $first_head[$x][5], order => $first_head[$x][6]};
    my $power_len = (scalar @{$first_dis[$x]})/4;
    for(my $y=0;$y<$power_len;$y++)
      {
      $retro{solo}{$p}{coords}{$y} = {SEQ_start => $first_dis[$x][($y*4)],
                                      SEQ_end => $first_dis[$x][($y*4)+1],
                                      TE_start => $first_dis[$x][($y*4)+2],
                                      TE_end => $first_dis[$x][($y*4)+3]};
      }
    }
  }
 elsif($first_dis eq 'PAIR')
  {
  my $count = @first_head;
  for(my $x=0;$x<$count;$x++)
    {
    my $p = "p".$x;
    $retro{pair}{$p} = {type => $first_head[$x][1], dir => $first_head[$x][2], bsr => $first_head[$x][3],
                        group => $first_head[$x][4], level => $first_head[$x][5], order => $first_head[$x][6]};
    my $power_len = (scalar @{$first_dis[$x*3]})/4;
    for(my $y=0;$y<$power_len;$y++)
      {
      $retro{pair}{$p}{coords}{L}{$y} = {SEQ_start => $first_dis[$x*3][($y*4)],
                                         SEQ_end => $first_dis[$x*3][($y*4)+1],
                                         TE_start => $first_dis[$x*3][($y*4)+2],
                                         TE_end => $first_dis[$x*3][($y*4)+3]};
      }
    $power_len = (scalar @{$first_dis[($x*3)+1]})/4;
    for(my $y=0;$y<$power_len;$y++)
      {
      $retro{pair}{$p}{coords}{R}{$y} = {SEQ_start => $first_dis[($x*3)+1][($y*4)],
                                         SEQ_end => $first_dis[($x*3)+1][($y*4)+1],
                                         TE_start => $first_dis[($x*3)+1][($y*4)+2],
                                         TE_end => $first_dis[($x*3)+1][($y*4)+3]};
      }
    if($first_dis[($x*3)+2][0] ne 'blank')
      {
      $power_len = (scalar @{$first_dis[($x*3)+2]})/4;
      for(my $y=0;$y<$power_len;$y++)
        {
        $retro{pair}{$p}{coords}{M}{$y} = {SEQ_start => $first_dis[($x*3)+2][($y*4)],
                                           SEQ_end => $first_dis[($x*3)+2][($y*4)+1],
                                           TE_start => $first_dis[($x*3)+2][($y*4)+2],
                                           TE_end => $first_dis[($x*3)+2][($y*4)+3]};
        }
      }
    }
  }
 elsif($first_dis eq 'FRAG-NLTR')
  {
  my $count = @first_head;
  my $z=0;
  my $w=0;
  for(my $x=0;$x<$count;$x++)
    {
    my $type;
    my $p;
    if($first_head[$x][0] =~ /^f/)
      {
      $type = 'frag';
      $p = "f".$w;
      $w++;
      }
     elsif($first_head[$x][0] =~ /^n/)
      {
      $type = 'nltr';
      $p = "n".$z;
      $z++;
      }
    $retro{$type}{$p} = {type => $first_head[$x][1], dir => $first_head[$x][2], bsr => $first_head[$x][3],
                        group => $first_head[$x][4], level => $first_head[$x][5], order => $first_head[$x][6]};
    my $power_len = (scalar @{$first_dis[$x]})/4;
    for(my $y=0;$y<$power_len;$y++)
      {
      $retro{$type}{$p}{coords}{$y} = {SEQ_start => $first_dis[$x][($y*4)],
                                       SEQ_end => $first_dis[$x][($y*4)+1],
                                       TE_start => $first_dis[$x][($y*4)+2],
                                       TE_end => $first_dis[$x][($y*4)+3]};
      }
    }
  }
if($_[1])
  {
  my $second_dis = $_[1];
  if($second_dis eq 'SOLO')
    {
    my $count = @second_head;
    for(my $x=0;$x<$count;$x++)
      {
      my $p = "s".$x;
      $retro{solo}{$p} = {type => $second_head[$x][1], dir => $second_head[$x][2], bsr => $second_head[$x][3],
                          group => $second_head[$x][4], level => $second_head[$x][5], order => $second_head[$x][6]};
      my $power_len = (scalar @{$second_dis[$x]})/4;
      for(my $y=0;$y<$power_len;$y++)
        {
        $retro{solo}{$p}{coords}{$y} = {SEQ_start => $second_dis[$x][($y*4)],
                                        SEQ_end => $second_dis[$x][($y*4)+1],
                                        TE_start => $second_dis[$x][($y*4)+2],
                                        TE_end => $second_dis[$x][($y*4)+3]};
        }
      }
    }
   elsif($second_dis eq 'PAIR')
    {
    my $count = @second_head;
    for(my $x=0;$x<$count;$x++)
      {
      my $p = "p".$x;
      $retro{pair}{$p} = {type => $second_head[$x][1], dir => $second_head[$x][2], bsr => $second_head[$x][3],
                          group => $second_head[$x][4], level => $second_head[$x][5], order => $second_head[$x][6]};
      my $power_len = (scalar @{$second_dis[$x*3]})/4;
      for(my $y=0;$y<$power_len;$y++)
        {
        $retro{pair}{$p}{coords}{L}{$y} = {SEQ_start => $second_dis[$x*3][($y*4)],
                                           SEQ_end => $second_dis[$x*3][($y*4)+1],
                                           TE_start => $second_dis[$x*3][($y*4)+2],
                                           TE_end => $second_dis[$x*3][($y*4)+3]};
        }
      $power_len = (scalar @{$second_dis[($x*3)+1]})/4;
      for(my $y=0;$y<$power_len;$y++)
        {
        $retro{pair}{$p}{coords}{R}{$y} = {SEQ_start => $second_dis[($x*3)+1][($y*4)],
                                           SEQ_end => $second_dis[($x*3)+1][($y*4)+1],
                                           TE_start => $second_dis[($x*3)+1][($y*4)+2],
                                           TE_end => $second_dis[($x*3)+1][($y*4)+3]};
        }
      if($second_dis[($x*3)+2][0] ne 'blank')
        {
        $power_len = (scalar @{$second_dis[($x*3)+2]})/4;
        for(my $y=0;$y<$power_len;$y++)
          {
          $retro{pair}{$p}{coords}{M}{$y} = {SEQ_start => $second_dis[($x*3)+2][($y*4)],
                                             SEQ_end => $second_dis[($x*3)+2][($y*4)+1],
                                             TE_start => $second_dis[($x*3)+2][($y*4)+2],
                                             TE_end => $second_dis[($x*3)+2][($y*4)+3]};
          }
        }
      }
    }
   elsif($second_dis eq 'FRAG-NLTR')
    {
    my $count = @second_head;
    my $z=0;
    my $w=0;
    for(my $x=0;$x<$count;$x++)
      {
      my $type;
      my $p;
      if($second_head[$x][0] =~ /^f/)
        {
        $type = 'frag';
        $p = "f".$w;
        $w++;
        }
       elsif($second_head[$x][0] =~ /^n/)
        {
        $type = 'nltr';
        $p = "n".$z;
        $z++;
        }
      $retro{$type}{$p} = {type => $second_head[$x][1], dir => $second_head[$x][2], bsr => $second_head[$x][3],
                          group => $second_head[$x][4], level => $second_head[$x][5], order => $second_head[$x][6]};
      my $power_len = (scalar @{$second_dis[$x]})/4;
      for(my $y=0;$y<$power_len;$y++)
        {
        $retro{$type}{$p}{coords}{$y} = {SEQ_start => $second_dis[$x][($y*4)],
                                        SEQ_end => $second_dis[$x][($y*4)+1],
                                        TE_start => $second_dis[$x][($y*4)+2],
                                        TE_end => $second_dis[$x][($y*4)+3]};
        }
      }
    }
  }
}
###########################################################################

sub OVERLAP_S_FN
{
my $first_dis = $_[0];
my $second_dis = $_[1];
LOAD_DISCREPANCY($first_dis,$second_dis);
my @temp_dis_1 = ();
my @temp_dis_2 = ();
my @temp_head_1 = ();
my @temp_head_2 = ();
for(my $w=0;$w<@first_dis;$w++)
  {
  @temp_dis_2 = ();
  @temp_head_2 = ();
  my @remove_list_1 = ();
  for(my $x=0;$x<@second_dis;$x++)
    {
    my @remove_list_2 = ();
    for(my $y=0;$y<(scalar @{$first_dis[$w]})/4;$y++)
      {
      for(my $z=0;$z<(scalar @{$second_dis[$x]})/4;$z++)
        {
        if($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)+1])
          {
          my $remove = "$w $y";
          push @remove_list_1, [split(/ /,$remove)];
          }
         elsif($first_dis[$w][$y*4] <= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][($z*4)+1])
          {
          my $remove = "$x $z";
          push @remove_list_2, [split(/ /,$remove)];
          }
         elsif($first_dis[$w][$y*4] <= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)+1] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][$z*4])
          {
          if(($first_dis[$w][($y*4)+1] - $first_dis[$w][$y*4]) <= ($second_dis[$x][($z*4)+1] - $second_dis[$x][$z*4]))
            {$first_dis[$w][($y*4)+1] = $second_dis[$x][$z*4]-1;
}
           else
            {$second_dis[$x][$z*4] = $first_dis[$w][($y*4)+1]+1;
}
          }
         elsif($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][($z*4)+1] && $first_dis[$w][$y*4] <= $second_dis[$x][($z*4)+1])
          {
          if(($first_dis[$w][($y*4)+1] - $first_dis[$w][$y*4]) <= ($second_dis[$x][($z*4)+1] - $second_dis[$x][$z*4]))
            {$first_dis[$w][$y*4] = $second_dis[$x][($z*4)+1]+1;
}
           else
            {$second_dis[$x][($z*4)+1] = $first_dis[$w][$y*4]-1;
}
          }
        }
      }
    my $found = 0;
    my $found1 = 1;
    my $temp_dis = '';
    for(my $z=0;$z<(scalar @{$second_dis[$x]})/4;$z++)
      {
      for(my $a=0;$a<@remove_list_2;$a++)
        {
        if($remove_list_2[$a][1] eq $z)
          {
          $found1 = 0;
          }
        }
      if($found1 == 1 || @remove_list_2 == 0)
        {
        $temp_dis = $temp_dis . " $second_dis[$x][$z*4] $second_dis[$x][($z*4)+1] $second_dis[$x][($z*4)+2] $second_dis[$x][($z*4)+3]";
        $found = 1;
        }
      }
    if($found == 1)
      {
      $temp_dis =~ s/^\s//g;
      push @temp_dis_2, [split(/ /,$temp_dis)];
      my $temp_head = "$second_head[$x][0] $second_head[$x][1] $second_head[$x][2] $second_head[$x][3] $second_head[$x][4] $second_head[$x][5] $second_head[$x][6]";
      push @temp_head_2, [split(/ /,$temp_head)];
      }
    }
  @second_dis = @temp_dis_2;
  @second_head = @temp_head_2;
  my $found = 0;
  my $found1 = 1;
  my $temp_dis = '';
  for(my $y=0;$y<(scalar @{$first_dis[$w]})/4;$y++)
    {
    for(my $z=0;$z<@remove_list_1;$z++)
      {
      if($remove_list_1[$z][1] eq $y)
        {
        $found1 = 0;
        }
      }
    if($found1 == 1 || @remove_list_1 == 0)
      {      
      $temp_dis = $temp_dis . " $first_dis[$w][$y*4] $first_dis[$w][($y*4)+1] $first_dis[$w][($y*4)+2] $first_dis[$w][($y*4)+3]";
      $found = 1;
      }
    }
  if($found == 1)
    {
    $temp_dis =~ s/^\s//g;
    push @temp_dis_1, [split(/ /,$temp_dis)];
    my $temp_head = "$first_head[$w][0] $first_head[$w][1] $first_head[$w][2] $first_head[$w][3] $first_head[$w][4] $first_head[$w][5] $first_head[$w][6]";
    push @temp_head_1, [split(/ /,$temp_head)];
    }
  }
@first_head = @temp_head_1;
@first_dis = @temp_dis_1;
UNLOAD_DISCREPANCY($first_dis,$second_dis);
}
###########################################################################

sub OVERLAP_S_M
{
my $first_dis = $_[0];
my $second_dis = $_[1];
LOAD_DISCREPANCY($first_dis,$second_dis);
my @temp_dis_1 = ();
my @temp_head_1 = ();
for(my $w=0;$w<@first_dis;$w++)
  {
  my @remove_list_1 = ();
  for(my $x=2;$x<@second_dis;$x=$x+3)
    {
    if($second_dis[$x][0] ne 'blank')
      {
      for(my $y=0;$y<(scalar @{$first_dis[$w]})/4;$y++)
        {
        for(my $z=0;$z<(scalar @{$second_dis[$x]})/4;$z++)
          {
          if($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)+1])
            {
            my $remove = "$w $y";
            push @remove_list_1, [split(/ /,$remove)];
            }
           elsif($first_dis[$w][$y*4] <= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)+1] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][$z*4])
            {
            if(($first_dis[$w][($y*4)+1] - $first_dis[$w][$y*4]) <= ($second_dis[$x][($z*4)+1] - $second_dis[$x][$z*4]))
              {
              $first_dis[$w][($y*4)+1] = $second_dis[$x][$z*4]-1;
              }
             else
              {
              $second_dis[$x][$z*4] = $first_dis[$w][($y*4)+1]+1;
              }
            }
           elsif($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][($z*4)+1] && $first_dis[$w][$y*4] <= $second_dis[$x][($z*4)+1])
            {
            if(($first_dis[$w][($y*4)+1] - $first_dis[$w][$y*4]) <= ($second_dis[$x][($z*4)+1] - $second_dis[$x][$z*4]))
              {
              $first_dis[$w][$y*4] = $second_dis[$x][($z*4)+1]+1;
              }
             else
              {
              $second_dis[$x][($z*4)+1] = $first_dis[$w][$y*4]-1;
              }
            }
          }
        }
      }
    }
  my $found = 0;
  my $found1 = 1;
  my $temp_dis = '';
  for(my $y=0;$y<(scalar @{$first_dis[$w]})/4;$y++)
    {
    for(my $z=0;$z<@remove_list_1;$z++)
      {
      if($remove_list_1[$z][1] eq $y)
        {
        $found1 = 0;
        }
      }
    if($found1 == 1 || @remove_list_1 == 0)
      {
      $temp_dis = $temp_dis . " $first_dis[$w][$y*4] $first_dis[$w][($y*4)+1] $first_dis[$w][($y*4)+2] $first_dis[$w][($y*4)+3]";
      $found = 1;
      }
    }
  if($found == 1)
    {
    $temp_dis =~ s/^\s//g;
    push @temp_dis_1, [split(/ /,$temp_dis)];
    my $temp_head = "$first_head[$w][0] $first_head[$w][1] $first_head[$w][2] $first_head[$w][3] $first_head[$w][4] $first_head[$w][5] $first_head[$w][6]";
    push @temp_head_1, [split(/ /,$temp_head)];
    }
  }
@first_head = @temp_head_1;
@first_dis = @temp_dis_1;
UNLOAD_DISCREPANCY($first_dis,$second_dis);
}
###########################################################################

sub OVERLAP_P_M
{
my $first_dis = $_[0];
my $second_dis = $_[1];
LOAD_DISCREPANCY($first_dis,$second_dis);
my @temp_dis_1 = ();
my @temp_head_1 = ();
for(my $w=0;$w<@first_dis;$w=$w+3)
  {
  my @remove_list_1 = ();
  for(my $x=2;$x<@second_dis;$x=$x+3)
    {
    if($second_dis[$x][0] ne 'blank')
      {
      my $r_end = scalar @{$first_dis[$w+1]};
      for(my $z=0;$z<(scalar @{$second_dis[$x]})/4;$z++)
        {
        if($first_dis[$w][0] >= $second_dis[$x][$z*4] && $first_dis[$w+1][$r_end-3] <= $second_dis[$x][($z*4)+1])
          {
          for(my $y=0;$y<(scalar @{$first_dis[$w]})/4;$y++)
            {
            my $remove = "$w $y";
            push @remove_list_1, [split(/ /,$remove)];
            }
          for(my $y=0;$y<(scalar @{$first_dis[$w+1]})/4;$y++)
            {
            my $remove = "$w $y";
            push @remove_list_1, [split(/ /,$remove)];
            }
          for(my $y=0;$y<(scalar @{$first_dis[$w+2]})/4;$y++)
            {
            my $remove = "$w $y";
            push @remove_list_1, [split(/ /,$remove)];
            }
          }
        }
      }
    }
  my $found = 0;
  my $found1 = 1;
  my $temp_dis = '';
  for(my $x=0;$x<3;$x++)
    {
    $temp_dis = '';
    my $w_c = $w + $x;
    if($first_dis[$w_c][0] ne 'blank')
      {
      for(my $y=0;$y<(scalar @{$first_dis[$w_c]})/4;$y++)
        {
        for(my $z=0;$z<@remove_list_1;$z++)
          {
          if($remove_list_1[$z][1] eq $y)
            {
            $found1 = 0;
            }
          }
        if($found1 == 1 || @remove_list_1 == 0)
          {
          $temp_dis = $temp_dis . " $first_dis[$w_c][$y*4] $first_dis[$w_c][($y*4)+1] $first_dis[$w_c][($y*4)+2] $first_dis[$w_c][($y*4)+3]";
          $found = 1;
          }
        }
      }
    else
     {$temp_dis = 'blank';}
    $temp_dis =~ s/^\s//g;
    push @temp_dis_1, [split(/ /,$temp_dis)];
    }
  if($found == 1)
    {
    $temp_dis =~ s/^\s//g;
    my $temp_head = "$first_head[$w/3][0] $first_head[$w/3][1] $first_head[$w/3][2] $first_head[$w/3][3] $first_head[$w/3][4] $first_head[$w/3][5] $first_head[$w/3][6]";
    push @temp_head_1, [split(/ /,$temp_head)];
    }
  }
@first_head = @temp_head_1;
@first_dis = @temp_dis_1;
UNLOAD_DISCREPANCY($first_dis);
}
###########################################################################

sub DISCREPANCY
{
my $first_dis = $_[0];
my $second_dis = $_[1];
LOAD_DISCREPANCY($first_dis,$second_dis);
my $do_delete = 0;
my @temp_dis_1 = ();
my @temp_head_1 = ();
do{
  $do_delete = 0;
  my @coord_prob = ();
  my $first_dis_count = @first_dis;
  if($first_dis eq 'PAIR')
    {$first_dis_count = $first_dis_count / 3;}
  for(my $w=0;$w<$first_dis_count;$w++)
    {
    my $first_last;
    my $w_c = $w;
    my $set1_size = 0;
    if($first_dis eq 'PAIR')
      {
      $w_c = $w*3;
      $first_last = $first_dis[$w_c+1][(scalar @{$first_dis[$w_c+1]}) - 3];
      for(my $a=0;$a<scalar @{$first_dis[$w_c]};$a=$a+4)
        {$set1_size = $set1_size + ($first_dis[$w_c][$a+1] - $first_dis[$w_c][$a]);}
      for(my $a=0;$a<scalar @{$first_dis[$w_c+1]};$a=$a+4)
        {$set1_size = $set1_size + ($first_dis[$w_c+1][$a+1] - $first_dis[$w_c+1][$a]);}
      if($first_dis[$w_c+2][0] ne 'blank')
        {
        for(my $a=0;$a<scalar @{$first_dis[$w_c+2]};$a=$a+4)
          {$set1_size = $set1_size + ($first_dis[$w_c+2][$a+1] - $first_dis[$w_c+2][$a]);}
        }
      }
     else
      {
      $first_last = $first_dis[$w][(scalar @{$first_dis[$w]}) - 3];
      for(my $a=0;$a<scalar @{$first_dis[$w_c]};$a=$a+4)
        {$set1_size = $set1_size + ($first_dis[$w_c][$a+1] - $first_dis[$w_c][$a]);}
      }
    my $second_dis_count = @second_dis;
    if($second_dis eq 'PAIR')
      {$second_dis_count = $second_dis_count / 3;}
    for(my $x=0;$x<$second_dis_count;$x++)
      {
      my $second_last;
      my $x_c = $x;
      my $set2_size = 0;
      if($second_dis eq 'PAIR')
        {
        $x_c = $x*3;
        $second_last = $second_dis[$x_c+1][(scalar @{$second_dis[$x_c+1]}) - 3];
        for(my $a=0;$a<scalar @{$second_dis[$x_c]};$a=$a+4)
          {$set2_size = $set2_size + ($second_dis[$x_c][$a+1] - $second_dis[$x_c][$a]);}
        for(my $a=0;$a<scalar @{$second_dis[$x_c+1]};$a=$a+4)
          {$set2_size = $set2_size + ($second_dis[$x_c+1][$a+1] - $second_dis[$x_c+1][$a]);}
        if($second_dis[$x_c+2][0] ne 'blank')
          {
          for(my $a=0;$a<scalar @{$second_dis[$x_c+2]};$a=$a+4)
            {$set2_size = $set2_size + ($second_dis[$x_c+2][$a+1] - $second_dis[$x_c+2][$a]);}
          }
        }
       else
        {
        $second_last = $second_dis[$x][(scalar @{$second_dis[$x]}) - 3];
        for(my $a=0;$a<scalar @{$second_dis[$x_c]};$a=$a+4)
          {$set2_size = $set2_size + ($second_dis[$x_c][$a+1] - $second_dis[$x_c][$a]);}
        }
      if($first_dis[$w_c][0] < $second_dis[$x_c][0] && $first_last > $second_dis[$x_c][0] && $first_last < $second_last)
        {
        my $coord_prob_count = @coord_prob;
        my $coord_found = 0;
        if($coord_prob_count == 0)
          {
          my $coord_set = "$w 1 $set1_size";
          push @coord_prob, [split(/ /,$coord_set)];
          $coord_found = 1;
          }
        for(my $y=0;$y<$coord_prob_count;$y++)
          {
          if($w == $coord_prob[$y][0] && $coord_found == 0)
            {
            $coord_prob[$y][1]++;
            $coord_found = 1;
            }
          }
        if($coord_found == 0)
          {
          my $coord_set = "$w 1 $set1_size";
          push @coord_prob, [split(/ /,$coord_set)];
          $coord_found = 1;
          }
        }
       elsif($first_dis[$w_c][0] > $second_dis[$x_c][0] && $first_dis[$w_c][0] < $second_last && $first_last > $second_last)
        {
        my $coord_prob_count = @coord_prob;
        my $coord_found = 0;
        if($coord_prob_count == 0)
          {
          my $coord_set = "$w 1 $set2_size";
          push @coord_prob, [split(/ /,$coord_set)];
          $coord_found = 1;
          }
        for(my $y=0;$y<$coord_prob_count;$y++)
          {
          if($w == $coord_prob[$y][0] && $coord_found == 0)
            {
            $coord_prob[$y][1]++;
            $coord_found = 1;
            }
          }
        if($coord_found == 0)
          {
          my $coord_set = "$w 1 $set2_size";
          push @coord_prob, [split(/ /,$coord_set)];
          $coord_found = 1;
          }
        }
      }
    }
  my $coord_prob_count = @coord_prob;
  if($coord_prob_count > 0)
    {
    @temp_dis_1 = ();
    @temp_head_1 = ();
    $do_delete = 1;
    @coord_prob = sort {$a->[2] <=> $b->[2]} @coord_prob;
    @coord_prob = sort {$b->[1] <=> $a->[1]} @coord_prob;
    for(my $x=0;$x<$coord_prob[0][0];$x++)
      {
      if($first_dis eq 'PAIR')
        {
        for(my $z=0;$z<3;$z++)
          {
          my $temp_dis = $first_dis[($x*3)+$z][0];
          for(my $y=1;$y<(scalar @{$first_dis[($x*3)+$z]});$y++)
            {
            $temp_dis = $temp_dis . " " . $first_dis[($x*3)+$z][$y];
            }
          push @temp_dis_1, [split(/ /,$temp_dis)];
          }
        }
       else
        {
        my $temp_dis = $first_dis[$x][0];
        for(my $y=1;$y<(scalar @{$first_dis[$x]});$y++)
          {
          $temp_dis = $temp_dis . " " . $first_dis[$x][$y];
          }
        push @temp_dis_1, [split(/ /,$temp_dis)];
        }
      my $temp_head = "$first_head[$x][0] $first_head[$x][1] $first_head[$x][2] $first_head[$x][3] $first_head[$x][4] $first_head[$x][5] $first_head[$x][6]";
      push @temp_head_1, [split(/ /,$temp_head)];
      }
    for(my $x=$coord_prob[0][0]+1;$x<$first_dis_count;$x++)
      {
      if($first_dis eq 'PAIR')
        {
        for(my $z=0;$z<3;$z++)
          {
          my $temp_dis = $first_dis[($x*3)+$z][0];
          for(my $y=1;$y<(scalar @{$first_dis[($x*3)+$z]});$y++)
            {
            $temp_dis = $temp_dis . " " . $first_dis[($x*3)+$z][$y];
            }
          push @temp_dis_1, [split(/ /,$temp_dis)];
          }
        }
       else
        {
        my $temp_dis = $first_dis[$x][0];
        for(my $y=1;$y<(scalar @{$first_dis[$x]});$y++)
          {
          $temp_dis = $temp_dis . " " . $first_dis[$x][$y];
          }
        push @temp_dis_1, [split(/ /,$temp_dis)];
        }
      my $temp_head = "$first_head[$x][0] $first_head[$x][1] $first_head[$x][2] $first_head[$x][3] $first_head[$x][4] $first_head[$x][5] $first_head[$x][6]";
      push @temp_head_1, [split(/ /,$temp_head)];
      }
    my $seperate_size = scalar @{$first_dis[$coord_prob[0][0]]};
    if($first_dis ne 'PAIR' && $first_dis ne 'SOLO' && $seperate_size > 4)
      {
      for(my $a=0;$a<$seperate_size;$a=$a+4)
        {
        my $temp_dis = "$first_dis[$coord_prob[0][0]][$a] $first_dis[$coord_prob[0][0]][$a+1] $first_dis[$coord_prob[0][0]][$a+2] $first_dis[$coord_prob[0][0]][$a+3]";
        push @temp_dis_1, [split(/ /,$temp_dis)];
        my $temp_head = "$first_head[$coord_prob[0][0]][0] $first_head[$coord_prob[0][0]][1] $first_head[$coord_prob[0][0]][2] $first_head[$coord_prob[0][0]][3] $first_head[$coord_prob[0][0]][4] $first_head[$coord_prob[0][0]][5] $first_head[$coord_prob[0][0]][6]";
        push @temp_head_1, [split(/ /,$temp_head)];
        }
      }
    @first_head = @temp_head_1;
    @first_dis = @temp_dis_1;
    if($first_dis eq $second_dis)
      {
      @second_head = @temp_head_1;
      @second_dis = @temp_dis_1;
      }
    }
  }until($do_delete == 0);
if($first_dis ne $second_dis)
  {
  UNLOAD_DISCREPANCY($first_dis,$second_dis);
  }
 elsif($first_dis eq $second_dis)
  {
  UNLOAD_DISCREPANCY($first_dis);
  }
}

