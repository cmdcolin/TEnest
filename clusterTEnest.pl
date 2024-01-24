#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);

#########################################################################################
# A submission script for TEnest for use on clustered computers.  This version uses the qsub
# submission, clustered computers using different submissions will need slight alterations.
#
# This script cuts up your input sequence and runs TEnest on each section on a different node
# concurrently (assuming there are enough nodes).  The script will loop until all sections have
# run on a node.  The input sequence is cut into equal length segments if the length is greater
# than the REQUIRED option value of '--node_length'.  '--node_length' rather than number of
# nodes is used to allow more sequence cuts than number of nodes available. I.E. if your full 
# input sequence length is 1002Kb, and you set the --node_length to 200Kb, you will have 6 
# sequence sections each with length of 167Kb.  (It's set up this way to not get a small last
# cut - in the previous example, if we used 5 sections of 200Kb we'd have a last section of
# 2Kb and TEnest would get annoyed with you).  Once each node is complete, the sequences and 
# annotations are rejoined, the non-annotated sequence length is counted and if greater than
# the '--node_length' another loop of cutting the input sequence and sending it to nodes is
# run, if it is less than the '--node_length' a final round of TEnest is run on the whole
# sequence.  If --node_length calculates to a number of sequences more than the available
# number of nodes the remaining sequences will be queued until nodes open.
#
# In computational time this will take longer than running TEnest.pl, but when spread out to
# cluster nodes the actual time until completion is much quicker.
#
# Options used here will be copied over to the TEnest options, so use the command line options
# just as you would when running TEnest.pl, except you must input a '--node_length' value.
#
################################################################################################
#                                              MAIN                                            #
################################################################################################

#PRE-PROCESS OPTIONS
my @input = @ARGV;
my $node_remove;
my $input_sequence = $ARGV[0];
for(my $x=1;$x<@ARGV;$x++)
  {
  if($ARGV[$x] eq '--node_length')
    {$node_remove = $x;}
  }
if(!$ARGV[$node_remove+1])
  {print "You must enter a cut length value (--node_length X)\n";die;}
my $node_length = $ARGV[$node_remove+1];
splice(@ARGV,$node_remove,2);
my $options = '';
for(my $x=1;$x<@ARGV;$x++)
  {$options = $options . " " . $ARGV[$x];}
$options =~ s/^\s//g;

#COUNT INPUT SEQUENCE
my $date = DATE();
my $input_dir = $input_sequence . "_" . $date;
print "Beginning run, results will be found in $input_dir\n";
mkdir("$input_dir", 0755);
system("cp $input_sequence $input_dir");
open(INPUT, "$input_dir/$input_sequence");
my $seq = '';
while(<INPUT>)
  {
  my $line = $_;
  chomp($line);
  if($line !~ m/>/)
    {$seq = $seq . $line;}
  }
close(INPUT);
my @seq = split(//,$seq);
my $n_count = 0;
my @seq_n = ();
for(my $x=0;$x<@seq;$x++)
  {
  $seq_n[$x][0] = $seq[$x];
  $seq_n[$x][1] = 0;
  if($seq_n[$x][0] ne 'N' && $seq_n[$x][0] ne 'n')
    {
    $seq_n[$x][1] = $n_count;
    $n_count++;
    }
  }
my $seq_count = @seq;

#START PAIR QSUB INPUT LOOP
my $wrote_ltr = 1;
my $round = 0;
my @ltr = ();
my $pair_count = 0;
my $p = 'p0';
my $stop = 'continue';
while($n_count > $node_length && $wrote_ltr == 1 && $stop eq 'continue')
  {
#SPLIT INPUT SEQUENCE
print "@ the start... $n_count > $node_length && $wrote_ltr == 1\n";
  my $wrote_ltr = 0;
  my $start = 0;
  my $cut = 0;
  my $cut_size = ceil($n_count/ceil($n_count/$node_length));
  my $seq_length = $cut_size;
  my $real_cut_size = 0;
  my $cut_number = ceil($n_count/$node_length)-1;
  my $last_size = $n_count - ($cut_number * $cut_size);
  print "Round $round, $n_count non-N bases\n";
  print "  $cut_number sequence(s) of $cut_size bp, 1 with $last_size bp\n";
  for(my $y=0;$y<ceil($n_count/$node_length)-1;$y++)
    {
    for(my $x=0;$x<@seq_n;$x++)
      {
      if($seq_n[$x][1] == $seq_length)
        {
        $real_cut_size = $x;
        my $cut_seq = $input_sequence . "-" . $round . "-" . $cut;
        open(CUT,">$input_dir/$cut_seq");
        print CUT ">$cut_seq\n";
        my $seventy_count = 0;
        for(my $w=$start;$w<$x;$w++)
          {
          print CUT "$seq_n[$w][0]";
          $seventy_count++;
          if($seventy_count == 70)
            {
            print CUT "\n";
            $seventy_count = 0;
            }
          }
        print CUT "\n";
        close(CUT);
        $cut++;
        $start = $x;
        }
      }
    $seq_length = $seq_length + $cut_size;
    }
  my $cut_seq = $input_sequence . "-" . $round . "-" . $cut;
  open(CUT,">$input_dir/$cut_seq");
  print CUT ">$cut_seq\n";
  my $seventy_count = 0;
  for(my $w=$start;$w<@seq;$w++)
    {
    print CUT "$seq_n[$w][0]";
    $seventy_count++;
    if($seventy_count == 70)
      {
      print CUT "\n";
      $seventy_count = 0;
      }
    }
  print CUT "\n";
  $cut++;
  close(CUT);
print "wrote $wrote_ltr, $n_count\n";
#SEND TO QSUB WAIT FOR SEND CONFORMATION
  my @qsub_sent = ();
  my $script = 'q-script';
  system("touch $input_dir/$script");
  for(my $y=0;$y<$cut;$y++)
    {
    system("rm $input_dir/$script");
    open(SCRIPT, ">$input_dir/$script");
    print SCRIPT "./TE_nest.pl $input_dir/$input_sequence-$round-$y $options --output $input_dir --pair_only T";
    close(SCRIPT);
    my $out;
    do{
      $out = system("qsub -d . -k eo $input_dir/$script > $input_dir/q-out");
      sleep(2);
      open(OUT,"$input_dir/q-out");
      my $qout = <OUT>;
      chomp($qout);
      close(OUT);
#      system("rm $input_dir/q-out");
      push @qsub_sent, $qout;
      }until($out == 0);
    }
#  system("rm $input_dir/q-script");
#CHECK FOR COMPLETED FILES - check if file is off qstat, check if *.LTR exists, die if not found
  print "Waiting for PAIR node runs to complete";
  for(my $y=0;$y<$cut;$y++)
    {
    my $done = 0;
    do{
      my $grep = system("qstat | grep $qsub_sent[$y] > /dev/null");
      # $grep returns 0 if found, 256 if not found
      if($grep == 0)
        {
        sleep 5;
        print ".";
        }  #RAISE THIS
       else
        {$done = 1;}
      }until($done == 1);
    my $check_number = 0;
    my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-$round-$y > $input_dir/q-ls");
    open(OUT,"$input_dir/q-ls");
    my $qls = <OUT>;
    chomp($qls);
    close(OUT);
    system("rm $input_dir/q-ls");
    while(!(-e "$qls/$input_sequence-$round-$y.LTR") && $check_number < 3)
      {
      sleep 10;  #RAISE THIS
      $check_number++;
      }
    if(!(-e "$qls/$input_sequence-$round-$y.LTR"))
      {die "Round $round Cut $y has finished running, but *.LTR doesn't exist.  Exiting...\n";}
    }
  print "\n";
#REJOIN SEQUENCES AND *.LTR ANNOTATIONS
  @seq = ();
  my $seq = '';
print "wrote $wrote_ltr, $n_count\n";
  for(my $y=0;$y<$cut;$y++)
    {
    my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-$round-$y > $input_dir/q-ls");
    open(OUT,"$input_dir/q-ls");
    my $qls = <OUT>;
    chomp($qls);
    close(OUT);
    system("rm $input_dir/q-ls");
    open(SEQ,"$qls/$input_sequence-$round-$y.mask") or die "$input_sequence-$round-$y.mask not found, TE_nest on round $round cut $cut failed.\n";
    while(<SEQ>)
      {
      my $line = $_;
      chomp($line);
      if($line !~ m/>/)
        {$seq = $seq . $line;}
      }
    }
  close(SEQ);
  @seq = split(//,$seq);
  $n_count = 0;
  @seq_n = ();
print "wrote $wrote_ltr, $n_count\n";
  for(my $x=0;$x<@seq;$x++)
    {
    $seq_n[$x][0] = $seq[$x];
    $seq_n[$x][1] = 0;
    if($seq_n[$x][0] ne 'N' && $seq_n[$x][0] ne 'n')
      {
      $seq_n[$x][1] = $n_count;
      $n_count++;
      }
    }
print "wrote $wrote_ltr, $n_count\n";
  for(my $y=0;$y<$cut;$y++)
    {
    my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-$round-$y > $input_dir/q-ls");
    open(OUT,"$input_dir/q-ls");
    my $qls = <OUT>;
    chomp($qls);
    close(OUT);
    system("rm $input_dir/q-ls");
    open(LTR,"$qls/$input_sequence-$round-$y.LTR");
    while(<LTR>)
      {
      my $line = $_;
      chomp($line);
      my @ltr_tmp = split(/ /,$line);
      my $increase = $real_cut_size * $y;
      if($ltr_tmp[0] eq 'PAIR')
        {
        $p = 'p'.$pair_count;
        $ltr_tmp[1] = $p;
        $wrote_ltr = 1;
        my $ltr = "$ltr_tmp[0] $ltr_tmp[1] $ltr_tmp[2] $ltr_tmp[3] $ltr_tmp[4] $ltr_tmp[5] $ltr_tmp[6] $ltr_tmp[7]";
print "FOUND A PAIR WRITING writing $y $ltr -  wrote_ltr is $wrote_ltr\n";
        push @ltr, [split(/ /,$ltr)];
        $pair_count++;
        }
      if($ltr_tmp[0] =~ m/^p/)
        {
        $ltr_tmp[0] = $p;
        my $ltr = "$ltr_tmp[0] $ltr_tmp[1]";
        for(my $y=0;$y<(@ltr_tmp-2)/4;$y++)
          {
          $ltr_tmp[($y*4)+2] = $ltr_tmp[($y*4)+2] + $increase;
          $ltr_tmp[($y*4)+3] = $ltr_tmp[($y*4)+3] + $increase;
          $ltr = $ltr . " $ltr_tmp[($y*4)+2] $ltr_tmp[($y*4)+3] $ltr_tmp[($y*4)+4] $ltr_tmp[($y*4)+5]";
          }
        push @ltr, [split(/ /,$ltr)];
        }
      }
    close(LTR);
    }
  print "PAIRs found: $pair_count\n";
  $round++;
print "THE LAST - wrote $wrote_ltr, $n_count\n";
if($wrote_ltr == 0)
  {
  print "IT SHOULD STOP NOW\n";
  $stop = 'stop';
  }
print "What does the stop say? $stop\n";
  }
#  }until($n_count <= $node_length && $wrote_ltr == 0 && $stop eq 'stop');


#FINAL TEnest PAIR SUBMISSION
my $out_seq = $input_sequence . "-" . $round;
open(SEQ,">$input_dir/$out_seq");
print SEQ ">$out_seq\n";
my $seventy_count = 0;
for(my $w=0;$w<@seq_n;$w++)
  {
  print SEQ "$seq_n[$w][0]";
  $seventy_count++;
  if($seventy_count == 70)
    {
    print SEQ "\n";
    $seventy_count = 0;
    }
  }
print SEQ "\n";
close(SEQ);
#SEND LAST SEQUENCE TO QSUB TEnest
my @qsub_sent = ();
my $script = 'q-script';
open(SCRIPT, ">$input_dir/$script");
print SCRIPT "./TE_nest.pl $input_dir/$input_sequence-$round $options --output $input_dir --pair_only S";
close(SCRIPT);
my $out;
do{
  $out = system("qsub -d . -k eo $input_dir/$script > $input_dir/q-out");
  sleep(2);
  open(OUT,"$input_dir/q-out");
  my $qout = <OUT>;
  chomp($qout);
  system("rm $input_dir/q-out");
  push @qsub_sent, $qout;
  close(OUT);
  }until($out == 0);
system("rm $input_dir/q-script");
#CHECK FOR COMPLETED FILES - check if file is off qstat, check if *.LTR exists, die if not found
print "Waiting for final PAIR node run to complete";
my $done = 0;
do{
  my $grep = system("qstat | grep $qsub_sent[0] > /dev/null");
  # $grep returns 0 if found, 256 if not found
  if($grep == 0)
    {
    sleep 5;
    print ".";
    }  #RAISE THIS
   else
    {$done = 1;}
  }until($done == 1);
my $check_number = 0;
my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-$round > $input_dir/q-ls");
open(OUT,"$input_dir/q-ls");
my $qls = <OUT>;
chomp($qls);
close(OUT);
system("rm $input_dir/q-ls");
while(!(-e "$qls/$input_sequence-$round.LTR") && $check_number < 3)
  {
  sleep 10;  #RAISE THIS
  $check_number++;
  }
if(!(-e "$qls/$input_sequence-$round.LTR"))
  {die "Round $round has finished running, but *.LTR doesn't exist.  Exiting...\n";}
print "\n";
#FINAL *.LTR COMBO WRITE SOLO AND PAIR TO @ltr_final
open(LTR,"$qls/$input_sequence-$round.LTR");
my @ltr_final = ();
my $line = <LTR>;
$line = <LTR>;
while(<LTR>)
  {
  $line = $_;
  chomp($line);
  my @ltr_tmp = split(/ /,$line);
  if($ltr_tmp[0] eq 'PAIR')
    {
    $p = 'p'.$pair_count;
    $ltr_tmp[1] = $p;
    $wrote_ltr = 1;
    my $ltr = "$ltr_tmp[0] $ltr_tmp[1] $ltr_tmp[2] $ltr_tmp[3] $ltr_tmp[4] $ltr_tmp[5] $ltr_tmp[6] $ltr_tmp[7]";
    push @ltr_final, [split(/ /,$ltr)];
    $pair_count++;
    }
  elsif($ltr_tmp[0] =~ m/^p/)
    {
    $ltr_tmp[0] = $p;
    my $ltr = "$ltr_tmp[0] $ltr_tmp[1]";
    my $ltr_tmp_count = @ltr_tmp;
    for(my $y=0;$y<($ltr_tmp_count-2)/4;$y++)
      {$ltr = $ltr . " $ltr_tmp[($y*4)+2] $ltr_tmp[($y*4)+3] $ltr_tmp[($y*4)+4] $ltr_tmp[($y*4)+5]";}
    push @ltr_final, [split(/ /,$ltr)];
    }
  else
    {
    my $ltr = $ltr_tmp[0];
    for(my $y=1;$y<@ltr_tmp;$y++)
      {$ltr = $ltr . " $ltr_tmp[$y]";}
    push @ltr_final, [split(/ /,$ltr)];
    }
  }
close(LTR);
for(my $x=0;$x<@ltr;$x++)
  {
  my $ltr = $ltr[$x][0];
  for(my $y=1;$y<(scalar @{$ltr[$x]});$y++)
    {$ltr = $ltr . " $ltr[$x][$y]";}
  push @ltr_final, [split(/ /,$ltr)];
  }
#for(my $x=0;$x<@ltr_final;$x++)
#  {
#  for(my $y=0;$y<(scalar @{$ltr_final[$x]});$y++)
#    {
#    print "$ltr_final[$x][$y] ";
#    }
#  print "\n";
#  }

#LOAD SOLO-PAIR ARRAY TO HASH
our %retro = ();
LOAD_HASH(@ltr_final);

my @first_dis = ();
my @second_dis = ();
my @first_head = ();
my @second_head = ();

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
#MAKE SURE THERE ARE NO PAIRS IN MID LOCATIONS, IF SO REMOVE THE PAIR
OVERLAP_S_P('SOLO','PAIR');

#RELOAD INPUT SEQUENCE TO CORRECTLY N OUT PAIRS
open(INPUT, "$input_dir/$input_sequence");
$seq = '';
while(<INPUT>)
  {
  my $line = $_;
  chomp($line);
  if($line !~ m/>/)
    {$seq = $seq . $line;}
  }
close(INPUT);
@seq = split(//,$seq);
my $seq_count_11 = @seq;
print "first seq_count $seq_count_11\n";

my $seq_ref = PRINT_LTR('F',\@ltr_final,\@seq);
@seq = @$seq_ref;
$n_count = 0;
@seq_n = ();
my $seq_count_11 = @seq;
print "second seq_count $seq_count_11\n";
for(my $x=0;$x<@seq;$x++)
  {
  $seq_n[$x][0] = $seq[$x];
  $seq_n[$x][1] = 0;
print "$x$seq_n[$x][0]";
  if($seq_n[$x][0] ne 'N' && $seq_n[$x][0] ne 'n')
    {
    $seq_n[$x][1] = $n_count;
    $n_count++;
    }
  }
print "\n";

#SPLIT INPUT SEQUENCE FOR FRAG RUN
my $start = 0;
my $cut = 0;
my $cut_size = ceil($n_count/ceil($n_count/$node_length));
my $seq_length = $cut_size;
my $real_cut_size = 0;
my $cut_number = ceil($n_count/$node_length)-1;
my $last_size = $n_count - ($cut_number * $cut_size);
print "Fragment submission, $n_count non-N bases\n";
print "  $cut_number sequence(s) of $cut_size bp, 1 with $last_size bp\n";
for(my $y=0;$y<ceil($n_count/$node_length)-1;$y++)
  {
  for(my $x=0;$x<@seq_n;$x++)
    {
    if($seq_n[$x][1] == $seq_length)
      {
      $real_cut_size = $x;
      my $cut_seq = $input_sequence . "-f-" . $cut;
      open(CUT,">$input_dir/$cut_seq");
      print CUT ">$cut_seq\n";
      $seventy_count = 0;
      for(my $w=$start;$w<$x;$w++)
        {
        print CUT "$seq_n[$w][0]";
        $seventy_count++;
        if($seventy_count == 70)
          {
          print CUT "\n";
          $seventy_count = 0;
          }
        }
      print CUT "\n";
      close(CUT);
      $cut++;
      $start = $x;
      }
    }
  $seq_length = $seq_length + $cut_size;
  }
my $cut_seq = $input_sequence . "-f-" . $cut;
open(CUT,">$input_dir/$cut_seq");
print CUT ">$cut_seq\n";
$seventy_count = 0;
for(my $w=$start;$w<@seq;$w++)
  {
  print CUT "$seq_n[$w][0]";
  $seventy_count++;
  if($seventy_count == 70)
    {
    print CUT "\n";
    $seventy_count = 0;
    }
  }
print CUT "\n";
$cut++;
close(CUT);
print "wrote $wrote_ltr, $n_count\n";
#SEND TO QSUB WAIT FOR SEND CONFORMATION
@qsub_sent = ();
$script = 'q-script';
system("touch $input_dir/$script");
for(my $y=0;$y<$cut;$y++)
  {
  system("rm $input_dir/$script");
  open(SCRIPT, ">$input_dir/$script");
  print SCRIPT "./TE_nest.pl $input_dir/$input_sequence-f-$y $options --output $input_dir --frag_only T";
  close(SCRIPT);
  my $out;
  do{
    $out = system("qsub -d . -k eo $input_dir/$script > $input_dir/q-out");
    sleep(2);
    open(OUT,"$input_dir/q-out");
    my $qout = <OUT>;
    chomp($qout);
    close(OUT);
#    system("rm $input_dir/q-out");
    push @qsub_sent, $qout;
    }until($out == 0);
  }
#system("rm $input_dir/q-script");
#CHECK FOR COMPLETED FILES - check if file is off qstat, check if *.LTR exists, die if not found
print "Waiting for FRAG node runs to complete";
for(my $y=0;$y<$cut;$y++)
  {
  my $done = 0;
  do{
    my $grep = system("qstat | grep $qsub_sent[$y] > /dev/null");
    # $grep returns 0 if found, 256 if not found
    if($grep == 0)
      {
      sleep 5;
      print ".";
      }  #RAISE THIS
     else
      {$done = 1;}
    }until($done == 1);
  my $check_number = 0;
  my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-f-$y > $input_dir/q-ls");
  open(OUT,"$input_dir/q-ls");
  my $qls = <OUT>;
  chomp($qls);
  close(OUT);
  system("rm $input_dir/q-ls");
  while(!(-e "$qls/$input_sequence-f-$y.LTR") && $check_number < 3)
    {
    sleep 10;  #RAISE THIS
    $check_number++;
    }
  if(!(-e "$qls/$input_sequence-f-$y.LTR"))
    {die "Round FRAG Cut $y has finished running, but *.LTR doesn't exist.  Exiting...\n";}
  }
print "\n";
#REJOIN SEQUENCES AND *.LTR ANNOTATIONS
@seq = ();
$seq = '';
for(my $y=0;$y<$cut;$y++)
  {
  my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-f-$y > $input_dir/q-ls");
  open(OUT,"$input_dir/q-ls");
  my $qls = <OUT>;
  chomp($qls);
  close(OUT);
  system("rm $input_dir/q-ls");
  open(SEQ,"$qls/$input_sequence-f-$y.mask") or die "$input_sequence-f-$y.mask not found, TE_nest on round frag cut $cut failed.\n";
  while(<SEQ>)
    {
    my $line = $_;
    chomp($line);
    if($line !~ m/>/)
      {$seq = $seq . $line;}
    }
  }
close(SEQ);
@seq = split(//,$seq);
#$n_count = 0;
#@seq_n = ();
#for(my $x=0;$x<@seq;$x++)
#  {
#  $seq_n[$x][0] = $seq[$x];
#  $seq_n[$x][1] = 0;
#  if($seq_n[$x][0] ne 'N' && $seq_n[$x][0] ne 'n')
#    {
#    $seq_n[$x][1] = $n_count;
#    $n_count++;
#    }
#  }
@ltr_final = ();
my $frag_count = 0;
my $nltr_count = 0;
my $n = 'n0';
my $f = 'f0';
for(my $y=0;$y<$cut;$y++)
  {
  my $ls_directory = system("ls -1d $input_dir/\*_$input_sequence-f-$y > $input_dir/q-ls");
  open(OUT,"$input_dir/q-ls");
  my $qls = <OUT>;
  chomp($qls);
  close(OUT);
  system("rm $input_dir/q-ls");
  open(LTR,"$qls/$input_sequence-f-$y.LTR");
  while(<LTR>)
    {
    my $line = $_;
    chomp($line);
    my @ltr_tmp = split(/ /,$line);
    my $increase = $real_cut_size * $y;
    if($ltr_tmp[0] eq 'FRAG')
      {
      $f = 'f'.$frag_count;
      $ltr_tmp[1] = $f;
      $wrote_ltr = 1;
      my $ltr = "$ltr_tmp[0] $ltr_tmp[1] $ltr_tmp[2] $ltr_tmp[3] $ltr_tmp[4] $ltr_tmp[5] $ltr_tmp[6]";
print "writing $y $ltr\n";
      push @ltr_final, [split(/ /,$ltr)];
      $frag_count++;
      }
    if($ltr_tmp[0] =~ m/^f/)
      {
      $ltr_tmp[0] = $f;
      my $ltr = "$ltr_tmp[0]";
      for(my $y=0;$y<(@ltr_tmp-1)/4;$y++)
        {
        $ltr_tmp[($y*4)+1] = $ltr_tmp[($y*4)+1] + $increase;
        $ltr_tmp[($y*4)+2] = $ltr_tmp[($y*4)+2] + $increase;
        $ltr = $ltr . " $ltr_tmp[($y*4)+1] $ltr_tmp[($y*4)+2] $ltr_tmp[($y*4)+3] $ltr_tmp[($y*4)+4]";
        }
      push @ltr_final, [split(/ /,$ltr)];
      }
    if($ltr_tmp[0] eq 'NLTR')
      {
      $n = 'n'.$nltr_count;
      $ltr_tmp[1] = $n;
      $wrote_ltr = 1;
      my $ltr = "$ltr_tmp[0] $ltr_tmp[1] $ltr_tmp[2] $ltr_tmp[3] $ltr_tmp[4] $ltr_tmp[5] $ltr_tmp[6]";
print "writing $y $ltr\n";
      push @ltr_final, [split(/ /,$ltr)];
      $nltr_count++;
      }
    if($ltr_tmp[0] =~ m/^n/)
      {
      $ltr_tmp[0] = $n;
      my $ltr = "$ltr_tmp[0]";
      for(my $y=0;$y<(@ltr_tmp-1)/4;$y++)
        {
        $ltr_tmp[($y*4)+1] = $ltr_tmp[($y*4)+1] + $increase;
        $ltr_tmp[($y*4)+2] = $ltr_tmp[($y*4)+2] + $increase;
        $ltr = $ltr . " $ltr_tmp[($y*4)+1] $ltr_tmp[($y*4)+2] $ltr_tmp[($y*4)+3] $ltr_tmp[($y*4)+4]";
        }
      push @ltr_final, [split(/ /,$ltr)];
      }
    }
  close(LTR);
  }
LOAD_HASH(@ltr_final);
DISCREPANCY('FRAG-NLTR','FRAG-NLTR');
DISCREPANCY('FRAG-NLTR','PAIR');
DISCREPANCY('SOLO','FRAG-NLTR');
OVERLAP_S_FN('SOLO','FRAG-NLTR');
$seq_ref = PRINT_LTR('T',\@ltr_final,\@seq);

print "writing mask\n";
#PRINT MASKED SEQUENCE FILE
my $mask = $input_sequence;
$mask =~ s/.fasta//g;
$mask = $mask . '.mask';
my $seq_out = "$input_dir/$mask";
open(OUT_SEQ, ">$seq_out");
print OUT_SEQ ">$mask\n";
my $x = 0;
for($x=0;$x<(@seq/70)-1;$x++)
  {
  for(my $y=0;$y<70;$y++)
    {
    print OUT_SEQ "$seq[($x*70)+$y]";
    }
  print OUT_SEQ "\n";
  }
for(my $y=$x*70;$y<@seq;$y++)
  {
  print OUT_SEQ "$seq[$y]";
  }
print OUT_SEQ "\n";
close(OUT_SEQ);

################################################################################################
#                                       SUBROUTINES                                            #
################################################################################################

sub LOAD_HASH
{
my @ltr_final = @_;
my $total_inserts = 0;
my @line = ();
for(my $w=0;$w<@ltr_final;$w++)
  {
  my $line = $ltr_final[$w][0];
  for(my $x=1;$x<(scalar @{$ltr_final[$w]});$x++)
    {$line = $line . " $ltr_final[$w][$x]";}
  @line = split(/ /,$line);
  if($line[0] eq 'PAIR')
    {
    $total_inserts++;
    $retro{pair}{$line[1]} = {type => $line[2], dir => $line[3], bsr => $line[4], group => $line[5], order => $line[6], level => $line[7]};
    for(my $y=0;$y<3;$y++)
      {
      $w++;
      my $line = $ltr_final[$w][0];
      for(my $x=1;$x<(scalar @{$ltr_final[$w]});$x++)
        {$line = $line . " $ltr_final[$w][$x]";}
      @line = split(/ /,$line);
      my $pair_count = @line;
      $pair_count = ($pair_count-2)/4;
      my $section_size = 0;
      for(my $x=0;$x<$pair_count;$x++)
        {
        $retro{pair}{$line[0]}{coords}{$line[1]}{$x} = {SEQ_start => $line[($x*4)+2], SEQ_end => $line[($x*4)+3], TE_start => $line[($x*4)+4], TE_end => $line[($x*4)+5]};
        $section_size = $section_size + $retro{pair}{$line[0]}{coords}{$line[1]}{$x}{SEQ_end} - $retro{pair}{$line[0]}{coords}{$line[1]}{$x}{SEQ_start};
        }
      if($section_size != 0)
        {$retro{pair}{$line[0]}{size}{$line[1]} = $section_size;}
      if($y==0)
        {$retro{pair}{$line[0]}{start} = $retro{pair}{$line[0]}{coords}{$line[1]}{0}{SEQ_start};}
       elsif($y==1)
        {$retro{pair}{$line[0]}{end} = $retro{pair}{$line[0]}{coords}{$line[1]}{$pair_count-1}{SEQ_end};}
      }
    }
  }
my @seq_spots = ();
my $hash_pair = scalar keys %{$retro{pair}};
my $x = my $found = 0;
do{
  my $p = 'p'.$x;
  if($retro{pair}{$p})
    {
    my $amt = 0;
    my $hit_amt = scalar keys %{$retro{pair}{$p}{coords}{L}};
    for(my $z=0;$z<$hit_amt;$z++)
      {
      $retro{pair}{$p}{coords}{$amt}{SEQ_start} = $retro{pair}{$p}{coords}{L}{$z}{SEQ_start};
      $retro{pair}{$p}{coords}{$amt}{SEQ_end} = $retro{pair}{$p}{coords}{L}{$z}{SEQ_end};
      $retro{pair}{$p}{coords}{$amt}{TE_start} = $retro{pair}{$p}{coords}{L}{$z}{TE_start};
      $retro{pair}{$p}{coords}{$amt}{TE_end} = $retro{pair}{$p}{coords}{L}{$z}{TE_end};
      $amt++;
      }
    $hit_amt = scalar keys %{$retro{pair}{$p}{coords}{M}};
    for(my $z=0;$z<$hit_amt;$z++)
      {
      $retro{pair}{$p}{coords}{$amt}{SEQ_start} = $retro{pair}{$p}{coords}{M}{$z}{SEQ_start};
      $retro{pair}{$p}{coords}{$amt}{SEQ_end} = $retro{pair}{$p}{coords}{M}{$z}{SEQ_end};
      $retro{pair}{$p}{coords}{$amt}{TE_start} = $retro{pair}{$p}{coords}{M}{$z}{TE_start};
      $retro{pair}{$p}{coords}{$amt}{TE_end} = $retro{pair}{$p}{coords}{M}{$z}{TE_end};
      $amt++;
      }
    $hit_amt = scalar keys %{$retro{pair}{$p}{coords}{R}};
    for(my $z=0;$z<$hit_amt;$z++)
      {
      $retro{pair}{$p}{coords}{$amt}{SEQ_start} = $retro{pair}{$p}{coords}{R}{$z}{SEQ_start};
      $retro{pair}{$p}{coords}{$amt}{SEQ_end} = $retro{pair}{$p}{coords}{R}{$z}{SEQ_end};
      $retro{pair}{$p}{coords}{$amt}{TE_start} = $retro{pair}{$p}{coords}{R}{$z}{TE_start};
      $retro{pair}{$p}{coords}{$amt}{TE_end} = $retro{pair}{$p}{coords}{R}{$z}{TE_end};
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
my $solo_number_count = 0;
my $frag_number_count = 0;
my $nltr_number_count = 0;
my $seq_spots_count = @seq_spots;
for(my $w=0;$w<@ltr_final;$w++)
  {
  my $line = $ltr_final[$w][0];
  for(my $x=1;$x<(scalar @{$ltr_final[$w]});$x++)
    {$line = $line . " $ltr_final[$w][$x]";}
  @line = split(/ /,$line);
  if($line[0] eq 'SOLO')
    {
    my $unique = 0;
    my $number = 's'.$solo_number_count;
    my $type = $line[2];
    my $dir = $line[3];
    $w++;
    my $line = $ltr_final[$w][0];
    for(my $x=1;$x<(scalar @{$ltr_final[$w]});$x++)
      {$line = $line . " $ltr_final[$w][$x]";}
    @line = split(/ /,$line);
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
      $retro{solo}{$number}{group} = 0;
      $retro{solo}{$number}{order} = 0;
      $retro{solo}{$number}{level} = 0;
      $retro{solo}{$number}{bsr} = 0;
      $retro{solo}{$number}{start} = $retro{solo}{$number}{coords}{0}{SEQ_start};
      $retro{solo}{$number}{end} = $retro{solo}{$number}{coords}{$count-1}{SEQ_end};
      $solo_number_count++;
      }
    }
   elsif($line[0] eq 'FRAG')
    {
    my $seq_spots_count = @seq_spots;
    my $unique = 0;
    my $number = 'f'.$frag_number_count;
    my $type = $line[2];
    my $dir = $line[3];
    $w++;
    my $line = $ltr_final[$w][0];
    for(my $x=1;$x<(scalar @{$ltr_final[$w]});$x++)
      {$line = $line . " $ltr_final[$w][$x]";}
    @line = split(/ /,$line);
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
      $retro{frag}{$number}{group} = 0;
      $retro{frag}{$number}{order} = 0;
      $retro{frag}{$number}{level} = 0;
      $retro{frag}{$number}{bsr} = 0;
      $retro{frag}{$number}{start} = $retro{frag}{$number}{coords}{0}{SEQ_start};
      $retro{frag}{$number}{end} = $retro{frag}{$number}{coords}{$count-1}{SEQ_end};
      $frag_number_count++;
      }
    }
   elsif($line[0] eq 'NLTR')
    {
      my $seq_spots_count = @seq_spots;
      my $unique = 0;
      my $number = 'n'.$nltr_number_count;
      my $type = $line[2];
      my $dir = $line[3];
      $w++;
      my $line = $ltr_final[$w][0];
      for(my $x=1;$x<(scalar @{$ltr_final[$w]});$x++)
        {$line = $line . " $ltr_final[$w][$x]";}
      @line = split(/ /,$line);
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
      $retro{nltr}{$number}{group} = 0;
      $retro{nltr}{$number}{order} = 0;
      $retro{nltr}{$number}{level} = 0;
      $retro{nltr}{$number}{bsr} = 0;
      $retro{nltr}{$number}{start} = $retro{nltr}{$number}{coords}{0}{SEQ_start};
      $retro{nltr}{$number}{end} = $retro{nltr}{$number}{coords}{$count-1}{SEQ_end};
      $nltr_number_count++;
      }
    }
  }
#PRINT HASH
print "PRINT HASH\n";
#SOLO
my $hash_solo = scalar keys %{$retro{solo}};
print "SOLOs $hash_solo\n";
for(my $x=0;$x<$hash_solo;$x++)
  {
  my $s = 's'.$x;
  print "solo-$x type-$retro{solo}{$s}{type} dir-$retro{solo}{$s}{dir}";
  my $solo_count = scalar keys %{$retro{solo}{$s}{coords}};
  print " parts-$solo_count\n";
  }
#PAIR
$hash_pair = scalar keys %{$retro{pair}};
print "PAIRs\n";
for(my $x=0;$x<$hash_pair;$x++)
  {
  my $p = 'p' . $x;
  print "pair-$x type-$retro{pair}{$p}{type} dir-$retro{pair}{$p}{dir} rate-$retro{pair}{$p}{bsr}";
  print " group-$retro{pair}{$p}{group} level-$retro{pair}{$p}{level} order-$retro{pair}{$p}{order}\n";
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
    print " $part-parts-$pair_count";
    for(my $z=0;$z<$pair_count;$z++)
      {
      print "  $retro{pair}{$p}{coords}{$part}{$z}{SEQ_start} $retro{pair}{$p}{coords}{$part}{$z}{SEQ_end}";
      }
    print "\n";
    }
  }
#FRAG
my $hash_frag = scalar keys %{$retro{frag}};
print "FRAGs $hash_frag\n";
for(my $x=0;$x<$hash_frag;$x++)
  {
  my $f = 'f' . $x;
  print "frag-$f type-$retro{frag}{$f}{type} dir-$retro{frag}{$f}{dir}";
  my $frag_count = scalar keys %{$retro{frag}{$f}{coords}};
  print " parts-$frag_count\n";
  }
#NLTR
my $hash_nltr = scalar keys %{$retro{nltr}};
print "NLTRs $hash_nltr\n";
for(my $x=0;$x<$hash_nltr;$x++)
  {
  my $n = 'n'.$x;
  print "nltr-$x type-$retro{nltr}{$n}{type} dir-$retro{nltr}{$n}{dir}";
  my $nltr_count = scalar keys %{$retro{nltr}{$n}{coords}};
  print " parts-$nltr_count\n";
  }
}

################################################################################################

sub PRINT_LTR
{
my ($print,$ltr_final_ref,$seq_ref) = @_;
my @ltr_final = @$ltr_final_ref;
my @seq = @$seq_ref;
my $seq_count_11 = @seq;
print "in PRINT_LTR $seq_count_11\n";
#PRINT OUT TO FINAL FILE
if($print eq 'T')
  {
  my $LTR = $input_sequence;
  $LTR =~ s/.fasta//g;
  $LTR = $LTR.'.LTR';
  open(LTR, ">$input_dir/$LTR");
  print LTR "$seq_count\n";
  print LTR "$input_sequence\n";
  }
#SOLO
my $seq_count_11 = @seq;
print "in 1 PRINT_LTR $seq_count_11\n";
my $hash_solo = scalar keys %{$retro{solo}};
for(my $x=0;$x<$hash_solo;$x++)
  {
  my $s = 's'.$x;
  if($print eq 'T')
    {
    print LTR "SOLO $s $retro{solo}{$s}{type} $retro{solo}{$s}{dir} $retro{solo}{$s}{bsr} $retro{solo}{$s}{group} $retro{solo}{$s}{order} $retro{solo}{$s}{level}\n";
    print "??? SOLO $s $retro{solo}{$s}{type} $retro{solo}{$s}{dir} $retro{solo}{$s}{bsr} $retro{solo}{$s}{group} $retro{solo}{$s}{order} $retro{solo}{$s}{level}\n";
    print LTR "$s";
    }
  my $solo_count = scalar keys %{$retro{solo}{$s}{coords}};
  for(my $z=0;$z<$solo_count;$z++)
    {
    if($print eq 'T')
      {
      print LTR " $retro{solo}{$s}{coords}{$z}{SEQ_start} $retro{solo}{$s}{coords}{$z}{SEQ_end} $retro{solo}{$s}{coords}{$z}{TE_start} $retro{solo}{$s}{coords}{$z}{TE_end}";
      }
    for(my $y=$retro{solo}{$s}{coords}{$z}{SEQ_start};$y<$retro{solo}{$s}{coords}{$z}{SEQ_end}+1;$y++)
      {
      $seq[$y] = 'N';
      }
    }
  if($print eq 'T')
    {print LTR "\n";}
  }
my $seq_count_11 = @seq;
print "in 2 PRINT_LTR $seq_count_11\n";

#PAIR
my $hash_pair = scalar keys %{$retro{pair}};
for(my $x=0;$x<$hash_pair;$x++)
  {
  my $p = 'p' . $x;
  if($print eq 'T')
    {
    print LTR "PAIR $p $retro{pair}{$p}{type} $retro{pair}{$p}{dir} ";
    printf LTR ("%.3f", $retro{pair}{$p}{bsr});
    print LTR " $retro{pair}{$p}{group} $retro{pair}{$p}{order} $retro{pair}{$p}{level}\n";
    }
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
    if($print eq 'T')
      {print LTR "$p $part";}
    for(my $z=0;$z<$pair_count;$z++)
      {
      if($print eq 'T')
        {print LTR " $retro{pair}{$p}{coords}{$part}{$z}{SEQ_start} $retro{pair}{$p}{coords}{$part}{$z}{SEQ_end} $retro{pair}{$p}{coords}{$part}{$z}{TE_start} $retro{pair}{$p}{coords}{$part}{$z}{TE_end}";}
      for(my $y=$retro{pair}{$p}{coords}{$part}{$z}{SEQ_start};$y<$retro{pair}{$p}{coords}{$part}{$z}{SEQ_end}+1;$y++)
        {
        $seq[$z] = 'N';
        }
      }
    if($print eq 'T')
      {print LTR "\n";}
    }
  }
my $seq_count_11 = @seq;
print "in 3 PRINT_LTR $seq_count_11\n";

#FRAG
my $hash_frag = scalar keys %{$retro{frag}};
for(my $x=0;$x<$hash_frag;$x++)
  {
  my $f = 'f' . $x;
  if($print eq 'T')
    {
    print LTR "FRAG $f $retro{frag}{$f}{type} $retro{frag}{$f}{dir} $retro{frag}{$f}{bsr} $retro{frag}{$f}{group} $retro{frag}{$f}{order} $retro{frag}{$f}{level}\n";
    print LTR "$f";
    }
  my $frag_count = scalar keys %{$retro{frag}{$f}{coords}};
  for(my $z=0;$z<$frag_count;$z++)
    {
    if($print eq 'T')
      {print LTR " $retro{frag}{$f}{coords}{$z}{SEQ_start} $retro{frag}{$f}{coords}{$z}{SEQ_end} $retro{frag}{$f}{coords}{$z}{TE_start} $retro{frag}{$f}{coords}{$z}{TE_end}";}
    for(my $y=$retro{frag}{$f}{coords}{$z}{SEQ_start};$y<$retro{frag}{$f}{coords}{$z}{SEQ_end}+1;$y++)
      {
      $seq[$y] = 'N';
      }
    }
  if($print eq 'T')
    {print LTR "\n";}
  }
my $seq_count_11 = @seq;
print "in 4 PRINT_LTR $seq_count_11\n";

#NLTR
my $hash_nltr = scalar keys %{$retro{nltr}};
for(my $x=0;$x<$hash_nltr;$x++)
  {
  my $n = 'n'.$x;
  if($print eq 'T')
    {
    print LTR "NLTR $n $retro{nltr}{$n}{type} $retro{nltr}{$n}{dir} $retro{nltr}{$n}{bsr} $retro{nltr}{$n}{group} $retro{nltr}{$n}{order} $retro{nltr}{$n}{level}\n";
    print LTR "$n";
    }
  my $nltr_count = scalar keys %{$retro{nltr}{$n}{coords}};
  for(my $z=0;$z<$nltr_count;$z++)
    { 
    if($print eq 'T')
      {print LTR " $retro{nltr}{$n}{coords}{$z}{SEQ_start} $retro{nltr}{$n}{coords}{$z}{SEQ_end} $retro{nltr}{$n}{coords}{$z}{TE_start} $retro{nltr}{$n}{coords}{$z}{TE_end}";}
    for(my $y=$retro{nltr}{$n}{coords}{$z}{SEQ_start};$y<$retro{nltr}{$n}{coords}{$z}{SEQ_end}+1;$y++)
      {
      $seq[$y] = 'N';
      }
    }
  if($print eq 'T')
    {print LTR "\n";}
  }
if($print eq 'T')
  {close(LTR);}
my $seq_size = @seq;
my $seq_count_11 = @seq;
print "in PRINT_LTR last $seq_count_11\n";
return(\@seq);
}
################################################################################################

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
################################################################################################

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
      $retro{solo}{$p} = {type => $second_head[$x][1], dir => $second_head[$x][2], bsr => $second_head[$x][3],                          group => $second_head[$x][4], level => $second_head[$x][5], order => $second_head[$x][6]};
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
      $retro{pair}{$p} = {type => $second_head[$x][1], dir => $second_head[$x][2], bsr => $second_head[$x][3],                          group => $second_head[$x][4], level => $second_head[$x][5], order => $second_head[$x][6]};
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

sub OVERLAP_S_P
{
my $first_dis = $_[0];
my $second_dis = $_[1];
LOAD_DISCREPANCY($first_dis,$second_dis);
my @temp_dis_1 = ();
my @temp_head_1 = ();

#for(my $x=0;$x<@second_dis;$x=$x+3)
#  {
#  my $second_start = scalar @{$second_dis[$x]};
#  for(my $b=2;$b>0;$b--)
#    {
#    for(my $c=0;$c<scalar @{$second_dis[$x+$b]};$c++)
#      {
#      $second_dis[$x][$second_start] = $second_dis[$x+$b][$c];
#      $second_start++;
#      }
#    }
#for(my $z=0;$z<(scalar @{$second_dis[$x]});$z++)
#  {print "$second_dis[$x][$z] ";}
#print "\n";
#  }

for(my $w=0;$w<@first_dis;$w++)
  {
  my @remove_list_1 = ();
  for(my $x=0;$x<@second_dis;$x=$x+3)
    {
    if($second_dis[$x][0] ne 'blank')
      {
      for(my $y=0;$y<(scalar @{$first_dis[$w]})/4;$y++)
        {
        for(my $z=0;$z<(scalar @{$second_dis[$x]})/4;$z++)
          {
#print "$second_dis[$x][$z*4] $second_dis[$x][($z*4)+1] $second_dis[$x][($z*4)+2] $second_dis[$x][($z*4)+3]\n";
#          print "($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)+1])\n";
          if($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)+1])
            {
            my $remove = "$w $y";
            push @remove_list_1, [split(/ /,$remove)];
            }
           elsif($first_dis[$w][$y*4] <= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] <= $second_dis[$x][($z*4)
+1] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][$z*4])
            {
            if(($first_dis[$w][($y*4)+1] - $first_dis[$w][$y*4]) <= ($second_dis[$x][($z*4)+1] - $second_dis[$x][$z*
4]))
              {
              $first_dis[$w][($y*4)+1] = $second_dis[$x][$z*4]-1;
              }
             else
              {
              $second_dis[$x][$z*4] = $first_dis[$w][($y*4)+1]+1;
              }
            }
           elsif($first_dis[$w][$y*4] >= $second_dis[$x][$z*4] && $first_dis[$w][($y*4)+1] >= $second_dis[$x][($z*4)
+1] && $first_dis[$w][$y*4] <= $second_dis[$x][($z*4)+1])
            {
            if(($first_dis[$w][($y*4)+1] - $first_dis[$w][$y*4]) <= ($second_dis[$x][($z*4)+1] - $second_dis[$x][$z*
4]))
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
      $temp_dis = $temp_dis . " $first_dis[$w][$y*4] $first_dis[$w][($y*4)+1] $first_dis[$w][($y*4)+2] $first_dis[$w
][($y*4)+3]";
      $found = 1;
      }
    }
  if($found == 1)
    {
    $temp_dis =~ s/^\s//g;
    push @temp_dis_1, [split(/ /,$temp_dis)];
    my $temp_head = "$first_head[$w][0] $first_head[$w][1] $first_head[$w][2] $first_head[$w][3] $first_head[$w][4] $first_head[$w][5] $first_head[$w][6]";
    push @temp_head_1, [split(/ /,$temp_head)];
print "TEMP head $temp_head\n";
print "TEMP DIS $temp_dis\n";
    }
  }
@first_head = @temp_head_1;
@first_dis = @temp_dis_1;
UNLOAD_DISCREPANCY($first_dis,$second_dis);
}
################################################################################################

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
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
