#!/usr/bin/perl

# Add taxonomy to USEARCH Otu table
# 7.December 2016/ Pfeiffer; last modified 1/7/2018
# Contact: microbiawesome@gmail.com; 
# Copyright (C) Stefan Pfeiffer, 2016-2018, all rights reserved.# FASPA is a workflow for analysing Illumina paired-end sequence data. 
# This file is distributed without warranty 
# Cite as: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses, DOI:10.5281/yenodo.1302800

use strict;

if ($#ARGV != 2) {
  print "$ARGV[0] <first file> <second file> <target file>\n";
  exit 1;
}

my $F_OTUTABLE; my $F_TAXONOMY; my $F_OUTPUT;
my %OTUTABLE;
my %TAXONOMY;
my $comment='';

open($F_OTUTABLE, $ARGV[0]) || die "Cant open OTUTable $ARGV[0] \n";
while(<$F_OTUTABLE>) {
  my $OTUTEXT=$_;
  my @OTULINE=split(/\t/);
  if (! ( $OTULINE[0]=~ /^Otu/ || /^\s*\#/ || /^\s*$/) ) { die "No OTUTable  $_ \n"; }
  if (! (/^\s*\#/ || /^\s*$/ )) {
    $OTUTABLE{$OTULINE[0]}=$OTUTEXT; 
  } else {
    $comment.=$OTUTEXT;
  }
}

open($F_TAXONOMY, $ARGV[1]) || die "Cant open OTUTable $ARGV[1] \n";
while(<$F_TAXONOMY>) {
  my @OTULINE=split(/\t/);  my $Otu=$OTULINE[0];
  if (! ($Otu=~ /^Otu/ || /^\s*\#/ || /^\s*$/ )) { die "No Taxonomy File"; }
  if (! (/^\s*\#/ || /^\s*$/ )) {
    $TAXONOMY{$Otu}=$OTULINE[3];
    if ($OTULINE[3] eq 'd:Bacteria') {
      delete $OTUTABLE{$Otu};  
      delete $TAXONOMY{$Otu}; 
    }
  }
}

open($F_OUTPUT, ">", $ARGV[2]) || die "Cant open output file";

print $F_OUTPUT $comment;

my $line=0;
foreach my $Otu (keys %OTUTABLE) {
  my $i=0;
  my $taxonomy='';
  my @tax_split=split(',',$TAXONOMY{$Otu});
  foreach my $tax_field (@tax_split) {
    if ($i > 0) { $taxonomy.="; "; }
    # Split taxonomic level identifier from the information
    my @taxinfo=split(':',$tax_field);
    $taxonomy.="D_".$i."__".$taxinfo[1];
    $i++;
  }
  my $outputline=$OTUTABLE{$Otu}."\t".$taxonomy;
  $outputline=~s/[\"\n]//g;


  print $F_OUTPUT $outputline."\n";
  $line++;
}

close $F_OUTPUT;

print "File ".$ARGV[2]." created, contains $line Otu values\n";
