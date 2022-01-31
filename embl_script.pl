#!/usr/bin/env perl

use warnings;
use strict;

my $in_file=$ARGV[0];
my $out_file='';

if($in_file =~ /(\S+)(\.embl)/){
	my $base=$1;my $ext=$2;
	
	$out_file="${base}_formatted${ext}";
}else{
	die "input file is not in embl format\n";
}

open OUTPUT, ">$out_file";

print OUTPUT "Position\tnode\tcolour\ttaxa\tparent_base\treplace\n";

my $pos=0;
my $node='';
my $colour=0;
my $taxa='';
my $parent='';
my $replace='';

my $count=0;

open INPUT, "$in_file";
while(my $line=<INPUT>){
	chomp $line;
	
	$count++;
	
	if($line =~ /^FT\s+variation\s+(\d+)/){
		$pos=$1;
		
		$node='';
		$colour=0;
		$taxa='';
		$parent='';
		$replace='';
		
		$count=0;
		
	}elsif($line =~ /^FT\s+\/node\=\"([^\"]+)\"/){
		$node=$1;
	}elsif($line =~ /^FT\s+\/colour\=\"([^\"]+)\"/){
		$colour=$1;
	}elsif($line =~ /^FT\s+\/taxa\=\"([^\"]+)\"/){
		$taxa=$1;
	}elsif($line =~ /^FT\s+\/parent_base\=\"([^\"]+)\"/){
		$parent=$1;
	}elsif($line =~ /^FT\s+\/replace\=\"([^\"]+)\"/){
		$replace=$1;
		
		if($count == 5){
			print OUTPUT "$pos\t$node\t$colour\t$taxa\t$parent\t$replace\n";
		}else{
			die "variant $pos doesn't have enough attributes.\n";
		}
	}
}

