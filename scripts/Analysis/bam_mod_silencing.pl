#takes bam input from STDIN line by line and forces any specified mod calls below a threshold to 0
#Nicolas Altemose 2022
use warnings;
use strict;
my $tic = time;

my $usage = "perl ThresholdBAMModCalls_separateAC_v0.pl <mA Threshold 1-255 (set to >255 to eliminate A calls) [155]> <mC Threshold 1-255 (set to >255 to eliminate C calls) [155]>\n";

my $threshA=155;
my $threshC=155;

if(defined $ARGV[0]){
	$threshA = $ARGV[0];
	chomp($threshA);
}else{
	print STDERR "\nNo mA threshold specified. Using default thresholds of 155 for both mA and mC\n";
}
if(defined $ARGV[1]){
	$threshC = $ARGV[1];
	chomp($threshC);
}else{
	print STDERR "\nNo mC threshold specified. Using default threshold of 155 for mC\n";
}
die "You must specify a modification probability threshold from 1-255, or >255 to set all calls to 0\n\n$usage\n\n" unless($threshA>=1 && $threshC>=1);


while(my $line = <STDIN>){
	chomp($line);
	if($line=~/^\@/){
		print "$line\n";
	}elsif($line=~/^(.+)\s+Mm\:Z\:A\+a\,(.+)\;C\+m\,(.+)\;\s+Ml:B:C\,(\S+)$/){
		my $before=$1;
		my $stringC = $3;
		my $stringA = $2;
		my @valarray = split(',',$4);
		my $numC = 1 + $stringC=~tr/\,/\,/;
		my $numA = 1 + $stringA=~tr/\,/\,/;
		my $size = scalar @valarray;
		my @Carray = map($_<$threshC ? 0 : $_, @valarray[$numA..$size-1]);
		my @Aarray = map($_<$threshA ? 0 : $_, @valarray[0..$numA-1]);
		my $printstring = join(',',@Aarray).','.join(',',@Carray);
		#my $printstring = join(',',map($_<$threshA ? 0 : $_, @valarray[0..$numA-1])).','.join(',',map($_<$threshC ? 0 : $_, @valarray[$numA..$size-1]));
		print "$before\tMm:Z:A+a,$stringA;C+m,$stringC;\tMl:B:C,$printstring\n";
	}elsif($line=~/^(.+)\s+Mm\:Z\:C\+m\,(.+)\;A\+a\,(.+)\;\s+Ml:B:C\,(\S+)$/){
		my $before=$1;
		my $stringC = $2;
		my $stringA = $3;
		my @valarray = split(',',$4);
		my $numC = 1 + $stringC=~tr/\,/\,/;
		my $numA = 1 + $stringA=~tr/\,/\,/;
		my $size = scalar @valarray;
		my @Carray = map($_<$threshC ? 0 : $_, @valarray[0..$numC-1]);
		my @Aarray = map($_<$threshA ? 0 : $_, @valarray[$numC..$size-1]);
		my $printstring = join(',',@Carray).','.join(',',@Aarray);
		print "$before\tMm:Z:C+m,$stringC;A+a,$stringA;\tMl:B:C,$printstring\n";
	}elsif($line=~/^(.+)\s+Mm\:Z\:C\+m\,(.+)\;\s+Ml:B:C\,(\S+)$/){
		my $before=$1;
		my $stringC = $2;
		my @valarray = split(',',$3);
		my $numC = 1 + $stringC=~tr/\,/\,/;
		my $size = scalar(@valarray);
		my @Carray = map($_<$threshC ? 0 : $_, @valarray);
		my $printstring = join(',',@Carray);
		print "$before\tMm:Z:C+m,$stringC;\tMl:B:C,$printstring\n";
	}elsif($line=~/^(.+)\s+Mm\:Z\:A\+a\,(.+)\;\s+Ml:B:C\,(\S+)$/){
		my $before=$1;
		my $stringA = $2;
		my @valarray = split(',',$3);
		my $numA = 1 + $stringA=~tr/\,/\,/;
		my $size = scalar(@valarray);
		my @Aarray = map($_<$threshA ? 0 : $_, @valarray);
		my $printstring = join(',',@Aarray);
		print "$before\tMm:Z:A+a,$stringA;\tMl:B:C,$printstring\n";
	}else{
		die "ERROR could not find Ml tag in this entry:$line\n";
	}
}


#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
my $time1 = sprintf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
print STDERR "$time1\n";