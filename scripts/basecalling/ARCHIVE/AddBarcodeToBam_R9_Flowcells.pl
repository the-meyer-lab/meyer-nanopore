#takes bam input from STDIN line by line and adds barcode info from a barcoding summary file output by guppy_barcoder
#Nicolas Altemose 2022
use warnings;
use strict;

my $usage = "perl AddBarcodeToBam.pl <Path_to_barcoding_summary.txt> <space_separated_list_of_barcode_numbers_to_keep>\n";

my $infile1 = "barcoding_summary.txt";

my @barcodelist = qw(
01
02
03
04
05
06
07
08
09
10
11
12
13
14
15
16
17
18
19
20
24
);


if(defined $ARGV[0]){
	$infile1 = $ARGV[0];
	chomp($infile1);
}
if(defined $ARGV[1]){
	shift(@ARGV);
	@barcodelist = @ARGV;
}
unshift @barcodelist, '0';


my %barcodes;
my %barcodesMapped;
foreach my $barcode(@barcodelist){
	chomp($barcode);
	$barcodes{$barcode}=0;
	$barcodesMapped{$barcode}=0;
}

my %memory;
open(IN, $infile1);
my $header = <IN>;
my $ct=0;
while(my $line = <IN>){
	chomp($line);
	$ct++;
	if($line=~/^(\S+)\s+barcode(\d+)\s+/){
		if(exists $barcodes{$2}){
			$memory{$1}=$2;
			$barcodes{$2}++;
		}else{
			$memory{$1}=0;
			$barcodes{0}++;
		}
	}elsif($line=~/^(\S+)\s+unclassified/){
		$memory{$1}=0;
		$barcodes{0}++;
	}
}
close IN;



while(my $line = <STDIN>){
	chomp($line);
	if($line=~/^\@RG/){
		foreach my $barcode(@barcodelist){
			print "\@RG\tID:$barcode\tSM:barcode$barcode\n";
		}
	}elsif($line=~/^\@/){
		print "$line\n";
	}
	elsif($line=~/^((\S+)\s+.+)\s+RG:Z:\S+\s+(.+)$/){
		if(exists $memory{$2}){
			my $bc = $memory{$2};
			print "$1\tRG:Z:$bc\t$3\n";
			$barcodesMapped{$bc}++;
		}else{
			print "$1\tRG:Z:0\t$3\n";
			$barcodesMapped{0}++;
		}
	}else{
		die "ERROR:$line\n";
	}
}

my $sum1 = 0;
my $sum2 = 0;
print STDERR "Barcode\tNumber\tNumberMapped\n";
foreach my $barcode(@barcodelist){
	print STDERR "$barcode\t$barcodes{$barcode}\t$barcodesMapped{$barcode}\n";
	$sum1+=$barcodes{$barcode};
	$sum2+=$barcodesMapped{$barcode};
}
print STDERR "ALL\t$sum1\t$sum2\n";