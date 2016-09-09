#!/usr/bin/perl
use warnings;
use strict;

my %Viromes = (
"A0001" => "A0001",
"A0003" => "A0003", 
"A0007" => "A0007", 
"A0009" => "A0009", 
"A0011" => "A0011", 
"A0013" => "A0013", 
"A0019" => "A0019",
"A0021" => "A0021",
"A0023" => "A0023",
"A0027" => "A0027",
"A0035" => "A0035",
"A0037" => "A0037");

my $path1 = "/rsgrps1/mbsulli/agregory/ABOR_50_viromes/viral_engraphment";
chdir($path1) or die "Could not chdir to $path1 \n";

my %metadata;
my $metadata_file = "metadata_4_engraphment.txt";
open(META, $metadata_file) or die "Could not open $metadata_file \n";
while(my $line = <META>) {
	chomp $line;
	my @spl = split(/\t/, $line);
	my $virome = $spl[0];
	my $Ldonor = $spl[1];
	my $Sdonor = $spl[2];
	my $value = $Ldonor."|".$Sdonor;
	$metadata{$virome}=$value;
}
close META;

my %virome_totals;
my %pop;
my $pop_coverage = "transposed_all_ABOR_pop_coverage.csv";
my $counta = 0;
open(POP, $pop_coverage) or die "Could not open $pop_coverage \n";
while(my $line = <POP>) {
	chomp $line;
	my @spl = split(/,/, $line);
	my $virome = $spl[0]; 
	$pop{$virome}=$line;
	my $total;
	foreach my $j (1..$#spl) {
		my $number = $spl[$j];
		$total += $number;
	}
	$virome_totals{$virome}=$total;
}
close POP;

my $outfile = "engraphment_percentages_based_on_abundances.txt";
open(OUTFILE, '>', $outfile) or die "Could not open $outfile \n";

print OUTFILE "Virome"."\t"."T1_Ldonor_Sdonor"."\t"."T1_Ldonor"."\t"."T1_Sdonor"."\t"."Ldonor_Sdonor"."\t"."T1_only"."\t"."Ldonor_only"."\t"."Sonor_only"."\t"."other"."\n";

foreach my $vir (keys %Viromes) {
	my $T1_virome = $vir.".ST001";
	my $T3_virome = $vir.".ST003";
	my $donors = $metadata{$T1_virome};
	my @donor_split = split(/\|/, $donors);
	my $Ldonor = $donor_split[0];
	my $Sdonor1 = $donor_split[1];
	my $Sdonora;
	my $Sdonorb;
	if ($Sdonor1 =~ /and/) {
		my @Sdonor_split = split(/and/, $Sdonor1);
		$Sdonora = $Sdonor_split[0];
		$Sdonorb = "D.2014.".$Sdonor_split[1];
	}
	else {
		$Sdonora = $Sdonor1;
		$Sdonorb = "NA";
	}
	my $T3_total = $virome_totals{$T3_virome};
	my $Sdonorb_total;
	if ($Sdonorb ne "NA") {
		$Sdonorb_total = $virome_totals{$Sdonorb};
	}	
	my %T1_coverage;
	my %T3_coverage;
	my %Ldonor_coverage;
	my %Sdonor_coverage;
	my $T1_coverage_line = $pop{$T1_virome};
	my $T3_coverage_line = $pop{$T3_virome};
	my $Ldonor_coverage_line = $pop{$Ldonor};
	my $Sdonora_coverage_line = $pop{$Sdonora};
	my $Sdonorb_coverage_line;
	if ($Sdonorb ne "NA") {
		$Sdonorb_coverage_line = $pop{$Sdonorb};
	}
	my @T1_spl = split(/,/, $T1_coverage_line);
	my $countT1 = 0;
	for my $i (1..$#T1_spl){
		$countT1++;
		my $coverage = $T1_spl[$i];
		if ($coverage > 0) {
			$T1_coverage{$countT1}=$coverage;
		}
	}
	my @T3_spl = split(/,/, $T3_coverage_line);
	my $countT3 = 0;
	for my $j (1..$#T3_spl){
		$countT3++;
		my $coverage = $T3_spl[$j];
		if ($coverage > 0) {
			$T3_coverage{$countT3}=$coverage;
		}
	}
	my @Ldonor_spl = split(/,/, $Ldonor_coverage_line);
	my $countL = 0;
	for my $k (1..$#Ldonor_spl){
		$countL++;
		my $coverage = $Ldonor_spl[$k];
		if ($coverage > 0) {
			$Ldonor_coverage{$countL}=$coverage;
		}
	}
	my @Sdonora_spl = split(/,/, $Sdonora_coverage_line);
	my $countSa = 0;
	for my $l (1..$#Sdonora_spl){
		$countSa++;
		my $coverage = $Sdonora_spl[$l];
		if ($coverage > 0) {
			$Sdonor_coverage{$countSa}=$coverage;
		}
	}
	my @Sdonorb_spl;
	if ($Sdonorb ne "NA") {
		@Sdonorb_spl = split(/,/, $Sdonorb_coverage_line);
		my $countSb = 0;
		for my $m (1..$#Sdonorb_spl){
			$countSb++;
			my $coverage = $Sdonorb_spl[$m];
			if ($coverage > 0) {
				$Sdonor_coverage{$countSb}=$coverage;
			}
		}
	}
	my $T1_Ldonor_Sdonor_number = 0;
	my $T1_unique_number = 0;
	my $Ldonor_unique_number = 0;
	my $Sdonor_unique_number = 0;
	my $T1_n_Ldonor_number = 0;
	my $T1_n_Sdonor_number = 0;
	my $Ldonor_n_Sdonor_number = 0;
	my $other_number = 0;
	foreach my $T3_number (keys %T3_coverage) {
		if ((exists $T1_coverage{$T3_number}) && (exists $Ldonor_coverage{$T3_number}) && (exists $Sdonor_coverage{$T3_number})) {
			my $cov = $T3_coverage{$T3_number};
			$T1_Ldonor_Sdonor_number+=$cov;
		}
		elsif ((exists $T1_coverage{$T3_number}) && (exists $Ldonor_coverage{$T3_number})) {
			my $cov = $T3_coverage{$T3_number};
			$T1_n_Ldonor_number+=$cov;
		}
		elsif ((exists $T1_coverage{$T3_number}) && (exists $Sdonor_coverage{$T3_number})) {
			my $cov = $T3_coverage{$T3_number};
			$T1_n_Sdonor_number+=$cov;
		}
		elsif ((exists $Ldonor_coverage{$T3_number}) && (exists $Sdonor_coverage{$T3_number})) {
			my $cov = $T3_coverage{$T3_number};
			$Ldonor_n_Sdonor_number+=$cov;
		}
		elsif (exists $T1_coverage{$T3_number}) {
			my $cov = $T3_coverage{$T3_number};
			$T1_unique_number+=$cov;
		}
		elsif (exists $Ldonor_coverage{$T3_number}) {
			my $cov = $T3_coverage{$T3_number};
			$Ldonor_unique_number+=$cov;
		}
		elsif (exists $Sdonor_coverage{$T3_number}) {
			my $cov = $T3_coverage{$T3_number};
			$Sdonor_unique_number+=$cov;
		}
		else {
			my $cov = $T3_coverage{$T3_number};
			$other_number+=$cov;
		}
	}
	my $present_in_T1_Ldonor_Sdonor = $T1_Ldonor_Sdonor_number/$T3_total;
	my $present_in_T1_Ldonor = $T1_n_Ldonor_number/$T3_total;
	my $present_in_T1_Sdonor = $T1_n_Sdonor_number/$T3_total;
	my $present_in_Ldonor_Sdonor = $Ldonor_n_Sdonor_number/$T3_total;
	my $present_in_T1_only = $T1_unique_number/$T3_total;
	my $present_in_Ldonor_only = $Ldonor_unique_number/$T3_total;
	my $present_in_Sdonor_only = $Sdonor_unique_number/$T3_total;
	my $present_in_other = $other_number/$T3_total;
	print OUTFILE $T3_virome."\t".$present_in_T1_Ldonor_Sdonor."\t".$present_in_T1_Ldonor."\t".$present_in_T1_Sdonor."\t".$present_in_Ldonor_Sdonor."\t".$present_in_T1_only."\t".$present_in_Ldonor_only."\t".$present_in_Sdonor_only."\t".$present_in_other."\n";
}
close OUTFILE;
	
