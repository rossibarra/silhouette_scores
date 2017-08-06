#!/usr/bin/perl


use warnings;
use strict;
use Data::Dumper;
use Getopt::Std;
use List::Util qw(sum min max);
use List::MoreUtils qw(first_index);
use bignum ( p => -5 );
use Clone qw(clone);


my %flags;
getopts( "i:o:r:p:gd:e:l:f:", \%flags );
unless ( defined $flags{i} && $flags{o} )
{        print "\nThis script uses any distance matrix (without headers) with population assignation for each individual to calculate silhouette scores of their fit to their population. Matrix has to be symmetric (column 1: ID\@POP, col. 2-: IBS. This step is important before doing any admixture tests\n\nUSAGE: silhoutte.pl -i input -o output [options]\n\n";
        print "\t-i\t [FILE] input file\n";
        print "\t-o\t [FILE] output file prefix\n";
	print "\t-r\t [STRING] comma separated list of populations one would like to remove\n";
	print "\t-p\t [STRING] comma separated list of populations one would like to pool together\n";
	print "\t-g\t [BOOLEAN] switch if you would like to employ associated R script to produce graphics\n";
	print "\t-d\t [INTEGER] define fraction of indviduals knived in all populations efficiency\n";
	print "\t-e\t [INTEGER] define fraction of indviduals knived for single population at a time diagnostic\n";
	print "\t-l\t [INTEGER] define cutoff score for indviduals you would like to keep\n";
	print "\t-f\t [INTEGER] specify fraction of lowest score individuals per population to be removed\n\n\n";

        exit ;
}

	##### FLAGS INIT

my $file = $flags{i};
my $outfile = $flags{o};


my @remove;
if($flags{r})
{
	push @remove , split ",",$flags{r};	
}

my @pool;
if($flags{p})
{
        @pool = split ",",$flags{p};
}

if($flags{e})
{
	my $diagno = $flags{e};
	open (DIAOUT, ">${outfile}_diagnostics.sil") or die "Cannot open diagnostic output\n";
}

my $lcutoff = 0;

if($flags{l})
{
        open (LISTOUT, ">${outfile}_list.sil") or die "Cannot open diagnostic output\n";
	$lcutoff = $flags{l};
}

	##### R SCRIPT 

#	library(gridExtra)
#	arg <-commandArgs()
#	ind <- read.table(paste(arg[6],"_ind.sil",sep=""), header = F)
#	pop <- read.table(paste(arg[6],"_pop.sil",sep=""), header = F)
#	palette <- sample(colours(), length(table(ind$V1)))
#	cols <- c()
#	for (i in 1:length(table(ind$V1)))
#	{
#	        cols <- c(cols, rep(palette[i],table(ind$V1)[i]))
#	}
#	pdf(paste(arg[6],"_ind.pdf",sep=""))
#	barplot(ind$V2, col = cols, border = NA)
#	legend("topright", legend = names(table(ind$V1)), cex = 0.5, col = palette, pch = 15,  box.lwd = 0)
#	dev.off()
#	pdf(paste(arg[6],"_pop.pdf",sep=""))
#	barplot(pop$V2, col = palette, border = NA)
#	legend("topright", legend = pop$V1, cex = 0.5, col = palette, pch = 15,  box.lwd = 0)
#	plot(NA, ylab = "", xlab = "", ylim = c(0,1), xlim = c(0,1),xaxt = 'n', yaxt = 'n', frame = F)
#	grid.table(pop)
#	dev.off()


	##### FILE HANDLES INIT


open(INFILE, "$file") or die "Cannot open input";
open(INDOUT, "| sort -k 1,1 -k 2,2n >${outfile}_ind.sil") or die "Cannot open output";
open(POPOUT, ">${outfile}_pop.sil") or die "Cannot open output";


	##### GIVE INDIVIDUALS INDEX NUMBER SEGREGATE THEM INTO POPULATIONS #####


my $index1 = 0;
my %pop_list;	# [hash: keys - populations, values - [hash: keys - indexes of individuals, values - IDB matrix line and ]]

while(<INFILE>)
{
	chomp;
	my @line1 = split "\t", $_;
	my @pop1 = split "\@", $line1[0];
	${$pop_list{$pop1[1]}}{$index1} = $_;
	$index1++;
}
	# in this set-up indexing starts from 0, conveniently IBD values also start from position 0 (column 1) of @IBDs (dowstream the code)


	##### APPLY REMOVE AND POOL OPTIONS

if($flags{r})
{
	foreach my $rm (@remove)
	{
		if ($pop_list{$rm})
		{
			delete $pop_list{$rm};
		} else {
			die "\nERROR\tPopulation you want to remove doesn't exist\n\n";
		}
	}
}

if($flags{p})
{
	foreach my $pl (@pool[1 .. $#pool])
	{
		if ($pop_list{$pl} && $pop_list{$pool[0]})
		{
			foreach (keys %{$pop_list{$pl}})
			{
				${$pop_list{$pool[0]}}{$_} = ${$pop_list{$pl}}{$_};
			}
			delete $pop_list{$pl};
		} else {
			die "\nERROR\tPopulation you want to pool doesn't exist\n\n";
		}
	}
}


#print Dumper \%pop_list;

print "\nDONE\tanalyzing population assignation\n";




	##### CALCULATE INDIVIDUAL SCORES #####


my %pop_scores;	# [hash: keys - populations, values - [hash: keys - indexes of individuals, values - silhouette scores ]]
my %nemesis;	# [hash: keys - populations, values - [array of foreign population names that current individual is closes to]]

foreach my $pop1 (keys %pop_list)
{
	foreach my $ind (sort keys %{$pop_list{$pop1}})
	{
		my @line2 = split "\t",${$pop_list{$pop1}}{$ind};		
		# JRI: line changed here and below to allow all tab-delimited input files
		#my @IBDs = split "\t", $line2[1]; 
		my @IBDs = @line2[1..$#line2];
		
		### getting average distance to own population
		my $with_sum = sum @IBDs[keys %{$pop_list{$pop1}}];
		my $with_mean = $with_sum / (scalar (keys %{$pop_list{$pop1}}) - 1);

		### getting average distance to all other populations
		my %out_means;
		foreach my $pop2 (keys%pop_list)
		{
			next if($pop2 eq $pop1);
			my $sum = sum @IBDs[keys %{$pop_list{$pop2}}];
			my $mean = $sum / scalar keys %{$pop_list{$pop2}};
			$out_means{$pop2} = $mean;
		}
		
		### Save closest foreign population average distance
		my @sorted = sort { $a <=> $b } values %out_means;

		### Calculating Individual Silhouette score
		my $score;
		if ($with_mean > $sorted[0])
		{
			$score = (($sorted[0]-$with_mean)/$with_mean);
		} else {
			$score = (($sorted[0]-$with_mean)/$sorted[0]);
		}

		${$pop_scores{$pop1}}{$ind} = $score;	
		
		#### find which population is closest to this individual
		my @sorted_keys = sort { $out_means{$a} <=> $out_means{$b} } keys %out_means;
		push @{$nemesis{$pop1}},  $sorted_keys[0];

		if ($flags{d})
		{
			next;
		} else {	
			print INDOUT "$pop1\t$score\n" if ($pop1 ne "NA");
			print LISTOUT "$line2[0]\n" if ($flags{l} && $score > $lcutoff);
		}
	}
}

print "\nDONE\tproducing scores for individuals\n";


	##### OPTION: DELETE FRACTION OF INDVIDUALS FROM ALL POPULATIONS


if ($flags{d})
{
	my $cutoff = $flags{d};		
	# how does global cut-off of indviduals influences all population - MARKER FOR GLOBAL FILTERING
	foreach my $pop0 (keys%pop_scores)
	{
		my $tot = scalar keys %{$pop_scores{$pop0}};
		my $frac = $cutoff;
		my $no_rm = int ($tot * $frac +0.5) -1;
		
		my @pop_srt = sort { ${$pop_scores{$pop0}}{$a} <=> ${$pop_scores{$pop0}}{$b} } keys %{$pop_scores{$pop0}};	
		my @pop_rm = @pop_srt[0 .. $no_rm];
		foreach my $pr (@pop_rm)
		{
			delete ${$pop_list{$pop0}}{$pr};
		}
	}
	print "\nDONE\tdeleting $flags{d} weakest individuals\n";

	foreach my $pop1 (keys %pop_list)
	{
		foreach my $ind (sort keys %{$pop_list{$pop1}})
		{
			my @line2 = split "\t",${$pop_list{$pop1}}{$ind};		
			my @IBDs = @line2[1..$#line2];
	
			### getting average distance to own population
			my $with_sum = sum @IBDs[keys %{$pop_list{$pop1}}];
			my $with_mean = $with_sum / (scalar (keys %{$pop_list{$pop1}}) - 1);
	
			### getting average distance to all other populations
			my %out_means;
			foreach my $pop2 (keys%pop_list)
			{
				next if($pop2 eq $pop1);
				my $sum = sum @IBDs[keys %{$pop_list{$pop2}}];
				my $mean = $sum / scalar keys %{$pop_list{$pop2}};
				$out_means{$pop2} = $mean;
			}
			
			### Save closest foreign population average distance
			my @sorted = sort { $a <=> $b } values %out_means;
	
			### Calculating Individual Silhouette score
			my $score;
			if ($with_mean > $sorted[0])
			{
				$score = (($sorted[0]-$with_mean)/$with_mean);
			} else {
				$score = (($sorted[0]-$with_mean)/$sorted[0]);
			}
	
			${$pop_scores{$pop1}}{$ind} = $score;	
			
			#### find which population is closest to this individual
			my @sorted_keys = sort { $out_means{$a} <=> $out_means{$b} } keys %out_means;
			push @{$nemesis{$pop1}},  $sorted_keys[0];
	
			print INDOUT "$pop1\t$score\n" if ($pop1 ne "NA");
			print LISTOUT "$line2[0]\n" if ($flags{l} && $score > $lcutoff);
		}
	}

}


	##### CALCULATE POPULATION SILHOUETTE SCORES AND THEIR NEMSIS

#print Dumper \%nemesis;


foreach my $pop3 (sort keys%pop_scores)
{
	### calculate population silhouette scores
	my $out_sum = sum values %{$pop_scores{$pop3}};
	my $out_mean = $out_sum / (scalar keys %{$pop_scores{$pop3}});
	

	### calculate most frequent nemesis
	my %counts;
	foreach my $tt (@{$nemesis{$pop3}})
	{
		$counts{$tt}++;
	}
	my $freq = (max values %counts) / (sum values %counts);	
	my $nem = (sort {$counts{$a} <=> $counts{$b}} keys %counts)[0];
	if ($pop3 ne "NA")
	{
		print POPOUT "$pop3\t$out_mean\t${nem}\t${freq}\t";
		print POPOUT "\n";
	}
}


print "\nDONE\tproducing scores for populations\n";


	##### OPTION FOR DIAGNOSTICS

if ($flags{e})
{
	my $cutoff = $flags{e};

	print DIAOUT "ID";
	foreach my $pop4 (sort keys%pop_scores)
        {
		print DIAOUT "\t$pop4";
	}
	print DIAOUT "\n";
		# how does global cut-off of indviduals influences all population - MARKER FOR GLOBAL FILTERING
	foreach my $pop4 (sort keys%pop_scores)
	{
		print DIAOUT "$pop4\t";
		my $tot = scalar keys %{$pop_scores{$pop4}};
		my $frac = $cutoff;
		my $no_rm = int ($tot * $frac +0.5) -1;
		
		my @pop_srt = sort { ${$pop_scores{$pop4}}{$a} <=> ${$pop_scores{$pop4}}{$b} } keys %{$pop_scores{$pop4}};	
		my @pop_rm = @pop_srt[0 .. $no_rm];

		my %rm_list =  %{clone(\%pop_list)};
		my %rm_scores = %{ clone(\%pop_scores)};
		foreach my $pr (@pop_rm)
		{
			delete ${$rm_list{$pop4}}{$pr};
			delete ${$rm_scores{$pop4}}{$pr};
		}
		# re-calculate silhouette scores
 
		foreach my $pop1 (sort keys %pop_list)
		{
	       		foreach my $ind (keys %{$rm_list{$pop1}})
        		{
        		        my @line2 = split "\t",${$rm_list{$pop1}}{$ind};
				my @IBDs = @line2[1..$#line2];

        		        ### getting average distance to own population
        		        my $with_sum = sum @IBDs[keys %{$rm_list{$pop1}}];
        		        my $with_mean = $with_sum / (scalar (keys %{$rm_list{$pop1}}) - 1);

        		        ### getting average distance to all other populations
        		        my %out_means;
        		        foreach my $pop2 (keys%rm_list)
        		        {
        		                next if($pop2 eq $pop1);
        		                my $sum = sum @IBDs[keys %{$rm_list{$pop2}}];
        		                my $mean = $sum / scalar keys %{$rm_list{$pop2}};
        		                $out_means{$pop2} = $mean;
        		        }

        		        ### Save closest foreign population average distance
        		        my @sorted = sort { $a <=> $b } values %out_means;

        		        ### Calculating Individual Silhouette score
        		        my $score;
        		        if ($with_mean > $sorted[0])
        		        {
        		                $score = (($sorted[0]-$with_mean)/$with_mean);
        		        } else {
        		                $score = (($sorted[0]-$with_mean)/$sorted[0]);
        		        }

        		        ${$rm_scores{$pop1}}{$ind} = $score;		
			}
			my $rm_sum = sum values %{$rm_scores{$pop1}};
			my $rm_mean = $rm_sum / (scalar keys %{$rm_scores{$pop1}});
			my $out_sum = sum values %{$pop_scores{$pop1}};
			my $out_mean = $out_sum / (scalar keys %{$pop_scores{$pop1}});
			print DIAOUT $rm_mean - $out_mean, "\t";
		}
		print DIAOUT "\n";
	}
}
close INFILE;
close INDOUT;
close POPOUT;
close DIAOUT if ($flags{e});

if($flags{g})
{
	my $graph = "Rscript silhouette.r ${outfile}";
	my $warn = system $graph;
	if ($warn) {die "\nERROR\tDied: Rscript failed\n\n";}
	else { print "\nDONE\tproducing graphs\n\n";}
} 

##### DONE #####

__END__
