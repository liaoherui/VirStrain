#!/usr/bin/perl -w
use strict;

# given an MSA of n sequences, conduct progressive clustering based on column entropy
# the default input format is an aligned multi fasta file (not clustal format)

(@ARGV == 2)|| die "Usage: <aln file>  <number of sequences> \n";
my ($file, $seqn)=@ARGV;


my @seq;
my @seqname;
my $count=0;
my $string="";
open(IN,$file)||die "cannot open $file";
while(<IN>)
{
    chomp();
    if($_=~/^>/) # line with file name
    {
	my @items =split(/\|/,$_);
	#push(@seqname,$items[0]);

        if($string ne "")
        {
            #print $seqname[$#seqname]," ", length($string), "\n";
	    #push(@seq,$string);
	    my @seqchar=split(//,$string);
	    push(@seq,\@seqchar);
            $string="";
        }
	push(@seqname,$items[0]);
       
    }
    else
    {
	$_=~tr/ACGTN/acgtn/;
        $string.=$_;
    }

}
close(IN);

if($string ne "")
{
            #print $seqname[$#seqname]," ", length($string), "\n";                                                                                   
            #push(@seq,$string);                                                                                                                     
    my @seqchar=split(//,$string);
    push(@seq,\@seqchar);
    #$string="";
}


print "finish inputting ", $#seq+1, " sequences from the alignment\n";

my @mattercolumns;

my $len=length($string); #the last string
my @entropy=0;
my %count;
my @countarray;
for(my $i=0; $i<$len; $i++)
{
    $count{'a'}=0;
    $count{'c'}=0;
    $count{'g'}=0;
    $count{'t'}=0;
    $count{'-'}=0;
    #$count{'n'}=0;

    for(my $j=0; $j<$seqn; $j++)
    {
	my $char = $seq[$j][$i];
	#print $char, "\n";
	if($char=~/[acgt-]/)
	{
	    $count{$char}++;
	}
	#else
	#{
	#    $count{$char}=1;
	#}
    }
    my $entropy=0;
    my $total =0;
    foreach (keys  %count)
    {
	$total+=$count{$_};
    }
    foreach (keys %count)
    {
	if($count{$_}>0)
	{
	    my $percent=$count{$_}/$total;
	    $entropy-=$percent * log($percent);
	}
    }
    #push(@countarray, \%count);
    if($entropy>0)
    {
	print "column ",$i," ";
	foreach (keys %count)
	{
	    print  $_,":",$count{$_}," ";
	}
	print $entropy,"\n";
    }
}





