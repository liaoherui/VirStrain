#!/usr/bin/perl -w
use strict;

# given an MSA of n sequences, conduct progressive clustering based on column entropy
# the default input format is an aligned multi fasta file (not clustal format)

(@ARGV == 6)|| die "Usage: <aln file>  <number of sequences> <allow -? if yes, set an upperbound of -, 0=no> <start> <end>  <verbose (0,1,2 increasing order of output)>\n";
my ($file, $seqn, $dash, $start, $end, $verbose)=@ARGV;

print "input arguments: $file, $seqn, $dash, $start, $end\n";


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
	my @items =split(/\s+/,$_);
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

#handle the last sequence in the multi-fasta file
if($string ne "")
{
    #print $seqname[$#seqname]," ", length($string), "\n";                                                                                   
    #push(@seq,$string);                                                                                                                     
    my @seqchar=split(//,$string);
    push(@seq,\@seqchar);
    #$string="";
}


print "finish inputting ", $#seq+1, " sequences from the alignment\n";

my @columnset; # saving columns with entropy > 0. If needed, you can shrink the original alignment into this column set (save space)

my $len=length($string); #the last string
my @entropy=0;
my %count;
my @countarray;

#for(my $i=0; $i<$len; $i++)
for(my $i=$start; $i<=$end; $i++)
{
    $count{'a'}=0;
    $count{'c'}=0;
    $count{'g'}=0;
    $count{'t'}=0;
    $count{'-'}=0;
    #$count{'n'}=0; # ignore n/N

    for(my $j=0; $j<$seqn; $j++)
    {
	my $char = $seq[$j][$i];
	#print $char, "\n";
	if($char=~/[acgt-]/)
	{
	    $count{$char}++;
	}
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

# Only use columns without -. If - is allowed, you can say if $count{'-'}<10/100 etc. (limit the number of -)
    #if($entropy>0 and $count{'-'}==0)
    if($entropy>0 and $count{'-'}<=$dash)
    {
	my %newcount; # MUST use a new hash in order to build an array of hash!!!
	if($verbose) {print "column ",$i," ";}
	foreach (keys %count)
	{
	    if ($verbose) {print  $_,":",$count{$_}," ";}
	    $newcount{$_}=$count{$_};
	}
	if ($verbose) {print $entropy,"\n";}

	push(@columnset, $i);
	push(@countarray, \%newcount);
	#print "a is ", $countarray[$#countarray]{'a'},"\n";
    }
}

if($verbose==2)
{ 
for(my $i=0; $i<5; $i++)
{
    print $columnset[$i]," ";
    my $test=$countarray[$i];
    foreach (keys %$test)
    {
	print $_, ":", $test->{$_},":", $countarray[$i]{$_}, " ";
       
    }
    print "\n";
}
}
#exit(0); 


#clusterG is the number of sequences inside clusters
my $clusteredseq="";
my $clusterG=0;

my $columncount=@columnset;
my $min=100;
my $mincol=-1;

while ($columncount>0 && $clusterG<$seqn)
{
    for(my $i=0; $i<@columnset; $i++)
    {
	if($columnset[$i]>=0)
	{
	    my $c=$countarray[$i]; #reference
	    my $total =0;
	    my $entry=0;
	    foreach (keys  %$c)
	    {
		$total+=$c->{$_};
	    }
	    foreach (keys %$c)
	    {
		if($c->{$_}>0)
		{
		    my $percent=$c->{$_}/$total;
		    $entry-=$percent * log($percent);
		}
	    }
	    
	    if ($entry >0 && $entry < $min) {$min = $entry;  $mincol=$i;}
	    if ($entry <= 0 ) {$columnset[$i]=-1; $columncount--; }
	}
    }   
    
    if($mincol>=0)
    {
	print "column ", $columnset[$mincol]," ";

	my $markcol=$columnset[$mincol];
	#$columnset[$mincol]=-1; # NO soft removing this column from columnset, allow reusing this column 
	#$columncount--;

	my @lineset;
	my $mincharct=$seqn;
	my $minchar="";
	my $cnt=$countarray[$mincol];
	foreach (keys %$cnt)
	{
	    if($cnt->{$_} >0 && $cnt->{$_} < $mincharct) {$minchar=$_; $mincharct=$cnt->{$_};} 
	    print $_,"(",$cnt->{$_},") ";
	}	
	print "\n";

	if($verbose) {print " min char is ",$minchar,"\n";}

	for(my $cx=0; $cx<$seqn; $cx++)
	{
	    if($seq[$cx][$markcol] eq $minchar) 
	    {
		if($clusteredseq!~/$seqname[$cx]/)
		    {push(@lineset, $cx);} 
		print $seqname[$cx]," ";
		if($clusteredseq!~/$seqname[$cx]/) {$clusteredseq.=" $seqname[$cx]"; $clusterG++;}
	    }
	}
	if(@lineset>=0) {print $markcol," ", $minchar,"\n"};


	for(my $lx=0; $lx<@lineset; $lx++)
	{
	    for(my $ly=0; $ly<@columnset; $ly++)
	    {
		if($verbose==3) {print "column: ", $columnset[$ly]," ";}
	
		if($columnset[$ly]==-1) {next;}
		my $cur=$seq[$lineset[$lx]][$columnset[$ly]];
		if($verbose==3) {print $cur, " ";}
		if($cur=~/[acgt-]/){
		    if($verbose==3) {print "countarray[",$columnset[$ly],"]->",$cur," is ", $countarray[$ly]{$cur},"\n";}
		    $countarray[$ly]{$cur}--;
		    if($verbose==3) {print "countarray[",$columnset[$ly],"]->",$cur," is ", $countarray[$ly]{$cur},"\n";}
		}
	    }
	}

    }
    else #mincol<0 
    {
	$columncount=-1;

	if($verbose==1)
	{
	    for(my $i=0; $i<=$#countarray; $i++)
	    {
		print $columnset[$i]," ";
		my $test=$countarray[$i];
		foreach (keys %$test)
		{
		    print $_, ":", $countarray[$i]{$_}, " ";
		    
		}
		print "\n";
	    }
	}


    }
    $min=100;
    $mincol=-1;
	
}
    #CHECK whether the columnset is used up


print "The sequences in clusters are ", $clusterG, "\n";

if($verbose==2)
{
    print "The remaining columns with mutations: \n";
    for(my $i=0; $i<@columnset; $i++)
    {
	if($columnset[$i]>0) {print $columnset[$i]," ";}
    }
    print "\n";
    
    print $clusteredseq,"\n";
    my @seqlist=split(/\s+/,$clusteredseq);
    print "There are ", $#seqlist+1, " sequences in clusters\n";
}






