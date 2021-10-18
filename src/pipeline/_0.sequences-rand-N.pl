#!/mnt/lustre/genkvg/perl/bin/perl -w
use Math::Random::MT::Auto qw(rand shuffle);
#USAGE: ./0.sequences-rand-N.pl bigMultiFasta numberOfSequencesExtract > multiFastaSequesteredAndRandomized

$n=0;
$seq="";
$name="";
%sequences=();
open (S, "$ARGV[0]");
while (<S>)
{
if ($_=~/>([^\|]+)\|/)
    {
    $nm=$1;
    $nm=~s/\//-/g;
    $nm=~s/\,|\:|\(|\)/_/g;
    if ($name ne "")
	{
	$n++;
	$sequences{">$n.$name"}=$seq;
	print STDERR "$n -- $name\n";
	}
    $seq="";
    $name=$nm;
    }
else {$seq.=$_}
}
close S;

$n++;
print STDERR "$n -- $name\n";
$sequences{">$n.$name"}=$seq;

$it=0;
@sequences=keys %sequences;
for ($i=0;$i<1000;$i++)
{
$it++;
print STDERR "shuffling $it\r";
@sequencesshuffled = shuffle(@sequences);
@sequences=@sequencesshuffled;
}

for ($i=0;$i<$ARGV[1];$i++)
{
print "$sequences[$i]\n$sequences{$sequences[$i]}";
}
