#!/usr/bin/perl -w

$n=0;
$seq="";
$name="";
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
	print "$n -- sequences/$name.fasta\n";
	open (O, "> sequences/$name.fasta");
	print O ">$name\n";
	print O "$seq";
	close O;
	}
    $seq="";
    $name=$nm;
    }
else {$seq.=$_}
}
close S;

$n++;
print "$n -- sequences/$name.fasta\n";
open (O, "> sequences/$name.fasta");
print O ">$name\n";
print O "$seq";
close O;
