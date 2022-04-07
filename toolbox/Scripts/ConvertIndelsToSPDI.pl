#!/usr/local/bin/perl

use warnings;
use strict;

use SPDI;

if (@ARGV < 1) {
    print "Description:\n";
    print "  This script reads indels from the input file and convert them into those with SPDI format.\n";
    print "  Results will be saved the following three files:\n";
    print "  1. file with \"_spdi\" appended to the name: same as the input file with indels replaced with SPDI indels\n";
    print "  2. file with \"_sum\" appended to the name: summary of all the distinct indel replaced\n";
    print "  3. file with \"_err\" appended to the name: shows all indels that cannot be converted\n";
    print "Usage: ConvertIndelsToSPDI.pl <input_file>\n";
    print "Requirements:\n";
    print "  1. The input file is a plain text .csv file or .txt (tab-delimited) file.\n";
    print "  2. The input file is a table with the following three columns: pos (or position), ref, and alt.\n";
    print "     It is acceptable it includes other columns. But the above three columns should be unique.\n";
    print "  3. File SPDI.pm is found under the same directory where the script is located.\n";;
    print "\n";
    exit 0;
}

my $inFile = $ARGV[0];

ConvertSpdiForFile($inFile);

# Convert all indels in the input file to SPDI and saves results to output files
sub ConvertSpdiForFile
{
    my ($inFile) = @_;

    my $spdi = new SPDI();

    my ($baseName, $ext, $delim) = ("", "", "");
    if ($inFile =~ /(.+)\.(csv)$/i) {
	$baseName = $1;
	$ext = $2;
	$delim = "\,";
    }
    elsif ($inFile =~ /(.+)\.(txt)$/i) {
	$baseName = $1;
	$ext = $2;
	$delim = "\t";
    }
    else {
	die "\nERROR: input file should be either a .csv file or a .txt file.\n";
    }

    my $outFile = $baseName . "_spdi." . $ext;

    open FILE, $inFile or die "Couldn't open $inFile\n";
    my $header = <FILE>;
    chomp $header;
    $header =~ s/position/pos/ig;

    # Find the pos, ref and alt colummns
    my ($posCol, $refCol, $altCol) = (-1, -1, -1);
    my @vars = split /$delim/, $header;
    for my $col (0 .. $#vars) {
	my $var = $vars[$col];
	if ($var =~ /^pos/i) {
	    $posCol = $col;
	}
	if ($var =~ /^ref$/i) {
	    $refCol = $col;
	}
	if ($var =~ /^alt$/i) {
	    $altCol = $col;
	}
    }

    die "\nERROR: didn't find pos column in file $inFile\n" if ($posCol < 0);
    die "\nERROR: didn't find ref column in file $inFile\n" if ($refCol < 0);
    die "\nERROR: didn't find alt column in file $inFile\n" if ($altCol < 0);

    # Read all lines of file into memory. Save all mutations;
    my @allLines = ();
    my @hasIndels = ();
    my @linePoss = ();
    my @lineRefs = ();
    my @lineAlts = ();

    my %allMuts = ();
    my %allLineMuts = ();
    my %allIndels = ();
    my %allLineIndels = ();

    my $lineNo = 0;
    while(<FILE>) {
	chomp;
	my $line = $_;
	$line =~ s/\s*$//;

	my @vals = split /$delim/, $line;
	next if ($line !~ /\w/);

	my $hasIndel = 0;
	my ($pos, $ref, $alt) = ("", "", "");

	$pos = $vals[$posCol] if ($posCol <= $#vals);
	$ref = $vals[$refCol] if ($refCol <= $#vals);
	$alt = $vals[$altCol] if ($altCol <= $#vals);
	$ref = "" if ($ref eq "-");
	$alt = "" if ($alt eq "-");

	$pos = 0 if ($pos !~ /^\d+$/);

	my $mut = "$pos\t$ref\t$alt";

	$allMuts{"$pos\t$ref\t$alt"} = 1;
	$allLineMuts{"$lineNo\t$pos\t$ref\t$alt"} = 1;
	my $refLen = length($ref);
	my $altLen = length($alt);

	if ($refLen != 1 || $altLen != 1) {
	    $allIndels{"$pos\t$ref\t$alt"} = 1;
	    $allLineIndels{"$lineNo\t$pos\t$ref\t$alt"} = 1;
	    $hasIndel = 1;
	}

	push @allLines, $_;
	push @hasIndels, $hasIndel;
	push @linePoss, $pos;
	push @lineRefs, $ref;
	push @lineAlts, $alt;

	$lineNo++;
    }
    close FILE;

    my $numMuts = keys %allMuts;
    my $numLineMuts = keys %allLineMuts;
    my $numIndels = keys %allIndels;
    my $numLineIndels = keys %allLineIndels;

    print "\nFile has $lineNo rows.\n";
    print "\tFound $numMuts distinct mutations from $numLineMuts lines\n";
    print "\tFound $numIndels distinct indels from $numLineIndels lines with indels\n";

    my @saveLines = (); # lines to be saved to output file

    my $flankLen = 2;
    my %allNewIndels = ();
    my %wrongMuts = ();

    for my $lineNo (0 .. $#allLines) {
	if ($hasIndels[$lineNo]) {
	    my $oPos = $linePoss[$lineNo];
	    my $oRef = $lineRefs[$lineNo];
	    my $oAlt = $lineAlts[$lineNo];

	    if (!$oPos || (!$oRef && !$oAlt)) {
		push @saveLines, $allLines[$lineNo];
		next;
	    }

	    my $flag = 0;
	    my ($pos, $ref, $alt, $err) = $spdi->ConvertIndelToSPDI($oPos, $oRef, $oAlt, $flag);

	    $wrongMuts{"$oPos\t$oRef\t$oAlt"} = $err if ($err);
	    unless ($pos) {
		push @saveLines, $allLines[$lineNo];
		next;
	    }

	    my @vals = split /$delim/, $allLines[$lineNo];
	    my $newLine = "";
	    for my $col (0 .. $#vals) {
		my $val = $vals[$col];
		if ($col == $posCol) {
		    $val = $pos;
		}
		elsif ($col == $refCol) {
		    $val = $ref ? $ref : "-";
		}
		elsif ($col == $altCol) {
		    $val = $alt ? $alt : "-";
		}
		$newLine .= $val;
		$newLine .= "$delim" if ($col < $#vals);
	    }
	    push @saveLines, $newLine;

	    # Show flanking sequences for debugging and testing
	    my $lPos = $oPos < $pos ? $oPos : $pos;;
	    my $fLen = length($oRef);
	    $fLen = length($oAlt) if (length($oAlt) > $fLen);
	    $fLen = length($ref)  if (length($ref)  > $fLen);
	    $fLen = length($alt)  if (length($alt)  > $fLen);

	    my $rPos = $lPos + $fLen + $flankLen;
	    $lPos -= $flankLen;

	    my $flankSeq = "";
	    for my $i ($lPos .. $rPos) {
		my $nt = " ";
		if ($i > 0 && $i+1 <= $spdi->{_genome_len}) {
		    $nt = $spdi->{_ref_nucs}->[$i-1];
		}
		$flankSeq .= $nt;
	    }

	    # Save all converted indels
	    if ($allNewIndels{"$pos\t$ref\t$alt"}) {
		$allNewIndels{"$pos\t$ref\t$alt"}->{"$oPos\t$oRef\t$oAlt"} = 1;
	    }
	    else {
		$allNewIndels{"$pos\t$ref\t$alt"} = {"$oPos\t$oRef\t$oAlt" => 1};
	    }
	}
	else {
	    push @saveLines, $allLines[$lineNo];
	}

	print "\tChecked $lineNo lines\n" if ($lineNo > 0 && $lineNo % 10000 == 0);
    }

    my $numNewIndels = keys %allNewIndels;
    print "\nConverted $numIndels indels to $numNewIndels SPDI indels\n";

    # Save indels that cannot be converted to "_err" file
    my $numWrongMuts = keys %wrongMuts;
    if ($numWrongMuts > 0) {
	my $errFile = $baseName . "_err.txt";
	open FILE, ">$errFile" or die "Couldn't open $errFile for writing\n";
	print FILE "pos\tref\talt\terror\n";
	foreach my $mut (keys %wrongMuts) {
	    my $err = $wrongMuts{$mut};
	    my $showErr = $err eq "SNP" ? "=> $err" : "";
	    $showErr = $err;
	    print FILE "$mut\t$showErr\n";
	}
    }
    close FILE;

    # Save summary to "_sum" file
    my $sumFile = $baseName . "_sum.txt";

    open FILE, ">$sumFile" or die "Couldn't open $sumFile for writing\n";
    print FILE "pos\tref\talt\tnew_pos\tnew_ref\tnew_alt\n";
    my $numNewInds = keys %allNewIndels;
    print "\nTotal $numNewInds indels are converted\n\n";
    foreach my $ind (sort keys %allNewIndels) {
	my ($pos, $ref, $alt) = split /\t/, $ind;
	$ref = "-" unless ($ref);
	$alt = "-" unless ($alt);
	my $lineTail = "$pos\t$ref\t$alt";
	my %oldInds = %{$allNewIndels{$ind}};
	foreach my $oldInd (keys %oldInds) {
	    my ($opos, $oref, $oalt) = split /\t/, $oldInd;
	    my $line = "$opos\t$oref\t$oalt\t$lineTail";
	    print FILE "$line\n" if ($pos != $opos || $ref ne $oref || $alt ne $oalt);
	}
    }
    close FILE;
    print "Summary saved to $sumFile\n";

    open FILE, ">$outFile" or die "Couldn't open $outFile for writing\n";
    print FILE "$header\n";
    for my $line (@saveLines) {
	print FILE "$line\n";
    }
    close FILE;

    print "Converted file saved to $outFile\n\n";
}
