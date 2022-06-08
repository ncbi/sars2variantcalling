#!/usr/local/bin/perl

use strict;
use warnings 'all';
use Test::More 'no_plan';

use_ok('SPDI');
my $spdi = new SPDI();

ok($spdi, 'Can create instance of SPDI');

is( $spdi->{_genome_len} => 29903, 'Instance is correct');

# testing TrimLeft
Test_TrimLeft("AGCT", "AG",   "CT", "",   2, "TrimLeft test 1");
Test_TrimLeft("AGCT", "AGTT", "CT", "TT", 2, "TrimLeft test 2");
Test_TrimLeft("AGCT", "AGCA", "T",  "A",  3, "TrimLeft test 3");
Test_TrimLeft("CAT",  "CAT",  "",   "",   3, "TrimLeft test 4");
Test_TrimLeft("CA",   "GT",   "CA", "GT", 0, "TrimLeft test 5");
Test_TrimLeft("",     "AG",   "",   "AG", 0, "TrimLeft test 6");
Test_TrimLeft("AA",   "",     "AA", "",   0, "TrimLeft test 7");

# testing TrimRight
Test_TrimRight("CT",   "AGCT",  "",  "AG", 2, "TrimRight test 1");
Test_TrimRight("GACT", "AGCT", "GA", "AG", 2, "TrimRight test 2");
Test_TrimRight("AT",   "G",    "AT", "G",  0, "TrimRight test 3");
Test_TrimRight("CAT",  "CAT",  "",   "",   3, "TrimRight test 4");
Test_TrimRight("",     "AG",   "",   "AG", 0, "TrimRight test 5");
Test_TrimRight("AA",   "",     "AA", "",   0, "TrimRight test 6");

# testing Blossom
Test_Blossom(14847, "ACT",  14844, "ACTACTA",     "ACTA",     "Blossom test 1");
Test_Blossom(21884, "ACT",  21883, "TACTACT",     "TACT",     "Blossom test 2");
Test_Blossom(4741 , "TTC",  4737,  "CTTCTTCTTCT", "CTTCTTCT", "Blossom test 3");
Test_Blossom(9431 , "GTA",  9431,  "GTAGTA",      "GTA",      "Blossom test 4");
Test_Blossom(3377 , "GTA",  3377,  "GTA",         "",         "Blossom test 5");
Test_Blossom(3918 , "GTA",  3917,  "AGTA",        "A",        "Blossom test 6");
Test_Blossom(4868,  "GTA",  4868,  "GTA",         "",         "Blossom test 7");
Test_Blossom(4021 , "C",    4021,  "C",           "",         "Blossom test 8");
Test_Blossom(4021,  "A",    4017,  "AAAAA",       "AAAA",     "Blossom test 9");
Test_Blossom(5226 , "A",    5225,  "AAAAA",       "AAAA",     "Blossom test 10");
Test_Blossom(3481 , "A",    3481,  "A",           "",         "Blossom test 11");
Test_Blossom(3481 , "G",    3481,  "GGG",         "GG",       "Blossom test 12");

# testing ConvertIndelToSPDI
Test_ConvertIndelToSPDI(3301,  "AGA",   "AG",   3302,  "GA",       "G",       "",    "ConvertIndel test 1");
Test_ConvertIndelToSPDI(3481,  "A",     "AG",   3482,  "G",        "GG",      "",    "ConvertIndel test 2");
Test_ConvertIndelToSPDI(1001,  "",      "A",    1001,  "AAAA",     "AAAAA",   "",    "ConvertIndel test 3");
Test_ConvertIndelToSPDI(1001,  "-",     "A",    1001,  "AAAA",     "AAAAA",   "",    "ConvertIndel test 4");
Test_ConvertIndelToSPDI(11074, "CT",    "C",    11075, "TTTTTTTT", "TTTTTTT", "",    "ConvertIndel test 5");
Test_ConvertIndelToSPDI(11082, "T",     "",     11075, "TTTTTTTT", "TTTTTTT", "",    "ConvertIndel test 6");
Test_ConvertIndelToSPDI(14747, "A",     "",     14747, "AA",       "A",       "",    "ConvertIndel test 7");
Test_ConvertIndelToSPDI(14746, "GAATT", "GATT", 14747, "AA",       "A",       "",    "ConvertIndel test 8");
Test_ConvertIndelToSPDI(14746, "GAA",   "GA",   14747, "AA",       "A",       "",    "ConvertIndel test 9");
Test_ConvertIndelToSPDI(14843, "AACT",  "A",    14844, "ACTACTA",  "ACTA",    "",    "ConvertIndel test 10");
Test_ConvertIndelToSPDI(21881, "GGTAC", "GG",   21883, "TACTACT",  "TACT",    "",    "ConvertIndel test 11");
Test_ConvertIndelToSPDI(9430,  "C",     "CGTA", 9431,  "GTA",      "GTAGTA",  "",    "ConvertIndel test 12");
Test_ConvertIndelToSPDI(12297, "",      "A",    12297, "AA",       "AAA",     "",    "ConvertIndel test 13");
Test_ConvertIndelToSPDI(11696, "TATGA", "TA",   11697, "ATGAT",    "AT",      "",    "ConvertIndel test 14");
Test_ConvertIndelToSPDI(11696, "TATGAT","TA",   11698, "TGATT",    "T",       "",    "ConvertIndel test 15");
Test_ConvertIndelToSPDI(18668, "GT",    "G",    18668, "GT",       "G",       "",    "ConvertIndel test 16");
Test_ConvertIndelToSPDI(1,     "ATT",   "CATT", 1,     "",         "C",       "",    "ConvertIndel test 17");
Test_ConvertIndelToSPDI(12297, "-",     "A",    12297, "AA",       "AAA",     "",    "ConvertIndel test 18");
Test_ConvertIndelToSPDI(18669, "T",     "-",    18668, "GT",       "G",       "",    "ConvertIndel test 19");
Test_ConvertIndelToSPDI(13200, "CGG",   "TGG",  13200, "C",        "T",       "SNP", "ConvertIndel test 20");
Test_ConvertIndelToSPDI(13200, "CGG",   "CAG",  13201, "G",        "A",       "SNP", "ConvertIndel test 21");
Test_ConvertIndelToSPDI(13201, "G",     "A",    13201, "G",        "A",       "SNP", "ConvertIndel test 22");
Test_ConvertIndelToSPDI(13201, "G",     "G",    13201, "G",        "G",       "SNP", "ConvertIndel test 23");

# Two cases from https://jira.ncbi.nlm.nih.gov/browse/SRAB-396
Test_ConvertIndelToSPDI(11287, "GTCTGGTTTT",   "G", 11287, "GTCTGGTTTT", "G", "", "ConvertIndel test 24");
Test_ConvertIndelToSPDI(21764, "ATACATG",      "A", 21765, "TACATGT",    "T", "", "ConvertIndel test 25");

# testing ConvertIndelToSPDI: unconvertable inputs
Test_ConvertIndelToSPDI(1001, "G",   "B",   0, "",     "",   "Invalid alt 'B'.", "ConvertIndel test 25");
Test_ConvertIndelToSPDI(1001, "9",   "A",   0, "",     "",   "Invalid ref '9'.", "ConvertIndel test 26");
Test_ConvertIndelToSPDI(1001, "G",   "AaT", 0, "",     "",   "Invalid alt 'AaT'.", "ConvertIndel test 27");
Test_ConvertIndelToSPDI(1001, "T",   "A",   0, "",     "",   "Wrong ref 'T' at pos 1001, should be 'G'.", "ConvertIndel test 28");
Test_ConvertIndelToSPDI(1001, "",    "",    0, "",     "",   "Empty mutation", "ConvertIndel test 29");
Test_ConvertIndelToSPDI(1001, "GAA", "CTT", 0, "GAA", "CTT", "Invalid mutation", "ConvertIndel test 30");
Test_ConvertIndelToSPDI(1001, "GAA", "CA",  0, "GA",  "C",   "Invalid mutation", "ConvertIndel test 31");
Test_ConvertIndelToSPDI(1001, "GA",  "CAT", 0, "GA",  "CAT", "Invalid mutation", "ConvertIndel test 32");
Test_ConvertIndelToSPDI(99999,"G",   "A",   0, "",    "",    "Invalid pos 99999", "ConvertIndel test 33");


sub Test_TrimLeft
{
    my ($oldRef, $oldAlt, $expRef, $expAlt, $expPos, $testName) = @_;

    my ($newRef, $newAlt, $newPos) = $spdi->TrimLeft($oldRef, $oldAlt);
    cmp_ok($newRef, 'eq', $expRef, $testName);
    cmp_ok($newAlt, 'eq', $expAlt, $testName);
    cmp_ok($newPos, '==', $expPos, $testName);
}

sub Test_TrimRight
{
    my ($oldRef, $oldAlt, $expRef, $expAlt, $expPos, $testName) = @_;

    my ($newRef, $newAlt, $newPos) = $spdi->TrimRight($oldRef, $oldAlt);
    cmp_ok($newRef, 'eq', $expRef, $testName);
    cmp_ok($newAlt, 'eq', $expAlt, $testName);
    cmp_ok($newPos, '==', $expPos, $testName);
}

sub Test_Blossom
{
    my ($oldPos, $bud, $expPos, $expRef, $expAlt, $testName) = @_;

    my ($newRef, $newAlt, $newPos) = $spdi->Blossom($oldPos, $bud);
    cmp_ok($newPos, '==', $expPos, $testName);
    cmp_ok($newRef, 'eq', $expRef, $testName);
    cmp_ok($newAlt, 'eq', $expAlt, $testName);
}

sub Test_ConvertIndelToSPDI
{
    my ($oldPos, $oldRef, $oldAlt, $expPos, $expRef, $expAlt, $expErr, $testName) = @_;

    my ($newPos, $newRef, $newAlt, $error) = $spdi->ConvertIndelToSPDI($oldPos, $oldRef, $oldAlt);
    cmp_ok($newPos, '==', $expPos, $testName);
    cmp_ok($newRef, 'eq', $expRef, $testName);
    cmp_ok($newAlt, 'eq', $expAlt, $testName);
    cmp_ok($error, 'eq', $expErr, $testName) if ($error || $expErr);
}

