#!/usr/bin/perl
use strict;

my $CORRFITIN = 100;
my $CORRFITOUT = 101;

my $fitterlist;
my $modellist;

my @fitterarray;
my @modelarray;
my $res;

$fitterlist = `ls CFFitter?*.cxx`;
$modellist  = `ls SourceModel?*.cxx`;

$res = $fitterlist;
$res =~ s/.cxx//g;
@fitterarray = split(/\n/,$res);

$res = $modellist;
$res =~ s/.cxx//g;
@modelarray = split(/\n/,$res);

# print $#fitterarray;
# print $#modelarray;

$fitterlist =~ s/.cxx/\" << endl/g;
$modellist  =~ s/.cxx/\" << endl/g;

$fitterlist =~ s/CFFitter/ << \"\# CFFitter/g;
$modellist  =~ s/SourceModel/ << \"\# SourceModel/g;

# print $fitterlist;
# print $modellist;

open(CORRFITIN, "CorrFit.cxx.in");
open(CORRFITOUT, ">CorrFit.cxx");

while (<CORRFITIN>)
  {
    if (m/CFList/) { print CORRFITOUT $fitterlist; }
    elsif (m/SMList/) { print CORRFITOUT $modellist; }
    elsif (m/SMCreate/) 
      {
	foreach (@modelarray) {
	  print CORRFITOUT "if (sSourceModelName == \"".$_."\") \{ tModel = (".$_." *) new ".$_."(); \}\n";
	}
      }
    elsif (m/CFCreate/) 
      {
	foreach (@fitterarray) {
	  print CORRFITOUT "if (sCFFitterName == \"".$_."\") \{ tFitter = (".$_." *) new ".$_."(); \}\n";
	}
      }
    elsif (m/SMIncludeList/) 
      {
	foreach (@modelarray) {
	  print CORRFITOUT "#include \"".$_.".h\"\n";
	}
      }
    elsif (m/CFIncludeList/) 
      {
	foreach (@fitterarray) {
	  print CORRFITOUT "#include \"".$_.".h\"\n";
	}
      }
    else { print CORRFITOUT $_; }
  }

close CORRFITIN;
close CORRFITOUT;


