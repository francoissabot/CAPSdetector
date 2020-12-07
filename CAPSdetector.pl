#!/usr/bin/perl 

###################################################################################################################################
#
# Copyright 2017-2018 IRD
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD
# Version 1 and latter written by Francois Sabot and Camille Carrette
####################################################################################################################################



use strict;
use Getopt::Long;
use Data::Dumper;


my ($file, $ezList, $outfile, $outFormat, $startType, $refFolder, $qual, $mq, $mqZero, $del, $readPosRS, $help,$outfileprim,$PrimThermoPath);

my $courriel="francois.sabot-at-ird.fr";
my ($nomprog) = $0 =~/([^\/]+)$/;
my $MessAbruti ="\nUsage:
\t$nomprog -i inputTabulatedFile -o outFile -e enzymeList -r referenceFolder [-f outFormat -s startValue -q QUAL -m MQ -n MQ0 -d Dels, -p ReadPosRankSum] -pr OutputPrimerFile


From a tabulated file, formatted with Chromosome - Position (VCF of BED-like), will extract the different positions and detected using EMBOSS restrict the putative restriction sites associated with a provided enzyme list.

The referenceFolder contains the chromosomes/scaffold files (one sequence per file) in FASTA, with EXACTLY the same file name starting as the name of the reference sequence (e.g. the file for Chr1 is 'Chr1' or 'Chr1.fasta', but not 'seqChr1.fasta').

The output will be a basic tabulated file or a VCF (only if input is VCF).

The outFormat can be VCF/vcf or tabulated. In standard it is tabulated Chrom - Position - Enzymes.

The starting point can be 1 (VCF) or 0 (BED-like). The 1-based value (VCF) is used as default.

If the file is a VCF, the data will be considered only if PASS Filter (automatic),  Calling QUAL of 200 (-q parameter), MQ of 20 (-m parameter), MQ0 = 0 (-n parameter), Dels = 0.00 (-d parameter), ReadPosRankSum < -5 (-p parameter).

The output file for primer activate the display of the primer left right and internal of each positions (-pr shouldn't be equal to '0')

If the message : 'Unable to open file /sastack.ds' appeare, you have to change the path of your primer3_config with -prpath

No filter will be applied on non VCF files

        contact: $courriel\n\n";


#Standard values
$startType = 1;
$outFormat = "vcf";
$qual = 200;
$mq = 20;
$mqZero = 0;
$del = 0;
$readPosRS = -5;
$outfileprim = 0;
$PrimThermoPath ="/opt/Geneious/bundledPlugins/com.biomatters.plugins.primerDesign.PrimerDesignPlugin/com/biomatters/plugins/primerDesign/primer3_2_3_4/primer3_config/";

unless (@ARGV) 
        {
        print "\nType --help for more informations\n\n";
        exit;
        }


GetOptions("prout|help|?|h" => \$help,
	   "i|in=s"=> \$file,
	   "o|out=s"=> \$outfile,
	   "e|ez=s"=> \$ezList,
	   "f|format=f"=> \$outFormat,
	   "s|start=s" => \$startType,
	   "r|ref=s" => \$refFolder,
	   "q|qual=s" => \$qual,
	   "m|mq=s" => \$mq,
	   "n|mqzero=s" => \$mqZero,
	   "d|del=s" => \$del,
	   "p|readpos=s" => \$readPosRS,
	   "pr|prim=s" => \$outfileprim,
	   "prpath|primPath=s"=> \$PrimThermoPath		
	   );

if ($help)
	{
	print $MessAbruti,"\n";
	exit;
	}



#Opening outfile
open (my $fhOut, ">", $outfile) or die("\nCannot create the outfile $outfile:\n$!\n");


#Reading file and checking the format
my $vcfControl = &checkFormatVcf($file); # 1 if VCF, 0 if simple tabulated or BED

#Position control
my $positionCorrector = 0;
$positionCorrector = -1 if $startType == 0;



#generate the individual names in vcf input
my $indName = &indNames if $vcfControl == 1;


my $firstOutLine =  "#Chrom\tPosition\tREF\tEnzyme\tALT1\taltEnzyme1\tALT2\taltEnzyme2\tALT3\taltEnzyme3\tALT4\taltEnzyme4\tPrimerLeft0\tPrimerRight0\tPrimerInternal0\tlengthOfInternalPrimer0".$indName."\n";
print $fhOut $firstOutLine;

#Reading and extraction for each position
open (my $fhIn, "<", $file) or die ("\nCannot open $file:\n$!\n");

while (my $line = <$fhIn>)
{
	chomp $line;
	next if $line =~ m/^#/;
	
	my @fields = split /\t/, $line;
	
	
	
	#Filtering on VCF limits
	if ($vcfControl > 0)
	{
		#Filtering on calling quality and filter value
		next if $fields[5] < $qual;
		next if $fields[6] ne "PASS";
		
		#Splitting INFOS fields
		my @data = split /;/, $fields[7];
		my %infos;
		while (@data)
		{
			my $data = shift @data;
			my ($clef, $valeur) = split /=/, $data;
			$infos{$clef} = $valeur;
			
			
			
		}
		
		#filtering on infos values
		if (defined $infos{MQ})
		{		
			next if ($infos{MQ} < $mq);	
		}
		if (defined $infos{MQ0})
		{
			next if ($infos{MQ0} < $mqZero);
			
		}
		if (defined $infos{Dels})
		{
			next if $infos{Dels} < $del;
		}
		if (defined $infos{ReadPosRankSum})
		{
			next if $infos{ReadPosRankSum} < $readPosRS;
		}
		
	}
	
	#Starting analysis
	my $chromosome = $fields[0];
	my $position = $fields[1] + $positionCorrector;
	my $ref;
	my @alt;	
		if($vcfControl > 0)
		{	
			$ref=$fields[3];
			push @alt, $fields[4];
		}else
		{
			$ref=$fields[2];
			push @alt ,$fields[4];
			push @alt ,$fields[5];
			push @alt ,$fields[6];
		}
	
	#value for identify the statue of the individuals
	my $RestrictSite=0;
	my $RestrictSiteAlt=0;
	my @outline;
	my $sequence = &extractSequence($refFolder."/".$chromosome,$position,$ref,$chromosome);
	my $ezOk = &checkEz($sequence,$ezList);
	my $ezAltOk;
	my @primer;	
	
	#If -pr = 0 no display of primers 
	push @primer , "---";
	push @primer , "---";
	push @primer , "---";
	push @primer , "---";
	@primer= &primer($sequence) if $outfileprim ne 0;#content primer left-right-internal 0 and the length of the internal	
	
	#enzyme test on reference base
	if(ref $ezOk)
	{
		$RestrictSite=0;

	                while (@{$ezOk})
			{
	     			my $localEz = shift @{$ezOk};
	     			my $localStart = shift @{$ezOk};
	     			my $localStop = shift @{$ezOk};
	 			push @outline , $chromosome."\t".$position."\t".$ref."\t".$localEz."\t";
			}
	}else
	{
		$RestrictSite=1;
		push @outline ,  $chromosome."\t".$position."\t".$ref."\t---\t";
	}   
	
	#enzyme test on altenrative bases

	#Manage the case of multiple mutations 
	my @list = split /,/,$alt[0] if $vcfControl > 0;
	push @list , split("/",$alt[0]) if $vcfControl ==0;
	push @list , split("/",$alt[1]) if $vcfControl ==0;
	push @list , split("/",$alt[2]) if $vcfControl ==0;
	
	my $i = 0;
	while (@list or $i !=4)
	{
		my $alt = shift @list;
		my $sequencealt = &extractSequence($refFolder."/".$chromosome,$position,$alt,$chromosome);
		$ezAltOk=&checkEz($sequencealt,$ezList);
	       	$alt = "-" , $ezAltOk = 0 unless($alt);
			

		if(ref $ezAltOk)
		{
			$RestrictSiteAlt=1;	
			        while (@{$ezAltOk})
				{
		 			my $altEz = shift@{$ezAltOk};
					my $altStart = shift@{$ezAltOk};
					my $altStop = shift@{$ezAltOk};
					push @outline ,  $alt."\t".$altEz."\t";
					$i++;
				}
		}else
		{
			$RestrictSiteAlt=0;
			push @outline ,  $alt."\t---\t";
			$i++;
		}
	}
		
		#Primers's display
		push @outline ,  $primer[0]."\t".$primer[1]."\t".$primer[2]."\t".$primer[3]."\t";

		#state of each individual (Cut/notCut) if vcf Format
		my $n=9;
		while($fields[$n])
		{
			my @individu = split (":", $fields[$n]);
			@individu = split ("/", $individu[0]);
			my $Cut= "notCut";
				if($individu[0]==$RestrictSite and $individu[1]==$RestrictSiteAlt or $individu[0]==0 and $individu[1]==1)
				{		
					$Cut = "Cut";	
				}
				if(ref $ezOk and ref $ezAltOk)
				{
					$Cut = "Cut";	
				}

			push @outline ,  $Cut."\t";	
			$n++;	
		}
	

	#Line by line's display
	my $Outline;

		while(@outline)
		{
			$Outline = $Outline . shift @outline;
		}
	
	my @field;
	push @field , split /\t/ , $Outline;
	my $Enz = $field[3].$field[5].$field[7].$field[9];
	
	#avoid print lineS without any restrictions sites (before of after the muatations) 
	print $fhOut $Outline."\n" if $Enz ne "------------";	

}

close $fhIn;
close $fhOut;

print "\nFinished\n";

exit;

###################################
#
#	SUB PROGRAMS
#
###################################

#Primer3_core apply on each positions's sequences
sub primer{
	my ($sequence)=@_;
	
	#create the new file for the intput of primer3 for each positions
	open (my $fhOutPrim, ">",$outfileprim) or die("\nCannot create the outfile $outfileprim:\n$!\n");
 	
	#standard values for the input file for primer3_core
	my $SeqTarget = 37.21;
	my $PrimTask = "generic";
	my $PrimPickLeft = 1;
	my $primPickInt = 1;
	my $PrimPickRight = 1;
	my $PrimOptSize = 18;
	my $PRimMinSize =15;
	my $PrimMaxSize=21;
	my $PrimMaxNS = 1;
	my $PrimerProdSizeRange = "75-100";
	my $P3FileFlag = 1;
	my $SeqIntExcluedReg =37.21;
	my $PrimExplainFlag = 1;

my $PrimFile="SEQUENCE_TARGET=$SeqTarget
PRIMER_TASK=$PrimTask
PRIMER_PICK_LEFT_PRIMER=$PrimPickLeft
PRIMER_PICK_INTERNAL_OLIGO=$primPickInt
PRIMER_PICK_RIGHT_PRIMER=$PrimPickRight
PRIMER_OPT_SIZE=$PrimOptSize
PRIMER_MIN_SIZE=$PRimMinSize
PRIMER_MAX_SIZE=$PrimMaxSize
PRIMER_MAX_NS_ACCEPTED=$PrimMaxNS
PRIMER_PRODUCT_SIZE_RANGE=$PrimerProdSizeRange
P3_FILE_FLAG=$P3FileFlag
SEQUENCE_INTERNAL_EXCLUDED_REGION=$SeqIntExcluedReg
PRIMER_EXPLAIN_FLAG=$PrimExplainFlag
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$PrimThermoPath
=
SEQUENCE_ID=test1
SEQUENCE_TEMPLATE=$sequence\n=";
print $fhOutPrim $PrimFile;

	my $primCom = "primer3_core $outfileprim";
	my $prim =`$primCom` or die ("\nCannot execute the following command:\n$primCom\n$!\n");
	my @prim = split /\n/,$prim;
	my $i=0;
	my %info;	
		
		while (@prim !="")
		{
			my $prim = shift @prim;
			my ($clef, $valeur) = split /=/, $prim;
			$info{$clef} = $valeur;
		}
		#picking informations in the primer3 output
	my @list;
	push @list , $info{PRIMER_LEFT_0_SEQUENCE};
	push @list , $info{PRIMER_RIGHT_0_SEQUENCE};
	push @list , $info{PRIMER_INTERNAL_0_SEQUENCE};
	push @list , length($info{PRIMER_INTERNAL_0_SEQUENCE});
	
	return @list;	
}

#display of the individual's names
sub indNames{
	open (my $fhIn, "<", $file) or die ("\nCannot open $file:\n$!\n");
	
	while (my $line = <$fhIn>)
	{
		chomp $line;
		next if $line =~ m/^##/;
		next if $line =~ m/Chr/;
		my @fields = split /\t/, $line;
		my $Names=$fields[9];
		my $i=10; 
			while($fields[$i])
			{
				$Names=$Names."\t".$fields[$i];
				$i++;		
			}

	return $Names;
	}

}

#extraction sequence with samtools
sub extractSequence{
	my ($ref,$loc,$base,$Chr) = @_;
	
	
	#extraction of 400 bases arrond the ref or alt base
	my $start1=$loc-200;
	my $stop1=$loc-1;
	my $start2=$loc+1;
	my $stop2=$loc+200;	
	my $extractCom = "samtools faidx $ref '$Chr':$start1-$stop1 '$Chr':$start2-$stop2";
	my $sequence =`$extractCom` or die ("\nCannot execute the following command:\n$extractCom\n$!\n");	
	$sequence =~ s/>\w+\n//;
	$sequence =~ s/\n//g;
	my @sequences = split(">",$sequence);
	my $i=1;	
		while($sequences[$i])
		{	
			$sequences[$i]=~s/Chr\d+:\d*-\d*//;
			$sequences[$i]=~s/:\d*-\d*//;
			$i++;	
		}
	my $sequence1=$sequences[1].$base.$sequences[2];
	
	return $sequence1;
}


#check restriction site in the sequence
sub checkEz{
	my ($sequence, $ezFile) = @_;
	my @list;
	my $restrictCom = "echo $sequence | restrict -filter -auto -enzymes \@".$ezFile." -filter -stdout";
	my $restrictOut =`$restrictCom 2> /dev/null` or die ("\nCannot execute the following command:\n$restrictCom\n$!\n");
	my @lines = split /\n/, $restrictOut;
	my %infoDigest;

	while (@lines)
	{
		my $local = shift @lines;
		next if $local =~ m/^#/;
		next if $local =~ m/^$/;
		next if $local =~ m/Start/;
		my @data = split /\s+/, $local;
		my $localisationHere= ($data[6]-200)."..".($data[7]-200);
		push @{$infoDigest{$data[4]}}, $localisationHere if defined $localisationHere;		
	}
	return 0 if scalar (keys %infoDigest) < 1;
	
	foreach my $enzyme (keys %infoDigest)
	{
		next if scalar @{$infoDigest{$enzyme}} > 1; # multiple cuts in the sequence
		
		my ($start, $stop) = split /\.\./, ${$infoDigest{$enzyme}}[0];
		next if ($stop > 5);
		next if ($start < -5);
		next if ($stop < -5);
		next if ($start >5);
		
		
		push @list, $enzyme;
		push @list, $start;
		push @list, $stop;	
	}
	
	return 0 if scalar @list == 0;
	return \@list;
}

sub checkFormatVcf{
    #Inspired from The vcftools vcf-validator, from the TOGGLE adapted version
    my ($localFile) = @_;

    #Parsing the file
    my @header;#List to gather the header
    my @listOfFields;
    open(my $inputHandle, "<",$localFile) or die ("Cannot open the file $localFile\n$!\n");

    while (my $line=<$inputHandle>)
    {
        chomp $line;

        ##DEBUG print $line."\n";
        if ($line =~ m/^#/)
        {
           push (@header, $line);
	   next;
        }
        else
        {
            @listOfFields = split /\t/, $line;
	    last;
        }
    }

    #Check if the first line of the header is including the version
    my $versionLine=defined $header[0]? $header[0]:undef;
    return 0 unless (defined $versionLine); #Thrown an error if the $version cannot be obtained (misformatted line)

    my @version=split /VCFv/,$versionLine;
    
    return 0 if (scalar(@version)==0); #Thrown an error if the $version cannot be obtained (misformatted line)
    ##DEBUG print "DEBUG: $0: vcf version $versionLine : ".scalar(@version)." : $version[1] \n";
    
    eval ($version[1] == $version[1]) or return 0; #Verifying if the value obtained is numerical.

    # Check the first line format as recommanded
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Individual1
    #Chr1    228069  .       A       .       31.23   LowQual AN=2;DP=1;MQ=37.00;MQ0=0        GT:DP   0/0:1
    if (scalar @listOfFields < 10) #Less than the 10 minimum fields
    {
        return 0;
    }

    #Check if the second field (the position) is numerical
    eval ($listOfFields[1] == $listOfFields [1]) or return 0; #Verifying if numerical.

    close $inputHandle;

    return 1; #Return correct if all check are Ok
}
