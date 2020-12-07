CAPSdetector README
 -----------------------------------------------------------------------------------------------
 Copyright 2017-2018 IRD

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/> or
 write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston,
 MA 02110-1301, USA.

 You should have received a copy of the CeCILL-C license with this program.
 If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>

 Intellectual property belongs to IRD
 Version 1 and latter written by Francois Sabot and Camille Carrette
 -----------------------------------------------------------------------------------------------

Installations you need

Samtools for the extraction of sequences in vcf of BED file.

restrict to search restriction site in each chromosom in reference file

Primer3_core if you want to have the primers markers with -pr.

 -----------------------------------------------------------------------------------------------
	Gitlab contain:
	-enzyme.txt is a list of enzymes that you can use to make test 
	-compeleTest.vcf is a vcf that you can use to make test
	-snp3k_S35.BED is a BED	that you can use to make test
	-CAPSdetector.pl the perl script of CAPSdetector
	-Chromosmes is the reference folder where we test the restriction site with the enzymes
 -----------------------------------------------------------------------------------------------

	Usage:

CAPSdetector.pl -i inputTabulatedFile -o outFile -e enzymeList -r referenceFolder -f outFormat -s startValue -q QUAL -m MQ -n MQ0 -d Dels, -p ReadPosRankSum -pr OutputPrimerFile


From a tabulated file, formatted with Chromosome - Position (VCF of BED-like), will extract the different positions and detected using EMBOSS restrict the putative restriction sites associated with a provided enzyme list.

The referenceFolder contains the chromosomes/scaffold files (one sequence per file) in FASTA, with EXACTLY the same file name starting as the name of the reference sequence (e.g. the file for Chr1 is 'Chr1' or 'Chr1.fasta', but not 'seqChr1.fasta').

The output will be a basic tabulated file or a VCF (only if input is VCF).

The outFormat can be VCF/vcf or tabulated. In standard it is tabulated Chrom - Position - Enzymes.

The starting point can be 1 (VCF) or 0 (BED-like). The 1-based value (VCF) is used as default.

If the file is a VCF, the data will be considered only if PASS Filter (automatic),  Calling QUAL of 200 (-q parameter), MQ of 20 (-m parameter), MQ0 = 0 (-n parameter), Dels = 0.00 (-d parameter), ReadPosRankSum < -5 (-p parameter).

The output file for primer activate the display of the primer left right and internal of each positions (-pr shouldn't be equal to '0')

If the message : 'Unable to open file /sastack.ds' appeare, you have to change the path of your primer3_config with -prpath

No filter will be applied on non VCF files

        contact:Francois.sabot@ird.fr
		Camille.Carrette@outlook.com

 -----------------------------------------------------------------------------------------------
