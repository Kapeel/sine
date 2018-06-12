## SINE background
SINEs are transcribed by RNA polymerase III, and are derived from one of three classes of Pol III–transcribed molecules (tRNA, 7SL, 5s rRNA).
While animal SINEs from all three classes are known, plant SINEs are exclusively derived from tRNA. 
To find SINEs, I apply the implementation of [Wenke et al. 2011](http://www.plantcell.org/content/early/2011/09/07/tpc.111.088682) in [SINE-finder](http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt), which searches for tRNA-derived SINEs containing RNA polymerase III A and B boxes near the polyA tail. 
The defaults are that A and B box consensus nucleotide sequences are RVTGG and GTTCRA, there is a 25–50 bp spacer between the A and B boxes, and there is a spacer of 20–500 bp between the B box and polyA tail.

## How to identify SINEs

The script ```run_sines.sh``` will download sine_finder, run sine_finder, parse results, assign to families, and output a GFF.
Each candidate SINE is clustered using VSEARCH and silix, to characterize families, and added to existing families if specified in the config file.
Each family is also matched to [Maize TE Consortium](http://www.maizetedb.org) exemplars (via 80-80-80 identity) to faciliatate comparison between genome versions.

## Output

 - `${GENOMENAME}.RST.gff3` is all SINEs identified. 
 - `${GENOMENAME}.RST.tabout` contains extended information, like the TSD length and mismatch
 - `${GENOMENAME}.RST.fa` is all the SINEs identified, renamed with their TEID, in fasta format.
 - `post-${GENOMENAME}.existingRST.fa` includes all previously identified copies with the newly identified ones. These need to be deposited somewhere, so the next annotator can use them to add to existing families. Jeff and I have talked about maizeGDB hosting these via a link so somebody annotating can download existing and upload when finished.
