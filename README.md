# bamFilters
A utility tool to filter alignments from bam files, using identity percent, low complexity and read coverage.

```
bamFilters -h 
--------------------------------------------------------------------------------------------
Usage:   bamFilters -b <in.bam> -o <out.bma> [-u <out.uniq.bam>] [options] 

Options: -b          FILE  input BAM file
         -o          FILE  output BAM file
         -i          INT   identity percent level, default is 0
         -a          INT   alignment percent level, default is 0
         -r          INT   filter out sequences which contain more than r% of low-complexity bases, default is 100%
         -n          INT   filter out sequences which contain less than n high-complexity bases, default is 0
         -s          INT   don't count deletions longer than s as parts of alignment (useful for spliced mapping)
         -u          FILE  output BAM file for uniq filter which select reads that mapped only at one position
         -z                bam files no genereted, only list of sequences names of reference with one or more match. 
                           No compatible with -u and -o. print on stdout 
         -y          FILE  produce stats files with percent id, ali for each alignement
         -v                verbose mode
         -h                help
         example : bamFilters -b ./in.bam -i 95 -a 80 -r 75 -n 30 -o ./out.bam -u ./out.uniq.bam 

--------------------------------------------------------------------------------------------
```
