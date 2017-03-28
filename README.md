Tool kits for structural variation project, including
* simulation
* analysis
* validation
* tests
* miscs.

# add-an-indel-to-fasta
usage: add-an-indel-to-fasta [-h] [--n_bases N_BASES] [--insert | --delete]
                             in_fasta out_fasta out_bed

Modify a reference fasta by either inserting a substring to a reference
sequence or deleting a substring from a reference sequence, make a truth SV
BED file for mapping reads of the original reference to the modified reference
file.

positional arguments:
  in_fasta           Input FASTA
  out_fasta          Output FASTA with indels
  out_bed            Ground truth SV calls in BED if modified fasta were used
                     as reference

optional arguments:
  -h, --help         show this help message and exit
  --n_bases N_BASES  Number of bases to insert to or delete from a reference
                     sequence
  --insert           Insert a string to a sequence in in_fasta, and make a BED
                     file with an Deletion SV call
  --delete           Delete a substring from a sequence in in_fasta, and make
                     a BED file with an Insertion SV call


# make-subreads-bam-from-fasta
usage: make-subreads-bam-from-fasta [-h] [--dry_run]
                                    in_fasta in_bam_fofn out_prefix

Make a subreads.bam/xml containing all reads of a fasta file

positional arguments:
  in_fasta     Input FASTA file
  in_bam_fofn  Input BAM fofn
  out_prefix   Output prefix

optional arguments:
  -h, --help   show this help message and exit
  --dry_run    Dry run


# make-subreads-bam-from-zmws 
usage: make-subreads-bam-from-zmws [-h] [--dry_run]
                                   in_txt in_bam_fofn out_prefix

Make a subreads.bam/xml containing all subreads of zmws.txt

positional arguments:
  in_txt       Input txt file containg ZMW names or read names
  in_bam_fofn  Input BAM fofn
  out_prefix   Output prefix

optional arguments:
  -h, --help   show this help message and exit
  --dry_run    Dry run, display CMDs

