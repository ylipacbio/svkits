Test add-an-indel-to-fasta

Set up
  $ . $TESTDIR/setup.sh
  $ in_fa=/pbi/dept/secondary/siv/testdata/../yli/sv/ecoli/Ecoli/genome/genome.fa #ecoli

  $ out_fa=$OUTDIR/ecoli_del_50.fasta
  $ out_bed=$OUTDIR/ecoli_del_50.bed
  $ rm -rf $out_bed $out_fa
  $ add-an-indel-to-fasta $in_fa $out_fa $out_bed --delete >/dev/null && echo $?
  0
  $ ls $out_fa >/dev/null && echo $?
  0
  $ ls $out_bed >/dev/null && echo $?
  0


  $ out_fa=$OUTDIR/ecoli_ins_100.fasta
  $ out_bed=$OUTDIR/ecoli_ins_100.bed
  $ rm -rf $out_bed $out_fa
  $ add-an-indel-to-fasta $in_fa $out_fa $out_bed --insert --n_bases 100 >/dev/null && echo $?
  0
  $ ls $out_fa >/dev/null && echo $?
  0
  $ ls $out_bed >/dev/null && echo $?
  0


