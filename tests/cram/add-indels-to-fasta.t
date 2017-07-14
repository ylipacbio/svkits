Set up
  $ . $TESTDIR/setup.sh
  $ in_fa=/pbi/dept/secondary/siv/testdata/../yli/sv/ecoli/Ecoli/genome/genome.fa #ecoli

  $ out_fa=$OUTDIR/ecoli_del_2_2_100_200.fasta
  $ out_bed=$OUTDIR/ecoli_del_2_2_100_200.bed
  $ add-indels-to-fasta $in_fa $out_fa $out_bed --delete --numbers 2,2 --lens 100,200 1>/dev/null 2>/dev/null && echo $?
  0
  $ ls $out_fa >/dev/null && echo $?
  0
  $ ls $out_bed >/dev/null && echo $?
  0

  $ out_fa=$OUTDIR/ecoli_ins_2_2_100_200.fasta
  $ out_bed=$OUTDIR/ecoli_ins_2_2_100_200.bed
  $ add-indels-to-fasta $in_fa $out_fa $out_bed --insert --numbers 2,2 --lens 100,200 1>/dev/null 2>/dev/null && echo $?
  0
  $ ls $out_fa >/dev/null && echo $?
  0
  $ ls $out_bed >/dev/null && echo $?
  0
