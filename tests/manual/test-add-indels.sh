ecoli_k12=/pbi/dept/secondary/siv/yli/sv/ecoli/ecoliK12_pbi_March2013_DEL_268834_LEN_5647.fasta
sr=/pbi/dept/secondary/siv/yli/sv/ecoli/subreads/ecoli_10X.subreadset.xml

add-indels-to-fasta $ecoli_k12 modified_ecoli_k12.fasta out_truth.bed --delete --numbers 2,2 --lens 100,200
pbsv run modified_ecoli_k12.fasta $sr out.bed
diff out.bed out_truth.bed

add-indels-to-fasta $ecoli_k12 modified_ecoli_k12_2.fasta out_truth_2.bed --insert --numbers 2,2 --lens 100,200
pbsv run modified_ecoli_k12_2.fasta $sr out_2.bed
diff out.bed out_truth.bed
