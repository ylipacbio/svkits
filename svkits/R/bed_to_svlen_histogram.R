#!/usr/bin/env Rscript

filter=dplyr::filter
select=dplyr::select
unite=tidyr::unite

library("optparse")
 
option_list = list(make_option(c("-b", "--bed"), type="character", default=NULL, 
                               help="Input BED file name", metavar="character"),
                   make_option(c("-n", "--name"), type="character", default="Name",
                               help="Name [default= %default]", metavar="character"),
                   make_option(c("-o", "--out"), type="character", default="out.svlen.histogram.png", 
                               help="output png file name [default= %default]", metavar="character")
                   ); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bed)){
      print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

input_bed = opt$bed
name = opt$name
out_png = opt$out

print(input_bed)
print(name)
print(out_png)


## Program

t <- fread("20170109_HG00733.inversions.UW.formated.bed")
g = ggplot(t, aes(x=SVLENGTH, fill=SVTYPE)) + geom_histogram() + ggtitle(name) + xlab('Structural Variant Length') + ylab('Count')
png(out_png)
g
dev.off()
