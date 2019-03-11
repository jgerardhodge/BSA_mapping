#!/usr/bin/env Rscript

#Load required libraries
library(optparse)

option_list = list(
  	make_option(c("-i", "--input_nmu_variants"), type="character", default=NULL, help="Path to input mpileup csv containing all NMU SNPs in the mutant bulk", metavar="NMU VARIANTS"),
   	make_option(c("-f", "--filter_variants"), type="character", default=NULL, help="Path to input mpileup csv containing all high frequency mutant specific variants", metavar="FILTER VARIANTS"),
	make_option(c("-s", "--genome_size"), type="character", default=NULL, help="Path to a tab delimited file containing sizes of chromosomes", metavar="GENOME SIZE"),
   	make_option(c("-c", "--chr"), type="integer", default=NULL, help="Number of chromosome of interest (optional)", metavar="CHROMOSOME"),
 	make_option(c("-w", "--win_size"), type="integer", default=1000, help="Window size for calculating average variant frequencies (default 1000)", metavar="WINDOW SIZE"),
 	make_option(c("-o", "--output"), type="character", default=paste('Bulk_filter_map_', Sys.Date(), sep=''), help="output directory name [default= %default]", metavar="OUTHANDLE")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

nmu_variants=read.csv(opt$input, col.names=c('LG', 'pos', 'ref', 'A', 'C', 'G', 'T', 'cov', 'fvar'))
filt_variants=read.csv(file=opt$filter_variants, col.names=c('LG', 'pos', 'ref', 'A', 'C', 'G', 'T', 'cov', 'fvar'))
genome_size=read.csv(file=opt$genome_size)

win_size=opt$win_size/2

cum_genome_size<-c(0, cumsum(genome_size[1:9,2]))

if (!is.null(opt$chr)){
	tiff(file=paste(opt$output, '.tif', sep=''), height=1200, width=1800)
	layout(matrix(c(1 ,2),2,1))
} else {
	tiff(file=paste(opt$output, '.tif', sep=''), height=600, width=1800)
}

plot (cum_genome_size[c(1,10)], c(0,1.15), xlab='Physical Position (Mb)', ylab='Variant Frequencies', main='Whole Genome Variant Frequency Map', col='white', xaxt='n', cex.main=2, cex.axis=1.25, cex.lab=1.25)
axis(1, at=seq(0, cum_genome_size[10], by=10000000), , labels=seq(0, cum_genome_size[10], by=10000000)/1000000, cex.axis=1.25, cex.lab=1.25)
text((cum_genome_size[c(1:9)]+cum_genome_size[c(2:10)])/2, rep(1.1, 9), labels=paste('Chr',1:9), cex=1.25)

for (LG in 1:9){

	print(paste('Iterating through chromosome ', LG, '...', sep=''))
	if (LG %% 2){
		col = 'chartreuse4'
	} else {
		col = 'goldenrod'
	}

	LG_nmu_variants<-nmu_variants[nmu_variants[,1] %in% genome_size[LG,1],]

	points(LG_nmu_variants[,2]+cum_genome_size[LG], LG_nmu_variants[,9], col=col, pch=19, cex=0.2)

	LG_filt_variants<-filt_variants[filt_variants[,1] %in% genome_size[LG,1],]

	points(LG_filt_variants[,2]+cum_genome_size[LG], LG_filt_variants[,9], pch=19, cex=0.7)
	points(LG_filt_variants[,2]+cum_genome_size[LG], LG_filt_variants[,9], pch=19, cex=0.6, col='red')

	win_pos<-seq(win_size, genome_size[LG,2], by=win_size)
	var_pos<-c()
	
	for (win in win_pos){
		var_pos<-c(var_pos, mean(LG_nmu_variants[LG_nmu_variants[,2]>(win-win_size) & LG_nmu_variants[,2]<=(win+win_size),9]))
	}

	lines(win_pos+cum_genome_size[LG], var_pos, lwd=2.5, col='red')
}

if (!is.null(opt$chr)){
	print(paste('Plotting specified chromosome ', opt$chr, '...', sep=''))

	plot (c(0, genome_size[opt$chr,2]), c(0,1.15), xlab='Physical Position (Mb)', ylab='Variant Frequencies', main=paste('Chr', opt$chr, 'Variant Frequency Map'), col='white', xaxt='n', cex.main=2, cex.axis=1.25, cex.lab=1.25)
	axis(1, at=seq(0, genome_size[opt$chr,2], by=1000000), , labels=seq(0, genome_size[opt$chr,2], by=1000000)/1000000, cex.axis=1.25, cex.lab=1.25)
	
	if (opt$chr %% 2){
		col = 'chartreuse4'
	} else {
		col = 'goldenrod'
	}

	LG_nmu_variants<-nmu_variants[nmu_variants[,1] %in% genome_size[opt$chr,1],]

	points(LG_nmu_variants[,2], LG_nmu_variants[,9], col=col, pch=19, cex=0.2)

	LG_filt_variants<-filt_variants[filt_variants[,1] %in% genome_size[opt$chr,1],]

	points(LG_filt_variants[,2], LG_filt_variants[,9], pch=19, cex=0.7)
	points(LG_filt_variants[,2], LG_filt_variants[,9], pch=19, cex=0.6, col='red')

	win_pos<-seq(win_size, genome_size[LG,2], by=win_size)
	var_pos<-c()
	
	for (win in win_pos){
		var_pos<-c(var_pos, mean(LG_nmu_variants[LG_nmu_variants[,2]>(win-win_size) & LG_nmu_variants[,2]<=(win+win_size),9]))
	}

	lines(win_pos, var_pos, lwd=2.5, col='red')
	text((20105893+20085604)/2, 1.1, labels='SvSPL9', cex=1.25)
}
 
dev.off()


