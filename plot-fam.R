#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  cat("Usage: ./plot-fam.R [base name]\n")
  quit()
} else {
  cat("Generating pdfs with pedigree structures from ", args[1], ".fam...\n\n",
      sep="")
}

famfile <- args[1]

if (substr(famfile, nchar(famfile)+1-4, nchar(famfile)) != ".fam") {
  famfile <- paste(famfile, ".fam", sep="")
} else {
  args[1] <- substr(famfile, 0, nchar(famfile)-4)
}

if (!file.exists(famfile)) {
  cat("ERROR:", famfile, "does not exist\n")
  quit()
}

dat <- read.table(famfile)

#plot pedigree
library(kinship2)

# Replace 0 parent ids with NA:
no_dad <- dat[,3]==0
dat[no_dad,3] <- NA
no_mom <- dat[,4]==0
dat[no_mom,4] <- NA

fam_id <- dat[,1]
status <- rep(1,length(dat[,6]))
status[dat[,6] != -9] <- 2

for(i in 1:length(fam_id)){
  sel <- dat[,1]==fam_id[i]

  ped <- pedigree(dat[sel,2], dat[sel,3], dat[sel,4], dat[sel,5], status[sel])

  pdf(file =paste(args[1], "-", fam_id[i], ".pdf", sep=""), width = 20, height = 8, onefile = TRUE, paper = "special")
  par(xpd=TRUE)
  # extra margin on left, right and a bit on top so that the sample names fit
  # change cex larger to make the circles/squares and labels larger relative to
  # the lines connecting samples: (may want to increase the margin values
  # in proportion
  print(plot.pedigree(ped,id=dat[sel,2], mar=c(0,4,1,4), cex=1))
  dev.off()
}
