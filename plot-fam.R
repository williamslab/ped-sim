#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
rm_fam_id <- TRUE

if (length(args) < 1 || length(args) > 2) {
  cat("Usage: ./plot-fam.R [-keep_fam_id] <base name>\n")
  quit()
} else {
  if (length(args) == 2) {
    if (args[1] == "-keep_fam_id") {
      rm_fam_id <- FALSE
      famfile <- args[2]
    } else {
      cat("Usage: ./plot-fam.R [-keep_fam_id] <base name>\n")
      quit()
    }
  } else { # just one command line option -- the filename
    famfile <- args[1]
  }
  cat("Generating pdfs with pedigree structures from ", args[1], ".fam...\n\n",
      sep="")
}


if (substr(famfile, nchar(famfile)+1-4, nchar(famfile)) != ".fam") {
  famfile <- paste(famfile, ".fam", sep="")
} else {
  args[1] <- substr(famfile, 0, nchar(famfile)-4)
}

if (!file.exists(famfile)) {
  cat("ERROR:", famfile, "does not exist\n")
  quit()
}

fam_df <- read.table(famfile)

if (rm_fam_id) {
  # remove the family id from each individuals' id -- it's not very helpful for
  # the plot
  fam_df[,2] <- substr(fam_df[,2], nchar(fam_df[,1])+2, nchar(fam_df[,2]))
  fam_df[,3] <- substr(fam_df[,3], nchar(fam_df[,1])+2, nchar(fam_df[,3]))
  fam_df[,4] <- substr(fam_df[,4], nchar(fam_df[,1])+2, nchar(fam_df[,4]))
}

#plot pedigree
library(kinship2)

# Replace 0 parent ids with NA:
no_dad <- fam_df[,3]==0
fam_df[no_dad,3] <- NA
no_mom <- fam_df[,4]==0
fam_df[no_mom,4] <- NA

fam_id <- fam_df[,1]
status <- rep(1,length(fam_df[,6]))
status[fam_df[,6] != -9] <- 2

for(i in 1:length(fam_id)){
  sel <- fam_df[,1]==fam_id[i]

  ped <- pedigree(fam_df[sel,2], fam_df[sel,3], fam_df[sel,4], fam_df[sel,5], status[sel])

  pdf(file =paste(args[1], "-", fam_id[i], ".pdf", sep=""), width = 20, height = 8, onefile = TRUE, paper = "special")
  par(xpd=TRUE)
  # extra margin on left, right and a bit on top so that the sample names fit
  # change cex larger to make the circles/squares and labels larger relative to
  # the lines connecting samples: (may want to increase the margin values
  # in proportion
  print(plot.pedigree(ped,id=fam_df[sel,2], mar=c(0,4,1,4), cex=1))
  dev.off()
}
