#!/mnt/software/R/3.5.1/bin/R

library(seqinr)

#arguments
args=(commandArgs(TRUE))

for(i in 1:length(args)) {
    try(eval(parse(text=args[[i]])))
}
if (!exists("f") && !exists("o")){
    cat("\n>>>NOPE\nplease define file with f=file and o=out\nFull command should look like this : \n$ R < tetra_nucl.R --no-save \"--args f='file.fa o='outfile.tetranucl' \"\n\n")
    quit()
}

input <- f
output <- o

seq=read.fasta(input)[[1]]

print("counting")

tetra_count = count(seq,freq=TRUE,wordsize=4)

write.csv(tetra_count, file=output, row.names=FALSE, quote=FALSE)
