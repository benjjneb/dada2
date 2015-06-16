## ----init----------------------------------------------------------------
#library("dadac"); packageVersion("dadac")
library("ggplot2"); packageVersion("ggplot2")

## ---- eval=TRUE----------------------------------------------------------
# System-agnostic file path
fooFile = system.file("extdata", "kmer.csv.gz", package="dadac")
# Read the table
foo = read.csv(fooFile)
colnames(foo) <- c("Kmer", "Align")
# Plot the Alignment versus Kmer distance as density plot
ggplot(data=foo, aes(x=Kmer, y=Align)) + 
  geom_bin2d(binwidth = c(0.005, 0.005))

