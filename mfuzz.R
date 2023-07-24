parser = argparse::ArgumentParser(description = 'Script for running mfuzz')
parser$add_argument('-i', '--input', dest = 'input', help = 'input mtx filename')
parser$add_argument('-p', '--prefix', dest = 'prefix', help = 'prefix of output')
opts = parser$parse_args()

library(dplyr)
library(ggplot2)
library(Mfuzz)

ex <- read.csv(opts$input, sep = '\t', as.is = TRUE, header = TRUE, row.names = 1)
ex.m <- as.matrix(ex)
eset <- new('ExpressionSet', exprs = ex.m)
eset <- filter.std(eset,min.std=0)
eset <- standardise(eset)
c <- 3
m <- mestimate(eset)
cl <- mfuzz(eset, c = c, m = m)
cl$size
#cl$cluster[cl$cluster == 1]
#cl$membership
time.labels <- c('RGC2','AC1RGC','AC1BC','BC2','ROD','CONE')

#time.labels.keep <- colnames(ex)

pdf(paste0(opts$prefix, '_mfuzz.pdf'))
mfuzz.plot2(eset, cl, mfrow = c(3, 2), 
            #colo='fancy',
            bg='black',
            col='white',
            ax.col='white',
            col.axis='white',
            col.lab='white',
            col.main='white',
            col.sub='white',
            centre.col='white',
            centre=TRUE,
            min.mem=0.5,
            new.window = FALSE, 
            time.labels = time.labels,
            x11=FALSE,
            Xwidth=200,
            )

acore <- acore(eset, cl, min.acore=0.8)
#acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))

write.table(acore_list, file = paste0(opts$prefix, '.csv'), sep='\t', quote=F, row.names=F)

