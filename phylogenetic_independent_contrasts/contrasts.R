library(phytools)
library(paleotree)
library(caper)
library(bestNormalize)
t <- read.tree("silva.nwk")
# drop zero-length terminal branches and collapse internal zero-length branches
nt <- di2multi(dropZLB(t))
df <- read.csv("TableS1_IntronPerKb_Annotation.csv")
tip <- as.data.frame(nt$tip.label)
names(tip) <- 'genome'
# create new table
d <- merge(df, tip, by='genome')

# try different transformations
# intron_per_kb
log_obj <- log10(d$intron_per_kb)
arcsinh_obj <- arcsinh_x(d$intron_per_kb)
boxcox_obj <- boxcox(d$intron_per_kb)
yeojohnson_obj <- yeojohnson(d$intron_per_kb)
orderNorm_obj <- orderNorm(d$intron_per_kb)
par(mfrow = c(2,3))
MASS::truehist(d$intron_per_kb, main = "No transformation", nbins = 12)
MASS::truehist(log_obj, main = "log10", nbins = 12)
MASS::truehist(arcsinh_obj$x.t, main = "Arcsinh", nbins = 12)
MASS::truehist(boxcox_obj$x.t, main = "Box Cox", nbins = 12)
MASS::truehist(yeojohnson_obj$x.t, main = "Yeo-Johnson", nbins = 12)
MASS::truehist(orderNorm_obj$x.t, main = "orderNorm", nbins = 12)
# genome_size
log_obj <- log10(d$genome_size)
arcsinh_obj <- arcsinh_x(d$genome_size)
boxcox_obj <- boxcox(d$genome_size)
yeojohnson_obj <- yeojohnson(d$genome_size)
orderNorm_obj <- orderNorm(d$genome_size)
par(mfrow = c(2,3))
MASS::truehist(d$genome_size, main = "No transformation", nbins = 12)
MASS::truehist(log_obj, main = "log10", nbins = 12)
MASS::truehist(arcsinh_obj$x.t, main = "Arcsinh", nbins = 12)
MASS::truehist(boxcox_obj$x.t, main = "Box Cox", nbins = 12)
MASS::truehist(yeojohnson_obj$x.t, main = "Yeo-Johnson", nbins = 12)
MASS::truehist(orderNorm_obj$x.t, main = "orderNorm", nbins = 12)
# num_CDS
log_obj <- log10(d$num_CDS)
arcsinh_obj <- arcsinh_x(d$num_CDS)
boxcox_obj <- boxcox(d$num_CDS)
yeojohnson_obj <- yeojohnson(d$num_CDS)
orderNorm_obj <- orderNorm(d$num_CDS)
par(mfrow = c(2,3))
MASS::truehist(d$num_CDS, main = "No transformation", nbins = 12)
MASS::truehist(log_obj, main = "log10", nbins = 12)
MASS::truehist(arcsinh_obj$x.t, main = "Arcsinh", nbins = 12)
MASS::truehist(boxcox_obj$x.t, main = "Box Cox", nbins = 12)
MASS::truehist(yeojohnson_obj$x.t, main = "Yeo-Johnson", nbins = 12)
MASS::truehist(orderNorm_obj$x.t, main = "orderNorm", nbins = 12)

# orderNorm Transformation
d$on_genome_size <- (orderNorm(d$genome_size))$x.t
d$on_intron_per_kb <- (orderNorm(d$intron_per_kb))$x.t
d$on_num_CDS <- (orderNorm(d$num_CDS))$x.t
intron <- comparative.data(nt, d, genome)

# phylogenetic independent contrasts
par(mfrow=c(2,2))
fit3 <- crunch(on_genome_size ~ on_intron_per_kb + on_num_CDS, data=intron)
fit3 <- caic.robust(fit3, robust=1)
caic.diagnostics(fit3)
crunchTab <- caic.table(fit3)
plot(on_genome_size ~ on_intron_per_kb + on_num_CDS, crunchTab)
summary(fit3)
# diagnosics
par(mfrow = c(2,3))
hist(fit3$mod$residuals)
plot(fit3)

# correlation tests
cor.test(crunchTab$on_genome_size, crunchTab$on_intron_per_kb, method='spearman')
cor.test(crunchTab$on_genome_size, crunchTab$on_num_CDS, method='spearman')
