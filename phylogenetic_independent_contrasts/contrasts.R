library(phylotools)
library(caper)
library(bestNormalize)
t <- read.tree("~/fungi/ref/silva/phylo.io.tre")
d <- read.table("~/fungi/ref/ensembl/FFP/malin_out/intronPerKb.silva.pic", header=F, sep="\t")
names(d) <- c("species","taxid","genome_size","num_CDS","num_CDS_introns","num_intron_retaining_CDS","intron_per_kb","num_CDS_introns_MCMC","intron_per_kb_MCMC")

#try different transformations
log_obj <- log10(d$intron_per_kb_MCMC)
arcsinh_obj <- arcsinh_x(d$intron_per_kb_MCMC)
boxcox_obj <- boxcox(d$intron_per_kb_MCMC)
yeojohnson_obj <- yeojohnson(d$intron_per_kb_MCMC)
orderNorm_obj <- orderNorm(d$intron_per_kb_MCMC)
par(mfrow = c(2,3))
MASS::truehist(d$intron_per_kb_MCMC, main = "No transformation", nbins = 12)
MASS::truehist(log_obj, main = "log10", nbins = 12)
MASS::truehist(arcsinh_obj$x.t, main = "Arcsinh", nbins = 12)
MASS::truehist(boxcox_obj$x.t, main = "Box Cox", nbins = 12)
MASS::truehist(yeojohnson_obj$x.t, main = "Yeo-Johnson", nbins = 12)
MASS::truehist(orderNorm_obj$x.t, main = "orderNorm", nbins = 12)

# orderNorm Transformation
d$on_genome_size <- (orderNorm(d$genome_size))$x.t
d$on_intron_per_kb_MCMC <- (orderNorm(d$intron_per_kb_MCMC))$x.t
d$on_num_CDS <- (orderNorm(d$num_CDS))$x.t
intron <- comparative.data(t, d, species)
fit1 <- crunch(on_genome_size ~ on_intron_per_kb_MCMC, data=intron, equal.branch.length=T)
fit2 <- crunch(on_genome_size ~ on_num_CDS, data=intron, equal.branch.length=T)
fit3 <- crunch(on_genome_size ~ on_intron_per_kb_MCMC + on_num_CDS, data=intron, equal.branch.length=T)
fit4 <- crunch(on_genome_size ~ 1, data=intron, equal.branch.length=T)
anova(fit1,fit3)
anova(fit2,fit3)
AIC(fit1,fit2,fit3,fit4)
# df       AIC
# fit1  2 76.967852
# fit2  2  5.983495
# fit3  3 -3.703498
# fit4  1 86.481392

par(mfrow=c(3,3))
crunchTab <- caic.table(fit3)
plot(on_genome_size ~ on_intron_per_kb_MCMC + on_num_CDS, crunchTab)
hist(fit3$mod$residuals)
plot(fit3)

summary(fit3)

#correlation analysis
cor.test(crunchTab$on_genome_size, crunchTab$on_intron_per_kb_MCMC, method='spearman')
cor.test(crunchTab$on_genome_size, crunchTab$on_num_CDS, method='spearman')
