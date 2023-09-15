library(here)
library(coda)
library(rcarbon)
source(here('src','utility.R'))
load(here('results','post_jp_abot.RData'))
load(here('results','post_gb_abot.RData'))

hpd.r.jp <- post.sample.core.jp.abot[,'r'] |> as.mcmc() |> HPDinterval() |> as.numeric() |> round(4) |> paste0(collapse='~')
median.r.jp <- post.sample.core.jp.abot[,'r'] |> as.mcmc() |> median() |> as.numeric() |> round(4)
hpd.k.jp <- post.sample.core.jp.abot[,'mu_k'] |> as.mcmc() |> HPDinterval() |> as.numeric() |> round(3) |> paste0(collapse='~')
median.k.jp <- post.sample.core.jp.abot[,'mu_k'] |> as.mcmc() |> median() |> as.numeric() |> round(3)
hpd.phi.jp  <- post.sample.core.jp.abot[,'phi'] |> as.mcmc() |> HPDinterval() |> round(2) |> paste0(collapse='~')
median.phi.jp  <- post.sample.core.jp.abot[,'phi'] |> as.mcmc() |> median() |> round(2)
hpd.m.jp <- post.sample.core.jp.abot[,'m'] |> as.mcmc() |> HPDinterval() |> as.numeric() |> round() |> BPtoBCAD() |> sort()
hpd.m.jp <- ifelse(hpd.m.jp<0,paste0('BC',abs(hpd.m.jp)),paste0('AD',hpd.m.jp)) |> paste0(collapse='~')
median.m.jp <- post.sample.core.jp.abot[,'m'] |> as.mcmc() |> median() |> as.numeric() |> round() |> BPtoBCAD()
median.m.jp  <- ifelse(median.m.jp<0,paste0('BC',abs(median.m.jp)),paste0('AD',median.m.jp))
rhats.jp  <- rhats.jp.abot[[1]][c('r','m','mu_k','phi'),1]

hpd.r.gb <- post.sample.core.gb.abot[,'r'] |> as.mcmc() |> HPDinterval() |> as.numeric() |> round(4) |> paste0(collapse='~')
median.r.gb <- post.sample.core.gb.abot[,'r'] |> as.mcmc() |> median() |> as.numeric() |> round(4)
hpd.k.gb <- post.sample.core.gb.abot[,'mu_k'] |> as.mcmc() |> HPDinterval() |> as.numeric() |> round(3) |> paste0(collapse='~')
median.k.gb <- post.sample.core.gb.abot[,'mu_k'] |> as.mcmc() |> median() |> as.numeric() |> round(3)
hpd.phi.gb  <- post.sample.core.gb.abot[,'phi'] |> as.mcmc() |> HPDinterval() |> round(2) |> paste0(collapse='~')
median.phi.gb  <- post.sample.core.gb.abot[,'phi'] |> as.mcmc() |> median() |> round(2)
hpd.m.gb <- post.sample.core.gb.abot[,'m'] |> as.mcmc() |> HPDinterval() |> as.numeric() |> round() |> BPtoBCAD() |> sort()
hpd.m.gb <- ifelse(hpd.m.gb<0,paste0('BC',abs(hpd.m.gb)),paste0('AD',hpd.m.gb)) |> paste0(collapse='~')
median.m.gb <- post.sample.core.gb.abot[,'m'] |> as.mcmc() |> median() |> as.numeric() |> round() |> BPtoBCAD()
median.m.gb  <- ifelse(median.m.gb<0,paste0('BC',abs(median.m.gb)),paste0('AD',median.m.gb))
rhats.gb  <- rhats.gb.abot[[1]][c('r','m','mu_k','phi'),1]



table1  <- data.frame(Region=c(rep('Japan',4),rep('Britain',4)),Parameters=rep(c('r','m','mu','phi'),2))
table1$Median <- c(median.r.jp,median.m.jp,median.k.jp,median.phi.jp,median.r.gb,median.m.gb,median.k.gb,median.phi.gb)
table1$HPD <- c(hpd.r.jp,hpd.m.jp,hpd.k.jp,hpd.phi.jp,hpd.r.gb,hpd.m.gb,hpd.k.gb,hpd.phi.gb)
table1$Rhats  <- c(rhats.jp,rhats.gb) |> round(4)


write.csv(table1,row.names=FALSE,here('figures_and_tables','table1.csv'))

