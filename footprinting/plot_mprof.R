rm(list=ls())
gc(verbose = T, reset = T, full = T)

Packages <- c("ggplot2")

# Load packages -----------------------------------------------
lapply(Packages, library, character.only = TRUE)

# IO ----------------------------------------------------------------------
filepath <- "~/dev/polyoligo/footprinting/footprinting_data/mprof.dat"
out <- "~/dev/polyoligo/footprinting/mprof.pdf"
  
# Paths -------------------------------------------------------------------
df <- read.csv(filepath, header=F, sep=" ")
colnames(df) <- c("CPU", "MEM")

# Main --------------------------------------------------------------------
p <- ggplot(data=df, aes(x=df$CPU, y=df$MEM)) +
  geom_point(size=3) +
  geom_smooth(method='lm',formula=y~x)+
  xlab("Number of CPUs") +
  ylab("RAM [MB]") +
  ylim(c(0, 2500))+
  theme_bw() +
  theme(
    panel.grid.minor=element_blank(),
    text = element_text(size=15, family="sans"),
  ) +
  geom_hline(yintercept=405.947, linetype=2)

ggsave(out, p, width=6.5, height=5)
