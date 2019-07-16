rm(list=ls())
gc(verbose = T, reset = T, full = T)

Packages <- c("ggplot2", "dplyr")

# Load packages -----------------------------------------------
lapply(Packages, library, character.only = TRUE)

# IO ----------------------------------------------------------------------
filepath <- "~/dev/polyoligo/footprinting/runtime.dat"
out <- "~/dev/polyoligo/footprinting/runtime.pdf"
out2 <- "~/dev/polyoligo/footprinting/runtime2.pdf"

# Paths -------------------------------------------------------------------
df <- read.csv(filepath, header=F, sep=" ")
colnames(df) <- c("N_MARKER", "N_CPU", "time")
df <- df[df$N_CPU<50, ]

# Functions ---------------------------------------------------------------
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# Main --------------------------------------------------------------------
df2 <- data_summary(df, varname="time", groupnames = c("N_MARKER","N_CPU"))
df2$N_CPU <- factor(df2$N_CPU)

p <- ggplot(data=df2, aes(x=N_MARKER, y=time/60, col=N_CPU)) + 
  geom_line(size=1) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=(time-sd)/60, ymax=(time+sd)/60), width=.1, size=1) +
  scale_x_continuous(trans='log10') +
  scale_colour_grey() +
  xlab("Number of markers") +
  ylab("Runtime [min]") +
  theme_bw() +
  theme(
    panel.grid.minor=element_blank(),
    text = element_text(size=15, family="sans"),
    legend.position="right"
  ) +
  labs(col="# CPUs")

ggsave(out, p, width=8, height=4.3)

p <- ggplot(data=df2, aes(x=N_MARKER, y=time, col=N_CPU)) + 
  geom_line(size=1) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=time-sd, ymax=time+sd), width=.1, size=1) +
  scale_x_continuous(trans='log10', limits=c(1, 100)) +
  scale_colour_grey() +
  ylab("Seconds") +
  ylim(c(0, 300)) +
  theme_bw() +
  theme(
    panel.grid.minor=element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size=15, family="sans"),
    legend.position="none"
  )

ggsave(out2, p, width=5, height=1.6)
