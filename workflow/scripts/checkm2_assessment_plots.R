###
### Visualization of bins quality
###

library(ggplot2)
library(vroom)
library(ggpubr)

# retrieving path to the CheckM2 merged report
args <- commandArgs(trailingOnly = TRUE)
path_to_report <- args[1]

# loading the report
report <- vroom(path_to_report, delim="\t")

theme_set(theme_pubclean())

# plotting the results
completeness <- ggplot(report, aes(x=binning, y=Completeness, color=binning)) + 
                geom_boxplot() + 
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
                labs(x = "Tool", y = "Predicted completeness") +
                ylim(0, 100) +
                facet_grid(. ~ assembly)
contamination <- ggplot(report, aes(x=binning, y=Contamination, color=binning)) + 
                 geom_boxplot() + 
                 geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
                 labs(x = "Tool", y = "Predicted contamination") +
                 ylim(0, 100)  +
                 facet_grid(. ~ assembly)                
n50 <- ggplot(report, aes(x=binning, y=Contig_N50, color=binning)) + 
       geom_boxplot() + 
       geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
       labs(x = "Tool", y = "N50") +
       facet_grid(. ~ assembly)
nb_bins <- ggplot(report, aes(x = binning, fill = binning)) +
           geom_bar() + 
           labs(x = "Tool", y = "Number of bins") +
           facet_grid(. ~ assembly)

# making the figure with the plots
arrange <- ggarrange(completeness, contamination, n50, nb_bins,
                     labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2,
                     common.legend = TRUE)

# saving the plot
path_to_save <- args[2]
ggsave(path_to_save, plot = arrange, scale = 2)
