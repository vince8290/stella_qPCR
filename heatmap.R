# Clear global environment
rm(list = ls())


##############################################################################################
# 1.Configuration
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set current dir as WD

xlabel = "putative MAX effector"
ylabel = "dpi"

rawdatafolder = "example3"         # rawfolder: name of the folder containing the raw data files (ex: "example1")
outputfolder = file.path(rawdatafolder, format(Sys.time(), "Results_heatmap_%Y%m%d_%H%M%S"))  # outfolder: name of the folder to which graphs are exported
raw_ext = "csv"  # raw_ext: Raw datafiles extension ("csv", "tsv", "txt") WITHOUT the dot (.)
sep = ";"   # column separator for raw data files (";" or "," or "\t")
dec = ","   # decimal character for raw data files("." or ",")
# ctmax = 33 # ctmax : maximal Ct value to consider. beyond this value, expression will be 0. no cutoff applied if ctmax=0
required_packages = c("tidyverse","tools")
columns_to_keep = c("id","cultivar","gene","rep","dpi", "Cp_gene","Cp_ref","Expression","exclude","E_gene","E_ref","file")
categories = c("< 0.008", "0.008 - 0.04", "0.04 - 0.2", "0.2 - 1", "1 - 5", "> 5")

#---------------------------------------------------------------------------------------------
# Graphics parameters
graph_unit = "cm"   # unit for graphs sizes ("cm", "in" or "mm" )
graph_width = 50    # graph width
graph_height = 20   # graph height
graph_dpi = 600     # graph resolution (typically between 100 and 1200)
plot_device = "png" # graph file extension

#---------------------------------------------------------------------------------------------
# Data Transformation (pick one)
# Transform the data for statistical analyses (choose below which transformation is best suited)
# trans <- function(x) x # No transformation (raw data)
# trans <- function(x) x^0.5 # Power Transformation
# trans <- log10 # Log10 Transformation
trans <- log2 # Log2 Transformation
#---------------------------------------------------------------------------------------------





##############################################################################################
# Run script

for (pkg in required_packages) {
  if(pkg %in% rownames(installed.packages())) {
    print(paste0(c(pkg, " - installed")))
    library(pkg, character.only = TRUE)
  } else {
    print(paste0(c(pkg, " - not installed")))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

##############################################################################################
# 2.Compile data
rawfolder = file.path(getwd(), rawdatafolder)
outfolder = file.path(getwd(), outputfolder)

# create results folder
dir.create(outfolder, showWarnings = FALSE)

# Get list of rawdata files
filelist = dir(path=rawfolder, pattern=paste0("*",".",raw_ext))

# Compile all rawdata files
rawdf <- data.frame()

# append each rawdata file to dataframe
for (i in 1:length(filelist)){
  tmp = read.csv2(file.path(rawfolder, filelist[i]), sep=sep, dec=dec)
  tmp$file = filelist[i]
  rawdf <- rbind(rawdf, tmp)
}


# data cleaning
rawdf$id = row.names(rawdf)

df <- rawdf[,columns_to_keep]
names(df)[names(df) == 'Expression'] <- 'old.expression'

# create gene.dpi
df$gene.dpi <- paste0(df$gene,".",df$dpi)
# unique(df[c("dpi","gene","gene.dpi")]) # Check if df is ok

# Convert columns format to factors
df$dpi <- factor(df$dpi, levels=unique(df$dpi))
df$cultivar <- factor(df$cultivar, levels=unique(df$cultivar))
df$gene <- factor(df$gene, levels=unique(df$gene))
df$gene.dpi <- factor(df$gene.dpi, levels=unique(df$gene.dpi))

# Recalculate expression level in R
df = df[is.na(df$exclude),]
df$exp_gene = (df$E_gene^(-df$Cp_gene))
df$exp_ref = (df$E_ref^(-df$Cp_ref))

df$expression = df$exp_gene / df$exp_ref

# convert Cp > ctmax to expression=0
# df[df$Cp_gene >=ctmax, "expression"] = 0

# check if calculated expression from R correspond to expression from excel
df$diff = df$old.expression - df$expression
sum(df$diff)  # must be 0

# Export final dataframe to .RData file
# save(df, file=file.path(outfolder,"df.RData"))
# write_excel_csv2(df, path=file.path(outfolder,"df.csv"))

##############################################################################################
# Exclude aberrant data
##############################################################################################

# calculate maximum median within each cultivar and gene and sort by descending max median
df2 = df %>%
  group_by(cultivar, gene, dpi) %>%
  summarise(med=median(expression, na.rm=TRUE), count=n()) %>%
  arrange(cultivar, gene, dpi) %>%
  group_by(cultivar, gene) %>%
  summarise(max_med=max(med)) %>%
  arrange(cultivar, gene, max_med)

df3 = merge(df, df2) %>%
  arrange(max_med, gene, rep, dpi)

cult="Maratelli"


# Compute means and medians for heatmap
df_heatmap = df3 %>%
  group_by(cultivar, gene, dpi, rep) %>%
  summarise(mean_exp=mean(expression),
            median_exp=median(expression), 
            count=n(),
            max_med=max_med) %>%
  arrange(cultivar, gene, dpi) %>%
  group_by(cultivar, gene) %>%
  ungroup

# create groups based on median
df_heatmap$heatmap=paste0("> 5")
df_heatmap[df_heatmap$median_exp<5,"heatmap"]="1 - 5"
df_heatmap[df_heatmap$median_exp<1,"heatmap"]="0.2 - 1"
df_heatmap[df_heatmap$median_exp<0.2,"heatmap"]="0.04 - 0.2"
df_heatmap[df_heatmap$median_exp<0.04,"heatmap"]="0.008 - 0.04"
df_heatmap[df_heatmap$median_exp<0.008,"heatmap"]="< 0.008"

# Order groups correctly
df_heatmap$heatmap = factor(df_heatmap$heatmap)
df_heatmap <- df_heatmap %>%
  arrange(max_med) %>%
  mutate(heatmap = factor(heatmap, levels=categories))



# Export final dataframe to .RData file
save(df_heatmap, file=file.path(outfolder,"df_heatmap.RData"))
write_excel_csv2(df_heatmap, path=file.path(outfolder,"df_heatmap.csv"))


#draw heatmap based on medians
ggplot(df_heatmap[(df_heatmap$cultivar==cult),],
       aes(x=reorder(gene, max_med), 
           y=dpi,
           fill=heatmap)) + 
  geom_tile() + theme_bw()+
  theme(axis.text = element_text(face="bold",size=10), axis.text.x=element_text(angle=90)) +
  facet_grid(rep~.) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(categories), "YlOrRd"))+
  ggtitle(paste0(cult))+
  labs(x=xlabel, y=ylabel) +
  theme(axis.text.x = element_text(margin=margin(10,0,0,0)))

# save heatmap to file
fn = paste0(cult, "_heatmap_med")
ggsave(filename= paste0(fn, ".", plot_device),
       path = outfolder,
       device=plot_device,
       width=graph_width,
       height=graph_height,
       unit=graph_unit,
       dpi=graph_dpi)  