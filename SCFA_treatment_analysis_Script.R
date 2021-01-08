############################################
#Main Figures
butyrate = read.csv('C:/Users/Simon/Downloads/Butyrate_MSc.csv') #<- replace with acetate/popionate paths when needed  
butyrate = butyrate[,c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                       30, 32, 34)]
butyrate$Dose = as.character(butyrate$Dose)
butyrate$Sex = factor(butyrate$Sex, levels = c('Male', 'Female'))
colnames(butyrate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
butyrate$N = as.character(butyrate$N)
### Plotting Time
#Plot function

male = subset(butyrate, butyrate$Sex == 'Male')
male$Dose = as.numeric(male$Dose)
male$Dose
female = subset(butyrate, butyrate$Sex == 'Female')
female$Dose = as.numeric(female$Dose)
female$Dose

#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Butyrate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            #legend.position="topright", #Let's you pick a position for the legend 
            legend.text = element_text(size= 10),
            legend.title = element_text(size = 12),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    geom_point(aes(y = gene, x =  Dose, fill = N), 
               size = 10, shape = 21, show.legend = TRUE) +
    scale_fill_manual(values = c('1' = '#e41a1c', 
                                 '2' = '#377eb8',
                                 '3' = '#4daf4a')) 
    
}
# Aryl-Hydrocarbon Receptor Pathway
ahrM   <- plot_data(male, male$Ahr, 'Ahr')
ahrF   <- plot_data(female, female$Ahr, 'Ahr')
ahrM
ahrF
s100bm <- plot_data(male, male$S100ß, 'S100ß')
s100bf <- plot_data(female, female$S100ß, 'S100ß')
s100bm
s100bf
il_22m <- plot_data(male, male$Il_22, 'Il 22')
il_22f <- plot_data(female, female$Il_22, 'Il 22')
il_22m
il_22f
ifnarm <- plot_data(male, male$Ifnar, 'Ifnar1')
ifnarf <- plot_data(female, female$Ifnar, 'Ifnar1')
ifnarm
ifnarf
gfapm  <- plot_data(male, male$Gfap, 'Gfap')
gfapf  <- plot_data(female, female$Gfap, 'Gfap')
gfapm
gfapf
cyp1b1m<- plot_data(male, male$Cyp1b1, 'Cyp1b1')
cyp1b1f<- plot_data(female, female$Cyp1b1, 'Cyp1b1')
cyp1b1m
cyp1b1f

#GABA/Glutamate Pathway
glulm <- plot_data(male, male$Glul, 'Glul')
glulf <- plot_data(female, female$Glul, 'Glul')
glulm
glulf
gad67m  <- plot_data(male, male$Gad67, 'Gad67')
gad67f  <- plot_data(female, female$Gad67, 'Gad67')
gad67m
gad67f
glud1m <- plot_data(male, male$Glud1, 'Glud')
glud1f <- plot_data(female, female$Glud1, 'Glud')
glud1m
glud1f

#HDACi Pathway
bdnfm <- plot_data(male, male$Bdnf, 'Bdnf')
bdnff <- plot_data(female, female$Bdnf, 'Bdnf') 
bdnfm
bdnff + geom_smooth(method = lm, se = FALSE)   
nfg1m <- plot_data(male, male$Ngf1, 'Ngf1')
nfg1f <- plot_data(female, female$Ngf1, 'Ngf1')
nfg1f
gdnfm <- plot_data(male, male$Gdnf, 'Gdnf')
gdnff <- plot_data(female, female$Gdnf, 'Gdnf')
gdnfm
gdnff
pgc1m <- plot_data(male, male$`Pgc1-a`, 'Pgc1-??')
pgc1f <- plot_data(female, female$`Pgc1-a`, 'Pgc1-??')
pgc1m
pgc1f + geom_smooth(method = lm, se = FALSE)
sp1m  <- plot_data(male, male$Sp1, 'Sp1')
sp1f  <- plot_data(female, female$Sp1, 'Sp1')
sp1m
sp1f

pacman::p_load(car, pheatmap, gplots, RColorBrewer)
myCol <- c('#fff7fb',
           '#ece7f2',
           '#d0d1e6',
           '#a6bddb',
           '#74a9cf',
           '#3690c0',
           '#0570b0',
           '#045a8d',
           '#023858')

myCol
heatmap = read.csv('C:/Users/Simon/Downloads/heatmap_msc_ace.csv')
heatmap
rownames(heatmap) = heatmap$Sex
heatmap = heatmap[,-1]
colnames(heatmap)
myBreaks <- seq(-2, 2, length.out=10)
heatmap = t(heatmap)
is.numeric(heatmap)
heatmap = as.matrix(heatmap)
heatmap = heatmap[,c(2:4,6:8)]
colnames(heatmap)
head(heatmap)
dev.off()
heatmap.2(scale(heatmap),
          Colv = 'none',
          Rowv = 'none',
          key.title = 'none',
          key.xlab = "Scaled Change",
          keysize = 2,
          dendrogram = 'none',
          density = "none",
          col     = myCol,
          breaks = myBreaks,
          trace= "none",
          tracecol="cyan",
          cexRow = 0.6,
          cexCol = 0.1,
          scale = 'none',
          labRow = NULL)
### Acetate
acetate = read.csv('C:/Users/Simon/Downloads/Acetate_MSc.csv') #<- replace with acetate/popionate paths when needed  
acetate = acetate[,c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                       30, 32, 34)]
acetate$Sex = factor(acetate$Sex, levels = c('Male', 'Female'))
colnames(acetate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
acetate$N = as.character(acetate$N)
### Plotting Time
#Plot function

male = subset(acetate, acetate$Sex == 'Male')
male$Dose = as.numeric(male$Dose)
male$Dose
female = subset(acetate, acetate$Sex == 'Female')
female$Dose = as.numeric(female$Dose)
female$Dose

#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Acetate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            #legend.position="topright", #Let's you pick a position for the legend 
            legend.text = element_text(size= 10),
            legend.title = element_text(size = 12),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    geom_point(aes(y = gene, x =  Dose, fill = N), 
               size = 10, shape = 21, show.legend = F) +
    scale_fill_manual(values = c('1' = '#e41a1c', 
                                 '2' = '#377eb8',
                                 '3' = '#4daf4a')) 
  
}
# Aryl-Hydrocarbon Receptor Pathway
ahrM   <- plot_data(male, male$Ahr, 'Ahr')
ahrF   <- plot_data(female, female$Ahr, 'Ahr')
ahrM + geom_smooth(method = lm, se = F)
ahrF
s100bm <- plot_data(male, male$S100ß, 'S100ß')
s100bf <- plot_data(female, female$S100ß, 'S100ß')
s100bm
s100bf
il_22m <- plot_data(male, male$Il_22, 'Il 22')
il_22f <- plot_data(female, female$Il_22, 'Il 22')
il_22m
il_22f
ifnarm <- plot_data(male, male$Ifnar, 'Ifnar1')
ifnarf <- plot_data(female, female$Ifnar, 'Ifnar1')
ifnarm
ifnarf
gfapm  <- plot_data(male, male$Gfap, 'Gfap')
gfapf  <- plot_data(female, female$Gfap, 'Gfap')
gfapm + geom_smooth(method = lm, se = F)
gfapf

cyp1b1m<- plot_data(male, male$Cyp1b1, 'Cyp1b1')
cyp1b1f<- plot_data(female, female$Cyp1b1, 'Cyp1b1')
cyp1b1m
cyp1b1f

#GABA/Glutamate Pathway
glulm <- plot_data(male, male$Glul, 'Glul')
glulf <- plot_data(female, female$Glul, 'Glul')
glulm
glulf
gad67m  <- plot_data(male, male$Gad67, 'Gad67')
gad67f  <- plot_data(female, female$Gad67, 'Gad67')
gad67m
gad67f
glud1m <- plot_data(male, male$Glud1, 'Glud1')
glud1f <- plot_data(female, female$Glud1, 'Glud1')
glud1m
glud1f

#HDACi Pathway
bdnfm <- plot_data(male, male$Bdnf, 'Bdnf')
bdnff <- plot_data(female, female$Bdnf, 'Bdnf') 
bdnfm
bdnff
nfg1m <- plot_data(male, male$Ngf1, 'Ngf1')
nfg1f <- plot_data(female, female$Ngf1, 'Ngf1')
nfg1m
nfg1f
gdnfm <- plot_data(male, male$Gdnf, 'Gdnf')
gdnff <- plot_data(female, female$Gdnf, 'Gdnf')
gdnfm
gdnff
pgc1m <- plot_data(male, male$`Pgc1-a`, 'Pgc1-??')
pgc1f <- plot_data(female, female$`Pgc1-a`, 'Pgc1-??')
pgc1m
pgc1f 
sp1m  <- plot_data(male, male$Sp1, 'Sp1')
sp1f  <- plot_data(female, female$Sp1, 'Sp1')
sp1m
sp1f

pacman::p_load(car, pheatmap, gplots, RColorBrewer)
myCol <- c('#fff7fb',
           '#ece7f2',
           '#d0d1e6',
           '#a6bddb',
           '#74a9cf',
           '#3690c0',
           '#0570b0',
           '#045a8d',
           '#023858')

myCol
heatmap = read.csv('C:/Users/Simon/Downloads/heatmap_msc_ace.csv')
heatmap
rownames(heatmap) = heatmap$Sex
heatmap = heatmap[,-1]
colnames(heatmap)
myBreaks <- seq(-2, 2, length.out=10)
heatmap = t(heatmap)
is.numeric(heatmap)
heatmap = as.matrix(heatmap)
heatmap = heatmap[,c(2:4,6:8)]
colnames(heatmap)
head(heatmap)
dev.off()
heatmap.2(scale(heatmap),
          Colv = 'none',
          Rowv = 'none',
          key.title = 'none',
          key.xlab = "Scaled Change",
          keysize = 2,
          dendrogram = 'none',
          density = "none",
          col     = myCol,
          breaks = myBreaks,
          trace= "none",
          tracecol="cyan",
          cexRow = 0.6,
          cexCol = 0.1,
          scale = 'none',
          labRow = NULL)

### Propionate
propionate = read.csv('C:/Users/Simon/Downloads/Propionate_MSc.csv') #<- replace with acetate/popionate paths when needed  
propionate = propionate[,c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                     30, 32, 34)]
propionate$Sex = factor(propionate$Sex, levels = c('Male', 'Female'))
colnames(propionate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
propionate$N = as.character(propionate$N)
### Plotting Time
#Plot function

male = subset(propionate, propionate$Sex == 'Male')
male$Dose = as.numeric(male$Dose)
male$Dose
female = subset(propionate, propionate$Sex == 'Female')
female$Dose = as.numeric(female$Dose)
female$Dose

#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Propionate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            #legend.position="topright", #Let's you pick a position for the legend 
            legend.text = element_text(size= 10),
            legend.title = element_text(size = 12),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    geom_point(aes(y = gene, x =  Dose, fill = N), 
               size = 10, shape = 21, show.legend = F) +
    scale_fill_manual(values = c('1' = '#e41a1c', 
                                 '2' = '#377eb8',
                                 '3' = '#4daf4a')) 
  
}
# Aryl-Hydrocarbon Receptor Pathway
ahrM   <- plot_data(male, male$Ahr, 'Ahr')
ahrF   <- plot_data(female, female$Ahr, 'Ahr')
ahrM
ahrF
s100bm <- plot_data(male, male$S100ß, 'S100ß')
s100bf <- plot_data(female, female$S100ß, 'S100ß')
s100bm
s100bf
il_22m <- plot_data(male, male$Il_22, 'Il 22')
il_22f <- plot_data(female, female$Il_22, 'Il 22')
il_22m
il_22f
ifnarm <- plot_data(male, male$Ifnar, 'Ifnar1')
ifnarf <- plot_data(female, female$Ifnar, 'Ifnar1')
ifnarm
ifnarf
gfapm  <- plot_data(male, male$Gfap, 'Gfap')
gfapf  <- plot_data(female, female$Gfap, 'Gfap')
gfapm
gfapf
cyp1b1m<- plot_data(male, male$Cyp1b1, 'Cyp1b1')
cyp1b1f<- plot_data(female, female$Cyp1b1, 'Cyp1b1')
cyp1b1m
cyp1b1f

#GABA/Glutamate Pathway
glulm <- plot_data(male, male$Glul, 'Glul')
glulf <- plot_data(female, female$Glul, 'Glul')
glulm
glulf
gad67m  <- plot_data(male, male$Gad67, 'Gad67')
gad67f  <- plot_data(female, female$Gad67, 'Gad67')
gad67m
gad67f
glud1m <- plot_data(male, male$Glud1, 'Glud1')
glud1f <- plot_data(female, female$Glud1, 'Glud1')
glud1m
glud1f

#HDACi Pathway
bdnfm <- plot_data(male, male$Bdnf, 'Bdnf')
bdnff <- plot_data(female, female$Bdnf, 'Bdnf') 
bdnfm
bdnff
nfg1m <- plot_data(male, male$Ngf1, 'Ngf1')
nfg1f <- plot_data(female, female$Ngf1, 'Ngf1')
nfg1m
nfg1f
gdnfm <- plot_data(male, male$Gdnf, 'Gdnf')
gdnff <- plot_data(female, female$Gdnf, 'Gdnf')
gdnfm
gdnff
pgc1m <- plot_data(male, male$`Pgc1-a`, 'Pgc1-??')
pgc1f <- plot_data(female, female$`Pgc1-a`, 'Pgc1-??')
pgc1m
pgc1f 
sp1m  <- plot_data(male, male$Sp1, 'Sp1')
sp1f  <- plot_data(female, female$Sp1, 'Sp1')
sp1m
sp1f

pacman::p_load(car, pheatmap, gplots, RColorBrewer)
myCol <- c('#fff7fb',
           '#ece7f2',
           '#d0d1e6',
           '#a6bddb',
           '#74a9cf',
           '#3690c0',
           '#0570b0',
           '#045a8d',
           '#023858')

myCol
heatmap = read.csv('C:/Users/Simon/Downloads/heatmap_msc_ace.csv')
heatmap
rownames(heatmap) = heatmap$Sex
heatmap = heatmap[,-1]
colnames(heatmap)
myBreaks <- seq(-2, 2, length.out=10)
heatmap = t(heatmap)
is.numeric(heatmap)
heatmap = as.matrix(heatmap)
heatmap = heatmap[,c(2:4,6:8)]
colnames(heatmap)
head(heatmap)
dev.off()
heatmap.2(scale(heatmap),
          Colv = 'none',
          Rowv = 'none',
          key.title = 'none',
          key.xlab = "Scaled Change",
          keysize = 2,
          dendrogram = 'none',
          density = "none",
          col     = myCol,
          breaks = myBreaks,
          trace= "none",
          tracecol="cyan",
          cexRow = 0.6,
          cexCol = 0.1,
          scale = 'none',
          labRow = NULL) Main Figures
############################################ 
#Supplemental Figures
#Butyrate
butyrate = read.csv('C:/Users/Simon/Downloads/Butyrate_MSc.csv') #<- replace with acetate/popionate paths when needed  
butyrate = butyrate[,c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                       30, 32, 34)]
butyrate$Dose = as.character(butyrate$Dose)
butyrate$Sex = factor(butyrate$Sex, levels = c('Male', 'Female'))
colnames(butyrate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
butyrate$N = as.character(butyrate$N)
butyrate$Sex = factor(butyrate$Sex, levels = c('Male', 'Female'))
butyrate$Dose = factor(butyrate$Dose, levels = c('0', '2.5', '12.5', '25'))

### Plotting Time
#Plot function
n1 = subset(butyrate, butyrate$N == '1')
n2 = subset(butyrate, butyrate$N == '2')
n3 = subset(butyrate, butyrate$N == '3')
#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Butyrate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            #legend.position="topright", #Let's you pick a position for the legend 
            legend.text = element_text(size= 10),
            legend.title = element_text(size = 12),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    geom_bar(aes(y = gene, x =  Dose, fill = Sex), width = 0.8, 
             color = 'black',
             stat = "summary",  position = position_dodge2()) +
    scale_fill_manual(values = c('Male' = '#e41a1c', 
                                 'Female' = '#377eb8')) +
    scale_y_continuous(expand = c(0,0))
  
}
# Aryl-Hydrocarbon Receptor Pathway
ahr1   <- plot_data(n1, n1$Ahr, 'Ahr')
ahr2   <- plot_data(n2, n2$Ahr, 'Ahr')
ahr3   <- plot_data(n3, n3$Ahr, 'Ahr')
s100b1 <- plot_data(n1, n1$S100ß, 'S100ß')
s100b2 <- plot_data(n2, n2$S100ß, 'S100ß')
s100b3 <- plot_data(n3, n3$S100ß, 'S100ß')
il_221 <- plot_data(n1, n1$Il_22, 'Il 22')
il_222 <- plot_data(n2, n2$Il_22, 'Il 22')
il_223 <- plot_data(n3, n3$Il_22, 'Il 22')
ifnar1 <- plot_data(n1, n1$Ifnar, 'Ifnar1')
ifnar2 <- plot_data(n2, n2$Ifnar, 'Ifnar1')
ifnar3 <- plot_data(n3, n3$Ifnar, 'Ifnar1')
gfap1  <- plot_data(n1, n1$Gfap, 'Gfap')
gfap2  <- plot_data(n2, n2$Gfap, 'Gfap')
gfap3  <- plot_data(n3, n3$Gfap, 'Gfap')
cyp1b11<- plot_data(n1, n1$Cyp1b1, 'Cyp1b1')
cyp1b12<- plot_data(n2, n2$Cyp1b1, 'Cyp1b1')
cyp1b13<- plot_data(n3, n3$Cyp1b1, 'Cyp1b1')
#GABA/Glutamate Pathway
glul1  <- plot_data(n1, n1$Glul, 'Glul')
glul2  <- plot_data(n2, n2$Glul, 'Glul')
glul3  <- plot_data(n3, n3$Glul, 'Glul')
gad671 <- plot_data(n1, n1$Gad67, 'Gad67')
gad672 <- plot_data(n2, n2$Gad67, 'Gad67')
gad673 <- plot_data(n3, n3$Gad67, 'Gad67')
glud11 <- plot_data(n1, n1$Glud1, 'Glud1')
glud12 <- plot_data(n2, n2$Glud1, 'Glud1')
glud13 <- plot_data(n3, n3$Glud1, 'Glud1')
#HDACi Pathway
bdnf1 <- plot_data(n1, n1$Bdnf, 'Bdnf')
bdnf2 <- plot_data(n2, n2$Bdnf, 'Bdnf') 
bdnf3 <- plot_data(n3, n3$Bdnf, 'Bdnf') 
nfg11 <- plot_data(n1, n1$Ngf1, 'Ngf1')
nfg12 <- plot_data(n2, n2$Ngf1, 'Ngf1')
nfg13 <- plot_data(n3, n3$Ngf1, 'Ngf1')
gdnf1 <- plot_data(n1, n1$Gdnf, 'Gdnf')
gdnf2 <- plot_data(n2, n2$Gdnf, 'Gdnf')
gdnf3 <- plot_data(n3, n3$Gdnf, 'Gdnf')
pgc11 <- plot_data(n1, n1$`Pgc1-a`, 'Pgc1-??')
pgc12 <- plot_data(n2, n2$`Pgc1-a`, 'Pgc1-??')
pgc13 <- plot_data(n3, n3$`Pgc1-a`, 'Pgc1-??')
sp11  <- plot_data(n1, n1$Sp1, 'Sp1')
sp12  <- plot_data(n2, n2$Sp1, 'Sp1')
sp13  <- plot_data(n3, n3$Sp1, 'Sp1')

pacman::p_load('patchwork')
ahr1 + ahr2 + ahr3 + s100b1 + s100b2 + s100b3 + il_221 + il_222 + il_223 +
  ifnar1 + ifnar2 + ifnar3 + gfap1 + gfap2 + gfap3 + cyp1b11 + cyp1b12 + cyp1b13 + 
  glul1 + glul2 + glul3 + gad671 + gad672 + gad673 + glud11 + glud12 + glud13 +
  bdnf1 + bdnf2 + bdnf3 + nfg11 + nfg12 + nfg13 + gdnf1 + gdnf2 + gdnf3 +
  pgc11 + pgc12 + pgc13 + sp11 + sp12 + sp13 + plot_layout(ncol = 3)

#Acetate
acetate = read.csv('C:/Users/Simon/Downloads/Acetate_MSc.csv') #<- replace with acetate/popionate paths when needed  
acetate = acetate[,c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                       30, 32, 34)]
acetate$Dose = as.character(acetate$Dose)
acetate$Sex = factor(acetate$Sex, levels = c('Male', 'Female'))
colnames(acetate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
acetate$N = as.character(acetate$N)
acetate$Sex = factor(acetate$Sex, levels = c('Male', 'Female'))
acetate$Dose = factor(acetate$Dose, levels = c('0', '150', '750', '1500'))

### Plotting Time
#Plot function
n1 = subset(acetate, acetate$N == '1')
n2 = subset(acetate, acetate$N == '2')
n3 = subset(acetate, acetate$N == '3')
#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Acetate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            legend.position="none", #Let's you pick a position for the legend 
            legend.text = element_text(size= 10),
            legend.title = element_text(size = 12),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    geom_bar(aes(y = gene, x =  Dose, fill = Sex), width = 0.8, 
             color = 'black',
             stat = "summary",  position = position_dodge2()) +
    scale_fill_manual(values = c('Male' = '#e41a1c', 
                                 'Female' = '#377eb8')) +
    scale_y_continuous(expand = c(0,0))
  
}
# Aryl-Hydrocarbon Receptor Pathway
ahr1   <- plot_data(n1, n1$Ahr, 'Ahr')
ahr2   <- plot_data(n2, n2$Ahr, 'Ahr')
ahr3   <- plot_data(n3, n3$Ahr, 'Ahr')
s100b1 <- plot_data(n1, n1$S100ß, 'S100ß')
s100b2 <- plot_data(n2, n2$S100ß, 'S100ß')
s100b3 <- plot_data(n3, n3$S100ß, 'S100ß')
il_221 <- plot_data(n1, n1$Il_22, 'Il 22')
il_222 <- plot_data(n2, n2$Il_22, 'Il 22')
il_223 <- plot_data(n3, n3$Il_22, 'Il 22')
ifnar1 <- plot_data(n1, n1$Ifnar, 'Ifnar1')
ifnar2 <- plot_data(n2, n2$Ifnar, 'Ifnar1')
ifnar3 <- plot_data(n3, n3$Ifnar, 'Ifnar1')
gfap1  <- plot_data(n1, n1$Gfap, 'Gfap')
gfap2  <- plot_data(n2, n2$Gfap, 'Gfap')
gfap3  <- plot_data(n3, n3$Gfap, 'Gfap')
cyp1b11<- plot_data(n1, n1$Cyp1b1, 'Cyp1b1')
cyp1b12<- plot_data(n2, n2$Cyp1b1, 'Cyp1b1')
cyp1b13<- plot_data(n3, n3$Cyp1b1, 'Cyp1b1')
#GABA/Glutamate Pathway
glul1  <- plot_data(n1, n1$Glul, 'Glul')
glul2  <- plot_data(n2, n2$Glul, 'Glul')
glul3  <- plot_data(n3, n3$Glul, 'Glul')
gad671 <- plot_data(n1, n1$Gad67, 'Gad67')
gad672 <- plot_data(n2, n2$Gad67, 'Gad67')
gad673 <- plot_data(n3, n3$Gad67, 'Gad67')
glud11 <- plot_data(n1, n1$Glud1, 'Glud1')
glud12 <- plot_data(n2, n2$Glud1, 'Glud1')
glud13 <- plot_data(n3, n3$Glud1, 'Glud1')
#HDACi Pathway
bdnf1 <- plot_data(n1, n1$Bdnf, 'Bdnf')
bdnf2 <- plot_data(n2, n2$Bdnf, 'Bdnf') 
bdnf3 <- plot_data(n3, n3$Bdnf, 'Bdnf') 
nfg11 <- plot_data(n1, n1$Ngf1, 'Ngf1')
nfg12 <- plot_data(n2, n2$Ngf1, 'Ngf1')
nfg13 <- plot_data(n3, n3$Ngf1, 'Ngf1')
gdnf1 <- plot_data(n1, n1$Gdnf, 'Gdnf')
gdnf2 <- plot_data(n2, n2$Gdnf, 'Gdnf')
gdnf3 <- plot_data(n3, n3$Gdnf, 'Gdnf')
pgc11 <- plot_data(n1, n1$`Pgc1-a`, 'Pgc1-??')
pgc12 <- plot_data(n2, n2$`Pgc1-a`, 'Pgc1-??')
pgc13 <- plot_data(n3, n3$`Pgc1-a`, 'Pgc1-??')
sp11  <- plot_data(n1, n1$Sp1, 'Sp1')
sp12  <- plot_data(n2, n2$Sp1, 'Sp1')
sp13  <- plot_data(n3, n3$Sp1, 'Sp1')
ahr1 + ahr2 + ahr3 + s100b1 + s100b2 + s100b3 + il_221 + il_222 + il_223 +
  ifnar1 + ifnar2 + ifnar3 + gfap1 + gfap2 + gfap3 + cyp1b11 + cyp1b12 + cyp1b13 + 
  glul1 + glul2 + glul3 + gad671 + gad672 + gad673 + glud11 + glud12 + glud13 +
  bdnf1 + bdnf2 + bdnf3 + nfg11 + nfg12 + nfg13 + 
  pgc11 + pgc12 + pgc13 + sp11 + sp12 + sp13 + plot_layout(ncol = 3)


#Propionate
propionate = read.csv('C:/Users/Simon/Downloads/Propionate_MSc.csv') #<- replace with acetate/popionate paths when needed  
propionate = propionate[c(1:24),c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                     30, 32, 34)]
propionate$Dose = as.character(propionate$Dose)
propionate$Sex = factor(propionate$Sex, levels = c('Male', 'Female'))
colnames(propionate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
propionate$N = as.character(propionate$N)
propionate$Sex = factor(propionate$Sex, levels = c('Male', 'Female'))
propionate$Dose = factor(propionate$Dose, levels = c('0', '3.5', '17.5', '35'))
### Plotting Time
#Plot function
n1 = subset(propionate, propionate$N == '1')
n2 = subset(propionate, propionate$N == '2')
n3 = subset(propionate, propionate$N == '3')
#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Propionate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            legend.position="none", #Let's you pick a position for the legend 
            legend.text = element_text(size= 10),
            legend.title = element_text(size = 12),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    geom_bar(aes(y = gene, x =  Dose, fill = Sex), width = 0.8, 
             color = 'black',
             stat = "summary",  position = position_dodge2()) +
    scale_fill_manual(values = c('Male' = '#e41a1c', 
                                 'Female' = '#377eb8')) +
    scale_y_continuous(expand = c(0,0))
  
}
# Aryl-Hydrocarbon Receptor Pathway
ahr1   <- plot_data(n1, n1$Ahr, 'Ahr')
ahr2   <- plot_data(n2, n2$Ahr, 'Ahr')
ahr3   <- plot_data(n3, n3$Ahr, 'Ahr')
s100b1 <- plot_data(n1, n1$S100ß, 'S100ß')
s100b2 <- plot_data(n2, n2$S100ß, 'S100ß')
s100b3 <- plot_data(n3, n3$S100ß, 'S100ß')
il_221 <- plot_data(n1, n1$Il_22, 'Il 22')
il_222 <- plot_data(n2, n2$Il_22, 'Il 22')
il_223 <- plot_data(n3, n3$Il_22, 'Il 22')
ifnar1 <- plot_data(n1, n1$Ifnar, 'Ifnar1')
ifnar2 <- plot_data(n2, n2$Ifnar, 'Ifnar1')
ifnar3 <- plot_data(n3, n3$Ifnar, 'Ifnar1')
gfap1  <- plot_data(n1, n1$Gfap, 'Gfap')
gfap2  <- plot_data(n2, n2$Gfap, 'Gfap')
gfap3  <- plot_data(n3, n3$Gfap, 'Gfap')
cyp1b11<- plot_data(n1, n1$Cyp1b1, 'Cyp1b1')
cyp1b12<- plot_data(n2, n2$Cyp1b1, 'Cyp1b1')
cyp1b13<- plot_data(n3, n3$Cyp1b1, 'Cyp1b1')
#GABA/Glutamate Pathway
glul1  <- plot_data(n1, n1$Glul, 'Glul')
glul2  <- plot_data(n2, n2$Glul, 'Glul')
glul3  <- plot_data(n3, n3$Glul, 'Glul')
gad671 <- plot_data(n1, n1$Gad67, 'Gad67')
gad672 <- plot_data(n2, n2$Gad67, 'Gad67')
gad673 <- plot_data(n3, n3$Gad67, 'Gad67')
glud11 <- plot_data(n1, n1$Glud1, 'Glud1')
glud12 <- plot_data(n2, n2$Glud1, 'Glud1')
glud13 <- plot_data(n3, n3$Glud1, 'Glud1')
#HDACi Pathway
bdnf1 <- plot_data(n1, n1$Bdnf, 'Bdnf')
bdnf2 <- plot_data(n2, n2$Bdnf, 'Bdnf') 
bdnf3 <- plot_data(n3, n3$Bdnf, 'Bdnf') 
nfg11 <- plot_data(n1, n1$Ngf1, 'Ngf1')
nfg12 <- plot_data(n2, n2$Ngf1, 'Ngf1')
nfg13 <- plot_data(n3, n3$Ngf1, 'Ngf1')
gdnf1 <- plot_data(n1, n1$Gdnf, 'Gdnf')
gdnf2 <- plot_data(n2, n2$Gdnf, 'Gdnf')
gdnf3 <- plot_data(n3, n3$Gdnf, 'Gdnf')
pgc11 <- plot_data(n1, n1$`Pgc1-a`, 'Pgc1-??')
pgc12 <- plot_data(n2, n2$`Pgc1-a`, 'Pgc1-??')
pgc13 <- plot_data(n3, n3$`Pgc1-a`, 'Pgc1-??')
sp11  <- plot_data(n1, n1$Sp1, 'Sp1')
sp12  <- plot_data(n2, n2$Sp1, 'Sp1')
sp13  <- plot_data(n3, n3$Sp1, 'Sp1')
ahr1 + ahr2 + ahr3 + s100b1 + s100b2 + s100b3 + il_221 + il_222 + il_223 +
  ifnar1 + ifnar2 + ifnar3 + gfap1 + gfap2 + gfap3 + cyp1b11 + cyp1b12 + cyp1b13 + 
  glul1 + glul2 + glul3 + gad671 + gad672 + gad673 + glud11 + glud12 + glud13 +
  bdnf1 + bdnf2 + bdnf3 + nfg11 + nfg12 + nfg13 + 
  pgc11 + pgc12 + pgc13 + sp11 + sp12 + sp13 + plot_layout(ncol = 3)



############################################
#Repeat for butyrate, acetate, propionate
m = subset(butyrate, butyrate$Sex == 'Male')
f = subset(butyrate, butyrate$Sex == 'Female')
one   <- (lm(Ahr ~ Dose*N, f))
summary(lm(Ahr ~ Dose*N, f))
two   <- (lm(S100ß ~ Dose, f))
three <- (lm(Gfap ~ Dose, f))
frr   <- (lm(Ifnar ~ Dose, f))
fiv   <- (lm(Il_22 ~ Dose, f))
six   <- (lm(Cyp1b1 ~ Dose, f))
sev   <- (lm(Glul ~ Dose, f))
eig   <- (lm(Bdnf ~ Dose, f))
nin   <- (lm(Ngf1 ~ Dose, f))
ten   <- (lm(Gdnf ~ Dose, f))
elvn  <- (lm(Gad67 ~ Dose, f))
twlv  <- (lm(Glud1 ~ Dose, f))
thrt  <- (lm(`Pgc1-a` ~ Dose, f))
frtn  <- (lm(Sp1 ~ Dose, f))
tab_model(one, two, three, frr, fiv,
          six, sev, eig, nin, ten,
          elvn, twlv, thrt, frtn, 
          show.fstat = T, collapse.ci = T,
          show.intercept = F, file = 'female_but.html')

one   <- (lm(Ahr ~ Dose, m))
two   <- (lm(S100ß ~ Dose, m))
three <- (lm(Gfap ~ Dose, m))
frr   <- (lm(Ifnar ~ Dose, m))
fiv   <- (lm(Il_22 ~ Dose, m))
six   <- (lm(Cyp1b1 ~ Dose, m))
sev   <- (lm(Glul ~ Dose, m))
eig   <- (lm(Bdnf ~ Dose, m))
nin   <- (lm(Ngf1 ~ Dose, m))
ten   <- (lm(Gdnf ~ Dose, m))
elvn  <- (lm(Gad67 ~ Dose, m))
twlv  <- (lm(Glud1 ~ Dose, m))
thrt  <- (lm(`Pgc1-a` ~ Dose, m))
frtn  <- (lm(Sp1 ~ Dose, m))
tab_model(one, two, three, frr, fiv,
          six, sev, eig, nin, ten,
          elvn, twlv, thrt, frtn, 
          show.fstat = T, collapse.ci = T,
          show.intercept = F, file = 'male_but.html')





