butyrate = read.csv(PathToButyrateDataSet) #<- replace with acetate/popionate paths when needed  
butyrate = butyrate[,c(1:5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28,
                       30, 32, 34)]
butyrate$Dose = as.character(butyrate$Dose)
butyrate$Dose = factor(butyrate$Dose, levels = c('0', '2.5', '12.5', '25'))
butyrate$Sex = factor(butyrate$Sex, levels = c('Male', 'Female'))
colnames(butyrate) = c(
  'ID', 'Sex', 'N', 'SCFA', 'Dose',
  'Ahr', 'S100ß', 'Gfap', 'Ifnar',
  'Il_22', 'Cyp1b1', 'Glul', 'Bdnf', 'Ngf1',
  'Gdnf', 'Gad67', 'Glud1', 'Pgc1-a', 'Sp1'
)
butyrate = as.data.frame(butyrate)
butyrate$Sex = as.character(butyrate$Sex)
butyrate$Dose = as.numeric(butyrate$Dose)

### Plotting Time
a <- colnames(butyrate[c(2,5:19)])
#Plot function
#Alter as needed for acetate and propionate plotting
plot_data = function(data, gene, title){
  ggplot(data, aes(y = gene, x = Dose)) +
    stat_summary(aes(group = interaction(Dose, Sex)), fun.data = mean_se, geom = "errorbar",
                 color = 'black', width = 0.5, position = position_dodge(width = 0.8)) +
    geom_bar(aes(y = gene, x =  Dose, fill = Sex), width = 0.8, 
             color = 'black',
             stat = "summary",  position = position_dodge2()) +
    scale_fill_manual(values = c('Male' = 'dodgerblue2', 'Female' = 'firebrick')) + 
    theme_bw() +
    ggtitle(title) +
    labs(x = 'Butyrate µM', y = 'ddCT Fold Change') +
    theme(  panel.grid.major = element_blank(), #Let's get rid of the grid inside the plot
            panel.grid.minor = element_blank(), #Let's get rid of the grid lines
            axis.line = element_line(colour = "black"),  #Change the color of the axes from default grey to black
            panel.border = element_blank(), #get rid of the border 
            plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm"), #Set the empty space outside of the plot
            legend.position="topright", #Let's you pick a position for the legend 
            legend.text = element_text(size=6),
            legend.title = element_text(size = 0),
            axis.title.x= element_text(size=16, color = 'black', vjust = 1),
            axis.title.y= element_text(size=16, color = 'black', vjust = 1),
            axis.text.x = element_text(size=16, color = 'black'),
            axis.text.y = element_text(size=16, color = 'black'),
            plot.title  = element_text(hjust = 0.5, size=20, face = 'bold', vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) 
}
# Aryl-Hydrocarbon Receptor Pathway
ahr   <- plot_data(butyrate, butyrate$Ahr_ddct, 'Ahr')
s100b <- plot_data(butyrate, butyrate$Sb100_ddct, 'S100B')
il_22 <- plot_data(butyrate, butyrate$IL_22_ddct, 'IL 22')
ifnar <- plot_data(butyrate, butyrate$IFNAR_ddct, 'IFNAR')
gfap  <- plot_data(butyrate, butyrate$GFAP_ddct, 'GFAP')
cyp1b1<- plot_data(butyrate, butyrate$Cyp1b1_ddct, 'Cyp1b1')
grid.arrange(ahr, s100b, il_22, ifnar, gfap, cyp1b1, nrow = 2)
#GABA/Glutamate Pathway
glal <- plot_data(butyrate, butyrate$Glal_ddct, 'Glal')
gad  <- plot_data(butyrate, butyrate$Gad_ddct, 'Gad')
glud <- plot_data(butyrate, butyrate$Glud_ddct, 'Glud')
grid.arrange(glal, glud, gad, nrow = 1)

#HDACi Pathway
bdnf <- plot_data(butyrate, butyrate$Bdnf_ddct, 'Bdnf')
nfg1 <- plot_data(butyrate, butyrate$Ngf1_ddct, 'Ng1')
gdnf <- plot_data(butyrate, butyrate$Gdnf_ddct, 'Gdnf')
pgc1 <- plot_data(butyrate, butyrate$Pgc1_ddct, 'Pgc1')
sp1  <- plot_data(butyrate, butyrate$Sp1_ddct, 'Sp1')
grid.arrange(bdnf, gdnf, nfg1, pgc1, sp1, nrow = 2)

######## Linear Regression
f = subset(butyrate, butyrate$Sex == 'Male')
m = subset(butyrate, butyrate$Sex == 'Female')
one   <- (lm(Ahr ~ Dose, f))
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

pgcm  <- (lm(`Pgc1-a` ~ Dose, m))
pgcf  <- (lm(`Pgc1-a` ~ Dose, f))
bdnfm  <- (lm(Bdnf ~ Dose, m))
bdnf  <- (lm(Bdnf ~ Dose, f))
bdnf_plot <- plot(y = f$Bdnf, x = f$Dose,
                  pch = 1, cex = 1.3, col = "black", cex.axis = 0.7, cex.lab = 0.7, main = "Dose-Response for Bdnf", xlab = "Butyrate (µM)", ylab = "Fold Change") +
  abline(bdnf, col = 'black')
d <- plot(y = m$`Pgc1-a`, x = m$Dose,
          pch = 1, cex = 1.3, col = "black", main = "Dose-Response for Pgc1", xlab = "Butyrate (µM)", ylab = "Fold Change",) +
  abline(pgcf, col = 'black')
### Load for heatmaps
pacman::p_load(car, pheatmap, gplots, RColorBrewer)

myCol <- c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858')

heatmap = read.csv('C:/Users/Simon/Downloads/heatmap_msc.csv')
rownames(heatmap) = heatmap$Sex
heatmap = heatmap[,2:15]
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
