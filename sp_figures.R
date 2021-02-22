# Define publication theme ####
plot_theme <- theme(plot.tag = element_text(color="black", face="bold"),
                    #strip.background=element_rect(color="black", size=1),
                    strip.text = element_text(size=9),
                    panel.grid = element_blank(),
                    panel.border = element_rect(color="black", fill=NA, size=1),
                    axis.text = element_text(color= "black", size=8),
                    axis.title = element_text(color="black", size=9),
                    legend.title = element_text(color="black", size=9),
                    plot.title=element_text(color="black", size=10),
                    legend.text = element_text(color="black", size=8))

#### FIGURE 1: MAP (made using GIS) ####
#### FIGURE 2: ALPHA DIVERSITY BOXPLOTS + BETA DIVERSITY ORDINATION (separated by type and location only) ####
Fig2.data <- reshape2::melt(sample_data)

levels(Fig2.data$variable)[levels(Fig2.data$variable)=="Observed_Extrap"] <- "ASV richness"
levels(Fig2.data$variable)[levels(Fig2.data$variable)=="Shannon_Extrap"] <- "Shannon diversity"

Fig2.datab <- Fig2.data
Fig2.datab$Source <- rep("All samples")
Fig2.data <- rbind(Fig2.datab, Fig2.data)
Fig2.data$Source <- factor(Fig2.data$Source, levels=c("All samples",
                                                      "Sooke",
                                                      "Nanaimo",
                                                      "Cowichan"))
rm(Fig2.datab)

Fig2a <- ggplot(subset(Fig2.data, variable %in% c("ASV richness")),
                 aes(Source,value)) +
  geom_boxplot(aes(fill=Type)) + 
  scale_fill_manual(values=c("darkorange1","cornflowerblue","chartreuse3")) +
  facet_wrap(~variable) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
  geom_vline(xintercept=c(1.5), color="black") +
  geom_vline(xintercept=c(2.5,3.5), color="grey") +
  theme_bw() + plot_theme +
  labs(x=NULL, y="ASV richness", tag="a") +
  theme(legend.position="none",
        legend.title=element_blank(),
        axis.title.y=element_blank())

Fig2b <- ggplot(subset(Fig2.data, variable %in% c("Shannon diversity")),
                 aes(Source,value)) +
  geom_boxplot(aes(fill=Type)) + 
  scale_fill_manual(values=c("darkorange1","cornflowerblue","chartreuse3")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
  facet_wrap(~variable, scales="free") +
  geom_vline(xintercept=c(1.5), color="black") +
  geom_vline(xintercept=c(2.5,3.5), color="grey") +
  theme_bw() + plot_theme +
  labs(x=NULL, y="Shannon diversity", tag="b") +
  theme(legend.position="none",
        legend.title=element_blank(),
        axis.title.y=element_blank())

# Aitchison distance
Fig2c <- ggplot(data =  bdiv.ordinations.all[["points"]][[3]], aes(PC1, PC2)) + 
  geom_point(aes(color = Type, shape=Source), size=2) +
  geom_path(data= bdiv.ordinations.all[["ellipses"]][[3]], aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  geom_path(data= bdiv.ordinations.all[["subellipses"]][[3]], aes(x=PC1, y=PC2, group=Merge), color="dimgrey", size=0.5, linetype=2) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x="PC1 (20.0%)", y="PC2 (15.4%)", tag="c") +
  #annotate("text", Inf, -Inf, label=paste0("Stress: ", 
  #                                         round(ordinations[["models"]][[1]]$stress, 3)), hjust=1.1, vjust=-0.4, size=3.125) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

grid.arrange(Fig2a, Fig2b, Fig2c,
             layout_matrix=rbind(c(1,3),c(2,3)))

Fig2a <- ggplotGrob(Fig2a)
Fig2b <- ggplotGrob(Fig2b)
maxWidth = grid::unit.pmax(Fig2a$widths[2:5], Fig2b$widths[2:5])
Fig2a$widths[2:5] <- as.list(maxWidth)
Fig2b$widths[2:5] <- as.list(maxWidth)

grid.arrange(Fig2a, Fig2b, Fig2c,
             layout_matrix=rbind(c(1,3),c(2,3)))

Fig2_grob <- arrangeGrob(Fig2a, Fig2b, Fig2c,
                         layout_matrix=rbind(c(1,3),c(2,3)))

ggsave("Fig2.adiv.bdiv.jpg", Fig2_grob, width=6.25, height=4.5, units="in", dpi=300)

rm(Fig2.data, Fig2_grob, Fig2a, Fig2b, Fig2c, maxWidth)

#### FIGURE 3: RELATIVE ABUNDANCE BAR CHARTS ####
y1 <- tax_glom(pseq.clr, taxrank="Class", NArm = FALSE)

y5 <- as.data.frame(cbind(tax_table(y1)))
y5$Phylum <- as.character(y5$Phylum)
y5$Phylum[y5$Phylum %in% c("Acidobacteria",
                           "Armatimonadetes",
                           "Chlamydiae",
                           "Chloroflexi",
                           "Deferribacteres",
                           "Deinococcus-Thermus",
                           "Euryarchaeota",
                           "Firmicutes",
                           "Fusobacteria",
                           "Gemmatimonadetes",
                           "Ignavibacteriae",
                           "Nitrospirae",
                           "Spirochaetes",
                           "Tenericutes",
                           "Thaumarchaeota")] <- "Other"
y5$Class <- as.character(y5$Class)
y5$Class[y5$Phylum=="Other"] <- "Taxa <1% Abundance"

y5$Class[y5$Phylum=="Verrucomicrobia"] <- "Verrucomicrobia"

y5$Class[y5$Phylum=="Planctomycetes"] <- "Planctomycetes"

y5$Class[y5$Class %in% c("Oligoflexia", "Epsilonproteobacteria",
                         "Deltaproteobacteria","Gammaproteobacteria")] <- "Other Proteobacteria"

y5$Phylum <- as.factor(y5$Phylum)
y5$Class <- as.factor(y5$Class)
seq <- rownames(y5)
y5 <- tax_table(y5)
rownames(y5) <- seq
colnames(y5) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax_table(y1) <- y5

y2 <- tax_glom(y1, taxrank="Class")
y2 <- prune_taxa(taxa_sums(y2) > 0, y2)
y2 <- transform_sample_counts(y2, function(x) 100*x/sum(x))

y4 <- psmelt(y2)
y4 <- y4 %>% dplyr::group_by(Type, Source, Phylum, Class) %>% dplyr::summarise(Abundance=mean(Abundance))


desired_order <- c("Other","Verrucomicrobia","Planctomycetes",
                   "Cyanobacteria","Actinobacteria","Bacteroidetes",
                   "Proteobacteria")
y4$Phylum <- factor(as.character(y4$Phylum), levels=desired_order)
y4 <- y4[order(y4$Phylum, y4$Class),]
colorCount <- length(unique(y4$Class))

desired_order <- c("Taxa <1% Abundance",
                   "Verrucomicrobia",
                   "Planctomycetes",
                   "Cyanobacteria",
                   "Actinobacteria",
                   "Bacteroidia",
                   "Cytophagia",
                   "Flavobacteriia",
                   "Sphingobacteriia",
                   "Alphaproteobacteria",
                   "Betaproteobacteria",
                   "Other Proteobacteria")
y4$Class <- factor(as.character(y4$Class), levels=desired_order)
y4 <- y4[order(y4$Class),]

myColors <- c("#999999", #Other#
              "#99CC99", #Verruco
              "#CC00FF", #Plancto
              "#006600", #Cyano#
              "#CC0033", #Actino
              "#FFCC00", ##Bacteroidia
              "#FF9900", ##Cytophagia
              "#FFCC66", ##Flavo
              "#CC9933", ##Sphingo
              "#3399FF", ##Alpha
              "#99CCFF", #Beta
              "#6699CC") #Other
names(myColors) <- levels(y4$Class)

desired_order <- c("sponge","water", "biofilm")
y4$Type <- factor(y4$Type, levels=desired_order)

desired_order <- c("Sooke","Nanaimo", "Cowichan")
y4$Source <- factor(y4$Source, levels=desired_order)

draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 3)
  
  grid::rectGrob(
    width = grid::unit(0.75, "npc"),
    height = grid::unit(0.75, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3

Fig3 <- ggplot(y4, aes(x=Type, y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  lemon::facet_rep_grid(~Source, repeat.tick.labels='left') +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  #scale_x_discrete(labels=function(Type) str_wrap(Location.Tag, width=3)) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() + plot_theme +
  theme(panel.spacing=unit(2, "lines"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color="black", size = 8, hjust=1, vjust=1, angle=45),
        axis.text.y=element_text(color="black", size = 8, vjust=1),
        strip.background=element_blank(), 
        legend.title=element_blank(),
        legend.position="right",
        legend.key.size=unit(1, "line")) +
  ylab("Mean relative abundance (%)") +
  xlab(NULL)

ggsave("Fig3.jpg", plot = Fig3, width=6.5, height=4, units="in")

rm(y1, y2, y3, y5, GeomBar, variable1, variable2, seq, myColors, desired_order, colorCount, Fig3)

#### FIGURE 4: ASV ABUNDANCE CORRELATIONS #### 
Fig4a <- ggplot(correlation.summary, aes(x=water, y=sponge, color=)) +
  geom_point(data=subset(correlation.summary, sponge !=0 & water!=0), color="darkgrey", size=0.75) +
  geom_point(data=subset(correlation.summary, sponge==0 | water==0), color="darkgrey", size=3, shape=3) +
  geom_point(data=subset(correlation.summary, sponge > 0 & water==0 & spwat.we.eBH < 0.05), color="darkorange1", size=3, shape=3) +
  geom_point(data=subset(correlation.summary, spwat.we.eBH < 0.05 & sponge > water & water > 0), aes(water,sponge), color="darkorange1", size=1) +
  geom_point(data=subset(correlation.summary, spwat.we.eBH < 0.05 & sponge < water), aes(water,sponge), color="cornflowerblue", size=1) +
  #geom_text_repel(data=subset(correlation.summary, spwat.we.eBH < 0.05 & sponge > water), aes(water,sponge,label=Genus), size=2.5) +
  geom_text_repel(data=subset(correlation.summary, spwat.we.eBH < 0.3 & water > 0.1 & sponge < water), aes(water,sponge,label=Genus), size=2.5) +
  scale_y_continuous(trans="log10", expand=c(0.035,0.035), limits=c(NA,NA), labels=comma) +
  scale_x_continuous(trans="log10", expand=c(0.035,0.035), limits=c(NA,NA), labels=comma) +
  theme_bw() + plot_theme + 
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title=element_blank(),
        legend.position="right") +
  labs(x="Water - Relative abundance (%)", 
       y="Sponge - Relative abundance (%)",
       tag = "a")

Fig4b <- ggplot(correlation.summary, aes(x=biofilm, y=sponge)) +
  geom_point(data=subset(correlation.summary, sponge !=0 & biofilm!=0), color="darkgrey", size=0.75) +
  geom_point(data=subset(correlation.summary, sponge==0 | biofilm==0), color="darkgrey", size=3, shape=3) +
  
  geom_point(data=subset(correlation.summary, sponge > 0 & biofilm==0 & spwat.we.eBH < 0.05), color="darkorange1", size=3, shape=3) +
  geom_point(data=subset(correlation.summary, sponge > biofilm & biofilm > 0 & spbio.we.eBH < 0.05), aes(biofilm,sponge), color="darkorange1", size=1) +
  
  geom_point(data=subset(correlation.summary, sponge < biofilm & biofilm > 0 & spbio.we.eBH < 0.05), aes(biofilm,sponge), color="chartreuse3", size=1) +
  #geom_point(data=subset(correlation.summary, (Genus %in% c("Haliea") | Class=="Cyanobacteria") & sponge < biofilm), aes(biofilm,sponge), color="chartreuse3", size=1) +
  #geom_text_repel(data=subset(correlation.summary, spbio.we.eBH < 0.05 & sponge > biofilm), aes(biofilm,sponge,label=Genus), size=2.5) +
  geom_text_repel(data=subset(correlation.summary, spbio.we.eBH < 0.05 & sponge < biofilm), aes(biofilm,sponge,label=Genus), size=2.5) +
  #geom_text_repel(data=subset(correlation.summary, sponge > 0.2 | biofilm > 0.2), aes(biofilm,sponge,label=Genus), size=2.5) +
  scale_y_continuous(trans="log10", expand=c(0.035,0.035), limits=c(NA,NA), labels=comma) +
  scale_x_continuous(trans="log10", expand=c(0.035,0.035), limits=c(NA,NA), labels=comma) +
  theme_bw() + plot_theme + 
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title=element_blank(),
        legend.position="right") +
  labs(x="Biofilm - Relative abundance (%)",
       y="Sponge - Relative abundance (%)",
       tag="b")

temp <- reshape2::melt(correlation.summary)
temp$yaxis <- seq(from=1, to=10, length.out=27995)
temp <- subset(temp, variable %in% c("sponge", "water", "biofilm"))
temp$variable <- factor(temp$variable, levels=c("sponge", "water", "biofilm"))

Fig4.legend <- cowplot::get_legend(ggplot(temp, aes(x=value, y=yaxis, color=variable)) +
                                     geom_point(size=3) +
                                     scale_color_manual(values=c("darkorange1", "cornflowerblue", "chartreuse3")) +
                                     theme_bw() + plot_theme +
                                     theme(legend.title=element_blank(),
                                           legend.position="bottom"))

grid.arrange(Fig4a, Fig4b, Fig4.legend,
             layout_matrix=rbind(c(1,2),c(3,3)), heights=c(0.9, 0.1))

gs <- arrangeGrob(Fig4a, Fig4b, Fig4.legend,
                  layout_matrix=rbind(c(1,2),c(3,3)), heights=c(0.9, 0.1))

ggsave(filename="enrichment.plot.signif.new.jpg", plot = gs,
       width = 6, height=3.5, units="in", dpi=500)

rm(temp, gs, Fig4a, Fig4b, Fig4.legend)

#### FIGURE 5: METAGENOME COGS HEAT MAP #####
names(sorting_list) <- colnames(sorting_data)
for(i in 1:length(sorting_list)){
  sorting_list[[i]]$Group <- rep(names(sorting_list)[i])
}
sorting_plot <- rbind(
  sorting_list$CRISPR,
  sorting_list$Transposase,
  sorting_list$FOGs,
  sorting_list$Cobalamin,
  sorting_list$Secretion,
  sorting_list$Restriction.Modification,
  sorting_list$Spore.Capsule
)

sorting_plot <- subset(sorting_plot, p.adj < 0.05)
sorting_plot <- subset(sorting_plot, abs(logFC) > 1)
sorting_plot$Group <- factor(sorting_plot$Group)
sorting_plot <- reshape2::melt(sorting_plot)

sorting_plot$Full.Name <- paste0(sorting_plot$COG, " ", sorting_plot$Name)
sorting_plot$Full.Name <- stringr::str_trunc(sorting_plot$Full.Name, 50)
sorting_plot$Full.Name <- as.factor(sorting_plot$Full.Name)
sorting_plot$value <- log(sorting_plot$value+1)

sorting_plot <- sorting_plot[order(sorting_plot$Group, sorting_plot$COG), ]
sorting_plot <- subset(sorting_plot, variable %in% c("SKE1", "SKE4", "SKE5",
                                                     "SKEW1", "SKEW4", "SKEW5"))

write.csv(sorting_plot, "trim.for.figure.csv")
sorting_plot <- read.csv("trimmed.for.figure.csv")


sorting_plot$Group <- factor(sorting_plot$Group)
sorting_plot$Full.Name <- factor(sorting_plot$Full.Name)

sorting_plot_left <- subset(sorting_plot, Group %in% c("Cobalamin", "CRISPR", "FOGs", "Restriction.Modification"))
sorting_plot_right <- subset(sorting_plot, !Group %in% c("Cobalamin", "CRISPR", "FOGs", "Restriction.Modification"))
sorting_plot_left$Group <- factor(sorting_plot_left$Group)
sorting_plot_right$Group <- factor(sorting_plot_right$Group)

# Create the left-hand side plot
rectangles_left <- data.frame(
  Group = levels(sorting_plot_left$Group),
  x1 = rep(0),
  x2 = rep(2.5),
  Number = rep(1),
  y1 = rep(0),
  y2 = rep(1)
)
for(j in c(1:nrow(rectangles_left))){
  rectangles_left[j,"Number"] <- (nrow(subset(sorting_plot_left, Group==levels(sorting_plot_left$Group)[j]))/6)
}
rectangles_left$Group <- factor(rectangles_left$Group)
rectangles_left <- rectangles_left %>% purrr::map_df(rev)
for(j in 1:nrow(rectangles_left)){
  rectangles_left[j,"y2"] <- sum(rectangles_left[c(1:j),"Number"])+0.5
  if(j > 1){
    rectangles_left[j,"y1"] <- rectangles_left[j-1, "y2"]}
}
rectangles_left$Number <- NULL

sorting_plot_left$Full.Name <- factor(sorting_plot_left$Full.Name)
sorting_plot_left$Full.Name <- forcats::fct_inorder(sorting_plot_left$Full.Name)
sorting_plot_left$variable <- factor(sorting_plot_left$variable)
sorting_plot_left$variable <- as.numeric(sorting_plot_left$variable)

# Create the right-hand side plot
rectangles_right <- data.frame(
  Group = levels(sorting_plot_right$Group),
  x1 = rep(6.5),
  x2 = rep(9),
  Number = rep(1),
  y1 = rep(0),
  y2 = rep(1)
)
for(j in c(1:nrow(rectangles_right))){
  rectangles_right[j,"Number"] <- (nrow(subset(sorting_plot_right, Group==levels(sorting_plot_right$Group)[j]))/6)
}
rectangles_right$Group <- factor(rectangles_right$Group)
rectangles_right <- rectangles_right %>% purrr::map_df(rev)
for(j in 1:nrow(rectangles_right)){
  rectangles_right[j,"y2"] <- sum(rectangles_right[c(1:j),"Number"])+0.5
  if(j > 1){
    rectangles_right[j,"y1"] <- rectangles_right[j-1, "y2"]}
}
rectangles_right$Number <- NULL

sorting_plot_right$Full.Name <- factor(sorting_plot_right$Full.Name)
sorting_plot_right$Full.Name <- forcats::fct_inorder(sorting_plot_right$Full.Name)
sorting_plot_right$variable <- factor(sorting_plot_right$variable)
sorting_plot_right$variable <- as.numeric(sorting_plot_right$variable)

# Generate the plots
left <- ggplot(sorting_plot_left) +
  geom_tile(aes(x=variable+2, y=Full.Name, fill=value)) +
  scale_fill_gradient2(guide = guide_colorbar(title.position="top",
                                              title.hjust=0.5,
                                              title="Log RPKM")) +
  scale_y_discrete(limits = rev(levels(sorting_plot_left$Full.Name))) +
  scale_x_continuous(expand=c(0,0), breaks=c(3,4,5,6,7,8),
                     labels=c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")) +
  geom_rect(data=rectangles_left,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="white", alpha=0.5) +
  geom_hline(yintercept=c(rectangles_left$y1)) +
  
  geom_text(data=subset(rectangles_left, Group =="Cobalamin"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Vitamin B12\nbiosynthesis", size=3, angle=90) +
  geom_text(data=subset(rectangles_left, Group =="CRISPR"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="CRISPR", size=3, angle=90) +
  geom_text(data=subset(rectangles_left, Group =="FOGs"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Structural\nmotifs", size=3, angle=90) +
  geom_text(data=subset(rectangles_left, Group =="Restriction.Modification"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Restriction/modification\nsystems", size=3, angle=90) +
  theme_bw() + plot_theme +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black", fill=NA)) +
  labs(x="Sample", y="COG")

right <- ggplot(sorting_plot_right) +
  geom_tile(aes(x=variable, y=Full.Name, fill=value)) +
  scale_fill_gradient2(guide = guide_colorbar(title.position="top",
                                              title.hjust=0.5,
                                              title="Log RPKM")) +
  scale_y_discrete(limits = rev(levels(sorting_plot_right$Full.Name)),
                   position="right") +
  scale_x_continuous(expand=c(0,0), breaks=c(1,2,3,4,5,6),
                     labels=c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")) +
  geom_rect(data=rectangles_right,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="white", alpha=0.5) +
  geom_hline(yintercept=c(rectangles_right$y1)) +
  geom_hline(yintercept=max(rectangles_right$y2)) +
  
  geom_text(data=subset(rectangles_right, Group =="Secretion"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Secretion\nsystems", size=3, angle=270) +
  geom_text(data=subset(rectangles_right, Group =="Spore.Capsule"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Sporulation", size=3, angle=270) +
  geom_text(data=subset(rectangles_right, Group =="Transposase"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Transposases", size=3, angle=270)+
  
  theme_bw() + plot_theme +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black", fill=NA)) +
  labs(x="Sample", y="COG")

grid.arrange(left, right, ncol=2)

master <- arrangeGrob(left, right, ncol=2)
ggsave("master.map.jpg", master, width=9.5, height=8, units="in", dpi=300)

rm(master, left, right, sorting_plot_left, rectangles_left, sorting_plot_right, rectangles_right, sorting_plot)

#### FIGURE 6: MAGs AND IMPORTANT FUNCTIONS ####
# (a) Bar plots for abundance
mag.sample_data$Sponge <- -1*rowMeans(mag.sample_data[,c("SKE1", "SKE5", "SKE5")])
mag.sample_data$Water <- rowMeans(mag.sample_data[,c("SKEW1", "SKEW5", "SKEW5")])

mag.sample_data$Spongesd <- rowSds(as.matrix(mag.sample_data[,c("SKE1", "SKE5", "SKE5")]))
mag.sample_data$Watersd <- rowSds(as.matrix(mag.sample_data[,c("SKEW1", "SKEW5", "SKEW5")]))

plota.data <- reshape2::melt(mag.sample_data)
plota.data.sd <- subset(plota.data, variable %in% c("Spongesd", "Watersd"))
plota.data <- subset(plota.data, variable %in% c("Sponge", "Water"))
plota.data$value[plota.data$value > 2] <- 2

colnames(plota.data.sd)[ncol(plota.data.sd)] <- "sd"
plota.data$sd <- plota.data.sd$sd
rm(plota.data.sd)

plota.data$MAG <- factor(plota.data$MAG, levels=c(" ", "  ", levels(plota.data$MAG)))
plota.data$Epithey <- paste0(plota.data$MAG, " ", plota.data$Species)
plota.data$smallname <- substr(plota.data$Epithey,1,nchar(plota.data$Epithey)-4)
plota.data$FinalPlotName <- paste0(plota.data$MAG, " ", plota.data$FinalPlotName)

relabund <- ggplot(plota.data, aes(x=MAG, y=value, fill=variable)) + geom_bar(stat="identity") +
  scale_fill_manual(values=c("darkorange1", "cornflowerblue"),
                    guide=guide_legend(title.position="top", title.hjust=0.5, title="Sample type", ncol=1)) + 
  labs(x="MAG", y="Relative\nabundance (%)") +
  coord_flip() +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, x=MAG), width=0.2) +
  scale_x_discrete(limits=rev(levels(plota.data$MAG)),
                   labels=c(as.character(rev(plota.data$FinalPlotName)), " ", "  ")) +
  #scale_y_log10() +
  scale_y_continuous(limits=c(-3, 3),
                     breaks=c(-2, 0, 2),
                     labels=c(2, 0, 2)) +
  theme_bw() + plot_theme +
  geom_hline(yintercept=0) +
  geom_rect(aes(ymin=0, ymax=Inf, xmin=25.5, xmax=28),
            color="white", fill="white", alpha=0.5) +
  theme(panel.border=element_blank(),
        legend.position="bottom",
        axis.text.y=element_text(color=c(rev(colors), "white", "white"))) +
  geom_segment(aes(x=-Inf, xend=25.5, y=-Inf, yend=-Inf), color="black") +
  geom_segment(aes(x=-Inf, xend=25.5, y=Inf, yend=Inf), color="black") +
  geom_segment(aes(x=-Inf, xend=-Inf, y=-Inf, yend=Inf), color="black") +
  geom_segment(aes(x=25.5, xend=25.5, y=-Inf, yend=Inf), color="black")

# (b) Heat map of genes
sorting_plot <- rbind(
  sorting_list$CRISPR,
  sorting_list$Transposase,
  sorting_list$FOGs,
  sorting_list$Cobalamin,
  sorting_list$Secretion,
  sorting_list$Restriction.Modification,
  sorting_list$Spore.Capsule
)

sorting_plot <- subset(sorting_plot, p.adj < 0.05)
sorting_plot <- subset(sorting_plot, abs(logFC) > 1)

mag.heat.data <- subset(mag.cog_counts,
                        rownames(mag.cog_counts) %in% sorting_plot$COG)

sorting_plot <- read.csv("trimmed.for.figure.csv")
mergedata <- sorting_plot[,c("COG", "Group", "Full.Name")]
mergedata <- dplyr::distinct(mergedata)

mag.heat.data$COG <- rownames(mag.heat.data)

mag.heat.data <- merge(mergedata,
                       mag.heat.data,
                       by="COG", all=FALSE)
mag.heat.data$Group <- factor(mag.heat.data$Group,
                              levels=c("Cobalamin", "CRISPR", "FOGs", "Restriction.Modification", "Secretion", "Spore.Capsule", "Transposase"))
mag.heat.data <- mag.heat.data[order(mag.heat.data$Group, mag.heat.data$COG), ]
mag.heat.data$Full.Name <- factor(mag.heat.data$Full.Name)
mag.heat.data$Full.Name <- forcats::fct_inorder(mag.heat.data$Full.Name)
mag.heat.data$COG <- factor(mag.heat.data$COG)
mag.heat.data$COG <- forcats::fct_inorder(mag.heat.data$COG)
substr(df$data,1,nchar(df$data)-3)

rectangles <- data.frame(
  Group = levels(mag.heat.data$Group),
  x1 = rep(0),
  x2 = rep(2.5),
  Number = rep(1),
  y1 = rep(25.5),
  y2 = rep(28)
)
for(j in c(1:nlevels(mag.heat.data$Group))){
  rectangles[j,"Number"] <- (nrow(subset(mag.heat.data, Group==levels(mag.heat.data$Group)[j])))
}
rectangles$Group <- factor(rectangles$Group)
for(j in 1:nrow(rectangles)){
  rectangles[j,"x2"] <- sum(rectangles[c(1:j),"Number"])+0.5
  if(j > 1){
    rectangles[j,"x1"] <- rectangles[j-1, "x2"]}
}
rectangles$Number <- NULL

mag.heat.data <- reshape2::melt(mag.heat.data)
mag.heat.data$variable <- factor(mag.heat.data$variable)
mag.heat.data$variable <- factor(mag.heat.data$variable,
                                 levels=rev(levels(mag.heat.data$variable)))
labels <- levels(mag.heat.data$variable)
mag.heat.data$variable <- as.numeric(mag.heat.data$variable)

mag.heat.data$value[mag.heat.data$value > 7] <- 7

#colors <- c(rep("#CC0033", 2), # Actinobacteria
#            rep("#999999", 1), # Other
#            rep("#6699CC", 10), # Proteobacteria
#            rep("#99CC99", 4), #Verruco
#            rep("#FF9900", 8)) #Bacteroidetes

colors <- c(rep("deeppink3", 2),
            rep("gray27", 1),
            rep("blue4", 10),
            rep("darkgreen", 4),
            rep("darkorange4", 8))

col2rgb("deeppink3")

magheat <- ggplot(mag.heat.data) + 
  geom_tile(aes(x=COG, y=variable, fill=value)) + 
  scale_fill_gradient(low="white", high="red",
                      guide=guide_colorbar(title.position="top",
                                            title.hjust=0.5,
                                            title="Number of copies")) +
  theme_bw() + plot_theme +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size=5),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position="bottom") +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(from=1, to=25, by=1)),
                     labels=labels) +
  geom_vline(xintercept=c(rectangles$x2), color="darkgrey") +
    geom_rect(data=rectangles,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="white", alpha=0.5) +

  
  geom_text(data=subset(rectangles, Group =="Cobalamin"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Vitamin B12\nbiosynthesis", size=2.75) +
  geom_text(data=subset(rectangles, Group =="CRISPR"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="CRISPR", size=2.75) +
  geom_text(data=subset(rectangles, Group =="FOGs"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Str.\nm.", size=2.75) +
  geom_text(data=subset(rectangles, Group =="Restriction.Modification"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="R/M\nsystems", size=2.75) +
  geom_text(data=subset(rectangles, Group =="Secretion"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Secretion\nsystems", size=2.75) +
  geom_text(data=subset(rectangles, Group =="Spore.Capsule"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Sporulation", size=2.75) +
  geom_text(data=subset(rectangles, Group =="Transposase"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Transposases", size=2.75)

grid.arrange(relabund, magheat, ncol=2)

relabund <- ggplotGrob(relabund)
magheat <- ggplotGrob(magheat)
relabund$heights <- magheat$heights

grid.arrange(relabund, magheat, ncol=2, widths=c(0.3, 0.8))

mag.summary.plot <- arrangeGrob(relabund, magheat, ncol=2, widths=c(0.35, 0.8))

ggsave("mag.summary.plot.jpg", mag.summary.plot, width=8.5, height=5, units="in", dpi=300)

rm(mag.summary.plot, magheat, relabund, mag.heat.data, mergedata, sorting_plot)

#### FIGURE S1: NEGATIVE CONTROLS ####
# (a) Read count histogram
data <- as.data.frame(sample_sums(pseq.original))
data <- subset(data, !rownames(data) %in% c("K1", "K2", "NCFILT1",
                                            "K0SKE2", "K0SKE3", "K0SKE5"))
colnames(data) <- "Counts"
fig.a <- ggplot(data, aes(x=Counts)) + geom_histogram(binwidth=1000) +
  theme_bw() + plot_theme +
  geom_vline(xintercept=mean(data$Counts), linetype="dashed") +
  labs(x="Number of Reads", y="Number of Samples", title="Read count distribution", tag="a") +
  scale_y_continuous(expand=c(0,0), limits=c(0,5))

# (b) Contaminant prevalences
data <- ctrl.contam.data$prevalence
data$contaminant <- factor(data$contaminant, levels=c("TRUE","FALSE"))

fig.b <- ggplot(data=subset(data, pa.pos > 0), aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  labs(x = "Prevalence (negative controls)", y = "Prevalence (true samples)", tag="b", 
       title = "Contaminant prevalence") +
  theme_bw() +
  guides(color=guide_legend(title="Contaminant")) +
  scale_color_manual(values=c("red", "darkgrey")) +
  plot_theme +
  theme(legend.position=c(0.95,0.95),
        legend.justification=c("right", "top"))

# (c) Contaminant abundances
data <- ctrl.contam.data$abundance
levels(data$Genus)[levels(data$Genus)=="Clostridium_sensu_stricto"] <- "Clostridium"
levels(data$Genus)[levels(data$Genus)=="Ruminococcus2"] <- "Ruminococcus"
levels(data$Genus)[levels(data$Genus)=="Ruminococcus2"] <- "Ruminococcus"
data$Genus <- factor(data$Genus, levels=c(levels(data$Genus), c("Uncl_Firmicutes", "Uncl_Clostridiales")))
data$Genus[data$Phylum=="Firmicutes" & is.na(data$Genus)=="TRUE" & is.na(data$Order)=="TRUE"] <- "Uncl_Firmicutes"
data$Genus[data$Phylum=="Firmicutes" & is.na(data$Genus)=="TRUE" & is.na(data$Order)=="FALSE"] <- "Uncl_Clostridiales"
data <- data[order(data$Phylum, data$Genus), ]
data$Genus <- forcats::fct_inorder(data$Genus)

fig.c <- ggplot(data, aes(x=ASV, y=means)) + geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  labs(x="ASV", y="Mean relative abundance (%)", title = "Contaminant abundances", tag="c") +
  scale_x_discrete(labels=as.character(data$Genus)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
  plot_theme +
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.position="bottom")

# (c2) Correlation between negative controls and experimental samples
data <- ctrl.nc.data$abundance.corr
data <- merge(data,
              as.data.frame(cbind(tax_table(pseq.filter))),
              by=0, all.y=FALSE)

fig.c2 <- ggplot(data, aes(x=experiment, y=`negative control`)) +
  geom_point(data=subset(data, experiment !=0 & `negative control` !=0), color="darkgrey", size=2) +
  geom_point(data=subset(data, experiment==0 | `negative control`==0), color="darkgrey", size=3, shape=3)+
  ggrepel::geom_text_repel(data=subset(data, experiment !=0 & `negative control` !=0),
                           aes(label=Genus)) +
  scale_y_continuous(trans="log10", expand=c(0.035,0.035), limits=c(NA,NA), labels=comma) +
  scale_x_continuous(trans="log10", expand=c(0.035,0.035), limits=c(NA,NA), labels=comma) +
  labs(x="ASV abundance in experiment (%)", y="ASV abundance in\nnegative controls (%)",
       title="Relative abundance: samples vs. controls", tag="c") +
  theme_bw() +
  plot_theme

# (d) ASV overlap among replicates
overlap <- data.frame(
  Sponge=c("SKE2", "SKE3", "SKE5"),
  Overlap=c(92.919, 95.223, 93.6997)
)

fig.d <- ggplot(overlap, aes(x=Sponge, y=Overlap)) +
  geom_bar(stat="identity") + theme_bw() + plot_theme +
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  labs(x="Sponge", y="ASV overlap (%)", title="Abundance-weighted ASV overlap", tag="d")

# (e) Alpha diversity among paired replicates
temp <- cbind(sample_data(ctrl.replicates), estimate_richness(ctrl.replicates))
temp$Sponge <- c("SKE2", "SKE3", "SKE5", "SKE2", "SKE3", "SKE5")
temp$Replicate <- c("1", "1", "1", "2", "2", "2")
temp$Replicate <- as.factor(as.character(temp$Replicate))
temp <- reshape2::melt(temp)
temp <- subset(temp, variable %in% c("Observed", "Shannon"))

levels(temp$variable)[levels(temp$variable)=="Observed"] <- "ASV richness"

fig.e <- ggplot(temp, aes(x=Sponge, y=value, fill=Replicate)) + 
  geom_bar(stat="identity", position="dodge") + theme_bw() + plot_theme +
  scale_fill_manual(values=c("black","grey")) +
  facet_wrap(~variable, scales="free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(panel.border=element_rect(color="black", fill=NA),
        legend.position="none") +
  labs(x="Sponge", y="Value", title="Alpha diversity in repeated samples", tag="e")

# (f) Ordination of replicates
temp.otu <- as.data.frame(otu_table(transform_sample_counts(ctrl.replicates, function(x) 100*x/sum(x))))
temp.dist <- vegdist(temp.otu, method="bray")
temp.pca <- capscale(temp.dist ~ 1)
temp.pca.scores <- cbind(sample_data(ctrl.replicates), vegan::scores(temp.pca, display="sites", choices=c(1:2)))

levels(temp.pca.scores$Type)[levels(temp.pca.scores$Type)=="sponge"] <- "original"
temp.pca.scores$Type <- factor(temp.pca.scores$Type, levels=c("original", "replicate"))
temp.pca.scores$SampleID <- c("SKE2_1", "SKE3_1", "SKE5_1", "SKE2_2", "SKE3_2", "SKE5_2")

fig.f <- ggplot(temp.pca.scores, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Type), size=0) +
  geom_text(aes(label=SampleID, color=Type)) +
  scale_color_manual(values=c("red", "deepskyblue3")) +
  theme_bw() + plot_theme +
  theme(panel.border=element_rect(color="black", fill=NA),
        legend.position="none") +
  scale_x_continuous(expand=c(0.1,0.1)) +
  labs(x=paste0("PCoA 1 (", round(100*temp.pca[["CA"]][["eig"]][1]/
                                    temp.pca$tot.chi, 1), "%)"),
       y=paste0("PCoA 2 (", round(100*temp.pca[["CA"]][["eig"]][2]/
                                    temp.pca$tot.chi, 1), "%)"),
       title = "Ordination of repeated samples", tag="f")


grid.arrange(fig.a, fig.b, fig.c2, fig.d, fig.e, fig.f, ncol=2, nrow=3)

fig.a <- ggplotGrob(fig.a)
fig.b <- ggplotGrob(fig.b)
fig.c2 <- ggplotGrob(fig.c2)
fig.d <- ggplotGrob(fig.d)
fig.e <- ggplotGrob(fig.e)
fig.f <- ggplotGrob(fig.f)

fig.a$widths <- fig.c2$widths
fig.f$widths <- fig.d$widths
fig.b$widths <- fig.d$widths

read.info.fig <- arrangeGrob(fig.a, fig.b, fig.c2, fig.d, fig.e, fig.f, ncol=2, nrow=3)

ggsave("read.info.fig.jpg", read.info.fig, width=6.75, height=8, units="in", dpi=300)

rm(data, fig.a, fig.b, fig.c, fig.c2, overlap, fig.d, temp, fig.e, temp.otu, temp.dist, temp.pca, temp.pca.scores, fig.f, read.info.fig)

#### FIGURE S2: RAREFACTION CURVES ####
rarefaction_curve_data$Type <- factor(rarefaction_curve_data$Type, levels=c("sponge", "water", "biofilm"))

a <- ggplot(subset(rarefaction_curve_data, Measure=="Observed"), aes(x=Depth, y=Alpha_diversity_mean, color=Type)) + 
  geom_line(aes(group=SampleID), size=1) + 
  facet_wrap(~Type) +
  theme_bw() + plot_theme +
  scale_color_manual(values=c("darkorange1", "cornflowerblue", "chartreuse3")) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none") +
  labs(x="Number of Reads", y="ASV richness")

b <- ggplot(subset(rarefaction_curve_data, Measure=="Shannon"), aes(x=Depth, y=Alpha_diversity_mean, color=Type)) + 
  geom_line(aes(group=SampleID), size=1) + 
  facet_wrap(~Type) +
  theme_bw() + plot_theme +
  scale_color_manual(values=c("darkorange1", "cornflowerblue", "chartreuse3")) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none") +
  labs(x="Number of Reads", y="Shannon diversity")

grid.arrange(a, b, nrow=2)

a <- ggplotGrob(a)
b <- ggplotGrob(b)
b$widths <- a$widths

grid.arrange(a, b, nrow=2, heights=c(0.45, 0.55))
rarefaction <- arrangeGrob(a, b, nrow=2, heights=c(0.45, 0.55))

ggsave("rarefaction.curves.jpg", rarefaction, width=5, height=4, units="in", dpi=300)

rm(a, b, rarefaction)

#### FIGURE S4: ALTERNATIVE DISTANCE MEASURES ####
bdiv.ordinations.all[["points"]][[1]]$Metric <- rep(c("Bray-Curtis"))
fig.a <- ggplot(data =  bdiv.ordinations.all[["points"]][[1]], aes(MDS1, MDS2)) + 
  geom_point(aes(color = Type, shape=Source), size=2) +
  geom_path(data= bdiv.ordinations.all[["ellipses"]][[1]], aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  #geom_path(data= bdiv.ordinations.all[["subellipses"]][[1]], aes(x=NMDS1, y=NMDS2, group=Merge), color="dimgrey", size=0.5, linetype=2) +
  facet_wrap(~Metric) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PCoA 1 (", round(100*bdiv.ordinations.all[["models"]][[1]][["CA"]][["eig"]][1]/
                  bdiv.ordinations.all[["models"]][[1]]$tot.chi, 1), "%)"),
       y=paste0("PCoA 2 (", round(100*bdiv.ordinations.all[["models"]][[1]][["CA"]][["eig"]][2]/
                  bdiv.ordinations.all[["models"]][[1]]$tot.chi, 1), "%)")
       ) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.title=element_text(size=8),
        plot.title=element_text(size=10),
        legend.position="none") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))
  
bdiv.ordinations.all[["points"]][[2]]$Metric <- rep(c("Jaccard"))
fig.b <- ggplot(data =  bdiv.ordinations.all[["points"]][[2]], aes(MDS1, MDS2)) + 
  geom_point(aes(color = Type, shape=Source), size=2) +
  geom_path(data= bdiv.ordinations.all[["ellipses"]][[2]], aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  #geom_path(data= bdiv.ordinations.all[["subellipses"]][[2]], aes(x=NMDS1, y=NMDS2, group=Merge), color="dimgrey", size=0.5, linetype=2) +
  facet_wrap(~Metric) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PCoA 1 (", round(100*bdiv.ordinations.all[["models"]][[2]][["CA"]][["eig"]][1]/
                                    bdiv.ordinations.all[["models"]][[2]]$tot.chi, 1), "%)"),
       y=paste0("PCoA 2 (", round(100*bdiv.ordinations.all[["models"]][[2]][["CA"]][["eig"]][2]/
                                    bdiv.ordinations.all[["models"]][[2]]$tot.chi, 1), "%)")
       ) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.title=element_text(size=8),
        plot.title=element_text(size=10),
        legend.position="none") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

bdiv.ordinations.all[["points"]][[4]]$Metric <- rep(c("Weighted UniFrac"))
fig.c <- ggplot(data =  bdiv.ordinations.all[["points"]][[4]], aes(MDS1, MDS2)) + 
  geom_point(aes(color = Type, shape=Source), size=2) +
  geom_path(data= bdiv.ordinations.all[["ellipses"]][[4]], aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  #geom_path(data= bdiv.ordinations.all[["subellipses"]][[4]], aes(x=NMDS1, y=NMDS2, group=Merge), color="dimgrey", size=0.5, linetype=2) +
  facet_wrap(~Metric) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PCoA 1 (", round(100*bdiv.ordinations.all[["models"]][[4]][["CA"]][["eig"]][1]/
                                    bdiv.ordinations.all[["models"]][[4]]$tot.chi, 1), "%)"),
       y=paste0("PCoA 2 (", round(100*bdiv.ordinations.all[["models"]][[4]][["CA"]][["eig"]][2]/
                                    bdiv.ordinations.all[["models"]][[4]]$tot.chi, 1), "%)")
       ) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.title=element_text(size=8),
        plot.title=element_text(size=10),
        legend.position="none") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

bdiv.ordinations.all[["points"]][[5]]$Metric <- rep(c("Unweighted UniFrac"))
fig.d <- ggplot(data =  bdiv.ordinations.all[["points"]][[5]], aes(MDS1, MDS2)) + 
  geom_point(aes(color = Type, shape=Source), size=2) +
  geom_path(data= bdiv.ordinations.all[["ellipses"]][[5]], aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  #geom_path(data= bdiv.ordinations.all[["subellipses"]][[5]], aes(x=NMDS1, y=NMDS2, group=Merge), color="dimgrey", size=0.5, linetype=2) +
  facet_wrap(~Metric) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PCoA 1 (", round(100*bdiv.ordinations.all[["models"]][[5]][["CA"]][["eig"]][1]/
                                    bdiv.ordinations.all[["models"]][[5]]$tot.chi, 1), "%)"),
       y=paste0("PCoA 2 (", round(100*bdiv.ordinations.all[["models"]][[5]][["CA"]][["eig"]][2]/
                                    bdiv.ordinations.all[["models"]][[5]]$tot.chi, 1), "%)")
       ) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.title=element_text(size=8),
        plot.title=element_text(size=10),
        legend.position="none") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

fig.legend <- cowplot::get_legend(
  ggplot(data =  bdiv.ordinations.all[["points"]][[5]], aes(MDS1, MDS2)) + 
  geom_point(aes(color = Type, shape=Source), size=2) +
  geom_path(data= bdiv.ordinations.all[["ellipses"]][[5]], aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        axis.title=element_text(size=8),
        plot.title=element_text(size=10),
        legend.title=element_blank(),
        legend.position="bottom") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3")))

grid.arrange(fig.a, fig.b, fig.c, fig.d, fig.legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
             heights=c(0.4,0.4,0.15))

fig.a <- ggplotGrob(fig.a)
fig.b <- ggplotGrob(fig.b)
fig.c <- ggplotGrob(fig.c)
fig.d <- ggplotGrob(fig.d)

fig.a$widths <- fig.c$widths
fig.b$widths <- fig.d$widths

fig.supp <- arrangeGrob(fig.a, fig.b, fig.c, fig.d, fig.legend,
                                     layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
                                     heights=c(0.4,0.4,0.15))

ggsave("other.distances.jpg", fig.supp, width=6.5, height=6.5, units="in", dpi=300)

rm(fig.a, fig.b, fig.c, fig.d, fig.legend, fig.supp)

#### FIGURE S5: AITCHISON ORDINATIONS FACETED BY RIVER ####
# (a) Facets by river ####
# Sooke
temp.clr <- subset_samples(pseq.clr, Source=="Sooke")
temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)
Sooke.dist <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")

Sooke.model <- rda(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))))
Sooke.points <- merge(vegan::scores(Sooke.model, display="sites", choices=c(1:3)),
                      subset(sample_data, Source=="Sooke")[,c(1:3)],
                      by=0, all=FALSE)

plot.new()
temp <- ordiellipse(Sooke.model, 
                    Sooke.points$Type, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)

Sooke.ellipse <- data.frame()

for(g in levels(Sooke.points$Type)){
  Sooke.ellipse <- rbind(Sooke.ellipse,
                         cbind(as.data.frame(with(Sooke.points[Sooke.points$Type==g,],
                                                  veganCovEllipse(temp[[g]]$cov,
                                                                  temp[[g]]$center,
                                                                  temp[[g]]$scale)))
                               ,Type=g))
}

rm(temp)

Sooke.model[["CA"]][["eig"]][[1]]/Sooke.model$tot.chi
Sooke.model[["CA"]][["eig"]][[2]]/Sooke.model$tot.chi

Sooke.points$Metric <- rep("Sooke")
Sooke <- ggplot(data =  Sooke.points, aes(PC1, PC2)) + 
  geom_point(aes(color = Type, shape = Type), size=2) +
  geom_path(data= Sooke.ellipse, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  facet_wrap(~Metric) +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*Sooke.model[["CA"]][["eig"]][[1]]/Sooke.model$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*Sooke.model[["CA"]][["eig"]][[2]]/Sooke.model$tot.chi, 1), "%)")) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

# Nanaimo
temp.clr <- subset_samples(pseq.clr, Source=="Nanaimo")
temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)
Nanaimo.dist <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")

Nanaimo.model <- rda(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))))
Nanaimo.points <- merge(vegan::scores(Nanaimo.model, display="sites", choices=c(1:3)),
                      subset(sample_data, Source=="Nanaimo")[,c(1:3)],
                      by=0, all=FALSE)

plot.new()
temp <- ordiellipse(Nanaimo.model, 
                    Nanaimo.points$Type, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)

Nanaimo.ellipse <- data.frame()

for(g in levels(Nanaimo.points$Type)){
  Nanaimo.ellipse <- rbind(Nanaimo.ellipse,
                         cbind(as.data.frame(with(Nanaimo.points[Nanaimo.points$Type==g,],
                                                  veganCovEllipse(temp[[g]]$cov,
                                                                  temp[[g]]$center,
                                                                  temp[[g]]$scale)))
                               ,Type=g))
}

rm(temp)

Nanaimo.model[["CA"]][["eig"]][[1]]/Nanaimo.model$tot.chi
Nanaimo.model[["CA"]][["eig"]][[2]]/Nanaimo.model$tot.chi

Nanaimo.points$Metric <- rep("Nanaimo")
Nanaimo <- ggplot(data =  Nanaimo.points, aes(PC1, PC2)) + 
  geom_point(aes(color = Type, shape = Type), size=2) +
  geom_path(data= Nanaimo.ellipse, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  facet_wrap(~Metric) +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*Nanaimo.model[["CA"]][["eig"]][[1]]/Nanaimo.model$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*Nanaimo.model[["CA"]][["eig"]][[2]]/Nanaimo.model$tot.chi, 1), "%)")) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

# Cowichan
temp.clr <- subset_samples(pseq.clr, Source=="Cowichan")
temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)
Cowichan.dist <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")

Cowichan.model <- rda(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))))
Cowichan.points <- merge(vegan::scores(Cowichan.model, display="sites", choices=c(1:3)),
                        subset(sample_data, Source=="Cowichan")[,c(1:3)],
                        by=0, all=FALSE)

plot.new()
temp <- ordiellipse(Cowichan.model, 
                    Cowichan.points$Type, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)

Cowichan.ellipse <- data.frame()

for(g in levels(Cowichan.points$Type)){
  Cowichan.ellipse <- rbind(Cowichan.ellipse,
                           cbind(as.data.frame(with(Cowichan.points[Cowichan.points$Type==g,],
                                                    veganCovEllipse(temp[[g]]$cov,
                                                                    temp[[g]]$center,
                                                                    temp[[g]]$scale)))
                                 ,Type=g))
}

rm(temp)

Cowichan.model[["CA"]][["eig"]][[1]]/Cowichan.model$tot.chi
Cowichan.model[["CA"]][["eig"]][[2]]/Cowichan.model$tot.chi

Cowichan.points$Metric <- rep("Cowichan")
Cowichan <- ggplot(data =  Cowichan.points, aes(PC1, PC2)) + 
  geom_point(aes(color = Type, shape = Type), size=2) +
  geom_path(data= Cowichan.ellipse, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  facet_wrap(~Metric) +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*Cowichan.model[["CA"]][["eig"]][[1]]/Cowichan.model$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*Cowichan.model[["CA"]][["eig"]][[2]]/Cowichan.model$tot.chi, 1), "%)")) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

# Assemble top facet
legend <- cowplot::get_legend(ggplot(data =  Cowichan.points, aes(PC1, PC2)) + 
                                geom_point(aes(color = Type, shape = Type), size=2) +
                                geom_path(data= Cowichan.ellipse, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
                                theme_bw() + plot_theme +
                                facet_wrap(~Metric) +
                                guides(color=guide_legend(ncol=3), shape=guide_legend(ncol=3)) +
                                labs(x="PC1 (41.2%)", y="PC2 (12.5%)", tag="b") +
                                theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
                                      panel.border = element_rect(color="black", fill=NA),
                                      axis.line = element_line(colour = "black", size=0.5),
                                      legend.title=element_blank(),
                                      legend.position="bottom",
                                      legend.box="horizontal") +
                                scale_color_manual(values=c("darkorange1","cornflowerblue","chartreuse3")))

grid.arrange(Sooke, Nanaimo, Cowichan, legend,
             layout_matrix=rbind(c(1,2,3),c(4,4,4)),
             heights=c(0.8, 0.2))

ords.by.source <- arrangeGrob(Sooke, Nanaimo, Cowichan, legend,
                            layout_matrix=rbind(c(1,2,3),c(4,4,4)),
                            heights=c(0.8, 0.2))

ggsave("ords.facet.by.source.jpg", ords.by.source,
       width=6.5, height=3.5, units="in", dpi=300)

# (b) Facets by source ####
# Sponge plot
bdiv.ordinations.sp[["points"]][[3]]$Metric <- rep("sponge")
sponge <- ggplot(data =  bdiv.ordinations.sp[["points"]][[3]], aes(PC1, PC2)) + 
  geom_point(aes(color = Source, shape = Source), size=2) +
  geom_path(data= bdiv.ordinations.sp[["ellipses"]][[3]], aes(x=PC1, y=PC2, colour=Source), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  facet_wrap(~Metric) +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*bdiv.ordinations.sp[["models"]][[3]][["CA"]][["eig"]][[1]]/
                                 bdiv.ordinations.sp[["models"]][[3]]$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*bdiv.ordinations.sp[["models"]][[3]][["CA"]][["eig"]][[2]]/
                                 bdiv.ordinations.sp[["models"]][[3]]$tot.chi, 1), "%)")) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkgoldenrod1","brown1","mediumorchid2"))

# Water plot
temp.clr <- subset_samples(pseq.clr, Type=="water")
temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)
water.dist <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")

water.model <- rda(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))))
water.points <- merge(vegan::scores(water.model, display="sites", choices=c(1:3)),
                      subset(sample_data, Type=="water")[,c(1:3)],
                      by=0, all=FALSE)

plot.new()
temp <- ordiellipse(water.model, 
                    water.points$Source, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)
  
water.ellipse <- data.frame()

for(g in levels(water.points$Source)){
  water.ellipse <- rbind(water.ellipse,
                         cbind(as.data.frame(with(water.points[water.points$Source==g,],
                                                  veganCovEllipse(temp[[g]]$cov,
                                                                  temp[[g]]$center,
                                                                  temp[[g]]$scale)))
                               ,Source=g))
}

rm(temp)

water.points$Metric <- rep("water")
water <- ggplot(data =  water.points, aes(PC1, PC2)) + 
  geom_point(aes(color = Source, shape = Source), size=2) +
  geom_path(data= water.ellipse, aes(x=PC1, y=PC2, colour=Source), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  facet_wrap(~Metric) +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*water.model[["CA"]][["eig"]][[1]]/
                                 water.model$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*water.model[["CA"]][["eig"]][[2]]/
                                 water.model$tot.chi, 1), "%)")) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkgoldenrod1","brown1","mediumorchid2"))

# Biofilm plot
temp.clr <- subset_samples(pseq.clr, Type=="biofilm")
temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)
biofilm.dist <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")

biofilm.model <- rda(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))))
biofilm.points <- merge(vegan::scores(biofilm.model, display="sites", choices=c(1:3)),
                      subset(sample_data, Type=="biofilm")[,c(1:3)],
                      by=0, all=FALSE)

plot.new()
temp <- ordiellipse(biofilm.model, 
                    biofilm.points$Source, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)

biofilm.ellipse <- data.frame()

for(g in levels(biofilm.points$Source)){
  biofilm.ellipse <- rbind(biofilm.ellipse,
                         cbind(as.data.frame(with(biofilm.points[biofilm.points$Source==g,],
                                                  veganCovEllipse(temp[[g]]$cov,
                                                                  temp[[g]]$center,
                                                                  temp[[g]]$scale)))
                               ,Source=g))
}

rm(temp)

biofilm.points$Metric <- rep("biofilm")
biofilm <- ggplot(data =  biofilm.points, aes(PC1, PC2)) + 
  geom_point(aes(color = Source, shape = Source), size=2) +
  geom_path(data= biofilm.ellipse, aes(x=PC1, y=PC2, colour=Source), size=0.5, linetype=1) +
  theme_bw() + plot_theme +
  facet_wrap(~Metric) +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*biofilm.model[["CA"]][["eig"]][[1]]/
                                 biofilm.model$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*water.model[["CA"]][["eig"]][[2]]/
                                 biofilm.model$tot.chi, 1), "%)")) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        panel.border = element_rect(color="black", fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(values=c("darkgoldenrod1","brown1","mediumorchid2"))

# Assemble plots
legend <- cowplot::get_legend(ggplot(data =  biofilm.points, aes(PC1, PC2)) + 
                                geom_point(aes(color = Source, shape = Source), size=2) +
                                geom_path(data= biofilm.ellipse, aes(x=PC1, y=PC2, colour=Source), size=0.5, linetype=1) +
                                theme_bw() + plot_theme +
                                facet_wrap(~Metric) +
                                guides(color=guide_legend(ncol=3), shape=guide_legend(ncol=3)) +
                                labs(x="PC1 (27.7%)", y="PC2 (22.3%)", tag="c") +
                                theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
                                      panel.border = element_rect(color="black", fill=NA),
                                      axis.line = element_line(colour = "black", size=0.5),
                                      legend.title=element_blank(),
                                      legend.position="bottom",
                                      legend.box="horizontal") +
                                scale_color_manual(values=c("darkgoldenrod1","brown1","mediumorchid2")))

grid.arrange(sponge, water, biofilm, legend,
             layout_matrix=rbind(c(1,2,3),c(4,4,4)),
             heights=c(0.8, 0.2))

ords.by.type <- arrangeGrob(sponge, water, biofilm, legend,
                            layout_matrix=rbind(c(1,2,3),c(4,4,4)),
                            heights=c(0.8, 0.2))

ggsave("ords.facet.by.type.jpg", ords.by.type,
       width=6.5, height=3.5, units="in", dpi=300)

grid.arrange(ords.by.source, ords.by.type, nrow=2)
ords.Aitchison.facets <- arrangeGrob(ords.by.source, ords.by.type, nrow=2)

ggsave("ords.Aitchison.facets.jpg", ords.Aitchison.facets,
       width=6.5, height=5.5, units="in", dpi=300)

rm(sponge, water, biofilm, temp.clr, legend, ords.by.type, 
   water.dist, water.model, water.points, water.ellipse,
   biofilm.dist, biofilm.model, biofilm.points, biofilm.ellipse)


rm(Sooke, Nanaimo, Cowichan, temp.clr, legend, ords.by.source, 
   Sooke.dist, Sooke.model, Sooke.points, Sooke.ellipse,
   Nanaimo.dist, Nanaimo.model, Nanaimo.points, Nanaimo.ellipse,
   Cowichan.dist, Cowichan.model, Cowichan.points, Cowichan.ellipse)

#### FIGURE S6: RANDOM FOREST MODELS ####
rf.graph.data.type <- rf.results[["Gini_scores"]][["type"]][order(-rf.results[["Gini_scores"]][["type"]]$MeanDecreaseGini),]
rf.graph.data.source <- rf.results[["Gini_scores"]][["sponges_by_location"]][order(-rf.results[["Gini_scores"]][["sponges_by_location"]]$MeanDecreaseGini),]
rf.graph.data.all <- rf.results[["Gini_scores"]][["location_and_type"]][order(-rf.results[["Gini_scores"]][["location_and_type"]]$MeanDecreaseGini),]

# Insert ASVID as a row name.
rf.graph.data.type$ASVID <- rownames(rf.graph.data.type)
rf.graph.data.source$ASVID <- rownames(rf.graph.data.source)
rf.graph.data.all$ASVID <- rownames(rf.graph.data.all)

rf.graph.data.type <- rf.graph.data.type[c(1:20), ]
rf.graph.data.source <- rf.graph.data.source[c(1:20), ]
rf.graph.data.all <- rf.graph.data.all[c(1:20), ]

# Order plots by ASVID.
rf.graph.data.type$ASVID <- forcats::fct_inorder(rf.graph.data.type$ASVID)
rf.graph.data.source$ASVID <- forcats::fct_inorder(rf.graph.data.source$ASVID)
rf.graph.data.all$ASVID <- forcats::fct_inorder(rf.graph.data.all$ASVID)

# Determine which grouping (segment: small/large, individual: identity, site: intestinal site) each taxon is most abundant in.
rf.graph.data.type$Max <- as.factor(colnames(rf.graph.data.type[,10:12])[max.col(rf.graph.data.type[,10:12], ties.method="first")])
rf.graph.data.source$Max <- as.factor(colnames(rf.graph.data.source[,10:12])[max.col(rf.graph.data.source[,10:12], ties.method="first")])
rf.graph.data.all$Max <- as.factor(colnames(rf.graph.data.all[,10:18])[max.col(rf.graph.data.all[,10:18], ties.method="random")])

rf.graph.data.type$Max <- factor(rf.graph.data.type$Max, levels=c("sponge", "water", "biofilm"))
rf.graph.data.source$Max <- factor(rf.graph.data.source$Max, levels=c("Sooke_sponge","Nanaimo_sponge", "Cowichan_sponge"))
levels(rf.graph.data.source$Max)[levels(rf.graph.data.source$Max)=="Sooke_sponge"] <- "Sooke"
levels(rf.graph.data.source$Max)[levels(rf.graph.data.source$Max)=="Nanaimo_sponge"] <- "Nanaimo"
levels(rf.graph.data.source$Max)[levels(rf.graph.data.source$Max)=="Cowichan_sponge"] <- "Cowichan"

# Plot for segment.
rf.plot.type <- ggplot(rf.graph.data.type, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Abundance)) + 
  scale_y_discrete(limits=rev(levels(rf.graph.data.type$ASVID)), labels=rev(rf.graph.data.type$Genus)) +
  scale_color_manual(name="Sample type", values=c("darkorange1","cornflowerblue","chartreuse3")) +
  scale_size_continuous(name="Abundance (%)") +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n97.30%", size=3) +
  theme_bw() + plot_theme +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=1)) +
  guides(size=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=0)) +
  labs(tag="a", y=NULL, x="Mean decrease Gini coefficient", title="All samples,\ndiscriminated by type")

rf.plot.source <- ggplot(rf.graph.data.source, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Abundance)) + 
  scale_y_discrete(limits=rev(levels(rf.graph.data.source$ASVID)), labels=rev(rf.graph.data.source$Genus)) +
  scale_color_manual(name="River", values=c("darkgoldenrod1","brown1","mediumorchid2")) +
  scale_size_continuous(name="Abundance (%)") +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n100%", size=3) +
  theme_bw() + plot_theme +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=1)) +
  guides(size=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=0)) +
  labs(tag="b", y=NULL, x="Mean decrease Gini coefficient", title="Sponges,\ndiscriminated by river")

grid.arrange(rf.plot.type, rf.plot.source, ncol=2)

rf.plot.type <- ggplotGrob(rf.plot.type)
rf.plot.source <- ggplotGrob(rf.plot.source)
rf.plot.source$heights <- rf.plot.type$heights

rf.models <- arrangeGrob(rf.plot.type, rf.plot.source, ncol=2)
ggsave("rf.model.results.jpg", rf.models, width=7, height=6, units="in", dpi=300)

rm(rf.models, rf.plot.type, rf.plot.source, rf.graph.data.all, rf.graph.data.source, rf.graph.data.type)

#### FIGURE S7: VENN DIAGRAMS + HISTOGRAMS SHOWING ASV OVERLAP ####
# Triple Venn - sponge/water/biofilm
plot.new()
Venn1 <- grid::grobTree(VennDiagram::draw.triple.venn(1959,  # total number of ASVs in sponges
                                   2073,  # total number of ASVs in water
                                   2042,  # total number of ASVs in biofilm
                                   1669,  # total sponge + water overlap
                                   1704,  # total water + biofilm overlap
                                   1563,  # total sponge + biofilm overlap
                                   1407, # total all three overlap
                                   #category=c("sponge","water","biofilm"),
                                   fill=c("darkorange1","cornflowerblue","chartreuse3"),
                                   cex=0.75,
                                   cat.cex=1,
                                   cat.dist=0.07))

# Triple Venn - SKE/NMO/COW
plot.new()
Venn2 <- grid::grobTree(VennDiagram::draw.triple.venn(1246,  # SKE (all ASVs)
                                   1143,  # NMO (all ASVs)
                                   1238,  # COW (all ASVs)
                                  745,  # SKE+NMO
                                   698,  # NMO+COW
                                   757,  # SKE+COW
                                   532, # all
                                   #category=c("Sooke","Nanaimo","Cowichan"),
                                   fill=c("darkgoldenrod1","brown1","mediumorchid2"),
                                   cex=0.75,
                                   cat.cex=1,
                                   cat.dist=0.07))

VennPlot <- arrangeGrob(Venn1, Venn2, ncol=2)
ggsave("spongeVennPlot.jpg", plot=VennPlot, width = 4, height=2, units="in", dpi=300)

theme.histogram <- theme(legend.position="right",
                         legend.title=element_blank(),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black"),
                         axis.text.x = element_text(color="black", size = 8, vjust=1),
                         axis.text.y = element_text(color="black", size = 8, vjust=1),
                         axis.title=element_text(color="black", size=8),
                         plot.margin=unit(c(5,10,5,5), "points"),
                         strip.background=element_blank(),
                         strip.text=element_blank())

plot.sponge <- ggplot() +
    geom_point(data=subset(Venn.histograms[["sponge"]], prevalence < 4 | mean < 0.01),
               aes(x=100*percent, y=mean), size=1) +
  geom_point(data=subset(Venn.histograms[["sponge"]], prevalence >= 4 & mean >= 0.01),
             aes(x=100*percent, y=mean), size=2, shape=18, color="red") +
    scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
    scale_y_continuous(limits=c(0,NA)) +
    theme_bw() +
    theme.histogram +
    labs(x="prevalence (%)", y="rel. abundance (%)")

plot.water <-  ggplot() +
  geom_point(data=subset(Venn.histograms[["water"]], prevalence < 4 | mean < 0.01),
             aes(x=100*percent, y=mean), size=1) +
  geom_point(data=subset(Venn.histograms[["water"]], prevalence >= 4 & mean >= 0.01),
             aes(x=100*percent, y=mean), size=2, shape=18, color="red") +
    scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
    scale_y_continuous(limits=c(0,NA)) +
    theme_bw() +
    theme.histogram +
    labs(x="prevalence (%)", y="rel. abundance (%)")

plot.biofilm <-  ggplot() +
  geom_point(data=subset(Venn.histograms[["biofilm"]], prevalence < 4 | mean < 0.01),
             aes(x=100*percent, y=mean), size=1) +
  geom_point(data=subset(Venn.histograms[["biofilm"]], prevalence >= 4 & mean >= 0.01),
             aes(x=100*percent, y=mean), size=2, shape=18, color="red") +
    scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
    scale_y_continuous(limits=c(0,NA)) +
    theme_bw() +
    theme.histogram +
    labs(x="prevalence (%)", y="rel. abundance (%)")

plot.SKE <- ggplot() +
  geom_point(data=subset(Venn.histograms[["SKE_sponge"]], prevalence < 4 | mean < 0.01),
             aes(x=100*percent, y=mean), size=1) +
  geom_point(data=subset(Venn.histograms[["SKE_sponge"]], prevalence >= 4 & mean >= 0.01),
             aes(x=100*percent, y=mean), size=2, shape=18, color="red") +
  scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
  scale_y_continuous(limits=c(0,NA)) +
  theme_bw() +
  theme.histogram +
  labs(x="prevalence (%)", y="rel. abundance (%)")

plot.NMO <- ggplot() +
  geom_point(data=subset(Venn.histograms[["NMO_sponge"]], prevalence < 4 | mean < 0.01),
             aes(x=100*percent, y=mean), size=1) +
  geom_point(data=subset(Venn.histograms[["NMO_sponge"]], prevalence >= 4 & mean >= 0.01),
             aes(x=100*percent, y=mean), size=2, shape=18, color="red") +
  scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
  scale_y_continuous(limits=c(0,NA)) +
  theme_bw() +
  theme.histogram +
  labs(x="prevalence (%)", y="rel. abundance (%)")

plot.COW <-  ggplot() +
  geom_point(data=subset(Venn.histograms[["COW_sponge"]], prevalence < 4 | mean < 0.01),
             aes(x=100*percent, y=mean), size=1) +
  geom_point(data=subset(Venn.histograms[["COW_sponge"]], prevalence >= 4 & mean >= 0.01),
             aes(x=100*percent, y=mean), size=2, shape=18, color="red") +
  scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
  scale_y_continuous(limits=c(0,NA)) +
  theme_bw() +
  theme.histogram +
  labs(x="prevalence (%)", y="rel. abundance (%)")

grid.arrange(plot.sponge, plot.water, plot.biofilm,
             plot.SKE, plot.NMO, plot.COW,
             ncol=2, nrow=3)

unique.prev.abund <- arrangeGrob(plot.sponge, plot.water, plot.biofilm,
                                 plot.SKE, plot.NMO, plot.COW,
                                 ncol=2, nrow=3)

ggsave(filename="unique.prev.abund.jpg", plot = unique.prev.abund, width = 4, height=6, units="in", dpi=300)

legend.data <- data.frame(
  Type=c(rep("other ASV", 5), rep("ASV unique to\nsample type/location", 5)),
  mean=rnorm(10, mean=5, sd=2),
  percent=rnorm(10, mean=5, sd=2)
)

legend <- cowplot::get_legend(ggplot(legend.data, aes(x=percent, y=mean)) +
                                geom_point(aes(color=Type, shape=Type)) +
                                scale_shape_manual(values=c(18, 16)) +
                                scale_color_manual(values=c("red", "black")) +
                                scale_x_continuous(limits=c(NA,100), breaks = c(0,20,40,60,80,100)) + 
                                scale_y_continuous(limits=c(0,NA)) +
                                theme_bw() +
                                theme.histogram +
                                labs(x="prevalence (%)", y="rel. abundance (%)"))

ggsave(filename="histogram.legend.jpg", plot = legend, width = 4, height=6, units="in", dpi=300)

rm(plot.sponge, plot.water, plot.biofilm,
   plot.SKE, plot.NMO, plot.COW, unique.prev.abund,
   Venn1, Venn2, VennPlot)

#### FIGURE S8: DIFFERENTIALLY ABUNDANT TAXA AMONG TYPES & SITES ####
temp.physeq <- tax_glom(pseq.clr.sub, taxrank="Genus")
temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
temp.physeq <- subset_taxa(temp.physeq, Genus %in% c("Sediminibacterium",
                                                     "Uncl_Rhodospirillales",
                                                     "Comamonas",
                                                     # "Fluviicola",
                                                     "Pseudarcicella",
                                                     "Flavobacterium",
                                                     "Polynucleobacter",
                                                     # "Rhodoluna",
                                                     "Acidovorax",
                                                     "Turneriella"))
temp.data <- psmelt(temp.physeq)
temp.data$PlotName <- temp.data$Genus

temp.physeq <- tax_glom(pseq.clr.sub, taxrank="Phylum")
temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
temp.physeq <- subset_taxa(temp.physeq, Phylum %in% c("Actinobacteria",
                                                      "Planctomycetes",
                                                      "Parcubacteria",
                                                      "Verrucomicrobia"))
temp.data2 <- psmelt(temp.physeq)
temp.data2$Class <- rep(NA)
temp.data2$Order <- rep(NA)
temp.data2$Family <- rep(NA)
temp.data2$Genus <- rep(NA)
temp.data2$PlotName <- temp.data2$Phylum

temp.physeq <- tax_glom(pseq.clr.sub, taxrank="Class")
temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
temp.physeq <- subset_taxa(temp.physeq, Class %in% c("Alphaproteobacteria"))
temp.data3 <- psmelt(temp.physeq)
temp.data3$Order <- rep(NA)
temp.data3$Family <- rep(NA)
temp.data3$Genus <- rep(NA)
temp.data3$PlotName <- temp.data3$Class

temp.data <- rbind(temp.data, temp.data2)
rm(temp.data2, temp.data3, temp.physeq)

desired_order <- c("Actinobacteria", "Parcubacteria", "Planctomycetes", "Verrucomicrobia",
                   "Acidovorax", "Comamonas", "Flavobacterium", "Polynucleobacter", "Pseudarcicella",  
                   "Uncl_Rhodospirillales", "Sediminibacterium", "Turneriella")
temp.data$PlotName <- factor(temp.data$PlotName, levels=desired_order)
#levels(temp.data$PlotName)[levels(temp.data$PlotName)=="Cyanobacteria/Chloroplast"] <- "Cyanobacteria"
temp.data$Type <- factor(temp.data$Type, levels=c("sponge", "water", "biofilm"))
temp.data$Source <- factor(temp.data$Source, levels=c("Sooke", "Nanaimo", "Cowichan"))

temp.data <- temp.data %>% dplyr::group_by(PlotName, Type, Source) %>%
  dplyr::summarise(mean=mean(Abundance),
                   sd=sd(Abundance))

taxa.by.rivers <- ggplot(temp.data, aes(x=Source, y=mean, fill=Type)) + 
  geom_col(position="dodge") + 
  #geom_boxplot() +
  facet_wrap(~PlotName, scales="free_y") +
  geom_errorbar(aes(ymin=pmax(mean-sd, 0), ymax=mean+sd, x=Source), position=position_dodge(0.9), width=0.2) +
  scale_fill_manual(values=c("darkorange1","cornflowerblue","chartreuse3")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_bw() + plot_theme +
  labs(x="River", y="Mean relative abundance (%)") +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, vjust=1, size=9),
        strip.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        legend.title=element_blank(),
        strip.text=element_text(color="black", size=9, face="italic"))

ggsave(filename="taxa.by.rivers.jpg", plot = taxa.by.rivers, width = 6.5, height=6, units="in", dpi=300)

#### FIGURE S9: METAGENOME (1) TAXONOMy AND (2) ORDINATIONS #####
# (a) Relative abundance bar charts ####
# Amplicons
y1 <- subset_samples(pseq.clr, sample_names(pseq.clr) %in% c("SKE1", "SKE4",
                                                             "SKE5", "SKEW1",
                                                             "SKEW4", "SKEW5"))
y1 <- tax_glom(y1, taxrank="Class", NArm = FALSE)

y5 <- as.data.frame(cbind(tax_table(y1)))
y5$Phylum <- as.character(y5$Phylum)
y5$Phylum[y5$Phylum %in% c("Acidobacteria",
                           "Armatimonadetes",
                           "Chlamydiae",
                           "Chloroflexi",
                           "Deferribacteres",
                           "Deinococcus-Thermus",
                           "Euryarchaeota",
                           "Firmicutes",
                           "Fusobacteria",
                           "Gemmatimonadetes",
                           "Ignavibacteriae",
                           "Nitrospirae",
                           "Spirochaetes",
                           "Tenericutes",
                           "Thaumarchaeota")] <- "Other"
y5$Class <- as.character(y5$Class)
y5$Class[y5$Phylum=="Other"] <- "Taxa <1% Abundance"

y5$Class[y5$Phylum=="Verrucomicrobia"] <- "Verrucomicrobia"
y5$Class[y5$Phylum=="Planctomycetes"] <- "Planctomycetes"

y5$Class[y5$Class %in% c("Oligoflexia", "Epsilonproteobacteria",
                         "Deltaproteobacteria","Gammaproteobacteria")] <- "Other Proteobacteria"

y5$Phylum <- as.factor(y5$Phylum)
y5$Class <- as.factor(y5$Class)
seq <- rownames(y5)
y5 <- tax_table(y5)
rownames(y5) <- seq
colnames(y5) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax_table(y1) <- y5

y2 <- tax_glom(y1, taxrank="Class")
y2 <- prune_taxa(taxa_sums(y2) > 0, y2)
y2 <- transform_sample_counts(y2, function(x) 100*x/sum(x))

y4 <- psmelt(y2)
y4 <- y4 %>% dplyr::group_by(Type, Phylum, Class) %>% dplyr::summarise(Abundance=mean(Abundance))

desired_order <- c("Other","Verrucomicrobia","Planctomycetes",
                   "Cyanobacteria","Actinobacteria","Bacteroidetes",
                   "Proteobacteria")
y4$Phylum <- factor(as.character(y4$Phylum), levels=desired_order)
y4 <- y4[order(y4$Phylum, y4$Class),]
colorCount <- length(unique(y4$Class))

desired_order <- c("Taxa <1% Abundance",
                   "Verrucomicrobia",
                   "Planctomycetes",
                   "Cyanobacteria",
                   "Actinobacteria",
                   "Bacteroidia",
                   "Cytophagia",
                   "Flavobacteriia",
                   "Sphingobacteriia",
                   "Alphaproteobacteria",
                   "Betaproteobacteria",
                   "Other Proteobacteria")
y4$Class <- factor(as.character(y4$Class), levels=desired_order)
y4 <- y4[order(y4$Class),]

# Metagenomes
y1 <- tax_glom(metagenome.pseq, taxrank="Class", NArm = FALSE)

y5 <- as.data.frame(cbind(tax_table(y1)))
y5$Phylum <- as.character(y5$Phylum)
y5$Phylum[y5$Phylum %in% c("Acidobacteria",
                           "Armatimonadetes",
                           "Chlamydiae",
                           "Chloroflexi",
                           "Deferribacteres",
                           "Deinococcus-Thermus",
                           "Euryarchaeota",
                           "Firmicutes",
                           "Fusobacteria",
                           "Gemmatimonadetes",
                           "Ignavibacteriae",
                           "Nitrospirae",
                           "Spirochaetes",
                           "Tenericutes",
                           "Thaumarchaeota")] <- "Other"
y5$Class <- as.character(y5$Class)
y5$Class[y5$Phylum=="Other"] <- "Taxa <1% Abundance"

y5$Class[y5$Phylum=="Verrucomicrobia"] <- "Verrucomicrobia"
y5$Class[y5$Phylum=="Planctomycetes"] <- "Planctomycetes"
y5$Class[y5$Class %in% c("Oligoflexia", "Epsilonproteobacteria",
                         "Deltaproteobacteria","Gammaproteobacteria")] <- "Other Proteobacteria"
y5$Class[y5$Phylum=="Cyanobacteria"] <- "Cyanobacteria"
y5$Class[y5$Phylum=="Actinobacteria"] <- "Actinobacteria"

y5$Class[y5$Class %in% c("Oligoflexia", "Epsilonproteobacteria",
                         "Deltaproteobacteria","Gammaproteobacteria",
                         "Zetaproteobacteria", "Hydrogenophilalia",
                         "Acidithiobacillia")] <- "Other Proteobacteria"
y5$Class[y5$Phylum=="Protebacteria" & is.na(y5$Class)=="TRUE"] <- "Other Proteobacteria"

y5$Class[y5$Class %in% c("Aquificae",
                         "Chitinophagia",
                         "Saprospiria",
                         "Caldisericia",
                         "Calditrichae",
                         "Chlorobia",
                         "Chrysiogenetes",
                         "Coprothermobacteria",
                         "Dictyoglomia",
                         "Endomicrobia",
                         "Elusimicrobia",
                         "Fibrobacteria",
                         "Kiritimatiellae",
                         "Lentisphaeria",
                         "Synergistia",
                         "Thermodesulfobacteria",
                         "Thermotogae")] <- "Taxa <1% Abundance"
y5$Class <- as.character(y5$Class)
y5$Phylum[y5$Class=="Taxa <1% Abundance"] <- "Other"

y5$Phylum <- as.factor(y5$Phylum)
y5$Class <- as.factor(y5$Class)
seq <- rownames(y5)
y5 <- tax_table(y5)
rownames(y5) <- seq
colnames(y5) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table(y1) <- y5

y2 <- tax_glom(y1, taxrank="Class")
y2 <- prune_taxa(taxa_sums(y2) > 0, y2)
y2 <- transform_sample_counts(y2, function(x) 100*x/sum(x))

y7 <- psmelt(y2)
y7 <- y7 %>% dplyr::group_by(Type, Phylum, Class) %>% dplyr::summarise(Abundance=mean(Abundance))

desired_order <- c("Other","Verrucomicrobia","Planctomycetes",
                   "Cyanobacteria","Actinobacteria","Bacteroidetes",
                   "Proteobacteria")
y7$Phylum <- factor(as.character(y7$Phylum), levels=desired_order)
y7 <- y7[order(y7$Phylum, y7$Class),]
colorCount <- length(unique(y7$Class))

desired_order <- c("Taxa <1% Abundance",
                   "Verrucomicrobia",
                   "Planctomycetes",
                   "Cyanobacteria",
                   "Actinobacteria",
                   "Bacteroidia",
                   "Cytophagia",
                   "Flavobacteriia",
                   "Sphingobacteriia",
                   "Alphaproteobacteria",
                   "Betaproteobacteria",
                   "Other Proteobacteria")
y7$Class <- factor(as.character(y7$Class), levels=desired_order)
y7 <- y7[order(y7$Class),]

# Put the plots together
y4$Experiment <- rep("amplicons")
y7$Experiment <- rep("metagenomes")

yx <- rbind(y4, y7)

myColors <- c("#999999", #Other#
              "#99CC99", #Verruco
              "#CC00FF", #Plancto
              "#006600", #Cyano#
              "#CC0033", #Actino
              "#FFCC00", ##Bacteroidia
              "#FF9900", ##Cytophagia
              "#FFCC66", ##Flavo
              "#CC9933", ##Sphingo
              "#3399FF", ##Alpha
              "#99CCFF", #Beta
              "#6699CC") #Other
names(myColors) <- levels(yx$Class)

desired_order <- c("sponge","water")
yx$Type <- factor(yx$Type, levels=desired_order)

fig.a <- ggplot(yx, aes(x=Experiment, y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  lemon::facet_rep_grid(~Type, repeat.tick.labels='left') +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  #scale_x_discrete(labels=function(Type) str_wrap(Location.Tag, width=3)) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() + plot_theme +
  theme(panel.spacing=unit(2, "lines"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color="black", size = 8, hjust=1, vjust=1, angle=45),
        axis.text.y=element_text(color="black", size = 8, vjust=1),
        strip.background=element_blank(), 
        legend.title=element_blank(),
        legend.position="right",
        legend.key.size=unit(1, "line")) +
  ylab("Mean relative abundance (%)") +
  xlab(NULL) +
  labs(tag="a")

# (b) Taxonomy ordination ####
temp.pca <- as.data.frame(t(otu_table(metagenome.pseq)))
temp.dist <- vegdist(temp.pca, method="bray")
adonis(temp.dist ~ Type, metagenome.sample_data, permutations=1000)
permutest(betadisper(temp.dist, metagenome.sample_data$Type))
temp.pca <- capscale(temp.dist ~ 1)
temp.pca.scores <- merge(vegan::scores(temp.pca,
                                display="sites", choices=c(1:2)),
                         sample_data(metagenome.pseq),
                         by=0, all=TRUE)
plot.new()
temp <- ordiellipse(temp.pca, 
                    temp.pca.scores$Type, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)

df_ell <- data.frame()

for(g in levels(temp.pca.scores$Type)){
  df_ell <- rbind(df_ell,
                  cbind(as.data.frame(with(temp.pca.scores[temp.pca.scores$Type==g,],
                                           veganCovEllipse(temp[[g]]$cov,
                                                           temp[[g]]$center,
                                                           temp[[g]]$scale)))
                        ,Type=g))
}

temp.pca.scores$Metric <- rep("Taxonomy")
fig.b <- ggplot(data = temp.pca.scores, aes(MDS1, MDS2)) + 
  geom_point(aes(color = Type, shape=Type), size=2) +
  geom_path(data= df_ell, aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  facet_wrap(~Metric) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PCoA 1 (", round(100*temp.pca[["CA"]][["eig"]][1]/
                                    temp.pca$tot.chi, 1), "%)"),
       y=paste0("PCoA 2 (", round(100*temp.pca[["CA"]][["eig"]][2]/
                                    temp.pca$tot.chi, 1), "%)"),
       tag="b"
  ) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        plot.title=element_text(size=10),
        legend.position="none") +
  scale_color_manual(values=c("darkorange1","cornflowerblue"))

# (c) COG ordination ####
cog_comm_matrix <- as.data.frame(t(metagenome.cog_rpkm[,c("SKE1","SKE4","SKE5", "SKEW1", "SKEW4", "SKEW5")]))
cog_comm_matrix <- log(cog_comm_matrix+0.01)
cog_ord_PCA <- rda(cog_comm_matrix)
temp.dist <- vegdist(cog_comm_matrix, method="euclidean")
adonis(temp.dist ~ Type, metagenome.sample_data, permutations=1000)
permutest(betadisper(temp.dist, metagenome.sample_data$Type))

metagenome.sample_data <- cbind(metagenome.sample_data,
                                vegan::scores(cog_ord_PCA, display="sites", choices = c(1:3)))
colnames(metagenome.sample_data)[c(7:9)] <- c("COG_logPC1", "COG_logPC2", "COG_logPC3")

plot.new()
temp <- ordiellipse(cog_ord_PCA, 
                    metagenome.sample_data$Type, 
                    display="sites", 
                    kind="sd", 
                    conf=0.90, 
                    label=T)

df_ell <- data.frame()

for(g in levels(metagenome.sample_data$Type)){
  df_ell <- rbind(df_ell,
                  cbind(as.data.frame(with(metagenome.sample_data[metagenome.sample_data$Type==g,],
                                           veganCovEllipse(temp[[g]]$cov,
                                                           temp[[g]]$center,
                                                           temp[[g]]$scale)))
                        ,Type=g))
}

metagenome.sample_data$Metric <- rep("COGs")
fig.c <- ggplot(data = metagenome.sample_data, aes(COG_logPC1, COG_logPC2)) + 
  geom_point(aes(color = Type, shape=Type), size=2) +
  geom_path(data= df_ell, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  facet_wrap(~Metric) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
  labs(x=paste0("PC1 (", round(100*cog_ord_PCA[["CA"]][["eig"]][1]/
                                 cog_ord_PCA$tot.chi, 1), "%)"),
       y=paste0("PC2 (", round(100*cog_ord_PCA[["CA"]][["eig"]][2]/
                                 cog_ord_PCA$tot.chi, 1), "%)"),
       tag="c"
  ) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        plot.title=element_text(size=10),
        legend.position="none") +
  scale_color_manual(values=c("darkorange1","cornflowerblue"))

legend <- cowplot::get_legend(ggplot(data = metagenome.sample_data, aes(COG_logPC1, COG_logPC2)) + 
                                 geom_point(aes(color = Type, shape=Type), size=2) +
                                 geom_path(data= df_ell, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
                                 facet_wrap(~Metric) +
                                 theme_bw() + plot_theme +
                                 guides(color=guide_legend(ncol=1, title.position="top"),
                                        shape=guide_legend(ncol=1, title.position="top")) +
                                 labs(x=paste0("PC1 (", round(100*cog_ord_PCA[["CA"]][["eig"]][1]/
                                                                cog_ord_PCA$tot.chi, 1), "%)"),
                                      y=paste0("PC2 (", round(100*cog_ord_PCA[["CA"]][["eig"]][2]/
                                                                cog_ord_PCA$tot.chi, 1), "%)"),
                                      tag="c"
                                 ) +
                                 theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
                                       axis.line = element_line(colour = "black", size=0.5),
                                       axis.title=element_text(size=8),
                                       plot.title=element_text(size=10),
                                       legend.position="bottom",
                                       legend.title=element_text(hjust=0.5)) +
                                 scale_color_manual(values=c("darkorange1","cornflowerblue")))

grid.arrange(fig.a,
             fig.b, fig.c,
             legend,
             layout_matrix=rbind(c(NA,1,1,1,1,NA),
                                 c(2,2,2,3,3,3),
                                 c(4,4,4,4,4,4)),
             heights=c(0.5,0.4,0.1))

metagenome.summary <- arrangeGrob(fig.a,
                                  fig.b, fig.c,
                                  legend,
                                  layout_matrix=rbind(c(NA,1,1,1,1,NA),
                                                      c(2,2,2,3,3,3),
                                                      c(4,4,4,4,4,4)),
                                  heights=c(0.5,0.4,0.1))

ggsave("metagenome.summary.jpg", metagenome.summary,
       height=8, width=6.5, units="in", dpi=300)

#### FIGURE S10: METAGNOME COG CLASSES ####
plot1.data <- reshape2::melt(metagenome.class_rpkm)
plot1.data$value <- log(plot1.data$value)
plot1.data$Class <- factor(plot1.data$Class)

plot1 <- ggplot(plot1.data, aes(x=variable, y=Class, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high=muted("blue"),
                      guide=guide_colorbar(title.position="top", title.hjust=0.5, title="Log RPKM")) +
  scale_y_discrete(limits=rev(levels(plot1.data$Class))) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Sample", y="COG class", tag="a") +
  theme_bw() + plot_theme +
  theme(legend.position="bottom")

plot2.data <- metagenome.class_sums
plot2.data$All <- NULL
plot2.data$Sponge <- 100*plot2.data$Sponge / plot2.data$Total_Count
plot2.data$Water <- 100*plot2.data$Water / plot2.data$Total_Count
plot2.data$Gap <- 100 - plot2.data$Water - plot2.data$Sponge
plot2.data$Total_Count <- NULL

plot2.data <- reshape2::melt(plot2.data)
plot2.data$variable <- factor(plot2.data$variable, levels=c("Water", "Gap", "Sponge"))
plot2.data$Class <- factor(plot2.data$Class)

plot2 <- ggplot(plot2.data, aes(x=Class, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=rev(levels(plot2.data$Class))) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("cornflowerblue", "white", "darkorange1"),
                    guide=guide_legend(title.position="top",
                                       title.hjust=0.5,
                                       ncol=2,
                                       title="Sample type",
                                       reverse=TRUE)) +
  theme_bw() + plot_theme +
  coord_flip() +
  labs(x = "COG class", y = "Percent of Genes", tag="b") +
  theme(legend.position = "bottom")

grid.arrange(plot1, plot2, ncol=2)

plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)
plot1$heights <- plot2$heights

grid.arrange(plot1, plot2, ncol=2)

cog.class.summary.plot <- arrangeGrob(plot1, plot2, ncol=2)
ggsave("cog.class.summary.plot.jpg", cog.class.summary.plot, width=6.5, height=6, units="in", dpi=300)

rm(plot1, plot2, plot1.data, plot2.data, cog.class.summary.plot)

#### FIGURE S11: COGs ENRICHED IN WATER SAMPLES ####
suppl_plot <- rbind(
  sorting_list$Flagella,
  sorting_list$Pilus,
  sorting_list$ABC.Transport,
  sorting_list$Nitrogen,
  sorting_list$NaHK,
  sorting_list$Drugs
)

suppl_plot <- subset(suppl_plot, p.adj < 0.05)
suppl_plot <- subset(suppl_plot, abs(logFC) > 1)
suppl_plot$Group <- factor(suppl_plot$Group)
suppl_plot <- reshape2::melt(suppl_plot)

suppl_plot$Full.Name <- paste0(suppl_plot$COG, " ", suppl_plot$Name)
suppl_plot$Full.Name <- stringr::str_trunc(suppl_plot$Full.Name, 50)
suppl_plot$Full.Name <- as.factor(suppl_plot$Full.Name)
suppl_plot$value <- log(suppl_plot$value+1)

suppl_plot <- suppl_plot[order(suppl_plot$Group, suppl_plot$COG), ]
suppl_plot <- subset(suppl_plot, variable %in% c("SKE1", "SKE4", "SKE5",
                                                     "SKEW1", "SKEW4", "SKEW5"))

write.csv(suppl_plot, "trim.for.figure.csv")
suppl_plot <- read.csv("trimmed.for.figure.csv")
suppl_plot$Group <- factor(suppl_plot$Group)
suppl_plot$Full.Name <- factor(suppl_plot$Full.Name)

suppl_plot_left <- subset(suppl_plot, Group %in% c("NaHK", "Flagella", "Pilus"))
suppl_plot_right <- subset(suppl_plot, !Group %in% c("NaHK", "Flagella", "Pilus"))
suppl_plot_left$Group <- factor(suppl_plot_left$Group)
suppl_plot_right$Group <- factor(suppl_plot_right$Group)

# Create the left-hand side plot
rectangles_left <- data.frame(
  Group = levels(suppl_plot_left$Group),
  x1 = rep(0),
  x2 = rep(2.5),
  Number = rep(1),
  y1 = rep(0),
  y2 = rep(1)
)
for(j in c(1:nrow(rectangles_left))){
  rectangles_left[j,"Number"] <- (nrow(subset(suppl_plot_left, Group==levels(suppl_plot_left$Group)[j]))/6)
}
rectangles_left$Group <- factor(rectangles_left$Group)
rectangles_left <- rectangles_left %>% purrr::map_df(rev)
for(j in 1:nrow(rectangles_left)){
  rectangles_left[j,"y2"] <- sum(rectangles_left[c(1:j),"Number"])+0.5
  if(j > 1){
    rectangles_left[j,"y1"] <- rectangles_left[j-1, "y2"]}
}
rectangles_left$Number <- NULL

suppl_plot_left$Full.Name <- factor(suppl_plot_left$Full.Name)
suppl_plot_left$Full.Name <- forcats::fct_inorder(suppl_plot_left$Full.Name)
suppl_plot_left$variable <- factor(suppl_plot_left$variable)
suppl_plot_left$variable <- as.numeric(suppl_plot_left$variable)

# Create the right-hand side plot
rectangles_right <- data.frame(
  Group = levels(suppl_plot_right$Group),
  x1 = rep(6.5),
  x2 = rep(9),
  Number = rep(1),
  y1 = rep(0),
  y2 = rep(1)
)
for(j in c(1:nrow(rectangles_right))){
  rectangles_right[j,"Number"] <- (nrow(subset(suppl_plot_right, Group==levels(suppl_plot_right$Group)[j]))/6)
}
rectangles_right$Group <- factor(rectangles_right$Group)
rectangles_right <- rectangles_right %>% purrr::map_df(rev)
for(j in 1:nrow(rectangles_right)){
  rectangles_right[j,"y2"] <- sum(rectangles_right[c(1:j),"Number"])+0.5
  if(j > 1){
    rectangles_right[j,"y1"] <- rectangles_right[j-1, "y2"]}
}
rectangles_right$Number <- NULL

suppl_plot_right$Full.Name <- factor(suppl_plot_right$Full.Name)
suppl_plot_right$Full.Name <- forcats::fct_inorder(suppl_plot_right$Full.Name)
suppl_plot_right$variable <- factor(suppl_plot_right$variable)
suppl_plot_right$variable <- as.numeric(suppl_plot_right$variable)

# Generate the plots 
left <- ggplot(suppl_plot_left) +
  geom_tile(aes(x=variable+2, y=Full.Name, fill=value)) +
  scale_fill_gradient2(guide = guide_colorbar(title.position="top",
                                              title.hjust=0.5,
                                              title="Log RPKM")) +
  scale_y_discrete(limits = rev(levels(suppl_plot_left$Full.Name))) +
  scale_x_continuous(expand=c(0,0), breaks=c(3,4,5,6,7,8),
                     labels=c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")) +
  geom_rect(data=rectangles_left,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="white", alpha=0.5) +
  geom_hline(yintercept=c(rectangles_left$y1)) +
  
  geom_text(data=subset(rectangles_left, Group =="NaHK"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Sodium, hydrogen, and\npotassium transport", size=3, angle=90) +
  geom_text(data=subset(rectangles_left, Group =="Flagella"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Flagellar\nproteins", size=3, angle=90) +
  geom_text(data=subset(rectangles_left, Group =="Pilus"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Pilus\nproteins", size=3, angle=90) +
  theme_bw() + plot_theme +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black", fill=NA)) +
  labs(x="Sample", y="COG")

right <- ggplot(suppl_plot_right) +
  geom_tile(aes(x=variable, y=Full.Name, fill=value)) +
  scale_fill_gradient2(guide = guide_colorbar(title.position="top",
                                              title.hjust=0.5,
                                              title="Log RPKM")) +
  scale_y_discrete(limits = rev(levels(suppl_plot_right$Full.Name)),
                   position="right") +
  scale_x_continuous(expand=c(0,0), breaks=c(1,2,3,4,5,6),
                     labels=c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")) +
  geom_rect(data=rectangles_right,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="white", alpha=0.5) +
  geom_hline(yintercept=c(rectangles_right$y1)) +
  geom_hline(yintercept=max(rectangles_right$y2)) +
  
  geom_text(data=subset(rectangles_right, Group =="ABC.Transport"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="ABC transporters", size=3, angle=270) +
  geom_text(data=subset(rectangles_right, Group =="Nitrogen"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Nitrogen\ncycle", size=3, angle=270) +
  geom_text(data=subset(rectangles_right, Group =="Drugs"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
            label="Efflux", size=3)+
  
  theme_bw() + plot_theme +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black", fill=NA)) +
  labs(x="Sample", y="COG")

grid.arrange(left, right, ncol=2)

suppl.master <- arrangeGrob(left, right, ncol=2)
ggsave("suppl.master.map.jpg", suppl.master, width=9.5, height=11, units="in", dpi=300)

#### FIGURE S12: COOL NITROGEN CYCLE FIGURE ####
pathways <- lapply(openxlsx::getSheetNames("steroid_chitin_nitrogen.xlsx"),
                   openxlsx::read.xlsx, xlsxFile="steroid_chitin_nitrogen.xlsx")
names(pathways) <- c("steroid", "chitin", "nitrogen")

desired_order <- c("Unclassified",
                   "Unclassified bacteria",
                   
                   "Actinobacteria",
                   "Bacteroidetes",
                   "Chloroflexi",
                   "Cyanobacteria",                   
                   "Planctomycetes",
                   "Verrucomicrobia",

                   "Alphaproteobacteria",
                   "Betaproteobacteria",
                   "Other Proteobacteria")

myColors <- c("#999999", #Other#
              "#333333", #Other bacitera
              
              "#CC0033", #Actino        
              "#FF9900", ##Bacteroidetes
              "burlywood3", #Chloroflexi
              "#006600", #Cyano#              
              "#CC00FF", #Plancto              
              "#99CC99", #Verruco

              "#3399FF", ##Alpha
              "#99CCFF", #Beta
              "#6699CC") #Other

names(myColors) <- desired_order

nitrogen <- pathways$nitrogen
nitrogen <- reshape2::melt(nitrogen)
nitrogen <- subset(nitrogen, variable %in% c("SKE1", "SKE4", "SKE5",
                                           "SKEW1", "SKEW4", "SKEW5"))
nitrogen$variable <- factor(nitrogen$variable, levels = c("SKE1", "SKE4", "SKE5",
                                                        "SKEW1", "SKEW4", "SKEW5"))
nitrogen <- subset(nitrogen, domain=="Bacteria")

nitrogen$phylum <- factor(nitrogen$phylum)
levels(nitrogen$phylum)[levels(nitrogen$phylum)=="unclassified"] <- "Unclassified"
levels(nitrogen$phylum)[levels(nitrogen$phylum)=="Uncl_Bacteria"] <- "Unclassified bacteria"
levels(nitrogen$phylum)[levels(nitrogen$phylum)=="Uncl_Proteobacteria"] <- "Other Proteobacteria"
levels(nitrogen$phylum)[levels(nitrogen$phylum)=="Gammaproteobacteria"] <- "Other Proteobacteria"

#nitrogen$phylum[nitrogen$domain=="Eukarya"] <- "Unclassified"

nitrogen$phylum <- factor(nitrogen$phylum, levels=desired_order)
nitrogen$HMM_family <- factor(nitrogen$HMM_family,
                              levels=c("NapA", "NarG", "NirA", "NirB",
                                       "NirS", "NirK", "NorB", "NosZ", 
                                       "AmoA Bacteria", "NxrA", "NxrB",
                                       "NifH", "NifD"))

a <- ggplot(subset(nitrogen, HMM_family %in% c("NapA", "NarG")),
              aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~HMM_family, scales="free_y",
             labeller = label_wrap_gen(), ncol=2) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x=NULL, y="RPKM", title="Nitrate reduction") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        legend.position="none")

b <-  ggplot(subset(nitrogen, HMM_family %in% c("NirA", "NirB")),
             aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~HMM_family, scales="free_y",
             labeller = label_wrap_gen(), ncol=2) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x=NULL, y="RPKM", title="Nitrite reduction to ammonium") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title.y = element_blank(),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        legend.position="none")

c <-  ggplot(subset(nitrogen, HMM_family %in% c("NirS", "NirK", "NorB", "NosZ")),
             aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~HMM_family, scales="free_y",
             labeller = label_wrap_gen(), ncol=4) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x=NULL, y="RPKM", title="Denitrification") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        legend.position="none")

d <-  ggplot(subset(nitrogen, HMM_family %in% c("AmoA Bacteria", "NxrA", "NxrB")),
             aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~HMM_family, scales="free_y",
             labeller = label_wrap_gen(), ncol=3) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x=NULL, y="RPKM", title="Nitrification") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        legend.position="none")

e <-  ggplot(subset(nitrogen, HMM_family %in% c("NifH", "NifD")),
             aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~HMM_family, scales="free_y",
             labeller = label_wrap_gen(), ncol=2) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x=NULL, y="RPKM", title="Nitrogen fixation") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        legend.position="none")

legend <- cowplot::get_legend(ggplot(subset(nitrogen, HMM_family %in% c("NirA", "NirB")),
                                     aes(x=variable, y=value, fill=phylum)) +
                                geom_bar(stat="identity") +
                                facet_wrap(~HMM_family, scales="free_y",
                                           labeller = label_wrap_gen(), ncol=2) +
                                scale_fill_manual(values=myColors, name = "Taxon",
                                                  guide = guide_legend(ncol=2, title.position="top", title.hjust=0.5)) +
                                scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                                labs(y="RPKM", title="Nitrite reduction to ammonium") +
                                geom_vline(xintercept=3.5) +
                                theme_bw() + plot_theme +
                                theme(panel.border = element_rect(color="black", fill=NA),
                                      axis.text.x = element_text(angle=45, hjust=1),
                                      plot.title=element_text(size=10, face="bold", hjust=0.5),
        legend.position="bottom"))

grid.arrange(a, b,
             c,
             d,
             e, legend,
             layout_matrix=rbind(c(1,1,2,2),
                                 c(3,3,3,3),
                                 c(4,4,4,NA),
                                 c(5,5,6,6)))

nitrogen.bars <- arrangeGrob(a, b,
                             c,
                             d,
                             e, legend,
                             layout_matrix=rbind(c(1,1,2,2),
                                                 c(3,3,3,3),
                                                 c(4,4,4,NA),
                                                 c(5,5,6,6)))

ggsave("nitrogen.bars.jpg", nitrogen.bars, width=6, height=6.5, units="in", dpi=300)

nitrogen <- nitrogen %>%
  dplyr::group_by(HMM_family, variable) %>% dplyr::summarise(RPKM=sum(value))

nitrogen$Type <- rep("sponge")
nitrogen <- as.data.frame(nitrogen)
for(i in 1:nrow(nitrogen)){
  if(nitrogen[i,"variable"] %in% c("SKEW1", "SKEW4", "SKEW5")){
    nitrogen[i,"Type"] <- "water"
  }
}

nitrogen <- nitrogen %>%
  dplyr::group_by(HMM_family, Type) %>% dplyr::summarise(RPKM=mean(RPKM))

nitrogen.dots <- ggplot(nitrogen, aes(x=Type, y=HMM_family)) +
  geom_point(aes(size=RPKM, color=Type)) +
  scale_color_manual(values=c("darkorange1", "cornflowerblue")) +
  theme_bw() + plot_theme +
  scale_size_continuous(range=c(1,11))

ggsave("nitrogen.dots.jpg", nitrogen.dots, width=3, height=5, units="in", dpi=300)

rm(nitrogen, nitrogen.dots, nitrogen.bars, a, b, c, d, e, legend)

#### FIGURE S13: CHITIN AND STEROID DEGRADATION ENZYMES ####
steroid <- pathways$steroid
steroid <- reshape2::melt(steroid)
steroid <- subset(steroid, variable %in% c("SKE1", "SKE4", "SKE5",
                                           "SKEW1", "SKEW4", "SKEW5"))
steroid$variable <- factor(steroid$variable, levels = c("SKE1", "SKE4", "SKE5",
                                                       "SKEW1", "SKEW4", "SKEW5"))
steroid$phylum <- factor(steroid$phylum)
levels(steroid$phylum)[levels(steroid$phylum)=="unclassified"] <- "Unclassified"
levels(steroid$phylum)[levels(steroid$phylum)=="Uncl_Bacteria"] <- "Unclassified bacteria"
levels(steroid$phylum)[levels(steroid$phylum)=="Uncl_Proteobacteria"] <- "Other Proteobacteria"
steroid$phylum <- factor(steroid$phylum, levels=desired_order)

steroid.plot <- ggplot(steroid, aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~HMM_family, scales="free_y",
             labeller = label_wrap_gen()) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x="Sample", y="RPKM") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none")

chitin <- pathways$chitin
chitin <- reshape2::melt(chitin)
chitin <- subset(chitin, variable %in% c("SKE1", "SKE4", "SKE5",
                                           "SKEW1", "SKEW4", "SKEW5"))
chitin$variable <- factor(chitin$variable, levels = c("SKE1", "SKE4", "SKE5",
                                                        "SKEW1", "SKEW4", "SKEW5"))
chitin$phylum <- factor(chitin$phylum)
levels(chitin$phylum)[levels(chitin$phylum)=="unclassified"] <- "Unclassified"
levels(chitin$phylum)[levels(chitin$phylum)=="Uncl_Bacteria"] <- "Unclassified bacteria"
levels(chitin$phylum)[levels(chitin$phylum)=="Uncl_Proteobacteria"] <- "Other Proteobacteria"
levels(chitin$phylum)[levels(chitin$phylum)=="Gammaproteobacteria"] <- "Other Proteobacteria"
chitin$phylum <- factor(chitin$phylum, levels=desired_order)

chitin <- subset(chitin, query_annotation %in% c("Chitinase 2", "Chitinase 63", "Chitinase A",
                                                 "Chitinase A1", "Chitinase C", "Chitobiase"))
chitin.plot <- ggplot(chitin, aes(x=variable, y=value, fill=phylum)) +
  geom_bar(stat="identity") +
  facet_wrap(~query_annotation, scales="free_y",
             labeller = label_wrap_gen()) +
  scale_fill_manual(values=myColors, name = "Taxon",
                    guide = guide_legend(ncol=2, title.position="top")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x="Sample", y="RPKM") +
  geom_vline(xintercept=3.5) +
  theme_bw() + plot_theme +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none")

chitin.for.legend <- chitin[c(1:5),]
chitin.for.legend$phylum <- rep("Actinobacteria")
chitin.for.legend <- rbind(chitin, chitin.for.legend)
chitin.for.legend$phylum <- factor(chitin.for.legend$phylum, levels=desired_order)

legend <- cowplot::get_legend(ggplot(chitin.for.legend, aes(x=variable, y=value, fill=phylum)) +
                                geom_bar(stat="identity") +
                                facet_wrap(~query_annotation, scales="free_y",
                                           labeller = label_wrap_gen()) +
                                scale_fill_manual(values=myColors, name = "Taxon",
                                                  guide = guide_legend(ncol=1, title.position="top")) +
                                scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                                labs(x="Sample", y="RPKM") +
                                geom_vline(xintercept=3.5) +
                                theme_bw() + plot_theme +
                                theme(panel.border = element_rect(color="black", fill=NA),
                                      axis.text.x = element_text(angle=45, hjust=1),
                                      legend.position="bottom",
                                      legend.title=element_blank()))

grid.arrange(chitin.plot, steroid.plot, legend,
             layout_matrix=rbind(c(1,3),c(2,3)), widths=c(0.8,0.2))

chitin.steroid <- arrangeGrob(chitin.plot+labs(tag="a"), steroid.plot+labs(tag="b"), legend,
                              layout_matrix=rbind(c(1,3),c(2,3)), widths=c(0.75,0.25))

ggsave("chitin.steroid.jpg", chitin.steroid, width=6.5, height=7, units="in", dpi=300)

rm(chitin, steroid, chitin.plot, steroid.plot, chitin.steroid, legend)

#### FIGURE ___: MAGs AND IMPORTANT FUNCTIONS
ggplot(mag.sample_data, aes(x=PC1, y=PC2, color=Type)) + 
  geom_point(size=3)  #+
#geom_text(aes(label=MAG), color="black")


#### TABLES: DIFFERENTIALLY ABUNDANT TAXA AND COGS ####
# Table S4: Abundance of select ASVs
temp <- transform_sample_counts(pseq.clr.sub, function(x) 100*x/sum(x))
temp <- subset_taxa(temp, Genus %in% c("Sediminibacterium", "Comamonas", "Uncl_Rhodospirillales"))
temp <- as.data.frame(otu_table(temp))
temp <- merge(temp, sample_data[,c("Type","Source")], by=0, all=TRUE)
temp <- temp %>% dplyr::group_by(Type, Source) %>%
  summarise_if(is.numeric, mean) %>%
  t() %>%
  as.data.frame()
write.csv(temp, "TableS4_key.ASVs.csv")
rm(temp)

# Table S5: Sponge-specific ASVs
TabS5 <- Venn.histograms[["sponge"]]
tax <- as.data.frame(cbind(tax_table(pseq.clr.sub)))
TabS5 <- merge(TabS5, tax, by=0, all=FALSE)
TabS5 <- TabS5[order(-TabS5$prevalence, -TabS5$mean), ]
write.csv(TabS5, "TableS5_sponge.specific.ASVs.csv")
rm(tax, TabS5)

# Table S6: River-specific ASVs
write.csv(river.specific, "TableS6_river.specific.ASVs.csv")

# Table S6: COG class enrichment
write.csv(metagenome.class_fisher_test, "TableS7_fisher.exact.csv")
