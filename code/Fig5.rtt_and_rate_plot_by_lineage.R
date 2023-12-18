##########
#Load libraries
##########
library(dplyr)
library(lubridate) 

library(ape) 
library(treeio) #tree_subset

library(emmeans)

library(ggplot2)
library(ggtree)
library(cowplot) # plot_grid

##########
#Define paths
##########
#to working dir
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

#to data files
timed_phy_file <- "data/Sup-Dat-Fig-3_timed-calibrated_phy.nexus"
full_phy_file <- "data/Sup-Dat-Fig-S1_Full_ML_phy.treefile"
full_phy_metadata_file <- "data/Sup-Table-1_Full_ML_phy_metadata.txt"

#to output files 
Fig5_file_png <- "out/Fig5/Fig5.raw.png"
Fig5_file_svg <- "out/Fig5/Fig5.raw.svg"

##########
#Define some "useful" variables
##########
#lineage colours
Lineage_col <- c(
	"A.6" = "#30ACEC",

	"B.1.36.16" = "#80C34F",

	"B.1.1.7" = "#E29D3E",

	"AY.85" = "#b96749",
	"AY.30" = "#D64A3B",

	"BA.1" = "#d64787",

	"BA.2" = "#a666e1",

	"BA.4" = "#94A088",

	"BA.5" = "#b3b383"
)

##########
#Load data
##########
timed_phy <- read.nexus(timed_phy_file)
full_phy <- read.tree(full_phy_file)
full_phy_metadata <- read.table(full_phy_metadata_file, header = TRUE, sep = "\t", quote = "") %>%
	select(Accession.ID, Collection.date, Pango.lineage, Continent, Country, Province, Region) %>%
	mutate(
		Th = ifelse(Region %in% c("BMR", "non-BMR", "Unknown"), "Th", "non-Th"),
		outlier = ifelse(Accession.ID %in% timed_phy$tip.label, "No", "Yes")
	)

ML_phy <- drop.tip(full_phy, full_phy_metadata %>% filter(outlier == "Yes") %>% pull(Accession.ID))

#get an ML tree for each lineage 
A_6_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_12717909", "EPI_ISL_12717950")),
	levels_back = 0
) 
B_1_36_16_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_5596905", "EPI_ISL_1074022")),
	levels_back = 0
) 
B_1_1_7_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("OL629465", "EPI_ISL_2433427")),
	levels_back = 0
) 
AY_30_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_3407839", "EPI_ISL_6695446")),
	levels_back = 0
) 
AY_85_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_7921019", "EPI_ISL_8097517")),
	levels_back = 0
) 
BA_2_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_12980932", "EPI_ISL_12862252")),
	levels_back = 0
) 
BA_4_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_13766253", "EPI_ISL_13388442")),
	levels_back = 0
) 
BA_5_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("OX284045", "EPI_ISL_14589616")),
	levels_back = 0
) 
BA_1_phy <- tree_subset(
	tree = ML_phy,
	node = MRCA(ML_phy, c("EPI_ISL_11442009", "EPI_ISL_10511333")),
	levels_back = 0
) 

#construct "rtt_dat"
ML_phy_plot <- ggtree(ML_phy) %<+% full_phy_metadata
rtt_dat <- ML_phy_plot$data %>% 
filter(isTip) %>% 
mutate(
	Pango.lineage = ifelse(label %in% A_6_phy$tip.label, "A.6", Pango.lineage),
	Pango.lineage = ifelse(label %in% B_1_36_16_phy$tip.label, "B.1.36.16", Pango.lineage),
	Pango.lineage = ifelse(label %in% B_1_1_7_phy$tip.label, "B.1.1.7", Pango.lineage),
	Pango.lineage = ifelse(label %in% AY_30_phy$tip.label, "AY.30", Pango.lineage),
	Pango.lineage = ifelse(label %in% AY_85_phy$tip.label, "AY.85", Pango.lineage),
	Pango.lineage = ifelse(label %in% BA_1_phy$tip.label, "BA.1", Pango.lineage),
	Pango.lineage = ifelse(label %in% BA_2_phy$tip.label, "BA.2", Pango.lineage),
	Pango.lineage = ifelse(label %in% BA_4_phy$tip.label, "BA.4", Pango.lineage),
	Pango.lineage = ifelse(label %in% BA_5_phy$tip.label, "BA.5", Pango.lineage)
) %>%
filter(label %in% c(
	A_6_phy$tip.label, 
	B_1_36_16_phy$tip.label,
	B_1_1_7_phy$tip.label,
	AY_30_phy$tip.label,
	AY_85_phy$tip.label,
	BA_1_phy$tip.label,
	BA_2_phy$tip.label,
	BA_4_phy$tip.label,
	BA_5_phy$tip.label
	)
) %>% 
select(label, x, Collection.date, Pango.lineage) %>% 
rename(rtt_dist = x) %>%
mutate(
	Collection.date.med = ifelse(
		!is.na(ymd(Collection.date)),
		decimal_date(ymd(Collection.date)),
		ifelse(
			!is.na(ym(Collection.date)),
			(decimal_date(floor_date(ym(Collection.date), "month")) + decimal_date(ceiling_date(ym(Collection.date), "month") - days(1)))/2,
			as.numeric(Collection.date) + 0.5
		)
	),
	Collection.date.lb = ifelse(
		!is.na(ymd(Collection.date)),
		NA,
		ifelse(
			!is.na(ym(Collection.date)),
			decimal_date(floor_date(ym(Collection.date), "month")),
			as.numeric(Collection.date)
		)
	),
	Collection.date.ub = ifelse(
		!is.na(ymd(Collection.date)),
		NA,
		ifelse(
			!is.na(ym(Collection.date)),
			decimal_date(ceiling_date(ym(Collection.date), "month")),
			decimal_date(as.Date(ISOdate(Collection.date, 12, 31)))
		)			
	)
)

#Make a best-fit model
mdl <- lm(rtt_dist ~ Pango.lineage + Collection.date.med:Pango.lineage, rtt_dat)

#add predicted values to the data
rtt_dat <- cbind(rtt_dat, predict(mdl, interval = "confidence"))

##########
##########
#Plot Figure 4 - rtt plot by lineage
##########
##########

##########
#Plot rtt plot by lineage
##########
rtt_plot <- 
ggplot(data = rtt_dat, mapping = aes(x = Collection.date.med, y = rtt_dist, group = Pango.lineage, color = Pango.lineage)) +
geom_pointrange(aes(xmin = Collection.date.lb, xmax = Collection.date.ub), shape = 20, linewidth = 0.1, size = 0.25, stroke = 0) + 
geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25, color = "black", linewidth = 0.1) + 
geom_line(aes(y = fit), linewidth = 0.5, color = "black") + 
geom_line(aes(y = fit), linewidth = 0.25) + 
scale_color_manual(
	name = "PANGO lineage",
	breaks = names(Lineage_col),
	values = Lineage_col
) +
scale_y_continuous(name = "Root-to-tip distance") + 
scale_x_continuous(name = "Collection date") + 
ggtitle("Root-to-tip regression") +
theme_bw() +
theme(
	legend.justification = c(0, 1), legend.position = c(0, 1), 
	legend.background = element_blank(),
	legend.key.size = unit(0.1, 'cm'),
	legend.title = element_text(size = 8),
	legend.text = element_text(size = 6),

	title = element_text(size = 8),
	axis.text = element_text(size = 6),
)

##########
#Plot rate estimate by lineage
##########
#make rate df
rate_dat <- emtrends(mdl, ~ Pango.lineage, var = "Collection.date.med") %>% data.frame %>%
mutate(Pango.lineage = factor(Pango.lineage, levels = c("Overall", "A.6", "B.1.36.16", "B.1.1.7", "AY.30", "AY.85", "BA.1", "BA.2", "BA.4", "BA.5"))) %>%
mutate(rate_info = paste(
	formatC(Collection.date.med.trend, digits = 3, format = "E"), "\n",
	"(",
	formatC(lower.CL, digits = 3, format = "E"), 
	"-",
	formatC(upper.CL, digits = 3, format = "E"), 
	 ")", sep = ""
	),
	rate_ratio = Collection.date.med.trend/0.00119193
) 

#make rate plot
rate_plot <- 
ggplot(data = rate_dat, aes(x = Pango.lineage, y = Collection.date.med.trend, group = Pango.lineage, color = Pango.lineage)) +
geom_pointrange(mapping = aes(ymin = lower.CL, ymax = upper.CL), size = 0.05) +
geom_text(aes(y = upper.CL, label = rate_info), angle = 45, hjust = 0.5, vjust = -.25, size = 6*0.35) +
scale_color_manual(
			breaks = names(Lineage_col),
			values = Lineage_col
		) +
geom_hline(yintercept = 0, linetype = "solid") +
geom_hline(yintercept = 0.00119193, linetype = "dashed") +
scale_y_continuous(name = "Rate estimate (s/n/y)") + 
scale_x_discrete(name = "PANGO lineage") + 
ggtitle("Evolutionary rate estimates") +
theme_bw() +
theme(
	legend.position = "none",
	title = element_text(size = 8),
	axis.text = element_text(size = 6),
)

##########
#Combine the plots together
##########
Fig5 <- plot_grid(rtt_plot, rate_plot, ncol = 2, rel_widths = c(1/3,2/3))

##########
#Save the plot to file
##########
ggsave(Fig5_file_png, plot = Fig5, width = 16, height = 8, units = "cm", dpi = 600)
ggsave(Fig5_file_svg, plot = Fig5, width = 16, height = 8, units = "cm", dpi = 600)