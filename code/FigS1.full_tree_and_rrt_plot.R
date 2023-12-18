##########
#Load libraries
##########
library(dplyr)
library(lubridate) 

library(ape) 

library(ggplot2)
library(ggtree)
library(ggbreak) #scale_x/y_break
library(ggplotify) #as.ggplot
library(ggimage) #geom_subview
library(cowplot) #get_legend, plot_grid

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
FigS1_file_png <- "out/FigS1/FigS1.raw.png"
FigS1_file_svg <- "out/FigS1/FigS1.raw.svg"

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

##########
#Figure S1 - Full ML tree plot
##########
full_phy_plot <- ggtree(
	full_phy,
	aes(color = outlier),
) %<+% full_phy_metadata

annotated_full_phy_plot <- full_phy_plot + 
geom_tippoint(
	aes(subset = (outlier == "Yes")),
	size = 0.50, shape = 4,
	show.legend = FALSE) +
scale_color_manual(
	breaks = c("No", "Yes"),
	values = c("No" = "black", "Yes" = "red"), 
	na.translate = TRUE, na.value = "black",
	guide = "none"
) +
scale_x_continuous(
	name = "Subtitutions per site",
	breaks = c(seq(0, 0.05, 0.01), 0.1, 0.11),
	limits = c(0, 0.11),
	expand = c(0,0)) + scale_x_break(breaks = c(0.05, 0.1), expand = F) + #, space = 0.1) + 
scale_y_continuous(expand = c(0,0)) + ylab("Collection date") + #added so that graphs can be aligned
geom_vline(xintercept = c(0.05, 0.1), linetype = "dashed", linewidth = 0.1) +
ggtitle("Full ML tree") +
coord_flip() + 
theme(
	title = element_text(size = 8),

	axis.text.y = element_text(size = 6),
	axis.ticks.y = element_line(colour = "black"),
	axis.text.x = element_text(colour = "white"), #added so that graphs can be aligned
	axis.ticks.x = element_line(colour = "white"), #added so that graphs can be aligned
	axis.title.x = element_text(size = 8.5, colour = "white"),

	panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
)

annotated_full_phy_plot <- as.ggplot(print(annotated_full_phy_plot))

##########
#Figure S1, rrt plot
##########
#get rrt dists from the full ML tree ggtree object
#construct "rrt_dat"
rrt_dat <- full_phy_plot$data %>% 
	filter(isTip) %>% 
	select(label, x, Collection.date, Pango.lineage, Region, Th, outlier) %>% 
	rename(rrt_dist = x) %>%
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

#plot the rrt plot
rrt_plot <- ggplot() +
geom_pointrange(
	data = rrt_dat %>% filter(outlier == "No"),
	mapping = aes(
		x = Collection.date.med,
		y = rrt_dist,
		xmin = Collection.date.lb,
		xmax = Collection.date.ub,
		#color = "No"
	), color = "black", shape = 19, linewidth = 0.2, size = 0.005) + 
geom_pointrange(
	data = rrt_dat %>% filter(outlier == "Yes"),
	mapping = aes(
		x = Collection.date.med,
		y = rrt_dist,
		xmin = Collection.date.lb,
		xmax = Collection.date.ub,
		color = "Yes"
	), shape = 4, linewidth = 0.2, size = 0.005) + 
scale_color_manual(
	name = NULL,
	breaks = c("Yes"), values = c("Yes" = "red"), #depict only the red outlier dots
	labels = c("Outlier")
	) +
scale_y_continuous(
	name = "Root-to-tip distance",
	breaks = c(seq(0, 0.05, 0.01), 0.1, 0.11),
	limits = c(0, 0.11),
	expand = c(0,0)) + scale_y_break(c(0.05, 0.1), expand = F) + 
scale_x_continuous(expand = c(0,0)) + xlab("Collection date") +
ggtitle("Root-to-tip distance vs collection date") +
theme_bw() +
theme(
	legend.background = element_blank(),
	legend.key.size = unit(0.1, 'cm'),
	legend.title = element_text(size = 8),
	legend.text = element_text(size = 6),

	#panel.border = element_blank(),
	panel.grid = element_blank(),
	plot.background = element_blank(),

	title = element_text(size = 8),
	axis.text = element_text(size = 6),
)

# Place the legend 
annotated_rrt_plot <- 
as.ggplot(print(rrt_plot + theme(legend.position = "none"))) + 
geom_subview(x = 0.85, y = 0.05, subview = get_legend(rrt_plot))

##########
#Combine the plots together
##########
FigS1 <- plot_grid(annotated_full_phy_plot, annotated_rrt_plot)

##########
#Save the plot to file
##########
ggsave(FigS1_file_png, plot = FigS1, width = 16, height = 8, units = "cm", dpi = 300)
ggsave(FigS1_file_svg, plot = FigS1, width = 16, height = 8, units = "cm", dpi = 300)
