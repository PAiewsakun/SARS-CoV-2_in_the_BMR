##########
#Load libraries
##########
library(ape) 
library(treeio) #tree_subset

library(dplyr)
library(tibble) #column_to_rownames
library(lubridate) #decimal_date

library(ggplot2)
library(ggtree)
library(ggnewscale) #new_scale_fill
library(ggbreak) #scale_x/y_break

library(cowplot)

##########
#Define paths
##########
#to working dir
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

#to data files
timed_phy_file <- "data/Sup-Dat-Fig-3_timed-calibrated_phy.nexus"
full_phy_metadata_file <- "data/Sup-Table-1_Full_ML_phy_metadata.txt"
Ne_file <- "data/Sup-Dat-Fig-3_BMR_Ne.txt"

#to output files 
Fig3_file_png <- "out/Fig3/Fig3.raw.png"
Fig3_file_svg <- "out/Fig3/Fig3.raw.svg"

##########
#Define some "useful" variables
##########
#provices in the BMR
BMR_provinces <- c("Bangkok", "Nakhon Pathom", "Pathum Thani", "Nonthaburi", "Samut Prakan", "Samut Sakhon")

#region colours
Region_col <- c(
	"BMR" = "#78c5d6", 
	"non-BMR" = "#459ba8", 
	"Unknown" = "#175b65",

	"North America" = "#bf62a6", 
	"South America" = "#e868a2",
	"Europe" = "#d64e12",
	"Africa" = "#f28c33", 
	"Asia" = "#c5d647", 
	"Oceania" = "#efdf48"
)

#lineage colours
Lineage_col <- c(
	"A.6" = "#30ACEC",
	"B.1" = "#3a8ec5",
	"B.1.1" = "#44709d",

	"B.1.36.16" = "#80C34F",

	"B.1.1.7" = "#E29D3E",

	"B.1.617.2" = "#9b8357",
	"AY.85" = "#b96749",
	"AY.30" = "#D64A3B",

	"BA.1" = "#d64787",
	"BA.1.1" = "#ca4f9e",
	"BA.1.17" = "#be57b4",
	"BA.1.17.2" = "#b25ecb",

	"BA.2" = "#a666e1",
	"BA.2.3" = "#9e65bd",
	"BA.2.10" = "#956498",
	"BA.2.10.1" = "#8d6374",
	"BA.2.9" = "#8f777b",
	"BA.2.9.5" = "#928c81",

	"BA.4.1" = "#94A088",

	"BA.5.1" = "#a3a985",
	"BA.5.2" = "#b3b383",
	"BA.5.2.1" = "#C2BC80",

	"Others" = "#bebebe", #grey
	"Unassigned" = "#5A5A5A" #darkgrey
)

#COVID wave label
Wave_label <- rbind(
	c(x = "2020-03-01", xend = "2020-05-31", y = 35000, yend = 35000, label = "1st wave\n(01/03/2020-31/05/2020)"),
	c(x = "2020-12-18", xend = "2021-02-28", y = 35000, yend = 35000, label = "2nd wave\n(18/12/2020-28/02/2021)"),
	c(x = "2021-04-01", xend = "2021-07-05", y = 36000, yend = 36000, label = "3rd wave\n(01/04/2021-05/07/2021)"),
	c(x = "2021-07-06", xend = "2022-01-05", y = 35000, yend = 35000, label = "4th wave\n(06/07/2021-05/01/2022)"),
	c(x = "2022-01-06", xend = "2022-10-10", y = 35000, yend = 35000, label = "5th wave\n(06/01/2022-10/10/2022)")
	) %>% as.data.frame %>% 
	mutate(
		x = as.numeric(decimal_date(as.Date(x))),
		xend = as.numeric(decimal_date(as.Date(xend))),
		y = as.numeric(y),
		yend = as.numeric(yend)
	)

#function to plot the phylogeny
phy_plot_fun <- function(
	phy, 
	mrsd = NULL, #most recent sample date - collection date of the youngest / most recent sample in the tree
	phy_metadata, 

	PANGO_lineage_col,
	Geo_region_col,

	x_ticks, x_lim,
	y_ticks, y_lim,

	gheatmap.font.size = 3,
	PANGO_lineage.legend.nrow = 4,
	Geo_region.legend.nrow = 3
	){

	#tree plot
	phy_plot <- ggtree(
		phy, 
		mrsd = mrsd, 
		size = 0.2) %<+% phy_metadata
	
	phy_plot <- phy_plot + 
		geom_tippoint(
			mapping = aes(color = Pango.lineage),
			size = 0.5,
			show.legend = FALSE) +
		scale_color_manual(
			breaks = names(PANGO_lineage_col),
			values = PANGO_lineage_col
		)

	#thickness of the heatmap, defined wrt to the tree's height
	o = (max(phy_plot$data$x) - min(phy_plot$data$x))/12

	#add PANGO lineage heatmap to the plot
	phy_plot <- gheatmap(
		p = phy_plot,
		data = phy_metadata %>% column_to_rownames(var = "Accession.ID") %>% select(Pango.lineage),
		offset = o*0, 
		width = 0.05, 
		color = setNames(PANGO_lineage_col[phy_metadata$Pango.lineage], phy_metadata$Accession.ID)[phy_plot$data %>% filter(isTip) %>% pull(label)], #need tile border colours to make the heatmap "looks" right; otherwise the colours would appear faded!
		colnames = TRUE, custom_column_labels = c("PANGO\nlineage"), colnames_position = "top",
		colnames_angle = 90, hjust = 0, #colnames_offset_x = -o/2, 
		font.size = gheatmap.font.size) +
		scale_fill_manual(
			breaks = names(PANGO_lineage_col),
			values = PANGO_lineage_col, 
			guide = guide_legend(title = "PANGO\nlineage", nrow = PANGO_lineage.legend.nrow, order = 1)
		)

	#add geo heatmap to the plot
	phy_plot <- gheatmap(
		p = phy_plot + new_scale_fill(),
		data = phy_metadata %>% column_to_rownames(var = "Accession.ID") %>% select(Region), 
		offset = o*1,
		width = 0.05,
		color = setNames(Geo_region_col[phy_metadata$Region], phy_metadata$Accession.ID)[phy_plot$data %>% filter(isTip) %>% pull(label)],
		colnames = TRUE, custom_column_labels = c("Region"), colnames_position = "top",
		colnames_angle = 90, hjust = 0, #colnames_offset_x = -o/2, 
		font.size = gheatmap.font.size) + 
		scale_fill_manual(
			breaks = names(Geo_region_col),
			values = Geo_region_col,
			guide = guide_legend(title = "Region", nrow = Geo_region.legend.nrow, order = 2)
		)

	#format the plot
	phy_plot <- phy_plot + 
		scale_x_ggtree() + #set the x axis more reasonably.
		theme_tree2() +
		theme(
			legend.position = "top",
			legend.direction = "horizontal", legend.box = "horizontal",
			legend.background = element_blank(),
			legend.margin = margin(c(0, 0, 0, 0)),
			legend.key.size = unit(0.25, 'cm'),
			legend.key = element_blank(),
			legend.title = element_text(size = 8),
			legend.text = element_text(size = 6),

			axis.title = element_text(size = 8),
			axis.text = element_text(size = 6),
			axis.text.x = element_text(size = 6, vjust = 0),

			plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
		) + 
		scale_y_continuous(
			breaks = y_ticks,
			limits = y_lim,
			expand = c(0, 0)
		) +
		scale_x_continuous(
			breaks = x_ticks,
			limits = x_lim,
			expand = expansion(mult = c(0, 0))) 

	return(phy_plot)
}

##########
#Load data
##########
#time-calibrated tree
timed_phy <- read.nexus(timed_phy_file)

#Metadata
full_phy_metadata <- read.table(full_phy_metadata_file, header = TRUE, sep = "\t", quote = "") %>%
	select(Accession.ID, Collection.date, Pango.lineage, Continent, Country, Province, Region) %>%
	mutate(
		Th = ifelse(Region %in% c("BMR", "non-BMR", "Unknown"), "Th", "non-Th"),
		outlier = ifelse(Accession.ID %in% timed_phy$tip.label, "No", "Yes")
	)

full_phy_metadata$Pango.lineage[full_phy_metadata$Pango.lineage == "unclassifiable"] <- "Unassigned"
full_phy_metadata$Pango.lineage[!(full_phy_metadata$Pango.lineage %in% names(Lineage_col))] <- "Others"

#filter out "outliers" from the full_phy_metadata df
timed_phy_metadata <- full_phy_metadata %>% filter(outlier == "No") 
timed_phy_metadata %>% count(outlier)

#BMR Ne estimate
Ne_est <- read.table(Ne_file, header = TRUE, sep = "\t", quote = "")

##########
#Figure 3 - time-calibrated tree plot
##########
#make the plot
timed_phy_main_plot <- phy_plot_fun(
	phy = timed_phy, 
	phy_metadata = timed_phy_metadata, 
	mrsd = timed_phy_metadata %>% pull(Collection.date) %>% as.Date() %>% sort(decreasing = T) %>% head(1),
	PANGO_lineage_col = Lineage_col,
	Geo_region_col = Region_col,
	x_ticks = c(2020, 2021, 2022), x_lim = c(2019.90, 2023.35),
	y_ticks = seq(0, nrow(timed_phy_metadata), 1000), y_lim = c(0, nrow(timed_phy_metadata)*1.08)
	)

#make label df for major lineages
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
get_min_max_y_ind_fun <- function(phy, tips, phy_plot_dat){
	phy <- tree_subset(
		tree = phy,
		node = MRCA(phy, tips),
		levels_back = 0
	)
	return(
	phy_plot_dat %>% filter(label %in% phy$tip.label) %>% 
	summarize(min_y = min(y), max_y = max(y), Pango.lineage = Mode(Pango.lineage))
	)
}

lineage_lab_df <- rbind(
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_12717909", "EPI_ISL_12717950"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_5596905", "EPI_ISL_1074022"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("OL629465", "EPI_ISL_2433427"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_3407839", "EPI_ISL_6695446"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_7921019", "EPI_ISL_8097517"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_12980932", "EPI_ISL_12862252"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_13766253", "EPI_ISL_13388442"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("OX284045", "EPI_ISL_14589616"), timed_phy_main_plot$data),
	get_min_max_y_ind_fun(timed_phy, c("EPI_ISL_11442009", "EPI_ISL_10511333"), timed_phy_main_plot$data)
) %>% 
mutate(Pango.lineage = forcats::fct_recode(Pango.lineage, "BA.2" = "BA.2", "BA.4" = "BA.4.1", "BA.5" = "BA.5.2", "BA.1" = "BA.1.1"))

#annotate the plot
annotated_timed_phy_main_plot <- timed_phy_main_plot + 
geom_segment(
	data = Wave_label, 
	mapping = aes(x = x, xend = x, y = 0, yend = 35000),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label, 
	mapping = aes(x = xend, xend = xend, y = 0, yend = 35000),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label, 
	mapping = aes(x = x, xend = xend, y = 35000, yend = 35000),
	colour = "black", linewidth = 0.5
) +
geom_text(
	data = Wave_label, 
	aes(x = (x + xend)/2, y = y, label = label), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F
) +

geom_segment(
	data = lineage_lab_df, 
	mapping = aes(x = 2022.9, xend = 2022.9, y = min_y, yend = max_y, color = Pango.lineage),
	linewidth = 0.5, show.legend = F
) +
geom_text(
	data = lineage_lab_df, 
	aes(x = 2022.9, y = (min_y + max_y)/2, label = Pango.lineage, color = Pango.lineage), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F, angle = 90
)

##########
#Figure 3 - Ne estimate plot
##########
ylim_ub <- 1000
Ne_plot <- Ne_est %>% 
	mutate(across(diff_ne_median:diff_ne_CI_ub, ~ ifelse(.x < 1, 1, .x))) %>% 
	mutate(across(diff_ne_median:diff_ne_CI_ub, ~ ifelse(.x > ylim_ub, ylim_ub, .x))) %>% 
ggplot() +
geom_line(aes(x = t, y = diff_ne_median), linewidth = 1) +
geom_ribbon(aes(x = t, ymin = diff_ne_CI_lb, ymax = diff_ne_CI_ub), alpha = 0.2) +
xlab("Time") + ylab(bquote(atop(Ne["BMR"] == phantom(), Ne["all"] - Ne["non-BMR sequences"]))) + 
scale_x_continuous(
	lim = c(2019.90, 2023.35),
	breaks = c(2020, 2021, 2022),
	expand = expansion(mult = c(0, 0)),
	) +
scale_y_continuous(lim = c(0, ylim_ub), expand = c(0,0)) +
theme_classic() +
theme(
	plot.margin = unit(c(0, 0, 0, 0), "pt"),

	axis.title = element_text(size = 8),
	axis.text = element_text(size = 6),
	axis.text.x = element_text(size = 6, vjust = 0),
	axis.text.y = element_text(size = 6, vjust = 0),

	axis.text.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.ticks.y.right = element_blank()
)
#add the wave lines
Ne_plot <- Ne_plot +
geom_segment(
	data = Wave_label, 
	mapping = aes(x = x, xend = x, y = 0, yend = ylim_ub),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label, 
	mapping = aes(x = xend, xend = xend, y = 0, yend = ylim_ub),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) 

##########
#Combine the plots together
##########
Fig3 <- plot_grid(
	annotated_timed_phy_main_plot + 
		ylab("Time-calibrated phylogeny") + 
		theme(
			legend.justification = c(1.1, 1), 
			legend.box.spacing = unit(0, "pt"), 

			axis.text.x = element_blank(),
			axis.ticks.x = element_blank()
		),
	Ne_plot,
	align = "v", ncol = 1, rel_heights = c(0.80, 0.20))

##########
#Save the plot to file
##########
ggsave(Fig3_file_png, plot = Fig3, width = 16, height = 24, units = "cm", dpi = 600)
ggsave(Fig3_file_svg, plot = Fig3, width = 16, height = 24, units = "cm", dpi = 300)