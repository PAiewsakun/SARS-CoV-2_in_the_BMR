##########
#Load libraries
##########
library(dplyr)
library(tibble) #column_to_rownames
library(lubridate) #decimal_date

library(ape) 
library(adephylo)
library(tidytree) #parent, MRCA
library(treeio) #tree_subset

library(ggplot2)
library(ggtree)
library(ggbreak) #scale_x/y_break
library(ggnewscale) #new_scale_fill

library(cowplot)

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
FigS2_file_png <- "out/FigS2/FigS2.raw.png"
FigS2_file_svg <- "out/FigS2/FigS2.raw.svg"

FigS3_file_png <- "out/FigS3/FigS3.raw.png"
FigS3_file_svg <- "out/FigS3/FigS3.raw.svg"

FigS4_file_png <- "out/FigS4/FigS4.raw.png"
FigS4_file_svg <- "out/FigS4/FigS4.raw.svg"

FigS5_file_png <- "out/FigS5/FigS5.raw.png"
FigS5_file_svg <- "out/FigS5/FigS5.raw.svg"

FigS6_file_png <- "out/FigS6/FigS6.raw.png"
FigS6_file_svg <- "out/FigS6/FigS6.raw.svg"

##########
#define some "useful" variables
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
			legend.position = "bottom",
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

#Wave label df
Wave_label <- 
	rbind(
	c(x = "2020-03-01", xend = "2020-05-31", label = "1st wave\n(01/03/2020-31/05/2020)"),
	c(x = "2020-12-18", xend = "2021-02-28", label = "2nd wave\n(18/12/2020-28/02/2021)"),
	c(x = "2021-04-01", xend = "2021-07-05", label = "3rd wave\n(01/04/2021-05/07/2021)"),
	c(x = "2021-07-06", xend = "2022-01-05", label = "4th wave\n(06/07/2021-05/01/2022)"),
	c(x = "2022-01-06", xend = "2022-10-10", label = "5th wave\n(06/01/2022-10/10/2022)")
	) %>% as.data.frame %>% 
	mutate(x = as.numeric(decimal_date(as.Date(x))), xend = as.numeric(decimal_date(as.Date(xend))))

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

full_phy_metadata$Pango.lineage[full_phy_metadata$Pango.lineage == "unclassifiable"] <- "Unassigned"
full_phy_metadata$Pango.lineage[!(full_phy_metadata$Pango.lineage %in% names(Lineage_col))] <- "Others"

#make a full timed_phy_plot to get the data
timed_phy_metadata <- full_phy_metadata %>% filter(outlier == "No") 
timed_phy_plot <- ggtree(
	timed_phy,
	mrsd = timed_phy_metadata %>% pull(Collection.date) %>% as.Date() %>% sort(decreasing = T) %>% head(1)
) %<+% timed_phy_metadata

##########
#Figure S2 - 1st wave tree, "2020-03-01" - "2020-05-31"
##########
first_wave_seq <- 
	timed_phy_plot$data %>% 
	filter(isTip) %>% 
	filter(x <= 2020.5) %>% 
	pull(label)

first_wave_timed_phy <- keep.tip(timed_phy, tip = first_wave_seq)
first_wave_timed_phy_metadata <- timed_phy_metadata %>% filter(Accession.ID %in% first_wave_timed_phy$tip.label)

first_wave_timed_phy_plot <- phy_plot_fun(
	phy = first_wave_timed_phy, 
	phy_metadata = first_wave_timed_phy_metadata, 
	mrsd = date_decimal(timed_phy_plot$data %>% filter(label %in% first_wave_timed_phy$tip.label) %>% pull(x) %>% max),

	PANGO_lineage_col = Lineage_col,
	Geo_region_col = Region_col,

	x_ticks = seq(2020, 2020.5, 0.1), x_lim = c(2019.90, 2020.6),
	y_ticks = NULL, y_lim = c(0, nrow(first_wave_timed_phy_metadata)*1.08),

	gheatmap.font.size = 2,
	PANGO_lineage.legend.nrow = 2,
	Geo_region.legend.nrow = 2
	)
 
annotated_first_wave_timed_phy_plot <- first_wave_timed_phy_plot + 
geom_tiplab(aes(subset = (Pango.lineage == "A.6") & (Th == "non-Th"), label = Country),
	size = 2, 
	offset = 0.01, 
	angle = 45) + 
geom_segment(
	data = Wave_label[1,], 
	mapping = aes(x = x, xend = x, y = 0, yend = nrow(first_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[1,], 
	mapping = aes(x = xend, xend = xend, y = 0, yend = nrow(first_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[1,], 
	mapping = aes(x = x, xend = xend, y = nrow(first_wave_timed_phy_metadata)*1.02, yend = nrow(first_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.5
) +
geom_text(
	data = Wave_label[1,], 
	aes(x = (x + xend)/2, y = nrow(first_wave_timed_phy_metadata)*1.02, label = label), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F
); annotated_first_wave_timed_phy_plot

#earliest sequence in the A.6 cluster
annotated_first_wave_timed_phy_plot$data %>% filter(label %in% 
tree_subset(
	tree = first_wave_timed_phy,
	node = MRCA(annotated_first_wave_timed_phy_plot, c("EPI_ISL_918162", "EPI_ISL_512857")),
	levels_back = 0
)$tip.label) %>% arrange(x) %>% pull(Collection.date) %>% head(1) 

#tMRCA of A.6
annotated_first_wave_timed_phy_plot$data %>% filter(isTip & Pango.lineage == "A.6") %>% arrange(x) %>% data.frame %>%head
annotated_first_wave_timed_phy_plot$data %>% filter(node == MRCA(annotated_first_wave_timed_phy_plot, c("EPI_ISL_918162", "EPI_ISL_512857"))) %>% 
pull(x) -> MRCA_A_6
MRCA_A_6; date_decimal(MRCA_A_6)

#save to file
ggsave(FigS2_file_png, plot = annotated_first_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)
ggsave(FigS2_file_svg, plot = annotated_first_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)

##########
#Figure S3 - 2nd wave B.1.36.16 tree: "2020-12-18" - "2021-02-28"
##########
#subset B.1.36.16 tree
second_wave_timed_phy <- tree_subset(
	tree = timed_phy,
	node = MRCA(timed_phy, (timed_phy_plot$data %>% filter(Pango.lineage == "B.1.36.16") %>% pull(label))),
	levels_back = 2
)

#make metadata df
second_wave_timed_phy_metadata <- timed_phy_metadata %>% filter(Accession.ID %in% second_wave_timed_phy$tip.label)

#remove long branches
long_terminal_B.1.36.16_tips <- timed_phy_plot$data %>% 
filter(isTip, label %in% second_wave_timed_phy$tip.label) %>%
filter(x > 2021.5) %>% pull(label)

second_wave_timed_phy <- drop.tip(second_wave_timed_phy, long_terminal_B.1.36.16_tips)
second_wave_timed_phy_metadata <- second_wave_timed_phy_metadata %>% filter(!(Accession.ID %in% long_terminal_B.1.36.16_tips))

#plot
second_wave_timed_phy_plot <- phy_plot_fun(
	phy = second_wave_timed_phy, 
	phy_metadata = second_wave_timed_phy_metadata, 
	mrsd = date_decimal(timed_phy_plot$data %>% filter(label %in% second_wave_timed_phy$tip.label) %>% pull(x) %>% max),

	PANGO_lineage_col = Lineage_col,
	Geo_region_col = Region_col,

	x_ticks = seq(2020, 2021.4, 0.1), x_lim = c(2020.45, 2021.6),
	y_ticks = NULL, y_lim = c(0, nrow(second_wave_timed_phy_metadata)*1.08),

	gheatmap.font.size = 2,
	PANGO_lineage.legend.nrow = 1,
	Geo_region.legend.nrow = 1
	)

annotated_second_wave_timed_phy_plot <- second_wave_timed_phy_plot + 
geom_tiplab(aes(subset = (Pango.lineage == "B.1.36.16") & (Th == "non-Th"), label = Country),
	size = 2, 
	offset = 0.01, 
	angle = 45) + 
geom_segment(
	data = Wave_label[2,], 
	mapping = aes(x = x, xend = x, y = 0, yend = nrow(second_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[2,], 
	mapping = aes(x = xend, xend = xend, y = 0, yend = nrow(second_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[2,], 
	mapping = aes(x = x, xend = xend, y = nrow(second_wave_timed_phy_metadata)*1.02, yend = nrow(second_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.5
) +
geom_text(
	data = Wave_label[2,], 
	aes(x = (x + xend)/2, y = nrow(second_wave_timed_phy_metadata)*1.02, label = label), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F
); annotated_second_wave_timed_phy_plot

#earliest sequence in the 3 clusters
annotated_second_wave_timed_phy_plot$data %>% filter(Country == "Thailand") %>% arrange(x) %>% pull(Collection.date) %>% head(1)

#tMRCA of the 3 B.1.36.16 clusters
#cluster 1
annotated_second_wave_timed_phy_plot$data %>% 
filter(node == MRCA(annotated_second_wave_timed_phy_plot, c("EPI_ISL_1138846", "EPI_ISL_1073955"))) %>% 
pull(x) -> MRCA_B_1_36_16_cluster1
MRCA_B_1_36_16_cluster1; date_decimal(MRCA_B_1_36_16_cluster1)

#cluster 2
annotated_second_wave_timed_phy_plot$data %>% 
filter(node == MRCA(annotated_second_wave_timed_phy_plot, c("EPI_ISL_3892049", "EPI_ISL_7852825"))) %>% 
pull(x) -> MRCA_B_1_36_16_cluster2
MRCA_B_1_36_16_cluster2; date_decimal(MRCA_B_1_36_16_cluster2)

#cluster 3
annotated_second_wave_timed_phy_plot$data %>% 
filter(node == MRCA(annotated_second_wave_timed_phy_plot, c("EPI_ISL_768529", "EPI_ISL_707709"))) %>% 
pull(x) -> MRCA_B_1_36_16_cluster3
MRCA_B_1_36_16_cluster3; date_decimal(MRCA_B_1_36_16_cluster3)

#save to file
ggsave(FigS3_file_png, plot = annotated_second_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)
ggsave(FigS3_file_svg, plot = annotated_second_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)

##########
#Figure S4 - 3rd wave B.1.1.7 tree: "2021-04-01" - "2021-07-05"
##########
#subset B.1.1.7 tree
third_wave_timed_phy <- tree_subset(
	tree = timed_phy,
	node = MRCA(timed_phy, c("EPI_ISL_5142766", "EPI_ISL_5161358")),
	levels_back = 0, root_edge = FALSE
)

#make metadata df
third_wave_timed_phy_metadata <- timed_phy_metadata %>% filter(Accession.ID %in% third_wave_timed_phy$tip.label)

#plot
third_wave_timed_phy_plot <- phy_plot_fun(
	phy = third_wave_timed_phy, 
	phy_metadata = third_wave_timed_phy_metadata, 
	mrsd = date_decimal(timed_phy_plot$data %>% filter(label %in% third_wave_timed_phy$tip.label) %>% pull(x) %>% max),

	PANGO_lineage_col = Lineage_col,
	Geo_region_col = Region_col,

	x_ticks = seq(2020, 2021.8, 0.1), x_lim = c(2020.7, 2022.05),
	y_ticks = NULL, y_lim = c(0, nrow(third_wave_timed_phy_metadata)*1.08),

	gheatmap.font.size = 2,
	PANGO_lineage.legend.nrow = 2,
	Geo_region.legend.nrow = 2
	)

annotated_third_wave_timed_phy_plot <- third_wave_timed_phy_plot + 
geom_tiplab(aes(subset = label %in% c("EPI_ISL_1098606", "EPI_ISL_1098607", "EPI_ISL_1098608"), label = Country),
	size = 1.75, offset = 0.01, angle = 45) + 
geom_segment(
	data = Wave_label[3,], 
	mapping = aes(x = x, xend = x, y = 0, yend = nrow(third_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[3,], 
	mapping = aes(x = xend, xend = xend, y = 0, yend = nrow(third_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[3,], 
	mapping = aes(x = x, xend = xend, y = nrow(third_wave_timed_phy_metadata)*1.02, yend = nrow(third_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.5
) +
geom_text(
	data = Wave_label[3,], 
	aes(x = (x + xend)/2, y = nrow(third_wave_timed_phy_metadata)*1.02, label = label), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F
); annotated_third_wave_timed_phy_plot

#Cambodian sequence info
annotated_third_wave_timed_phy_plot$data %>% 
filter(Country == "Cambodia") %>% arrange(x)

#sister clade info
annotated_third_wave_timed_phy_plot$data %>% 
filter(label %in% c("EPI_ISL_5142843", "EPI_ISL_1497603", "EPI_ISL_5142771", "EPI_ISL_5142715", "EPI_ISL_5142722", "EPI_ISL_5142724", "EPI_ISL_860080", "EPI_ISL_5142823", "EPI_ISL_833375", "EPI_ISL_5142766", "EPI_ISL_5142694", "OL629465")) %>% 
arrange(x)

#count sequences in the cluster by Region
annotated_third_wave_timed_phy_plot$data %>% filter(label %in% 
tree_subset(
	tree = third_wave_timed_phy,
	node = MRCA(annotated_third_wave_timed_phy_plot, c("EPI_ISL_3014857", "EPI_ISL_5161358")),
	levels_back = 0
)$tip.label) %>% select(Region) %>% table #BMR: 1252; non-BMR: 805

#earliest sequence in the cluster
annotated_third_wave_timed_phy_plot$data %>% filter(Country == "Thailand", label %in% 
tree_subset(
	tree = third_wave_timed_phy,
	node = MRCA(annotated_third_wave_timed_phy_plot, c("EPI_ISL_3014857", "EPI_ISL_5161358")),
	levels_back = 0
)$tip.label) %>% arrange(x) %>% pull(Collection.date) %>% head(1) 

#tMRCA of B.1.1.7
annotated_third_wave_timed_phy_plot$data %>% 
filter(node == MRCA(annotated_third_wave_timed_phy_plot, c("EPI_ISL_3014857", "EPI_ISL_5161358"))) %>% 
pull(x) -> MRCA_B_1_1_7
MRCA_B_1_1_7; date_decimal(MRCA_B_1_1_7)

#make a plot zoom in to the base
x_coors <- c(2020.8, 2021.6)
y_coors <- c(0, 150)
annotated_third_wave_timed_phy_plot_zoom_in <- annotated_third_wave_timed_phy_plot + 
geom_tiplab(aes(
		subset = (Th == "non-Th")|(label %in% c("EPI_ISL_5142843", "EPI_ISL_1497603", "EPI_ISL_5142771", "EPI_ISL_5142715", "EPI_ISL_5142722", "EPI_ISL_5142724", "EPI_ISL_860080", "EPI_ISL_5142823", "EPI_ISL_833375", "EPI_ISL_5142766", "EPI_ISL_5142694", "OL629465")),
		label = Country
	),
	size = 1.75, offset = 0.01, angle = 0) + 
coord_cartesian(
	xlim = x_coors, ylim = y_coors,
	expand = FALSE, default = FALSE, clip = "on") +
scale_x_continuous(breaks = seq(2020.9, 2021.5, 0.1))+
theme(
	legend.position = "none",
	axis.text.x = element_text(angle = 90),
	panel.background = element_rect(colour = "blue", fill = "white", linewidth = 1)
	)

#add a rectangle to the zoom-in area and combine the plots together
annotated_third_wave_timed_phy_plot <- ggdraw() +
draw_plot(
	annotated_third_wave_timed_phy_plot + 
	geom_rect(aes(xmin = x_coors[1], xmax = x_coors[2], ymin = y_coors[1], ymax = y_coors[2]),
		colour = "blue", fill = NA, linewidth = 0.5
	)
) +
draw_plot(annotated_third_wave_timed_phy_plot_zoom_in, x = 0.01, y = 0.13, width = 0.3, height = .86)

annotated_third_wave_timed_phy_plot

#save figure to file
ggsave(FigS4_file_png, plot = annotated_third_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)
ggsave(FigS4_file_svg, plot = annotated_third_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)

##########
#Figure S5 - 4th wave delta tree: "2021-07-06" - "2022-01-05"
##########
#subset Delta tree
forth_wave_timed_phy <- tree_subset(
	tree = timed_phy,
	node = MRCA(timed_phy, 
		c(timed_phy_plot$data %>% filter(Pango.lineage == "AY.30") %>% pull(label) %>% first(), timed_phy_plot$data %>% filter(Pango.lineage == "AY.85") %>% pull(label) %>% nth(2))),
	levels_back = 0, root_edge = FALSE
)

#make metadata df
forth_wave_timed_phy_metadata <- timed_phy_metadata %>% filter(Accession.ID %in% forth_wave_timed_phy$tip.label)
forth_wave_timed_phy_metadata <- forth_wave_timed_phy_metadata %>% 
mutate(Pango.lineage = ifelse(Pango.lineage %in% c("B.1", "B.1.1", "BA.1", "BA.1.1", "BA.2.10"), "Others", Pango.lineage))

#remove long terminal branches
long_terminal_delta_tips <- timed_phy_plot$data %>% 
filter(isTip, label %in% forth_wave_timed_phy$tip.label) %>%
filter(x > 2022.25) %>% pull(label)

forth_wave_timed_phy <- drop.tip(forth_wave_timed_phy, long_terminal_delta_tips)
forth_wave_timed_phy_metadata <- forth_wave_timed_phy_metadata %>% filter(!(Accession.ID %in% long_terminal_delta_tips))

#plot
forth_wave_timed_phy_plot <- phy_plot_fun(
	phy = forth_wave_timed_phy, 
	phy_metadata = forth_wave_timed_phy_metadata, 
	mrsd = date_decimal(timed_phy_plot$data %>% filter(label %in% forth_wave_timed_phy$tip.label) %>% pull(x) %>% max),

	PANGO_lineage_col = Lineage_col,
	Geo_region_col = Region_col,

	x_ticks = seq(2020.7, 2022.2, 0.1), x_lim = c(2020.6, 2022.52),
	y_ticks = NULL, y_lim = c(0, nrow(forth_wave_timed_phy_metadata)*1.08),

	gheatmap.font.size = 2,
	PANGO_lineage.legend.nrow = 2,
	Geo_region.legend.nrow = 2
	)

annotated_forth_wave_timed_phy_plot <- forth_wave_timed_phy_plot + 
geom_segment(
	data = Wave_label[4,], 
	mapping = aes(x = x, xend = x, y = 0, yend = nrow(forth_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[4,], 
	mapping = aes(x = xend, xend = xend, y = 0, yend = nrow(forth_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[4,], 
	mapping = aes(x = x, xend = xend, y = nrow(forth_wave_timed_phy_metadata)*1.02, yend = nrow(forth_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.5
) +
geom_text(
	data = Wave_label[4,], 
	aes(x = (x + xend)/2, y = nrow(forth_wave_timed_phy_metadata)*1.02, label = label), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F
); annotated_forth_wave_timed_phy_plot

#count sequences by Region
forth_wave_timed_phy_metadata %>% 
filter(Pango.lineage %in% c("AY.30", "AY.85") & Country == "Thailand" & Region != "Unknown") %>% 
select(Region, Pango.lineage) %>% table

#test if the proportions of BMR VS non-BMR differ significantly between AY.30 and AY.85
forth_wave_timed_phy_metadata %>% 
filter(Pango.lineage %in% c("AY.30", "AY.85") & Country == "Thailand" & Region != "Unknown") %>% 
select(Pango.lineage, Region) %>%
mutate(Pango.lineage = factor(Pango.lineage, levels = c("AY.85", "AY.30")),
Region = factor(Region, levels = c("non-BMR", "BMR"))) %>%
glm(formula = Pango.lineage ~ Region, family = "binomial", data = .) -> mdl
car::Anova(mdl, type = 3) 
summary(mdl)
lsmeans::lsmeans(mdl, pairwise ~ Region, type = "response")

#tMRCA of AY.30
forth_wave_timed_phy_plot$data %>% filter(node == MRCA(forth_wave_timed_phy, c("EPI_ISL_6036426", "EPI_ISL_3407839"))) %>% 
pull(x) -> MRCA_AY_30
MRCA_AY_30; date_decimal(MRCA_AY_30)

#earliest AY.30 sequence in Thailand
forth_wave_timed_phy_plot$data %>% filter(Country == "Thailand", Pango.lineage == "AY.30") %>% arrange(x) %>% pull(Collection.date) %>% head(1)

#tMRCA of AY.85
forth_wave_timed_phy_plot$data %>% filter(node == MRCA(forth_wave_timed_phy, c("EPI_ISL_6036486", "EPI_ISL_7199993"))) %>% 
pull(x) -> MRCA_AY_85
MRCA_AY_85; date_decimal(MRCA_AY_85)

#earliest AY.85 sequence in Thailand
forth_wave_timed_phy_plot$data %>% filter(Country == "Thailand", Pango.lineage == "AY.85") %>% arrange(x) %>% pull(Collection.date) %>% head(1)

#make a plot zoom in to the base of the AY.30 cluster
x_coors <- c(2020.9, 2021.9)
y_coors <- c(710, 840)
annotated_forth_wave_timed_phy_plot_zoom_in <- annotated_forth_wave_timed_phy_plot + 
geom_tiplab(aes(
		subset = (Th == "non-Th"),
		label = Country
	),
	size = 1.75, offset = 0.01, angle = 0) + 
coord_cartesian(
	xlim = x_coors, ylim = y_coors,
	expand = FALSE, default = FALSE, clip = "on") +
scale_x_continuous(breaks = seq(x_coors[1], x_coors[2], 0.1))+
theme(
	legend.position = "none",
	axis.text.x = element_text(angle = 90),
	panel.background = element_rect(colour = "blue", fill = "white", linewidth = 1)
	)
annotated_forth_wave_timed_phy_plot_zoom_in

#add a rectangle to the zoom-in area and combine the plots together
annotated_forth_wave_timed_phy_plot <- 
ggdraw() +
draw_plot(
	annotated_forth_wave_timed_phy_plot +
	geom_rect(aes(xmin = x_coors[1], xmax = x_coors[2], ymin = y_coors[1], ymax = y_coors[2]),
		colour = "blue", fill = NA, linewidth = 0.5
	)
) +
draw_plot(annotated_forth_wave_timed_phy_plot_zoom_in, x = 0.03, y = 0.37, width = 0.2, height = .62)
annotated_forth_wave_timed_phy_plot

#save figure to file
ggsave(FigS5_file_png, plot = annotated_forth_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)
ggsave(FigS5_file_svg, plot = annotated_forth_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)

##########
#Figure S6 - 5th wave omicron tree: "2022-01-06" - "2022-10-10"
##########
#subset Omicron tree
fifth_wave_timed_phy <- tree_subset(
	tree = timed_phy,
	node = MRCA(timed_phy, 
		c(timed_phy_plot$data %>% filter(Pango.lineage == "BA.1.17.2") %>% pull(label) %>% head(1),
		timed_phy_plot$data %>% filter(Pango.lineage == "BA.5.2.1") %>% pull(label) %>% head(1))),
	levels_back = 0, root_edge = FALSE
)

#make metadata df
fifth_wave_timed_phy_metadata <- timed_phy_metadata %>% filter(Accession.ID %in% fifth_wave_timed_phy$tip.label)

#plot
fifth_wave_timed_phy_plot <- phy_plot_fun(
	phy = fifth_wave_timed_phy, 
	phy_metadata = fifth_wave_timed_phy_metadata, 
	mrsd = date_decimal(timed_phy_plot$data %>% filter(label %in% fifth_wave_timed_phy$tip.label) %>% pull(x) %>% max),

	PANGO_lineage_col = Lineage_col,
	Geo_region_col = Region_col,

	x_ticks = seq(2021.3, 2022.8, 0.1), x_lim = c(2021.2, 2023.15),
	y_ticks = NULL, y_lim = c(0, nrow(fifth_wave_timed_phy_metadata)*1.08),

	gheatmap.font.size = 2,
	PANGO_lineage.legend.nrow = 4,
	Geo_region.legend.nrow = 3
	)

annotated_fifth_wave_timed_phy_plot <- fifth_wave_timed_phy_plot + 
geom_segment(
	data = Wave_label[5,], 
	mapping = aes(x = x, xend = x, y = 0, yend = nrow(fifth_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[5,], 
	mapping = aes(x = xend, xend = xend, y = 0, yend = nrow(fifth_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.1, linetype = "dashed"
) +
geom_segment(
	data = Wave_label[5,], 
	mapping = aes(x = x, xend = xend, y = nrow(fifth_wave_timed_phy_metadata)*1.02, yend = nrow(fifth_wave_timed_phy_metadata)*1.02),
	colour = "black", linewidth = 0.5
) +
geom_text(
	data = Wave_label[5,], 
	aes(x = (x + xend)/2, y = nrow(fifth_wave_timed_phy_metadata)*1.02, label = label), 
	hjust = 0.5, vjust = -0.25, size = 2, show.legend = F
); annotated_fifth_wave_timed_phy_plot

#tMRCA of BA.1
fifth_wave_timed_phy_plot$data %>% filter(node == MRCA(fifth_wave_timed_phy, c("EPI_ISL_10762036", "EPI_ISL_7747672"))) %>% 
pull(x) -> MRCA_BA_1
MRCA_BA_1; date_decimal(MRCA_BA_1)

#tMRCA of BA.2
fifth_wave_timed_phy_plot$data %>% filter(node == MRCA(fifth_wave_timed_phy, c("EPI_ISL_14669038", "EPI_ISL_11211315"))) %>% 
pull(x) -> MRCA_BA_2
MRCA_BA_2; date_decimal(MRCA_BA_2)

#tMRCA of BA.4
fifth_wave_timed_phy_plot$data %>% filter(node == MRCA(fifth_wave_timed_phy, c("EPI_ISL_13766253", "EPI_ISL_13276608"))) %>% 
pull(x) -> MRCA_BA_4
MRCA_BA_4; date_decimal(MRCA_BA_4)

#tMRCA of BA.5
fifth_wave_timed_phy_plot$data %>% filter(node == MRCA(fifth_wave_timed_phy, c("EPI_ISL_14583334", "EPI_ISL_13853856"))) %>% 
pull(x) -> MRCA_BA_5
MRCA_BA_5; date_decimal(MRCA_BA_5)

#save figure to file
ggsave(FigS6_file_png, plot = annotated_fifth_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)
ggsave(FigS6_file_svg, plot = annotated_fifth_wave_timed_phy_plot, width = 16, height = 22, units = "cm", dpi = 300)
