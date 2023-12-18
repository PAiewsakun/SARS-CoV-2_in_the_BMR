##########
#Load libraries
##########
library(dplyr)

library(ggplot2)
library(cowplot)

##########
#Define paths
##########
#to working dir
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

#to data files
Supplementary_Data_Fig_2ab_file <- "data/Sup-Dat-Fig-2ab_weekly-sequence-count-by-lineage.txt"

#to output files 
Fig2_file_png <- "out/Fig2/Fig2.raw.png"
Fig2_file_svg <- "out/Fig2/Fig2.raw.svg"

##########
#Define some "useful" variables
##########
#lineage colours
Lineage_col <- c(
	"A.6" = "#30ACEC",
	"B.1" = "#3a8ec5",
	"B.1.1" = "#44709d",

	"B.1.36.16" = "#80C34F",

	"B.1.1.7" = "#E29D3E",

	"AY.30" = "#D64A3B",
	"AY.85" = "#b96749",
	"B.1.617.2" = "#9b8357",

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

#a dummy df of min and max epiweek
Epiweek_dummy_count <- rbind(
	c(Pango.lineage = "Unassigned", Epiyear = 2020, Epiweek = 1, n = 0),
	c(Pango.lineage = "Unassigned", Epiyear = 2020, Epiweek = 53, n = 0),
	c(Pango.lineage = "Unassigned", Epiyear = 2021, Epiweek = 1, n = 0),
	c(Pango.lineage = "Unassigned", Epiyear = 2021, Epiweek = 52, n = 0),
	c(Pango.lineage = "Unassigned", Epiyear = 2022, Epiweek = 1, n = 0),
	c(Pango.lineage = "Unassigned", Epiyear = 2022, Epiweek = 41, n = 0)
	) %>% as.data.frame %>% 
	mutate(
		Epiyear = as.numeric(Epiyear), 
		Epiweek = as.numeric(Epiweek), 
		n = as.numeric(n)
	)

##########
##########
#Plot Figure 2 - sequence lineage diversity, BMR only
##########
##########
#load data
##"prodominant" lineages = having more than 30 sequences and occupying more than 5% of the sequencing data at some points in time
##N.epiweek = estimated number of sequences within the epi week
##N.Pango.lineage = the total number of sequences in the data collection belonging to the lineage
##n.Pango.lineage = number of sequences within the epi week, by lineage
##pct.Pango.lineage = pct of the sequences of that lineage within the epi week

Epiweek_lineage_in_BMR_count <- read.table(Supplementary_Data_Fig_2ab_file, header = TRUE, sep = "\t")

Wave_label <- rbind(
	c(Epiyear = 2020, x = 10, xend = 23, y = "B.1.1.7", yend = "B.1.1.7", label = "1st wave\n(01/03/2020-31/05/2020)", Pango.lineage = NA),
	c(Epiyear = 2020, x = 51+5/7, xend = 54, y = "AY.85", yend = "AY.85", label = "", Pango.lineage = NA),
	c(Epiyear = 2021, x = 1, xend = 9, y = "AY.85", yend = "AY.85", label = "2nd wave\n(18/12/2020-28/02/2021)", Pango.lineage = NA),
	c(Epiyear = 2021, x = 13+4/7, xend = 27+1/7, y = "BA.1.1", yend = "BA.1.1", label = "3rd wave\n(01/04/2021-05/07/2021)", Pango.lineage = NA),
	c(Epiyear = 2021, x = 27+2/7, xend = 53, y = "BA.2", yend = "BA.2", label = "4th wave\n(06/07/2021-05/01/2022)", Pango.lineage = NA),
	c(Epiyear = 2022, x = 1, xend = 1+3/7, y = "BA.2", yend = "BA.2", label = "", Pango.lineage = NA),
	c(Epiyear = 2022, x = 1+4/7, xend = 41, y = "B.1.1.7", yend = "B.1.1.7", label = "5th wave\n(06/01/2022-10/10/2022)", Pango.lineage = NA)
	) %>% as.data.frame %>% 
	mutate(Epiyear = as.numeric(Epiyear), x = as.numeric(x), xend = as.numeric(xend), y = as.factor(y), yend = as.factor(yend))

##########
#Plot Figure 2a - numbers of sequences by lineage, predominant lineages only
##########
Epiweek_lineage_in_BMR_bead_plot <- Epiweek_lineage_in_BMR_count %>% 
	mutate(
		Pango.lineage = ifelse(
			is.na(Predominant),
			"Unassigned",
			ifelse(
				Predominant,
				Pango.lineage,
				"Others"
			)
		)
	) %>% #label "non-predominant" lineages as "Others", and unassigned sequences as "Unassigned"
	group_by(Pango.lineage, Epiyear, Epiweek) %>% #combine counts of non-predominant lineages together
	summarise(
		n = sum(n.Pango.lineage),
	) %>% ungroup() %>% 
	mutate(n = ifelse(n >= 1, round(n), NA)) %>% #remove weeks with counts less than 1
	bind_rows(., Epiweek_dummy_count) %>% #add dummy min/max dates
	mutate(Pango.lineage = factor(Pango.lineage, rev(names(Lineage_col)))) %>% #reorder the lineage factor 
	filter(n > 0) %>% 

	ggplot(aes(x = Epiweek, y = Pango.lineage, colour = Pango.lineage))+
	geom_point(aes(size = n), shape = 1) +

	#add horizontal lines to indicate COVID-19 waves
	geom_segment(
		data = Wave_label, 
		aes(x = x, xend = xend, y = y, yend = yend),
		colour = "black", lwd = 0.5
	) +
	geom_text(
		data = Wave_label, 
		aes(x = (x + xend)/2, y = y, label = label), 
		hjust = 0.5, 
		vjust = 0.50, 
		size = 2, show.legend = F
	) +

	facet_grid(.~ Epiyear, 
		space = "free_x", 
		scales = "free_x", 
		labeller = labeller(Epiyear = c("2020" = "Year 2020", "2021" = "Year 2021", "2022" = "Year 2022"))
	) +

	scale_x_continuous(
		name = "Epidemiological week", 
		expand = c(0, 0), 
		position = "top",
		labels = function(x) sprintf("w.%s", x),
	) +

	scale_y_discrete(
		name = "PANGO lineage assignment",
	) +

	scale_color_manual(
		breaks = names(Lineage_col),
		values = Lineage_col
	) +

	scale_size_continuous(
		name = "Number of\nsequences"
	) +

	theme_bw() +
	theme(
		legend.justification = c(0, 0), legend.position = c(0, 0.1),
		legend.direction = "vertical", legend.box = "horizontal",
		legend.background = element_blank(),
		legend.key.size = unit(0.25, 'cm'),
		legend.key = element_blank(),
		legend.title = element_text(size = 6, face = "bold"),
		legend.text = element_text(size = 6),

		axis.title = element_text(size = 8),
		axis.text = element_text(size = 6),
		axis.text.y = element_text(size = 6),

		plot.margin = unit(c(0.5, 5, 0.5, 0.5), "pt"),
		panel.spacing.x = unit(0,"line"),

		strip.placement = 'outside',
		strip.text.x = element_text(size = 6),
		strip.background.x = element_blank()
	) + 
	guides(
		size = guide_legend(order = 2),
		colour = guide_legend(title = "PANGO\nlineage", ncol = 3, order = 1)
	) + 
	coord_cartesian(clip = "off")

##########
#Plot Figure 2b - sequence proportion by lineage, predominant lineages only
##########
Epiweek_lineage_in_BMR_prop_plot <- Epiweek_lineage_in_BMR_count %>% 
	mutate(
		Pango.lineage = ifelse(
			is.na(Predominant),
			"Unassigned",
			ifelse(
				Predominant,
				Pango.lineage,
				"Others"
			)
		)
	) %>% #label "non-predominant" lineages as "Others", and unassigned sequences as "Unassigned"
	group_by(Pango.lineage, Epiyear, Epiweek) %>% #combine counts of non-predominant lineages together
	summarise(
		n = sum(n.Pango.lineage),
	) %>% ungroup() %>% 
	mutate(n = ifelse(n >= 1, n, NA)) %>%
	bind_rows(., Epiweek_dummy_count) %>%
	mutate(Pango.lineage = factor(Pango.lineage, names(Lineage_col))) %>%
	group_by(Epiyear, Epiweek) %>% mutate(N = sum(n, na.rm = T), prop = n/N) %>%

	ggplot(aes(x = Epiweek, y = n, fill = Pango.lineage, label = Pango.lineage)) +
	geom_col(
		position = "fill", col = "white", lwd = 0.1
	) +
	facet_grid(.~ Epiyear, 
		space = "free_x", 
		scales = "free_x", 
	) +
	scale_x_continuous(
		name = "", 
		expand = c(0, 0), 
		position = "top",
		labels = function(x) sprintf("w.%s", x),
	) +
	scale_y_continuous(
		name = "Sequence\nproportion",
		expand = c(0, 0)
	) +

	scale_fill_manual(values = Lineage_col) +
	theme_bw() +
	theme(
		legend.position = "none",

		axis.title = element_text(size = 8),
		axis.title.x = element_blank(),
		axis.text = element_text(size = 6),
		axis.text.y = element_text(size = 6, vjust = 0),

		plot.margin = unit(c(0.5, 5, 0.5, 0.5), "pt"),
		panel.spacing.x = unit(0,"line"),

		strip.background = element_blank(),
		strip.text.x = element_blank()
	)

##########
#Combine all plots together
##########
Fig2 <- plot_grid(
			Epiweek_lineage_in_BMR_bead_plot,
			Epiweek_lineage_in_BMR_prop_plot,
			ncol = 1, align = "v", rel_heights = c(1, 0.30)
		) +
	theme(plot.background = element_rect(fill = "white", colour = NA))

##########
#Save the plot to file
##########
ggsave(Fig2_file_png, plot = Fig2, width = 16, height = 10, units = "cm", dpi = 300)
ggsave(Fig2_file_svg, plot = Fig2, width = 16, height = 10, units = "cm", dpi = 300)

