##########
#Load libraries
##########
library(dplyr)
library(lubridate) #ymd

library(sf) #st_read

library(ggplot2)
library(ggpmisc) #stat_poly_eq
library(ggforce) #facet_wrap_paginate
library(cowplot) #plot_grid
library(wesanderson) #wes_palette

##########
#Define paths
##########
#to working dir
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

#to data files
Supplementary_Data_Fig_1a_file <- "data/Sup-Dat-Fig-1a_total-case-count-during-each-wave-by-province.txt"
Supplementary_Data_Fig_1b_file <- "data/Sup-Dat-Fig-1b_weekly-case-count-by-region.txt"
Supplementary_Data_Fig_1c_file <- "data/Sup-Dat-Fig-1c_weekly-sequence-count-by-region.txt"

Th_map_file <- "data/Th_map/tha_admbnda_adm1_rtsd_20220121.shp"

#to output files
Fig1_file_png <- "out/Fig1/Fig1.raw.png"
Fig1_file_svg <- "out/Fig1/Fig1.raw.svg"

##########
#Define some "useful" variables
##########
#provices in the BMR
BMR_provinces <- c("Bangkok", "Nakhon Pathom", "Pathum Thani", "Nonthaburi", "Samut Prakan", "Samut Sakhon")

#province colours
Provinces_col <- c(
	"Bangkok" = palette.colors(palette = "Set 2")[2], 
	"Nakhon Pathom" = palette.colors(palette = "Set 2")[3], 
	"Pathum Thani" = palette.colors(palette = "Set 2")[4], 
	"Nonthaburi" = palette.colors(palette = "Set 2")[5], 
	"Samut Prakan" = palette.colors(palette = "Set 2")[6], 
	"Samut Sakhon" = palette.colors(palette = "Set 2")[7],
	"non-BMR" = "#bebebe", #grey
	"Unknown" = "#5A5A5A"
	)

#a dummy df of min and max epiweek
Epiweek_dummy_count <- rbind(
	c(Province = NA, Epiyear = 2020, Epiweek = 1, n = 0, Region = "Unknown"),
	c(Province = NA, Epiyear = 2020, Epiweek = 53, n = 0, Region = "Unknown"),
	c(Province = NA, Epiyear = 2021, Epiweek = 1, n = 0, Region = "Unknown"),
	c(Province = NA, Epiyear = 2021, Epiweek = 52, n = 0, Region = "Unknown"),
	c(Province = NA, Epiyear = 2022, Epiweek = 1, n = 0, Region = "Unknown"),
	c(Province = NA, Epiyear = 2022, Epiweek = 41, n = 0, Region = "Unknown")
	) %>% as.data.frame %>% 
	mutate(
		Epiyear = as.numeric(Epiyear), 
		Epiweek = as.numeric(Epiweek), 
		n = as.numeric(n)
	)

##########
##########
#Plot Figure 1 - time series of weekly new-case and sequence numbers
##########
##########
#load data
Wave_case_count <- read.table(Supplementary_Data_Fig_1a_file, header = TRUE, sep = "\t")
Epiweek_case_count <- read.table(Supplementary_Data_Fig_1b_file, header = TRUE, sep = "\t") %>% 
	mutate(Region = factor(Region, c("Unknown", "non-BMR", BMR_provinces)))
Epiweek_seq_count <- read.table(Supplementary_Data_Fig_1c_file, header = TRUE, sep = "\t") %>% 
	mutate(Region = factor(Region, c("Unknown", "non-BMR", BMR_provinces)))

Wave_label <- rbind(
	c(Epiyear = 2020, x = 10, xend = 23, y = 25000, yend = 25000, label = "1st wave\n(01/03/2020-31/05/2020)"),
	c(Epiyear = 2020, x = 51+5/7, xend = 54, y = 25000, yend = 25000, label = ""),
	c(Epiyear = 2021, x = 1, xend = 9, y = 25000, yend = 25000, label = "2nd wave\n(18/12/2020-28/02/2021)"),
	c(Epiyear = 2021, x = 13+4/7, xend = 27+1/7, y = 60000, yend = 60000, label = "3rd wave\n(01/04/2021-05/07/2021)"),
	c(Epiyear = 2021, x = 27+2/7, xend = 53, y = 90000, yend = 90000, label = "4th wave\n(06/07/2021-05/01/2022)"),
	c(Epiyear = 2022, x = 1, xend = 1+3/7, y = 90000, yend = 90000, label = ""),
	c(Epiyear = 2022, x = 1+4/7, xend = 41, y = 125000, yend = 125000, label = "5th wave\n(06/01/2022-10/10/2022)")
	) %>% as.data.frame %>% 
	mutate(Epiyear = as.numeric(Epiyear), x = as.numeric(x), xend = as.numeric(xend), y = as.numeric(y), yend = as.numeric(yend))

##########
#Plot Figure 1a - case count per 100,000 maps during each wave
##########
#load Thailand map
Th_map <- st_read(dsn = Th_map_file) %>% as.data.frame %>% select(geometry, ADM1_EN) 

#plot Thailand maps with case counts during the 5 COVID-19 waves
Th_map_with_case_count <- ggplot(
	merge(x = Wave_case_count, 
		y = Th_map, 
		by.x = "Province", 
		by.y = "ADM1_EN", 
		all = TRUE) %>% st_as_sf(na.fail = FALSE) 
	) +
	geom_sf(aes(fill = n_per_100000), lwd = 0.1) + 
	scale_fill_gradientn(name = "Total number\nof cases during\neach wave per\n100,000", colours = wes_palette("Zissou1", 100, type = "continuous"), trans = "log10") +

	facet_grid(
		. ~ Wave, 
		labeller = labeller(Wave = c(	"1st wave" = "1st wave\nBMR: 1,887\nnon-BMR: 1,148\nUnknown: 4",
				"2nd wave" = "2nd wave\nBMR: 18,865\nnon-BMR: 2,696\nUnknown: 109",
				"3rd wave" = "3rd wave\nBMR: 168,044\nnon-BMR: 92,185\nUnknown: 141",
				"4th wave" = "4th wave\nBMR: 679,138\nnon-BMR: 1,267,354\nUnknown: 3,750",
				"5th wave" = "5th wave\nBMR: 838,108\nnon-BMR: 1,609,698\nUnknown: 0")
			)
	) + 

	theme_void() + 
	theme(
		legend.position = "right",
		legend.title = element_text(size = 8),
		legend.text = element_text(size = 6),

		strip.text.x = element_text(size = 8)
	) 

wave_count <- Wave_case_count %>% pull(Wave) %>% unique %>% length

#zoom-in plots of the BMR regions
BMR_plots <- lapply(seq_len(wave_count), function(i) {
	Th_map_with_case_count + 
	facet_wrap_paginate( ~ Wave, nrow = 1, ncol = 1, page = i) +
	coord_sf(xlim = c(99.7, 101), ylim = c(13.4, 14.3), expand = FALSE) +
	guides(colour = "none", x = "none", y = "none") +
	theme(
		legend.position = "none",
		strip.background = element_blank(),
		strip.text.x = element_blank(),

		axis.title = element_blank(),

		plot.background = element_blank()
		)
	})

BMR_plots <- tibble(
	x = rep(100, wave_count),
	y = rep(0, wave_count),
	plot = BMR_plots,
	Wave = unique(Wave_case_count$Wave))

#combine the main Thailand maps with the BMR zoom-in plots 
Th_map_with_case_count <- Th_map_with_case_count + 
	geom_plot_npc(
		data = BMR_plots, 
		aes(npcx = x, npcy = y, label = plot, vp.width = 0.55, vp.height = 0.55)
	) +
	annotate(
		geom = "rect", 
		xmin = 99.7, xmax = 101, ymin = 13.4, ymax = 14.3,
		fill = NA, colour = "black", size = 0.5
	) 

##########
#Plot Figure 1b - Weekly new cases 
##########
#Main plot
Epiweek_case_count_plot <- 
	rbind(Epiweek_case_count, Epiweek_dummy_count) %>%
	group_by(Region, Epiyear, Epiweek) %>%
	summarise(n = sum(n)) %>% ungroup() %>% 
	ggplot() +
	geom_col(
		mapping = aes(x = Epiweek, y = n, fill = Region), #alternatively, x = interaction(Epiyear, Epiweek, lex.order = TRUE)
		position = "stack", col = "white", lwd = 0.1
	) +
	scale_fill_manual(breaks = names(Provinces_col), values = Provinces_col) +

	#add horizontal lines to indicate COVID-19 waves
	geom_segment(
		data = Wave_label, 
		aes(x = x, xend = xend, y = y, yend = yend),
		colour = "black", lwd = 0.5
	) +
	geom_text(
		data = Wave_label, 
		aes(x = (x + xend)/2, y = y, label = label), 
		hjust = 0.5, vjust = -0.25, size = 2
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
	scale_y_continuous(
		name = "Weekly cases",
		expand = c(0, 0)
	) +

	theme_bw() +
	theme(
		legend.justification = c(0, 1), legend.position = c(0.01, 0.99),
		legend.background = element_blank(),
		legend.key.size = unit(0.25, 'cm'),
		legend.title = element_text(size = 8),
		legend.text = element_text(size = 6),

		axis.title = element_text(size = 8),
		axis.text = element_text(size = 6),
		axis.text.y = element_text(size = 6, vjust = 0),

		panel.spacing.x = unit(0,"line"),

		strip.placement = 'outside',
		strip.text.x = element_text(size = 6),
		strip.background.x = element_blank()
	) + guides(fill = guide_legend(title = "Province", ncol = 1)) + coord_cartesian(clip = "off")

#zoom-in to the 1st wave, Epiyear = 2020, 10<= Epiweek <= 23
Epiweek_case_count_1st_wave_plot <- 
	rbind(Epiweek_case_count, Epiweek_dummy_count) %>%
	group_by(Region, Epiyear, Epiweek) %>%
	summarise(n = sum(n)) %>% ungroup() %>% 
	filter(Epiyear == 2020, Epiweek >= 9, Epiweek <= 25) %>%
	ggplot() +
	geom_col(
		mapping = aes(x = Epiweek, y = n, fill = Region), #alternatively, x = interaction(Epiyear, Epiweek, lex.order = TRUE)
		position = "stack", col = "white", lwd = 0.1
	) +
	scale_fill_manual(breaks = names(Provinces_col), values = Provinces_col) +

	scale_x_continuous(
		name = NULL, #"Epidemiological week"
		expand = c(0, 0), 
		position = "top",
		labels = function(x) sprintf("w.%s", x),
	) +
	scale_y_continuous(
		name = NULL, #"Weekly cases"
		expand = c(0, 0)
	) +
	ggtitle("1st wave") +

	theme_bw() +
	theme(
		legend.position = "none",

		plot.background = element_rect(fill = "transparent", colour = NA),

		plot.title = element_text(size = 6),
		axis.title = element_text(size = 8),
		axis.text = element_text(size = 6),
		axis.text.y = element_text(size = 6, vjust = 0),
	)

#zoom-in to the 2nd wave, Epiyear = 2020, 51.71429 <= Epiweek + Epiyear = 2021, 1 <= Epiweek <= 9
Epiweek_case_count_2nd_wave_plot <- 
	rbind(Epiweek_case_count, Epiweek_dummy_count) %>%
	group_by(Region, Epiyear, Epiweek) %>%
	summarise(n = sum(n)) %>% ungroup() %>% 
	filter((Epiyear == 2020 & Epiweek >= 51)|(Epiyear == 2021 & Epiweek <= 10)) %>%
	ggplot() +
	geom_col(
		mapping = aes(x = Epiweek, y = n, fill = Region), #alternatively, x = interaction(Epiyear, Epiweek, lex.order = TRUE)
		position = "stack", col = "white", lwd = 0.1
	) +
	scale_fill_manual(breaks = names(Provinces_col), values = Provinces_col) +

	facet_grid(.~ Epiyear, 
		space = "free_x", 
		scales = "free_x", 
		labeller = labeller(Epiyear = c("2020" = "Year 2020", "2021" = "Year 2021", "2022" = "Year 2022"))
	) +
	scale_x_continuous(
		name = NULL, #"Epidemiological week"
		breaks = c(5, 10, 53),
		expand = c(0, 0), 
		position = "top",
		labels = function(x) sprintf("w.%s", x),
	) +
	scale_y_continuous(
		name = NULL, #"Weekly cases"
		expand = c(0, 0)
	) +
	ggtitle("2nd wave") +

	theme_bw() +
	theme(
		legend.position = "none",

		plot.background = element_rect(fill = "transparent", colour = NA),

		plot.title = element_text(size = 6),
		axis.title = element_text(size = 8),
		axis.text = element_text(size = 6),
		axis.text.y = element_text(size = 6, vjust = 0),

		panel.spacing.x = unit(0, "line"),

		strip.placement = 'outside',
		strip.text.x = element_blank(),
		strip.background.x = element_blank()
	)

#cannot combine the 3 plots at this stage as we want to align Figure 1b and 1c.

##########
#Plot Figure 1c - Weekly new sequences
##########
#main plot
Epiweek_seq_count_plot <- 
	rbind(Epiweek_seq_count, Epiweek_dummy_count) %>%
	group_by(Region, Epiyear, Epiweek) %>%
	summarise(n = sum(n)) %>% ungroup() %>% 
	ggplot() +
	geom_col(
		mapping = aes(x = Epiweek, y = n, fill = Region), #alternatively, x = interaction(Epiyear, Epiweek, lex.order = TRUE)
		position = "stack", col = "white", lwd = 0.1
	) +
	scale_fill_manual(values = Provinces_col) +
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
	scale_y_reverse(
		name = "Weekly sequences",
		expand = c(0, 0)
	) +
	theme_bw() +
	theme(
		legend.position = "none", 

		axis.title = element_text(size = 8),
		axis.title.x = element_blank(),
		axis.text = element_text(size = 6),
		axis.text.y = element_text(size = 6, vjust = 1),


		plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
		panel.spacing.x = unit(0,"line"),

		strip.background = element_blank(),
		strip.text.x = element_blank()
	)

#plot of weekly seq ~ weekly cases, through time
Seq_vs_case_num_dat <- merge(
	Epiweek_case_count %>% 
		filter(Region %in% BMR_provinces) %>% 
		group_by(Epiyear, Epiweek) %>% 
		summarise(n = sum(n)) %>% ungroup(), 
	Epiweek_seq_count %>% 
		filter(Region %in% BMR_provinces) %>% 
		group_by(Epiyear, Epiweek) %>%
		summarise(n = sum(n)) %>% ungroup(), 
	by = c("Epiyear", "Epiweek"), all = TRUE
	) %>% 
	rename(Epiweek_case_count = n.x, Epiweek_seq_count = n.y) %>%
	mutate(
		seq_rate = Epiweek_seq_count/Epiweek_case_count*100,
		rough_date = ymd( sprintf("%d-01-01", Epiyear) ) + weeks( Epiweek - 1 )
	)

Seq_vs_case_num_plot_by_time <- Seq_vs_case_num_dat %>% 
ggplot(aes(x = Epiweek_case_count, y = Epiweek_seq_count, col = rough_date)) + 
	geom_point(size = 1) + 
	geom_smooth(method = "lm", formula = y~x-1) +
	stat_poly_eq(use_label(c("eq")), formula = y ~ x + 0, label.y = "top", label.x = "right", size = 2) +
	stat_correlation(method = "spearman", label.y = 0.8, label.x = "right", size = 2) +
	scale_x_continuous(name = "Cases") +
	scale_y_continuous(name = "Sequences") +
	scale_color_date(name = "Epiweek") +
	
	theme(
		axis.text = element_text(size = 6),
		axis.title = element_text(size = 6),
		
		plot.background = element_rect(fill = "transparent", colour = NA),

		legend.position = "right",
		legend.background = element_blank(),
		legend.margin = margin(c(0, 0, 0, 0)),
		legend.key.size = unit(0.25, 'cm'),
		legend.title = element_text(size = 6),
		legend.text = element_text(size = 6))

##########
#Combine all of the plots
##########
Fig1 <- plot_grid(
	Th_map_with_case_count,
	ggdraw() +
		draw_plot(
			plot_grid(
				Epiweek_case_count_plot, 
				Epiweek_seq_count_plot, 
				ncol = 1, align = "v", rel_heights = c(1.5, 1), labels = c('b', 'c')
			)
		) +
		draw_plot(Seq_vs_case_num_plot_by_time, x = 0.10, y = 0.02, width = 0.315, height = .25) + 
		draw_plot(Epiweek_case_count_1st_wave_plot, x = 0.22, y = 0.5, width = 0.15, height = 0.25) +
		draw_plot(Epiweek_case_count_2nd_wave_plot, x = 0.41, y = 0.6, width = 0.15, height = 0.27),
	ncol = 1, rel_heights = c(1, 1.5), labels = c('a', '')) +
	theme(plot.background = element_rect(fill = "white", colour = NA))

##########
#Save the plot to file
##########
ggsave(Fig1_file_png, plot = Fig1, width = 16, height = 16, units = "cm", dpi = 300)
ggsave(Fig1_file_svg, plot = Fig1, width = 16, height = 16, units = "cm", dpi = 300)
