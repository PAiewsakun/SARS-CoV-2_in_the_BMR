##########
#Load libraries
##########
library(dplyr)
library(lubridate)

library(ggplot2)

library(emmeans)

##########
#Define paths
##########
#to working dir
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

#to data files
Ne_file <- "data/Sup-Dat-Fig-3_BMR_Ne.txt"
Weekly_case_count_file <- "data/Sup-Dat-Fig-1b_weekly-case-count-by-region.txt"

#to output files 
Fig4_file_png <- "out/Fig4/Fig4.raw.png"
Fig4_file_svg <- "out/Fig4/Fig4.raw.svg"

##########
#Define some "useful" variables
##########
#COVID wave label
Wave_label <- rbind(
	c(x = "2020-03-01", xend = "2020-05-31", wave = "1st wave"),
	c(x = "2020-12-18", xend = "2021-02-28", wave = "2nd wave"),
	c(x = "2021-04-01", xend = "2021-07-05", wave = "3rd wave"),
	c(x = "2021-07-06", xend = "2022-01-05", wave = "4th wave"),
	c(x = "2022-01-06", xend = "2022-10-10", wave = "5th wave")
	) %>% as.data.frame %>% 
	mutate(
		x = as.numeric(decimal_date(as.Date(x))),
		xend = as.numeric(decimal_date(as.Date(xend)))
	)

wave_label_from_date_fun <- function(t) {
	wave <- Wave_label$wave[(Wave_label$x < t) * (t < Wave_label$xend) == 1]
	wave <- ifelse(identical(wave, character(0)), "Between waves", wave)
	return(wave)
}

##########
#Load data
##########
#BMR Ne estimate
Ne_est <- read.table(Ne_file, header = TRUE, sep = "\t", quote = "")

#BMR weekly case counts
Weekly_case_count <- read.table(Weekly_case_count_file, header = TRUE, sep = "\t", quote = "") %>%
	filter(!(Region %in% c("non-BMR", "Unknown"))) %>% 
	group_by(Epiyear, Epiweek) %>% summarise(Weekly_case_count = sum(n)) %>% ungroup() 

#Join the two df together
Ne_VS_weekly_case <- inner_join(
	x = Weekly_case_count,
	y = Ne_est %>%
		mutate(date = date(date_decimal(t)), Epiyear = year(date), Epiweek = week(date)) %>% 
		select(Epiyear, Epiweek, t, diff_ne_median),
	by = c("Epiyear", "Epiweek")) %>% 
	arrange(Epiyear, Epiweek) %>% 
	group_by(t) %>% mutate(wave = wave_label_from_date_fun(t)) %>% ungroup() %>% select(-t)

##########
#regression analyses between "Weekly_case_count" and "diff_ne_median"
##########
get_mean_CI_fun <- function(pred){
	fit <- pred$fit
	se <- pred$se.fit
	lb <- fit - se*qnorm(.975)
	ub <- fit + se*qnorm(.975)

	return (data.frame(fit = fit, lb = lb, ub = ub))
}
#overall
overall_mdl <- lm(Weekly_case_count ~ diff_ne_median, data = Ne_VS_weekly_case)
summary(overall_mdl) #R2 = 0.7039, adj_R2 = 0.7018
car::Anova(overall_mdl, type = "3") #sig. slopes: slope = 137.532; df = 1; F = 335.1359; p = <2e-16

#prediction
overall_mdl_pred <- predict(overall_mdl, Ne_VS_weekly_case, se.fit = T)
overall_mdl_pred <- data.frame(
	diff_ne_median = Ne_VS_weekly_case %>% select(diff_ne_median),
	get_mean_CI_fun(overall_mdl_pred),
	wave = "Overall"
)

#separate lines for each wave
all_wave_mdl <- lm(Weekly_case_count ~ diff_ne_median*wave, data = Ne_VS_weekly_case)
summary(all_wave_mdl)
car::Anova(all_wave_mdl, type = "3") #sig. slope diff. among waves: df = 5; F = 2.8277; p = 0.01852

#pair-wise slope comp.
emtrends(all_wave_mdl , pairwise ~ wave, var = "diff_ne_median") #4th sig. diff from 3rd and 5th

#excluding the 4th wave data from the analysis
no_4th_wave_mdl <- lm(Weekly_case_count ~ diff_ne_median*wave, data = Ne_VS_weekly_case, subset = (wave != "4th wave"))
summary(no_4th_wave_mdl)#R2 = 0.9406, adj_R2 = 0.9356
car::Anova(no_4th_wave_mdl, type = "3") #No sig. slope diff. among waves: df = 4; F = 1.1085; p = 0.35643

#update the model -> common slopes between all waves
no_4th_wave_mdl <- update(no_4th_wave_mdl, . ~ . - diff_ne_median:wave)
summary(no_4th_wave_mdl)
car::Anova(no_4th_wave_mdl, type = "3") #sig. slope: slope = 128.597; df = 1; F = 760.1878; p = <2e-16

#prediction
no_4th_wave_mdl_pred <- predict(no_4th_wave_mdl, Ne_VS_weekly_case %>% filter(wave != "4th wave"), se.fit = T)
no_4th_wave_mdl_pred <- data.frame(
	diff_ne_median = Ne_VS_weekly_case %>% filter(wave != "4th wave") %>% select(diff_ne_median),
	get_mean_CI_fun(no_4th_wave_mdl_pred),
	wave = Ne_VS_weekly_case %>% filter(wave != "4th wave") %>% select(wave)
	)

#4th wave model
Delta_wave_mdl <- lm(Weekly_case_count ~ diff_ne_median, data = Ne_VS_weekly_case, subset = (wave == "4th wave"))
summary(Delta_wave_mdl)
car::Anova(Delta_wave_mdl, type = "3") #sig. slope: slope = 316.4; df = 1; F =  8.4850; p = 0.007625

#prediction
Delta_wave_mdl_pred <- predict(Delta_wave_mdl, Ne_VS_weekly_case %>% filter(wave == "4th wave"), se.fit = T)
Delta_wave_mdl_pred <- data.frame(
	diff_ne_median = Ne_VS_weekly_case %>% filter(wave == "4th wave") %>% select(diff_ne_median),
	get_mean_CI_fun(Delta_wave_mdl_pred ),
	wave = Ne_VS_weekly_case %>% filter(wave == "4th wave") %>% select(wave)
	)

by_wave_mdl_pred <- rbind(no_4th_wave_mdl_pred, Delta_wave_mdl_pred)

##########
#Figure 4 - Ne estimate VS weekly case count 
##########
Ne_VS_case_plot <- 
ggplot() +
geom_point(data = Ne_VS_weekly_case, aes(x = diff_ne_median, y = Weekly_case_count, col = wave), size = 0.5) +  

#overall model
geom_ribbon(data = overall_mdl_pred, aes(x = diff_ne_median, ymin = lb, ymax = ub, col = wave, fill = wave), alpha = 0.1, linewidth = 0.5) +
geom_line(data = overall_mdl_pred, aes(x = diff_ne_median, y = fit, col = wave), linewidth = 0.5) +

#by wave model
geom_ribbon(data = by_wave_mdl_pred, aes(x = diff_ne_median, ymin = lb, ymax = ub, col = wave, fill = wave), alpha = 0.1, linewidth = 0.5) +
geom_line(data = by_wave_mdl_pred, aes(x = diff_ne_median, y = fit, col = wave), linewidth = 0.5) +
scale_color_manual(
	name = "COVID-19 Wave",
	values = c(
		"1st wave" = "#30ACEC",
		"2nd wave" = "#80C34F",
		"3rd wave" = "#E29D3E",
		"4th wave" = "#D64A3B",
		"5th wave" = "#a666e1",
		"Between waves" = "#000000",
		"Overall" = "#bebebe"),
	breaks = c("1st wave", "2nd wave", "3rd wave", "4th wave", "5th wave", "Between waves", "Overall"),
	guide = guide_legend(ncol = 1)
) +
scale_fill_manual(
	name = "COVID-19 Wave",
	values = c(
		"1st wave" = "#30ACEC",
		"2nd wave" = "#80C34F",
		"3rd wave" = "#E29D3E",
		"4th wave" = "#D64A3B",
		"5th wave" = "#a666e1",
		"Between waves" = "#000000",
		"Overall" = "#bebebe"),
	breaks = c("1st wave", "2nd wave", "3rd wave", "4th wave", "5th wave", "Between waves", "Overall"),
	guide = guide_legend(ncol = 1)
) +
xlab(bquote(Ne["all sequences"] - Ne["non-BMR sequences"])) + ylab("Weekly case count in the BMR") + 
coord_cartesian(ylim = c(0, 70000), xlim = c(0, 420), expand = FALSE) +
theme_classic() + 
theme(
	legend.justification = c(0, 1), legend.position = c(0, 1),
	legend.background = element_blank(),
	legend.margin = margin(c(0, 0, 0, 0)),
	legend.key.size = unit(0.25, 'cm'),
	legend.key = element_blank(),
	legend.title = element_text(size = 8),
	legend.text = element_text(size = 6),

	axis.title = element_text(size = 8),
	axis.text = element_text(size = 6),
	axis.text.x = element_text(size = 6, vjust = 0),
	axis.text.y = element_text(size = 6, vjust = 0),

	plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
)

##########
#Save the plot to file
##########
ggsave(Fig4_file_png, plot = Ne_VS_case_plot, width = 8, height = 8, units = "cm", dpi = 300)
ggsave(Fig4_file_svg, plot = Ne_VS_case_plot, width = 8, height = 8, units = "cm", dpi = 300)

