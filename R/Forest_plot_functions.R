#############################
# New forest plot functions #
#############################

facet_var_gen <- function (dat, col_num, group_num = 1) {
	#Generates a vector that can be added to a dataframe and used to split plots up 
	#group_num is more than 1 if you have multiple data points per categorical variable
	if (nrow(dat) %% 2 != 0) {
		facet_var <- rep(
			c(rep(1:(col_num-1), each = ceiling(nrow(dat) / col_num)),
			rep(col_num, times = nrow(dat) - (ceiling(nrow(dat) / col_num) * (col_num - 1)))),
			times = group_num)
	} else {
		facet_var <- rep(rep(1:col_num, each = ceiling(nrow(dat) / col_num)), times = group_num)
	}
}


forest_plot <- function(dat, col_num, group = NA, y_axis_var, units = NULL, title = NULL) {
#Format of data for plot (doesn't matter where in the data frame these things are and can have extra columns)
# y_axis_var Estimate 2.5 % 97.5 % 
# ----------  ----     --    ---
# ----------  ----     --    --- 
# ----------  ----     --    ---

require("tidyverse")
require("gridExtra")
require("ggplot2")


group_num <- length(unique(dat[[group]]))
if (group_num == 0) {
	group_num <- 1
}

met_num <- nrow(dat)/group_num/col_num

#
	#Produce a data frame that describes how the variables will be shaded
	if ((nrow(dat) / group_num) %% 2 != 0) {
		shading <- data.frame(
				min = 
					rep(c(rep(seq(from = 0.5, to = (met_num * col_num + 1)/col_num, by = 1), times = col_num - 1),
					seq(from = 0.5, to = ((nrow(dat) - (ceiling(met_num) * group_num * (col_num - 1))) / group_num), by = 1)),
					times = group_num),
	            max = 
	            	rep(c(rep(seq(from = 1.5, to = (met_num * col_num + 1)/col_num + 1, by = 1), times = col_num - 1),
					seq(from = 1.5, to = ((nrow(dat) - (ceiling(met_num) * group_num * (col_num - 1))) / group_num) + 1, by = 1)),
					times = group_num),
	            col = rep(c(0,1), length.out = nrow(dat)),
	            facet_var = facet_var)
	} else {
		shading <- data.frame(min = rep(seq(from = 0.5, to = met_num, by = 1), times = col_num * group_num),
	           	max = rep(seq(from = 1.5, to = met_num + 1, by = 1), times = col_num * group_num),
	           	col = rep(c(0,1), length.out = nrow(dat)),
	            facet_var = facet_var)
	}

	#Minimum and maximum x values
	min_val <- floor(min(dat[["2.5 %"]]))
	max_val <- ceiling(max(dat[["97.5 %"]]))

	#Produce all the plots separately
	plots <- list()
	for (i in 1:col_num) {
		test_forest_dat <- filter(dat, facet_var == i)
		test_shading_dat <- filter(shading, facet_var == i)
		plots[[i]] <- ggplot(test_forest_dat) +
		geom_point(aes_string(x = y_axis_var, y = "Estimate", ymin = "`2.5 %`", ymax = "`97.5 %`", colour = group, group = group), position=position_dodge(width = 0.9)) +
		geom_rect(data = test_shading_dat, aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = factor(col)) ) +
		scale_fill_manual(values = c("white", "gray80")) +
		geom_pointrange(aes_string(x = y_axis_var, y = "Estimate", ymin = "`2.5 %`", ymax = "`97.5 %`", colour = group, group = group), position=position_dodge(width = 0.9)) +
		geom_hline(yintercept = 0) +
		theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none") +
		scale_y_continuous(limit = c(min_val, max_val)) +
		coord_flip()
	}

	#Produce the legend for the plot
	legend <- cowplot::get_legend(
		ggplot(dat) +
		geom_point(aes_string(x = y_axis_var, y = "Estimate", ymin = "`2.5 %`", ymax = "`97.5 %`", colour = group, group = group), position=position_dodge(width = 0.9)) +
		scale_fill_manual(values = c("white", "gray80")) +
		geom_pointrange(aes_string(x = y_axis_var, y = "Estimate", ymin = "`2.5 %`", ymax = "`97.5 %`", colour = group, group = group), position=position_dodge(width = 0.9)) +
		theme(legend.direction = "vertical", legend.justification = "center", legend.title = element_blank())
	)


	marrangeGrob(plots, right = legend, bottom = units, top = title, ncol = col_num, nrow = 1, newpage = FALSE)
}





