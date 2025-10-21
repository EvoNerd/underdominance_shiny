# underdominance shiny app dependency for app.R
#   simulation and plotting functions are defined here
#author: Ana-Hermina Ghenu
#date: 2025-10-08

library(ggplot2)   # for plotting the time series
library(patchwork) # for combining ggplots in the time series

# define colour scheme for genotypes
# use rgb to set transparency (alpha = 0.75) # note that max value for alpha is 255 so 0.75*255 ~ 191
colours_genotypes <- c(AA = rgb(123, 45, 38, maxColorValue = 255, alpha = 191),  # burnt umber
                       Aa = rgb(136, 124, 156, maxColorValue = 255, alpha = 250)) # dark lilac
colours_genotypes["aA"] <- colours_genotypes["Aa"]
colours_genotypes["aa"] <- rgb(206, 184, 75, maxColorValue = 255, alpha = 225)  # old Gold

colour_p <- "#788277"

###################################
# define simulation functions
###################################

# recursion equation: a function to get the allele frequency in the next generation
next_p <- function(w_AA, w_heterozygote, current_p, W_bar){
  output <- (current_p^2*w_AA + current_p*(1-current_p)*w_heterozygote) / W_bar
  return(unname(output))
}

# a function to get the fitness of the genotypes
fitness_genotypes <- function(sel_coef, dom_coef)
  return(c(AA = 1 + sel_coef, Aa = 1 + dom_coef, aA = 1 + dom_coef, aa = 1.0))

# a function to get mean fitness of the population
mean_fitness <- function(w_genotypes, freq){
  output <- freq^2*w_genotypes["AA"] + 2*freq*(1-freq)*w_genotypes["Aa"] + (1-freq)^2*w_genotypes["aa"]
  return(unname(output))
}

# a function to simulate genotypes over time
sim_forward_time <- function(num_gens, freq_init, w_genotypes){
  # initialize data.frame for storage of allele frequencies (p)
    #starting from initial allele frequencies at generation 0
  df_p <- data.frame(gens = 0:num_gens,
                     p = c(freq_init, rep(NA, num_gens)))
  # loop through each generation
  for(i in 2:nrow(df_p)){
    W_bar <- mean_fitness(w_genotypes, df_p$p[i-1])
    # use current mean fitness to get allele frequency in the next generation
    df_p$p[i] <- next_p(w_genotypes["AA"], w_genotypes["Aa"], df_p$p[i-1], W_bar)
  }
  
  # create a data.frame for the genotype frequencies over time
  df_genotypes <- data.frame(gens = rep(df_p$gens, times = 3),
                             genotyp = rep(c("AA", "Aa", "aa"), each = nrow(df_p)),
                             freq = c(df_p$p^2, 2*df_p$p*(1-df_p$p), (1-df_p$p)^2))
  
  # this is a hybridization event so there's only homozygotes at time=0
  HW_heterozygotes_time0 <- df_genotypes$freq[df_genotypes$gens==0 & df_genotypes$genotyp=="Aa"]
  df_genotypes$freq[df_genotypes$gens==0 & df_genotypes$genotyp=="Aa"] <- 0
  df_genotypes$freq[df_genotypes$gens==0 & df_genotypes$genotyp=="AA"] <- df_genotypes$freq[df_genotypes$gens==0 & df_genotypes$genotyp=="AA"] + HW_heterozygotes_time0/2
  df_genotypes$freq[df_genotypes$gens==0 & df_genotypes$genotyp=="aa"] <- df_genotypes$freq[df_genotypes$gens==0 & df_genotypes$genotyp=="aa"] + HW_heterozygotes_time0/2

  # return the two data.frames as a list
  return(list(allele = df_p,
              genotypes = df_genotypes))
}

# a function to get the delta p values
get_p_byp <- function(w_genotypes) {
  # Initialize a vector of allele frequencies from 0 to 1
  p_values <- seq(0, 1, length.out = 1000)
  
  # Compute mean fitness and change in allele frequency
  W_bar_values <- sapply(p_values, function(p) mean_fitness(w_genotypes, p))
  delta_p_values <- sapply(1:length(p_values), function(i) 
    next_p(w_genotypes["AA"], w_genotypes["Aa"], p_values[i], W_bar_values[i]) - p_values[i])
  
  return(data.frame(p=p_values, delta_p=delta_p_values))
}

# a function to get the non-trivial equilibrium point
get_equilibrium_pt <- function(w_genotypes)
  unname((w_genotypes["Aa"] - w_genotypes["aa"]) / (2*w_genotypes["Aa"] - w_genotypes["AA"] - w_genotypes["aa"]))

# a function indicating the long-term steady state
get_steady_state <- function(w_genotypes, deltap.df){
  # some presets for plotting the arrows
  min_arrow_length = 0.07
  arrow_horizontal_whitespace = 0.08
  
  # underdominance always results in bistability
  if(w_genotypes["Aa"] < w_genotypes["aa"]) {
    # the non-trivial equilibrium point is *UNSTABLE*
    unstable_pt <- get_equilibrium_pt(w_genotypes)
      
    # draw the arrows symmetrically about y=0 at a half of the smallest max delta_p values
    y_val = min(abs(min(deltap.df$delta_p)), abs(max(deltap.df$delta_p)))/2
    # x values for left arrow (1): 
    arrow1_length <- max(unstable_pt-2*arrow_horizontal_whitespace, min_arrow_length)
    # x values for right arrow (2): 
    arrow2_length <- max(1-unstable_pt-2*arrow_horizontal_whitespace, min_arrow_length)
    output <- list(outcome = "bistable",
                   stable_pts = c(0, 1),
                   unstable_pt = unstable_pt,
                   arrow1 = c(y0 = y_val, x0 = unstable_pt - 0.5*arrow_horizontal_whitespace,
                              y1 = y_val, x1 = unstable_pt - arrow_horizontal_whitespace - arrow1_length),
                   arrow2 = c(y0 = -y_val, x0 = unstable_pt + .5*arrow_horizontal_whitespace,
                              y1 = -y_val, x1 = unstable_pt + arrow_horizontal_whitespace + arrow2_length)
    )
    
    # additivity / partial dominance always results in A going to fixation (1 stable attractor)
  } else if(w_genotypes["Aa"] >= w_genotypes["aa"] & 
            w_genotypes["Aa"] < w_genotypes["AA"]) {
    y_val <- 0.005 # this x-value looks nice
    # increase the whitespace of the arrow
    arrow_horizontal_whitespace <- 3*arrow_horizontal_whitespace
    arrow1_length <- 1-2*arrow_horizontal_whitespace
    output <- list(outcome = "p=1",
                   stable_pts = 1,
                   arrow1 = c(y0 = y_val, x0 = arrow_horizontal_whitespace,
                              y1 = y_val, x1 = arrow_horizontal_whitespace + arrow1_length)
    )
    
    # overdominance always results in 1 stable attractor at intermediate allele freq  
  } else if(w_genotypes["Aa"] > w_genotypes["AA"]){
    # the non-trivial equilibrium point is *STABLE*
    stable_pt <- get_equilibrium_pt(w_genotypes)
  
    # draw the arrows symmetrically about y=0 at half of the smallest max delta_p values
    y_val = min(abs(min(deltap.df$delta_p)), abs(max(deltap.df$delta_p)))/2
    # x values for left arrow (1): 
    arrow1_length <- max(stable_pt-2*arrow_horizontal_whitespace, min_arrow_length)
    # x values for right arrow (2): 
    arrow2_length <- max(1-stable_pt-2*arrow_horizontal_whitespace, min_arrow_length)
    output <- list(outcome = "p intermediate",
                   stable_pts = stable_pt,
                   arrow1 = c(y0 = -y_val, x0 = stable_pt - arrow_horizontal_whitespace - arrow1_length,
                              y1 = -y_val, x1 = stable_pt - arrow_horizontal_whitespace),
                   arrow2 = c(y0 = y_val, x0 = stable_pt + arrow_horizontal_whitespace + arrow2_length,
                              y1 = y_val, x1 = stable_pt + arrow_horizontal_whitespace)
    )
  } 
  
  
  return(output)
}



###################################
# define visualization functions
###################################

# plot fitness of the genotypes
plot_schematic <- function(w_genotypes) {
  # change the graphing settings
  par(cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(2.0, 3.5, 0.01, 0.02), # decrease the borders
      mgp = c(1.0, 0.1, 0) # decrease vertical white space btw tick labels and plot
      ) 
  
  # make the barplot in graphics (i.e., base R)
  barplot(
    w_genotypes, 
    cex.lab = 1.35,
    col = colours_genotypes[names(w_genotypes)], 
    border = "black", 
    ylim = c(0, 1.6), 
    ylab = "",
    xlab = "Genotype",
    yaxt = "n"  # Suppress default y-axis tick labels
  )
  
  # Add custom y-axis tick labels
  axis(2, at = seq(from=0, to=1.5, by=0.5), labels = c("0", "", "1", ""), mgp = c(3, 1, 0))
  # add custom y-axis label
  mtext("Fitness", side = 2, line = 2.0, adj=0.35, cex = 1.35)


  # add a box around the whole plot
  box()
}

plot_dp_by_p <- function(deltap.df, steady_state) {
  # change the graphing settings
  par(mgp = c(2, 0.5, 0), # decrease the space between tick labels and axis 
      mar = c(2.5, 3.2, 1.5, 0.8)) # decrease the borders
  
  # create the empty plot area
  plot(1, type = "n",
       cex.lab = 1.25,
       xlab = "",
       ylab = expression("Change in p: " * Delta * "p"),
       xlim = c(-0.01, 1.01), xaxs = "i",  # set a fixed padding for x axis
       ylim = range(deltap.df$delta_p),
       xaxt = "n", yaxt = "n"  # Suppress default x & y axis tick labels
       )
  
  # Add custom x-axis label with controlled margin
  mtext("Allele frequency (p)", side = 1, line = 1.5, cex = 1.25)

  # Add a dotted horizontal line at y = 0
  abline(h = 0, lty = 3, lwd = 0.6, col = "black")
  
  # plot the points
  points(x = deltap.df$p, y = deltap.df$delta_p, pch = 20, col = colour_p,)
  
  # Add custom x-axis tick labels
  axis(1, at = seq(from=0, to=1, by=0.25), labels = c("0", "", "0.5", "", "1"))
  
  # custom y-axis tick labels
    # for this model the top side is always furthest from 0
  # add tick label at 0 and top
  axis(2, at = c(0, max(deltap.df$delta_p)),
       labels = signif(c(0, max(deltap.df$delta_p)), digits = 1))
  
  # plot the arrows and the attractors
  if(steady_state$outcome == "bistable"){
    arrows(x0 = steady_state$arrow1["x0"],
           y0 = steady_state$arrow1["y0"],
           x1 = steady_state$arrow1["x1"],
           y1 = steady_state$arrow1["y1"],
           length=0.2)
    arrows(x0 = steady_state$arrow2["x0"],
           y0 = steady_state$arrow2["y0"],
           x1 = steady_state$arrow2["x1"],
           y1 = steady_state$arrow2["y1"],
           length=0.2)
    #plot attractors
    points(x=1, y=0, pch=23, bg=substring(colours_genotypes["AA"], first=1, last=7),
           col="white", cex=1.2, lwd=2)
    points(x=0, y=0, pch=23, bg=substring(colours_genotypes["aa"], first=1, last=7),
           col="white", cex=1.2, lwd=2)
    # Add a legend with just the point for the aa attractor
    legend("bottomright", "Attractors", text.col="white", pch=18, col=colours_genotypes["aa"],
           cex=1.3, inset=c(0.12,0.9), xpd=TRUE, horiz=TRUE, bty="n"
    )
    # Add the legend for the AA attractor
    legend("bottomright", "Attractors", pch=18, col=colours_genotypes["AA"],
           cex=1.3, inset=c(0.03,0.9), xpd=TRUE, horiz=TRUE, bty="n"
    )
    
  } else if(steady_state$outcome == "p=1"){
    arrows(x0 = steady_state$arrow1["x0"],
           y0 = steady_state$arrow1["y0"],
           x1 = steady_state$arrow1["x1"],
           y1 = steady_state$arrow1["y1"],
           length=0.2)
    #plot attractor
    points(x=steady_state$stable_pts, y=0, pch=23, bg=substring(colours_genotypes["AA"], first=1, last=7),
           col="white", cex=1.2, lwd=2)
    # Add the legend for the attractor
    legend("bottomright", "Attractor", pch=18, col=colours_genotypes["AA"],
           cex=1.3, inset=c(0.03,0.9), xpd=TRUE, horiz=TRUE, bty="n"
    )
    
  } else if(steady_state$outcome == "p intermediate"){
    arrows(x0 = steady_state$arrow1["x0"],
           y0 = steady_state$arrow1["y0"],
           x1 = steady_state$arrow1["x1"],
           y1 = steady_state$arrow1["y1"],
           length=0.2)
    arrows(x0 = steady_state$arrow2["x0"],
           y0 = steady_state$arrow2["y0"],
           x1 = steady_state$arrow2["x1"],
           y1 = steady_state$arrow2["y1"],
           length=0.2)
    #plot attractor
    points(x=steady_state$stable_pts, y=0, pch=23, bg=substring(colours_genotypes["Aa"], first=1, last=7),
           col="white", cex=1.2, lwd=2)
    # Add the legend for the attractor
    legend("bottomright", "Attractor", pch=18, col=colours_genotypes["Aa"],
           cex=1.3, inset=c(0.03,0.9), xpd=TRUE, horiz=TRUE, bty="n"
    )
  }
}

# a GGPLOT function to plot genotype frequencies over finite time
ggplot_genotype_finite <- function(genotypes_df, genotype_levels){
  # get x-axis labels
  tickwidth_x <- floor(max(genotypes_df$gens)/3)
  labs_x <- c(0, seq(from=tickwidth_x, to=2*tickwidth_x, by=tickwidth_x))
  
  # define y-axis labels
  ticks_y <- seq(from = 0, to = 1, length.out = 5)
  
  # being very specific about the order in which the lines are plotted
  df <- genotypes_df
  df$genotyp <- factor(df$genotyp, levels = genotype_levels)
  df <- df[order(df$genotyp), ]
  
  # create the plot
  finit_plot <- ggplot(df,
                       aes(x = gens, y = freq, colour = genotyp)) +
    geom_line(linewidth = 2.5) +
    scale_colour_manual(values = colours_genotypes[genotype_levels]) +
    scale_x_continuous(limits = c(min(genotypes_df$gens), max(genotypes_df$gens)), expand = c(0, 0),
                       breaks = labs_x,
                       labels = labs_x) +
    scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0),
                       breaks = ticks_y,
                       labels = c(ticks_y[1], "", ticks_y[3], "", ticks_y[5])) +
    labs(x="Generation", y="Genotype Frequency") +
    theme_bw() + 
    theme(text = element_text(size=14),
          axis.text=element_text(size=12),
          axis.line.y=element_blank(),
          axis.title.y = element_text(hjust=0.9), # move y-axis label down
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "null",
          plot.margin = margin(t=3, r=2, b=0, l=0)
    )
  
  return(finit_plot)
}

# a GGPLOT function to plot the long-term genotype freq outcome
ggplot_genotype_inf <- function(allele_df, steady_state, genotype_levels){
  
  # if the system has 2 attractors
  if (steady_state$outcome == "bistable"){
    # we need to check the initial conditions as compared to the unstable equilibrium point
    index_stable_pt <- ifelse(allele_df$p[1] < steady_state$unstable_pt,
                              1, 2)
    # use the appropriate stable state value
    inf_df <- data.frame(genotype = c("AA", "Aa", "aa"),
                         freq = c(steady_state$stable_pts[index_stable_pt]^2,
                                  2*steady_state$stable_pts[index_stable_pt]*(1-steady_state$stable_pts[index_stable_pt]),
                                  (1-steady_state$stable_pts[index_stable_pt])^2))
    
    # if the system has only 1 attractor
  } else {
    # we can add the line directly from the steady state
    inf_df <- data.frame(genotype = c("AA", "Aa", "aa"),
                         freq = c(steady_state$stable_pts^2,
                                  2*steady_state$stable_pts*(1-steady_state$stable_pts),
                                  (1-steady_state$stable_pts)^2))
  }
  
  # being very specific about the order in which the lines are plotted
  df <- inf_df
  df$genotype <- factor(df$genotype, levels = genotype_levels)
  df <- df[order(df$genotype), ]
  
  # create the plot
  inf_plot <- ggplot() +
    geom_hline(data = df, aes(yintercept = freq, colour = genotype),
               linewidth = 2.5) +
    scale_colour_manual(values = colours_genotypes[genotype_levels]) +
    scale_x_continuous(limits=c(0.5, 1.5), breaks = c(1), labels = "long-\nterm") +
    scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
    labs(x="") +
    theme_bw() + theme(text = element_text(size=14),
                       axis.text=element_text(size=12),
                       axis.line.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.y=element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       plot.margin = margin(t=3, r=2, b=0, l=0),
                       legend.position = "null"
    )
  return(inf_plot)
}


# a GGPLOT function to combine the finite and infinite time plots for genotype frequencies
ggplot_genotype_finiteANDoutcome <- function(sim_list, steady_state){
  # get the genotype dataframe from the list of sims
  genotypes_df <- sim_list$genotypes
  # and the allele dataframe too
  allele_df <- sim_list$allele
  
  # this is needed for consistent layering of the 3 genotype lines
  genotype_levels <- c("aa", "Aa", "AA")
  
  # make the plot
  wrap_plots(ggplot_genotype_finite(genotypes_df, genotype_levels),
             ggplot_genotype_inf(allele_df, steady_state, genotype_levels),
             ncol = 2, widths = c(3,1))
}



# a GGPLOT function to plot allele freq over finite time
ggplot_p_finite <- function(p_df){
  # get x-axis labels
  tickwidth_x <- floor(max(p_df$gens)/3)
  labs_x <- c(0, seq(from=tickwidth_x, to=2*tickwidth_x, by=tickwidth_x))
  
  # define y-axis labels
  ticks_y <- seq(from = 0, to = 1, length.out = 5)
  
  # create the plot
  finit_plot <- ggplot(p_df,
                       aes(x = gens, y = p)) +
    geom_line(colour = colour_p, linewidth = 2.5) +
    scale_x_continuous(limits = c(min(p_df$gens), max(p_df$gens)), expand = c(0, 0),
                       breaks = labs_x,
                       labels = labs_x) +
    scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0),
                       breaks = ticks_y,
                       labels = c(ticks_y[1], "", ticks_y[3], "", ticks_y[5])) +
    labs(x="Generation", y="Allele Frequency") +
    theme_bw() + 
    theme(text = element_text(size=14),
          axis.text=element_text(size=12),
          axis.line.y=element_blank(),
          axis.title.y = element_text(hjust=0.7), # move y-axis label down
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "null",
          plot.margin = margin(t=3, r=2, b=0, l=0)
    )
  
  return(finit_plot)
}

# a GGPLOT function to plot the long-term allele freq outcome
ggplot_p_inf <- function(p_df, steady_state){
  
  # if the system has 2 attractors
  if (steady_state$outcome == "bistable"){
    # we need to check the initial conditions as compared to the unstable equilibrium point
    index_stable_pt <- ifelse(p_df$p[1] < steady_state$unstable_pt,
                              1, 2)
    # use the appropriate stable state value
    inf_df <- data.frame(Value = steady_state$stable_pts[index_stable_pt])
    
    # if the system has only 1 attractor
  } else {
    # we can add the line directly from the steady state
    inf_df <- data.frame(Value = steady_state$stable_pts)
  }
  # create the plot
  inf_plot <- ggplot() +
    geom_hline(data = inf_df, aes(yintercept = Value), color = colour_p,
               linewidth = 2.5) +
    scale_x_continuous(limits=c(0.5, 1.5), breaks = c(1), labels = "long-\nterm") +
    scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
    labs(x="") +
    theme_bw() + theme(text = element_text(size=14),
                       axis.text=element_text(size=12),
                       axis.line.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.y=element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       plot.margin = margin(t=3, r=2, b=0, l=0),
                       legend.position = "null"
    )
  return(inf_plot)
}

# a GGPLOT function to combine the finite and infinite time plots for allele frequency
ggplot_p_finiteANDoutcome <- function(sim_list, steady_state){
  # get just the allele dataframe from the list of sims
  allele_df <- sim_list$allele
  
  # make the plot
  wrap_plots(ggplot_p_finite(allele_df),
             ggplot_p_inf(allele_df, steady_state),
             ncol = 2, widths = c(3,1))
}


