#
# author: Ana-Hermina Ghenu
# contact: ana.hermina.ghenu at gmail dot com
# date created: 14 May 2025
#
#--------------------------------------------------
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/
#

# load the environment

library(tidyverse) # for ggplot and pipe syntax
library(shiny) # for app
#library(patchwork) # for combining multiple plots into 1
library(bslib) # for most recent recommended UI options

# set the default plot style
fave_theme <- theme_light() + # see other options at https://ggplot2.tidyverse.org/reference/ggtheme.html
  theme(text = element_text(size=19), # larger text size for titles & axes
        panel.grid.major = element_blank(), # remove major gridlines
        panel.grid.minor = element_blank()) # remove minor gridlines
theme_set(fave_theme)

# increase the default thickness of lines for all geom objects in ggplot2
update_geom_defaults("line", list(linewidth = 2.5, alpha=0.8))

# define colour scheme for genotypes
  # use rgb to set transparency (alpha = 0.75) # note that max value for alpha is 255 so 0.75*255 ~ 191
colours_genotypes <- c(AA = rgb(178, 34, 34, maxColorValue = 255, alpha = 191),  # Firebrick
                       Aa = rgb(255, 165, 0, maxColorValue = 255, alpha = 191),  # Orange
                       aA = rgb(255, 165, 0, maxColorValue = 255, alpha = 191),  # Orange
                       aa = rgb(255, 215, 0, maxColorValue = 255, alpha = 191))  # Gold

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
  # initialize data.frame for storage starting from initial allele frequencies at generation 0
  df <- data.frame(gens = 0:num_gens,
                   p = c(freq_init, rep(NA, num_gens)))
  # loop through each generation
  for(i in 2:nrow(df)){
    W_bar <- mean_fitness(w_genotypes, df$p[i-1])
    # use current mean fitness to get allele frequency in the next generation
    df$p[i] <- next_p(w_genotypes["AA"], w_genotypes["Aa"], df$p[i-1], W_bar)
  }
  # finally, use the allele frequencies for each generation to get the genotypes
  df <- df %>% mutate(aa = (1-p)^2,
                      hetero = 2*p*(1-p),
                      AA = p^2)
  return(df)
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

# a function indicating the long-term steady state
get_steady_state <- function(w_genotypes, deltap.df){
  # some presets for plotting the arrows
  min_arrow_length = 0.07
  arrow_vertical_whitespace = 0.1
  
  # initialize output variable with NA values:
    # depending on outcome, only some may be replaced with non-NA values
  output <- list(outcome = NA, stable_pts = NA, unstable_pt = NA, arrow1 = NA, arrow2 = NA)
  
  # underdominance always results in bistability
  if(w_genotypes["Aa"] < w_genotypes["aa"]) {
    # a hacky way to approximate where the unstable equilibrium point is:
    unstable_pt <- mean(deltap.df$p[which(near(deltap.df$delta_p[12:990], 0, tol = 10^-3))])
    # draw the arrows symmetrically about x=0 at halfway of the smallest max delta_p values
    x_val = min(abs(min(deltap.df$delta_p)), abs(max(deltap.df$delta_p)))/2
    # y values for bottom arrow (1): vertically in middle of the 2 equil pts
    arrow1_length <- max(unstable_pt-2*arrow_vertical_whitespace, min_arrow_length)
    # y values for top arrow (2): vertically in middle of the 2 equil pts
    arrow2_length <- max(1-unstable_pt-2*arrow_vertical_whitespace, min_arrow_length)
    output <- list(outcome = "bistable",
                stable_pts = c(0, 1),
                unstable_pt = unstable_pt,
                arrow1 = c(x0 = -x_val, y0 = unstable_pt - arrow_vertical_whitespace,
                           x1 = -x_val, y1 = unstable_pt - arrow_vertical_whitespace - arrow1_length),
                arrow2 = c(x0 = x_val, y0 = unstable_pt + arrow_vertical_whitespace,
                           x1 = x_val, y1 = unstable_pt + arrow_vertical_whitespace + arrow2_length)
                )
    
  # additivity / partial dominance always results in A going to fixation (1 stable attractor)
  } else if(w_genotypes["Aa"] >= w_genotypes["aa"] & 
            w_genotypes["Aa"] < w_genotypes["AA"]) {
    x_val <- 0.005 # this x-value looks nice
    # increase the whitespace of the arrow
    arrow_vertical_whitespace <- 3*arrow_vertical_whitespace
    arrow1_length <- 1-2*arrow_vertical_whitespace
    output <- list(outcome = "p=1",
                   stable_pts = 1,
                   arrow1 = c(x0 = x_val, y0 = arrow_vertical_whitespace,
                              x1 = x_val, y1 = arrow_vertical_whitespace + arrow1_length)
                   )
    
  # overdominance always results in 1 stable attractor at intermediate allele freq  
  } else if(w_genotypes["Aa"] > w_genotypes["AA"]){
    # a hacky way to approximate where the stable equilibrium point is:
    stable_pt <- mean(deltap.df$p[which(near(deltap.df$delta_p[5:980], 0, tol = 10^-3))])
    # draw the arrows symmetrically about x=0 at one third of the smallest max delta_p values
    x_val = min(abs(min(deltap.df$delta_p)), abs(max(deltap.df$delta_p)))/3
    # y values for bottom arrow (1): vertically in middle of the 2 equil pts
    arrow1_length <- max(stable_pt-2*arrow_vertical_whitespace, min_arrow_length)
    # y values for top arrow (2): vertically in middle of the 2 equil pts
    arrow2_length <- max(1-stable_pt-2*arrow_vertical_whitespace, min_arrow_length)
    output <- list(outcome = "p intermediate",
                   stable_pts = stable_pt,
                   arrow1 = c(x0 = x_val, y0 = stable_pt - arrow_vertical_whitespace - arrow1_length,
                              x1 = x_val, y1 = stable_pt - arrow_vertical_whitespace),
                   arrow2 = c(x0 = -x_val, y0 = stable_pt + arrow_vertical_whitespace + arrow2_length,
                              x1 = -x_val, y1 = stable_pt + arrow_vertical_whitespace)
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
  par(cex.lab = 2, # increase size of axis labels
      cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(4.2, 4.4, 1, 0.1)) # decrease the borders
  
  # make the barplot in graphics (i.e., base R)
  barplot(
    w_genotypes, 
    col = colours_genotypes[names(w_genotypes)], 
    border = "black", 
    ylim = c(0, 2), 
    ylab = "Fitness",
    xlab = "Genotype"
  )
  
  # add a box around the whole plot
  box()
}

# a function to plot the genotypes over time
plot_t_genotypes <- function(sim_df) {
  # Define colors for genotypes
  temp_colours <- colours_genotypes[c(1,2,4)]
  
  # change the graphing settings
  par(cex.lab = 2, # increase size of axis labels
      cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(4.2, 4.4, 0.7, 0.85)) # decrease the borders
  
  # Set up an empty plot
  plot(sim_df$gens, sim_df$AA, type = "n",
       ylim = c(-0.01, 1.01), yaxs = "i", # remove y-axis padding
       xlab = "Generation", ylab = "Genotype frequency")
  
  # Add lines for each genotype
  lines(sim_df$gens, sim_df$AA, col = temp_colours[1], lwd = 7)
  lines(sim_df$gens, sim_df$hetero, col = temp_colours[2], lwd = 7)
  lines(sim_df$gens, sim_df$aa, col = temp_colours[3], lwd = 7)
  
  # Optionally, add a legend
  #legend("topright", legend = c("AA", "aA & Aa", "aa"), col = temp_colours, 
   #      lty = 1, lwd = 2.5, bty = "n")
}


# a function to plot the allele frequency over time
plot_t_allele <- function(sim_df) {
  # change the graphing settings
  par(cex.lab = 2, # increase size of axis labels
      cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(4.2, 4.4, 0.7, 0.85)) # decrease the borders
  
  # Set up an empty plot
  plot(sim_df$gens, sim_df$p, type = "n",
       ylim = c(-0.01, 1.01), yaxs = "i", # remove y-axis padding
       xlab = "Generation", ylab = "Allele frequency")
  
  # Add the line
  lines(sim_df$gens, sim_df$p, col = "black", lwd = 7)
}

plot_p_byp <- function(deltap.df, steady_state) {
  # change the graphing settings
  par(cex.lab = 2, # increase size of axis labels
      cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(4.2, 4.4, 0.7, 0.85)) # decrease the borders
  
  # plot the points
  plot(x = deltap.df$delta_p, y = deltap.df$p,
       pch = 16, col = "black", 
       xlab = expression("Change in allele frequency (" * Delta * "p)"), 
       ylab = "Allele frequency", xlim = range(deltap.df$delta_p),
       ylim = c(-0.01, 1.01), yaxs = "i")
  
  # plot the arrows
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
  } else if(steady_state$outcome == "p=1"){
    arrows(x0 = steady_state$arrow1["x0"],
           y0 = steady_state$arrow1["y0"],
           x1 = steady_state$arrow1["x1"],
           y1 = steady_state$arrow1["y1"],
           length=0.2)
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
  }
  
  # Add a dotted vertical line at x = 0
  abline(v = 0, lty = 3, lwd = 0.6, col = "black")
}


###################################
# make the R Shiny app
###################################


# user interface with fluidPage
ui <- fluidPage(
  
  # title
  titlePanel("Diploid model of selection with underdominance"),
  
  # top row is split into 3 columns of equal size (with 3 sliders)
  fluidRow(
    column(4,
           sliderInput(inputId = "h",
                       label = "Underdominance coefficient (h)",
                       value = -0.2, 
                       min = -1, max = 1,
                       step = 0.1)
           ),
    column(4,
           sliderInput(inputId = "p_init",
                       label = "Initial allele frequency of A (p)",
                       value = 0.32, 
                       min = 0, max = 1,
                       step = 0.02)
           ),
    column(4,
           sliderInput(inputId = "gen",
                       label = "Number of generations",
                       value = 50, 
                       min = 5, max = 950,
                       step = 5)
          )
    ),
  # middle row is 2 columns of equal size
  fluidRow(
    # interactively display the model parameters
    column(6,
           uiOutput("dynamEq"),
           plotOutput("fitLand")),
    # display the genotype frequencies over time
    column(6,
           plotOutput("genoTime")
           )
  ),
  # bottom row is 2 columns of equal size
  fluidRow(
    # display delta p by p
    column(6,
           plotOutput("delta_p")
           ),
    # display the allele frequency over time
    column(6,
           plotOutput("alleleTime")
           )
  )
  
)

# Server logic
server <- function(input, output) {
  # fix the selection coefficient
    # be careful to choose a value that does *NOT* result in neutral system
  s <- 0.206
  
  ######################
  # expressions
  ######################
  
  # a reactive expression to get the fitness of the genotypes
  w_gent <- reactive(fitness_genotypes(sel_coef = s,
                                       dom_coef = input$h))
  
  # a reactive expression to get delta p as a function of allele frequecy (p)
  delta_p <- reactive(get_p_byp(w_gent()))
  
  # a reactive expression to get the long-term model outcome
  outcome <- reactive(get_steady_state(w_genotypes = w_gent(), deltap.df = delta_p()))
  
  # a reactive expression to simulate evolution of the system over time
  sims <- reactive(sim_forward_time(num_gens = input$gen,
                                    freq_init = input$p_init,
                                    w_gent()))
  
  ######################
  # outputs
  ######################
  
  output$dynamEq <- renderUI({
    temp_genotypes <- w_gent()
    w_AA <- temp_genotypes["AA"]
    w_Aa <- temp_genotypes["Aa"]
    w_aa <- temp_genotypes["aa"]
    
    withMathJax(sprintf("$${\\tiny p_{t+1} = \\frac{p_t^2 \\cdot (1+s) + p_t(1-p_t) \\cdot (1+s+h)}{\\overline{w}} = \\frac{p_t^2 \\cdot %.01f + p_t(1-p_t) \\cdot %.01f}{p_t^2 \\cdot %.01f + 2p_t(1-p_t) \\cdot %.01f + (1-p_t)^2 \\cdot %.01f}}$$",
                        w_AA, w_Aa, w_AA, w_Aa, w_aa))
  })
  output$fitLand <- renderPlot({plot_schematic(w_gent())})
  output$genoTime <- renderPlot({plot_t_genotypes(sims())})
  output$delta_p <- renderPlot({plot_p_byp(delta_p(), outcome())})
  output$alleleTime <- renderPlot({plot_t_allele(sims())})
  
}

# Complete app with UI and server components
shinyApp(ui, server)

