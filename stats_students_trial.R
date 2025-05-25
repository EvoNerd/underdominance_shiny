#author: Ana-Hermina Ghenu
#date: 2025-05-25

# load the environment

library(tidyverse) # for ggplot and pipe syntax
library(shiny) # for app
library(bslib) # for most recent recommended UI options
library(thematic) # for converting R plots to have consistent theme as from bslib

# define colour scheme for genotypes
# use rgb to set transparency (alpha = 0.75) # note that max value for alpha is 255 so 0.75*255 ~ 191
colours_genotypes <- c(AA = rgb(178, 34, 34, maxColorValue = 255, alpha = 191),  # Firebrick
                       Aa = rgb(255, 165, 0, maxColorValue = 255, alpha = 191),  # Orange
                       aA = rgb(255, 165, 0, maxColorValue = 255, alpha = 191),  # Orange
                       aa = rgb(255, 215, 0, maxColorValue = 255, alpha = 191))  # Gold
colour_p <- "#2D1C56"

# set the theme for bslib objects
app_theme <- bs_theme(bootswatch = "simplex", # set simplex theme
                      version = 5,
                      primary = colour_p, # set the primary: this matches colour used for p
                      secondary = colours_genotypes["Aa"]) # secondary colour matches heterozygotes 

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
  arrow_vertical_whitespace = 0.08
  
  # initialize output variable with NA values:
  # depending on outcome, only some may be replaced with non-NA values
  output <- list(outcome = NA, stable_pts = NA, unstable_pt = NA, arrow1 = NA, arrow2 = NA)
  
  # underdominance always results in bistability
  if(w_genotypes["Aa"] < w_genotypes["aa"]) {
    # a hacky way to approximate where the unstable equilibrium point is:
    unstable_pt <- mean(deltap.df$p[which(near(deltap.df$delta_p[12:990], 0, tol = 10^-3))])
    # draw the arrows symmetrically about x=0 at a quarter of the smallest max delta_p values
    x_val = min(abs(min(deltap.df$delta_p)), abs(max(deltap.df$delta_p)))/4
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
      mar = c(4.2, 4.4, 0.02, 0.02)) # decrease the borders
  
  # make the barplot in graphics (i.e., base R)
  barplot(
    w_genotypes, 
    col = colours_genotypes[names(w_genotypes)], 
    border = "black", 
    ylim = c(0, 1.6), 
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
      mar = c(3.9, 4.4, 2.2, 0.85)) # decrease the borders
  
  # Set up an empty plot
  plot(sim_df$gens, sim_df$AA,
       type = "n",
       ylim = c(-0.01, 1.01), yaxs = "i", # remove y-axis padding
       xlab = "Generation",
       ylab = "Genotype frequency",
       yaxt = "n"  # Suppress default y-axis tick labels
       )
  
  # Add custom y-axis tick labels
  axis(2, at = seq(from=0, to=1, by=0.25), labels = c("0", "", "0.5", "", "1"))
  
  # Add lines for each genotype
  lines(sim_df$gens, sim_df$AA,     col = temp_colours[1], lwd = 7)
  lines(sim_df$gens, sim_df$hetero, col = temp_colours[2], lwd = 7)
  lines(sim_df$gens, sim_df$aa,     col = temp_colours[3], lwd = 7)
}


# a function to plot the allele frequency over time
plot_t_allele <- function(sim_df) {
  # change the graphing settings
  par(cex.lab = 2, # increase size of axis labels
      cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(3.9, 4.4, 2.0, 0.85)) # decrease the borders
  
  # Set up an empty plot
  plot(sim_df$gens, sim_df$p,
       type = "n",
       ylim = c(-0.01, 1.01), yaxs = "i", # remove y-axis padding
       xlab = "Generation",
       ylab = "Allele frequency (p)",
       yaxt = "n"  # Suppress default y-axis tick labels
       )
  
  # Add custom y-axis tick labels
  axis(2, at = seq(from=0, to=1, by=0.25), labels = c("0", "", "0.5", "", "1"))
  
  # Add the line
  lines(sim_df$gens, sim_df$p, col = colour_p, lwd = 7)
}

plot_p_byp <- function(deltap.df, steady_state) {
  # change the graphing settings
  par(cex.lab = 2, # increase size of axis labels
      cex.axis = 1.5, # increase size of tick mark labels 
      mar = c(3.9, 4.4, 2.0, 0.85)) # decrease the borders
  
  # plot the points
  plot(x = deltap.df$delta_p, y = deltap.df$p,
       pch = 16,
       col = colour_p, 
       xlab = expression("Change in allele frequency (" * Delta * " p)"), 
       ylab = "Allele frequency (p)",
       xlim = range(deltap.df$delta_p),
       ylim = c(-0.01, 1.01), yaxs = "i",
       yaxt = "n"  # Suppress default y-axis tick labels
       )
  
  # Add custom y-axis tick labels
  axis(2, at = seq(from=0, to=1, by=0.25), labels = c("0", "", "0.5", "", "1"))
  
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
# define the app
###################################

# call thematic before launching the shiny
  # this will integrate the theme from bslib to how the plots are displayed
thematic_shiny(font = font_spec("auto", scale = 1.8))

# user interface
ui <- fluidPage(

  theme = app_theme,
  
  # set all cards to have a white background
  tags$head(
    tags$style(HTML("
      .card {
        background-color: white !important;
      }
    "))
  ),
  
  # typeset math using latex code
  withMathJax(),
  # allow in-line LaTeX via $ in mathjax.
  tags$div(HTML("<script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$'], ['\\\\(','\\\\)']]}
            });
            </script >
            ")),

  # Title
  titlePanel("Diploid selection with (under)dominance"),
  
  # a subtitle / summary can be added here:
  markdown("Summary / goal of the tutorial."),
  
  navset_card_pill(# use pill-shaped buttons to navigate
    # explain what the tabs on this card do
    title = "Select a topic",
    
    nav_panel("Assumptions",
              markdown("Suppose you have two different alleles in a population, $A$ and $a$. On average, individuals with the $A$ allele are more fit than those with the $a$ allele: they are more likely to survive and reproduce. It would be reasonable to presume that the $A$ allele will eventually replace the $a$ allele in the population. In other words, if we use $p$ to denote the proportion of $A$ alleles in the population, it would be reasonable to expect that after a sufficiently long time $p=1$. 

However, for a diploid population, the outcome of selection depends on the fitness of heterozygous individuals.

This interactive tutorial will explain why the more fit allele, $A$, may not always replace the less fit allele, $a$. And it will allow you to explore how the fitness of heterozygous individuals impacts the outcome of selection.
                       
**Click the \"Dominance coefficient: $h$\" tab above to progress to the next instruction.**")),
    
    
    
    
    nav_panel("Dominance coefficient: $h$", markdown("Recall that for a diploid population there are three different genotypes:

- $aa$: homozygous for $a$. By convention these individuals have a fitness of one, $$w_{aa} = 1$$

- $AA$: homozygous for $A$. Let's suppose selection is very strong and that these individuals are 30% more likely to survive and reproduce than $aa$ homozygous individuals, so $$w_{AA}= w_{aa} + s = 1 + s = 1+0.3 =1.3$$

- $Aa$ or $aA$: heterozygous individuals have both $A$ and $a$ alleles. We assume that the two heterozygotes have equal fitness, $w_{Aa} = w_{aA}$. We will allow the fitness of these individuals to vary by defining their fitness as $$w_{Aa} = 1 + s - h = 1.3 - h$$ **The variable $h$ is called the \"dominance coefficient\".**

In the app below, the first slider allows you to manipulate the dominance coefficient. Intermediate values of the dominance coefficient ($h$ between 0 and 0.3) yield heterozygous individuals with a fitness value that is intermediate between the homozygous individuals:

[ADD FIGURE!!!: schematic for h=0.15]

Small values of the dominance coefficient ($h$ between -1 and 0) mean that heterozygous individuals have the lowest fitness value:

[ADD ANOTHER FIGURE!!!: schematic for h=-0.5]

**Try manipulating the first slider to see how the heterozygote fitness changes.**")),
    
    
    
    
    nav_panel("Recursion equation: $p_{t+1}$", markdown("The allele frequency of $A$ in the next generation ($p_{t+1}$) depends on how many $AA$ homozygous and heterozygous individuals there and the fitness of these individuals as compared to the mean population fitness in the current generation ($\\overline{w_t}$). This equation is called a recursion equation and it is written as follows:
$$ p_{t+1} = \\frac{p_t^2 \\cdot w_{AA} + p_t(1-p_t) \\cdot w_{Aa}}{\\overline{w_t}}$$

$$\\qquad \\qquad \\quad = \\frac{p_t^2 (1+s) + p_t (1-p_t) (1+s+h)}{\\overline{w_t}}$$
        
$$\\qquad \\quad = \\frac{p_t^2 (1.3) + p_t (1-p_t) (1.3+h)}{\\overline{w_t}}$$

A dynamic version of the recursion equation is shown in the app below. Changing the dominance coefficient ($h$) using the first slider will change the values in the recursion equation.

**Try manipulating the first slider to see how the recursion equation changes.**
                                                        
Check this box to see how the recursion equation can be coded in the R programming language:"),
              checkboxInput('show_recursion_eqn', " ", FALSE),
              uiOutput('recursion_eqn')),
    
    
    
    nav_panel("Rate of change: $\\Delta p$", markdown("")),
    
    
    
    nav_panel("Genotype frequencies", markdown(""))
  ),
  
  # Top row with sliders
  fluidRow(
    column(4,
           sliderInput(inputId = "h_slider",
                       label = "Dominance coefficient $(h)$",
                       value = 0.1, 
                       min = -1, max = 0.5,
                       step = 0.1)
    ),
    column(4,
           sliderInput(inputId = "p_init",
                       label = "Initial $A$ allele frequency $(p)$",
                       value = 0.02, 
                       min = 0.02, max = 0.98,
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
  
  # Middle and bottom rows replaced with cards in a responsive grid
  layout_columns(height="600px",
    card(
      #card_header("Dynamic Equilibrium & Fitness Landscape"),
      plotOutput("fitLand"),
      uiOutput("dynamEq"),
      card_style = "width: 3000;"
    ),
    card(
      #card_header("Genotype Frequencies Over Time"),
      plotOutput("genoTime")
    ),
    card(
      #card_header("Δp vs p"),
      plotOutput("delta_p")
    ),
    card(
      #card_header("Allele Frequency Over Time"),
      plotOutput("alleleTime")
    ),
    col_widths = c(6, 6, 6, 6)  # Optional: adjust widths for layout balance
  )
)


# Server logic
server <- function(input, output) {
  # fix the selection coefficient
  # be careful to choose a value that does *NOT* result in neutral system
  s <- 0.305
  
  ######################
  # expressions
  ######################
  
  # a reactive expression to get the fitness of the genotypes
  w_gent <- reactive(fitness_genotypes(sel_coef = s,
                                       dom_coef = input$h_slider))
  
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
  # typeset dynamic recursion equation
  output$dynamEq <- renderUI({
    temp_genotypes <- w_gent()
    w_AA <- temp_genotypes["AA"]
    w_Aa <- temp_genotypes["Aa"]
    # I would like to print the recursion equation with all variables as a reminder
      # but I can't get it to fit in the card width unless the font size is too small!!
    #withMathJax(sprintf("$${\\tiny p_{t+1} = \\frac{p_t^2 \\cdot (1+s) + p_t(1-p_t) \\cdot (1+s+h)}{\\overline{w}} = \\frac{p_t^2 \\cdot %.01f + p_t(1-p_t) \\cdot %.01f}{\\overline{w}}}$$",
    #                    w_AA, w_Aa))
    withMathJax(sprintf("$${ p_{t+1} = \\frac{p_t^2 \\cdot %.01f + p_t(1-p_t) \\cdot %.01f}{\\overline{w}}}$$",
                                            w_AA, w_Aa))
  })
  # render fitness landscape schematic:
  output$fitLand <- renderPlot({plot_schematic(w_gent())})
  # render genotypes over time plot
  output$genoTime <- renderPlot({plot_t_genotypes(sims())})
  # render p vs delta p plot
  output$delta_p <- renderPlot({plot_p_byp(delta_p(), outcome())})
  # render allele frequency over time plot
  output$alleleTime <- renderPlot({plot_t_allele(sims())})
  #
  output$recursion_eqn <- renderUI({
      if (!input$show_recursion_eqn) return()
        markdown("```
# next_p: recursion equation to get the allele frequency in the next generation
next_p <- function(w_AA, w_heterozygote, current_p, W_bar){
  output <- (current_p^2*w_AA + current_p*(1-current_p)*w_heterozygote) / W_bar
  return(output)
}
```")
    })
}

# Complete app with UI and server components
shinyApp(ui, server)
