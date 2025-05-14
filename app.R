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

#library(tidyverse) # for ggplot and pipe syntax
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
colours_genotypes <- c("AA" = "firebrick", "Aa" = "orange", "aA" = "orange", "aa" = "gold")

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
  return(c(aa = 1.0, aA = 1 + dom_coef, Aa = 1 + dom_coef, AA = 1 + sel_coef))

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

###################################
# define visualization functions
###################################

# plot fitness of the genotypes
plot_schematic <- function(w_genotypes) {
  # create df for nice plotting
  df <- data.frame(Fitness = w_genotypes,
                   Genotype = names(w_genotypes))
  # boxplot of fitness
  ggplot(df,
         aes(x = Genotype, y = Fitness, fill = Genotype)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colours_genotypes) +
    theme(legend.position="none") +
    scale_y_continuous(expand = c(0, 0))
}

# a function to plot the genotypes over time
plot_t_genotypes <- function(sim_df){
  # remove the column for allele frequency and change the df to long format
  plot_df <- sim_df %>% select(-p) %>%
    pivot_longer(!gens, names_to="genotype", values_to="freq")
  # make genotype a factor with levels for nice plotting
  plot_df$genotype <- factor(plot_df$genotype,
                             levels = c("AA", "hetero", "aa"))
  # rename the factor levels for nice plotting
  levels(plot_df$genotype) <- c("AA", "aA & Aa", "aa")
  ggplot(plot_df,
         aes(x=gens, y=freq, colour=genotype)) +
    geom_line() +
    scale_color_manual(values = unname(colours_genotypes)[-2]) +
    labs(x= "Generation", y="Genotype frequency", color = "Legend") +
    scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
    # can try removing the legend altogether
    theme(legend.position="none")
}

# a function to plot the allele frequency over time
plot_t_allele <- function(sim_df){
  # keep just the columns for time and allele frequency
  plot_df <- sim_df %>% select(gens, p)
  ggplot(plot_df,
         aes(x=gens, y=p)) +
    geom_line() +
    labs(x= "Generation", y="Allele frequency") +
    scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0))
}

# a function to plot delta p
plot_p_byp <- function(w_genotypes){
  # initialize a vector of allele frequencies from 0 to 1
  plot_df <- data.frame(p = seq(0, 1, length.out = 1000)) %>%
    # get the mean fitness at each allele frequency
    mutate(W_bar = mean_fitness(w_genotypes, p),
           # get delta_p by current allele frequency from that in the next generation
           delta_p = next_p(w_genotypes["AA"], w_genotypes["Aa"], p, W_bar) - p)
  # figure out how to draw the arrows
  # TO DO!!
  ggplot(plot_df,
         aes(x=delta_p, y=p)) +
    geom_vline(xintercept = 0, linewidth=0.6, linetype = "dotted") +
    geom_point() +
    labs(x= expression("Change in allele frequency (" * Delta * "p)"), y="Allele frequency") +
    scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0))
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
                       min = 5, max = 1000,
                       step = 5)
          )
    ),
  # middle row is 2 columns of equal size
  fluidRow(
    # interactively display the model parameters
    column(6,
           "some text here?",
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
           plotOutput("delta_p")),
    # display the allele frequency over time
    column(6,
           plotOutput("alleleTime"))
  )
  
)

# Server logic
server <- function(input, output) {
  # for now fix the selection coefficient to 0.2
  s <- 0.2
  
  ######################
  # expressions
  ######################
  
  # a reactive expression to get the fitness of the genotypes
  w_genotypes <- reactive(fitness_genotypes(sel_coef = s,
                                            dom_coef = input$h))
  
  # a reactive expression to simulate evolution of the system over time
  sims <- reactive(sim_forward_time(num_gens = input$gen,
                                    freq_init = input$p_init,
                                    w_genotypes()))
  
  ######################
  # outputs
  ######################
  
  output$fitLand <- renderPlot({plot_schematic(w_genotypes())})
  output$genoTime <- renderPlot({plot_t_genotypes(sims())})
  output$delta_p <- renderPlot({plot_p_byp(w_genotypes())})
  output$alleleTime <- renderPlot({plot_t_allele(sims())})
  
}

# Complete app with UI and server components
shinyApp(ui, server)

