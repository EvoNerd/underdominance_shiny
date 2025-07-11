# a shiny app to explore bistability in the context of underdominance
#   this is the main file. It loads dependencies then runs the app.
#author: Ana-Hermina Ghenu
#date: 2025-05-25

# load the environment
source("underdom_funs.R")
source("underdom_text.R")
library(shiny) # for app
library(bslib) # for most recent recommended UI options
library(thematic) # for converting R plots to have consistent theme as from bslib

# set the theme for bslib objects
app_theme <- bs_theme(bootswatch = "simplex", # set simplex theme
                      version = 5,
                      primary = substring(colours_genotypes["Aa"], first=1, last=7), # set the primary: this matches heterozygotes
                      secondary = colour_p) # secondary colour matches colour used for p

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
  markdown(text_subtitle),
  
  navset_card_pill(# use pill-shaped buttons to navigate
    # explain what the tabs on this card do
    title = "Select a topic:",
    
    nav_panel("Introduction",
              markdown(text_intro)),

    nav_panel("How to use the app", markdown(text_howto)),
    
    nav_panel("What's going on?", markdown(text_huh)),
    
    nav_panel("Key Insight", markdown(text_genotypes)),
    
    nav_panel("Things to try", markdown(text_try)),
    
    nav_panel("Optional: code", markdown(text_code)),
    
    nav_panel("Summary", markdown(text_summary))
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
                       label = "Initial $A$ allele frequency $(p_0)$",
                       value = 0.02, 
                       min = 0.02, max = 0.98,
                       step = 0.02)
    ),
    column(4,
           sliderInput(inputId = "gen",
                       label = "Number of generations",
                       value = 50, 
                       min = 5, max = 250,
                       step = 5)
    )
  ),
  
  # Middle and bottom rows replaced with cards in a responsive grid
  layout_columns(height="600px",
    card(
      card_header("1. Model schematic & Equation"),
      plotOutput("fitLand", height = "346px"),
      uiOutput("dynamEq"),
      card_style = "width: 3000;"
    ),
    card(
      card_header("3. Genotype Frequencies over Time"),
      plotOutput("genoTime")
    ),
    card(
      card_header("2. Change in Allele Frequency"),
      plotOutput("delta_p")
    ),
    card(
      card_header("4. Allele Frequency over Time"),
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
    # first define the equation by colouring the dynamic variables with their appropriate colours
    the_eqn <- paste0("$${ p_{t+1} = \\frac{p_t^2 \\cdot ",
                      "\\color{", substring(colours_genotypes["AA"], first=1, last=7),"}{%.01f}",
                      " + p_t(1-p_t) \\cdot ",
                      "\\color{", substring(colours_genotypes["Aa"], first=1, last=7),"}{%.01f}","}{\\overline{w}}}$$")
    # then typeset the equation using MathJax
    withMathJax(sprintf(the_eqn, w_AA, w_Aa))
  })
  # render fitness landscape schematic:
  output$fitLand <- renderPlot({plot_schematic(w_gent())})
  # render genotypes over time plot
  output$genoTime <- renderPlot({ggplot_genotype_finiteANDoutcome(sims(), outcome())})
  # render p vs delta p plot
  output$delta_p <- renderPlot({plot_dp_by_p(delta_p(), outcome())})
  # render allele frequency over time plot
  output$alleleTime <- renderPlot({ggplot_p_finiteANDoutcome(sims(), outcome())})
}

# Complete app with UI and server components
shinyApp(ui, server)
