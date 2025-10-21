# a shiny app whose format is adapted specifically for the experiment
#   this is the main file. It loads dependencies then runs the app.
#author: Ana-Hermina Ghenu
#date: 2025-10-08

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
thematic_shiny(font = font_spec("auto"))

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
  # reduce the font size of slider controls
  tags$style(HTML("
    .control-label {
      font-size: 14px;
    }
  ")),
  # I tried to reduce the size of the slider but I don't think it does anything
  tags$style(HTML("
    .js-range-slider {
      height: 0px !important;
    }
  ")),
  # reduce the font size of card headers
  tags$style(HTML("
    .card-header {
      font-size: 14px;
    }
  ")),

  
  # typeset math using latex code
  withMathJax(),
  # allow in-line LaTeX via $ in mathjax.
  tags$div(HTML("<script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$'], ['\\\\(','\\\\)']]}
            });
            </script >
            ")),

  HTML("<br>"),
  
  navset_card_pill(# use pill-shaped buttons to navigate
    # explain what the tabs on this card do
    title = "Select a topic:",
    
    nav_panel("Introduction",
              markdown(text_intro)),

    nav_panel("App How-to & Things to Try", markdown(text_howto)),
    
    nav_panel("What's going on?", markdown(text_huh)),
    
    nav_panel("Key Insight", markdown(text_genotypes)),
    
    nav_panel("Summary", markdown(text_summary))
  ),
  
  # Top row with sliders
  fluidRow(
    column(4,
           sliderInput(inputId = "h_slider",
                       label = "Dominance coefficient $(h)$",
                       value = 0.1, 
                       min = -0.4, max = 0.3,
                       step = 0.1)
    ),
    column(4,
           sliderInput(inputId = "p_init",
                       label = "Initial $A$ allele freq. $(p_0)$",
                       value = 0.02, 
                       min = 0.02, max = 0.5,
                       step = 0.04)
    ),
    column(4,
           sliderInput(inputId = "gen",
                       label = "Number of generations",
                       value = 60, 
                       min = 20, max = 200,
                       step = 20)
    )
  ),
  
  # Middle and bottom rows replaced with cards in a responsive grid
  layout_columns(height="440px", # sets the overall height of the 2 column layout
                 gap = "0rem", # removes the whitespace between cards
    card(
      card_header("1. Fitness landscape & Equation"),
      card_body(plotOutput("fitLand", height="125px"),
                uiOutput("dynamEq"),
                class = "align-items-center")
    ),
    card(
      card_header("3. Genotypes over Time"),
      plotOutput("genoTime", height="150px")
    ),
    card(
      card_header(HTML("2. Change in Allele Freq. <span style='white-space: nowrap;'>($\\Delta p = p_{t+1} - p_t$)</span>")),
      plotOutput("delta_p", height="150px")
    ),
    card(
      card_header("4. Allele Freq. ($p$) over Time"),
      plotOutput("alleleTime", height="150px")
    ),
    col_widths = breakpoints(
      sm = c(7, 5, 7, 5), # for portrait-mode phones, left column is wider than right column
      md = c(6, 6, 6, 6) # columns are equally sized when it gets wide enough
    )
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
    the_eqn <- paste0("$$\\scriptsize{ p_{t+1} = \\frac{p_t^2 \\cdot ",
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
