# underdominance shiny app dependency for app.R
#   tutorial text that appears on the top of the app is defined here
#author: Ana-Hermina Ghenu
#date: 2025-05-25

# use the same colours as used for the bslib theme
html_col1 <- sprintf("<span style='color: %s;'>", colour_p)
html_colAA <- sprintf("<span style='color: %s;'>", colours_genotypes["AA"])
html_colAa <- sprintf("<span style='color: %s;'>", colours_genotypes["Aa"])
html_colaa <- sprintf("<span style='color: %s;'>", colours_genotypes["aa"])


# subtitle text (appears underneath the app title)
text_subtitle <- paste0(html_col1, "When heterozygotes have the lowest fitness, this creates a system with two attractors where the final outcome depends on the initial allele frequency.</span>")




# Introduction
text_intro <- paste0("Imagine two alleles: an invasive, very beneficial allele, **$A$**, and a local, neutral allele, **$a$**.

You might expect $A$ to always fix in a diploid population, but this changes when we consider genotype fitness:

- ",html_colaa,"**$aa$ homozygote:**</span> fitness $w_{aa}= 1$

- ",html_colAA,"**$AA$ homozygote:**</span> has a 30% fitness advantage $w_{AA}= 1+s = 1.3$

- ",html_colAa,"**$Aa$ and $aA$ heterozygotes:**</span> fitness depends on a slider-controlled value, the **dominance coefficient ($h$)**.
<br>$w_{Aa}= w_{aA}= 1+s+h = 1.3 +h$

The fitness of heterozygotes determines whether the invasive allele $A$ can spread or not. That’s what you’ll explore here.

",html_col1,"*Play with the app using the sliders below or select another topic using the buttons above.*</span>")




# How to
text_howto <- paste0("#### Slider 1. Dominance coefficient ($h$)

Adjust the fitness of the heterozygotes ($Aa$ and $aA$) as follows:

- set $h < 0$: heterozygotes have lowest fitness. **This is a system with two attractors.**

- set $0 \\le h \\le 0.3$: heterozygote fitness is intermediate to homozygotes. **This is a system with one attractor.**

*",html_col1,"Slide $h$ to see how the heterozygotes' fitness changes on the fitness landscape (panel 1).</span>
<br>What other plots or equations change with $h$?*

#### Slider 2. Initial allele frequency ($p_0$)

Adjust the fraction of invasive $AA$ individuals that are introduced into the local $aa$ population at the start of the simulation.

*",html_col1,"Slide $p_0$ to see how the genotype frequencies change over time (panel 3).</span>
<br>What other plots change with $p_0$?*

",html_col1,"Set $h < 0$ and vary $p_0$:</span> *How does the starting point change the outcome?*

#### Slider 3. Number of generations

Choose how far into the near-future you want to see the population evolve.

",html_col1,"*How do panels 3 & 4 change when you adjust the third slider?*</span>")




# What's going on?
text_huh <- paste0("The allele frequency of $A$ changes over time based on how fit and frequent the genotypes are.

The difference equation ($p_{t+1}$ in the bottom of Panel 1) shows how $p$, the frequency of $A$, evolves each generation.
<br>We also visualize this in two plots:

- **Panel 2:** Shows in which direction and how quickly the frequency of the $A$ allele is changing.

- **Panel 4:** Shows how the $A$ allele spreads, or disappears, over generations.

*",html_col1,"Play with the sliders to explore how fitness ($h$) and initial conditions ($p_0$) shape long-term outcomes.</span>
<br>Can $A$ invade when rare? Which values of $h$ make this possible?*")




# Genotype frequencies
text_genotypes <- paste0("**Panel 3** shows how homozygote (",html_colAA,"**$AA$**</span> & ",html_colaa,"**$aa$**</span>) and heterozygote (",html_colAa,"**$Aa$ and $aA$**</span>) genotype frequencies change over time.

",html_col1,"**Key insight:**</span> When $A$ is rare, almost all $A$ alleles are in heterozygote genotypes. If heterozygotes are unfit ($h<0$), $A$ will not spread when rare -- even if ",html_colAA,"$AA$</span> is very fit!")




# Things to try
text_try <- paste0("- ",html_col1,"Set $h$ near 0.15:</span> $A$ always takes over.
<br>This is a system with only one attractor.

- ",html_col1,"Set $h$ below 0:</span> $A$ can spread only when it is introduced by a large invasion of $AA$ individuals (i.e., large $p_0$).
<br>This is a system with two attractors.

- ",html_col1,"Vary $A$'s initial frequency ($p_0$):</span> *How does the starting point change the outcome?*")




# A glimpse at the code
text_code <- "The difference equation in the app below describes how $p$, the frequency of $A$, evolves in the next generation ($p_{t+1}$). We can write it in math like this:
$$ p_{t+1} = \\frac{p_t^2 \\cdot w_{AA} + p_t(1-p_t) \\cdot w_{Aa}}{\\overline{w_t}} $$
We can also code it using the R programming language like this:
```
# next_p: difference equation to get the allele frequency in the next generation
next_p <- function(w_AA, w_heterozygote, current_p, W_bar){
  output <- (current_p^2*w_AA + current_p*(1-current_p)*w_heterozygote) / W_bar
  return(output)
}
```"




# Summary
text_summary <- "Even if $A$ is very beneficial, it might not fix. The fate of the invasive $A$ allele depends on:

- How **fit the heterozygote** is ($h$), and

- How **common** $A$ is to start with ($p_0$).

This is the essence of **underdominance**. When there is underdominance ($h < 0$), the system has two attractors and so, the outcome depends on the initial conditions."