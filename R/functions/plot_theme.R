## ----------------------------------------- ##
#          Publication Plot Theme ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#       Publication Plot Theme Helpers ----
#   Sets a consistent ggplot2 theme for all
#   manuscript figures (fonts, sizing, layout).
## ----------------------------------------- ##

set_plot_theme <- function() {
  
  # Register bundled font for cross-platform rendering
  fontfam <- initialize_font()
  
  # Set a global ggplot2 theme for the current R session
  theme_set(theme_ipsum_rc(base_family = fontfam) +
              theme(
                text = element_text(family = fontfam),
        
                plot.title    = element_text(size = 14),
                plot.subtitle = element_text(size = 12),
                axis.title.x  = element_text(size = 10),
                axis.title.y  = element_text(size = 10),
                axis.text.x   = element_text(size = 8),
                axis.text.y   = element_text(size = 8),
                strip.text.x  = element_text(size = 8),
                strip.text.y  = element_text(size = 8),
                
                legend.title  = element_text(size = 8),
                legend.text   = element_text(size = 6,
                                             margin = margin(l = 0.25, unit = 'mm'),
                                             vjust = 0.675),
                legend.key.size = unit(0.35, 'cm'),
                
                plot.background = element_blank(),
                plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      )
  )
  
  invisible(fontfam)
}