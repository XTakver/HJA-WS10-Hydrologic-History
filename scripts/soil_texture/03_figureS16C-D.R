# -----------------------------------------------------------------------------#
# Purpose
#   Generate the soil texture triangle figures (S16C & S16D) and their
#   shared legend by running run_triangle_figures() (defined in
#   R/functions/FigureS16C-D_rendering.R) inside a disposable callr
#   subprocess.
#
#   A subprocess is required because ggtern permanently patches ggplot2
#   internals when loaded, and this cannot be undone within a running R
#   session -- it would break any ggplot2 figure made afterward. Isolating
#   the ggtern-dependent code this way lets figure scripts be run
#   individually, in any order, without restarting R between them.
#
# Outputs (written to PD; see R/functions/FigureS16C-D_rendering.R)
#   - Fig.S16C_Texture_Triangle.pdf
#   - Fig.S16D_Texture_Triangle.pdf
#   - Fig.S16C-D_Legend.pdf
#
# Created: 01/23/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#

librarian::shelf(here, callr)

callr::r(function() {
  source(here::here('R', 'functions', 'FigureS16C-D_rendering.R'))
  run_triangle_figures()
})


