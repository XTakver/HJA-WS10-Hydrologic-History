# ------------------------------------------------------------ #
# initialize_font()
#   Register bundled Roboto Condensed fonts for reproducible
#   cross-platform plotting (e.g., CairoPDF) using showtext.
# ------------------------------------------------------------ #

initialize_font <- function() {
  
  # Dependency checks (fail fast with a clear message)
  if (!requireNamespace('sysfonts', quietly = TRUE)) stop('sysfonts is required.')
  if (!requireNamespace('showtext', quietly = TRUE)) stop('showtext is required.')
  
  # Font files are shipped with this repo for consistent rendering across machines.
  # Expected directory: assets/fonts/roboto_condensed/
  font_dir <- file.path('assets', 'fonts', 'roboto_condensed')
  
  # Register the font family name used by tinythemes/ggplot themes.
  # This prevents "font family not found" warnings on systems without Roboto installed.
  sysfonts::font_add(family = 'Roboto Condensed',
                     regular = file.path(font_dir, 'RobotoCondensed-Regular.ttf'),
                     bold = file.path(font_dir, 'RobotoCondensed-Bold.ttf'),
                     italic = file.path(font_dir, 'RobotoCondensed-Italic.ttf'),
                     bolditalic = file.path(font_dir, 'RobotoCondensed-BoldItalic.ttf'))
  
  # Optional alias: convenient if you prefer an underscore family name in your own code.
  sysfonts::font_add(family = 'roboto_condensed',
                     regular = file.path(font_dir, 'RobotoCondensed-Regular.ttf'),
                     bold = file.path(font_dir, 'RobotoCondensed-Bold.ttf'),
                     italic = file.path(font_dir, 'RobotoCondensed-Italic.ttf'),
                     bolditalic = file.path(font_dir, 'RobotoCondensed-BoldItalic.ttf'))
  
  # Enable showtext so devices (including CairoPDF) render from the bundled TTFs.
  # Note: for device-specific control during export, you can also use
  # showtext::showtext_begin() / showtext::showtext_end().
  showtext::showtext_auto(enable = TRUE)
  
  # Return the canonical family name for use in theme_set(), element_text(), etc.
  'Roboto Condensed'
}