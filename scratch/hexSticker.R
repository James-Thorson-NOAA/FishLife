
setwd( R'(C:\Users\james\OneDrive\Desktop\Git\FishLife\scratch)' )

#pak::pak("emilioxavier/hexSticker")
library(hexSticker)

library(magick)

img <- image_read("logo.png")

# Make white transparent
img <- image_transparent(img, "white")

image_write(img, "logo_transparent.png")

p = sticker(
  subplot = "logo_transparent.png",
  package = "",

  # package text
  p_family = "sans",
  p_size = 20,
  p_color = "white",

  # hexagon
  h_fill = "white",
  h_color = "black",
  h_size = 1,

  # image positioning
  s_x = 1,
  s_y = 1,
  #asp = 1.6,
  s_width = 0.67, # Only use s_width or s_height ... both results in only one used

  # text positioning
  p_x = 1,
  p_y = 0.6,
  width = 2500,
  height = 2500,

  filename = "FishLife-sticker.png"
) # + theme_sticker()
