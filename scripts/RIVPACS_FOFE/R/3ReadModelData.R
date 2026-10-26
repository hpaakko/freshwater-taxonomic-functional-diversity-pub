rivpacs_file <- switch(ecosystem,
  lakes  = "Strict_Henna_RIVPACS_data_lakes.RData",
  rivers = "Strict_Henna_RIVPACS_data_rivers.RData"
)

load(file.path("data", rivpacs_file))

print(ls())
