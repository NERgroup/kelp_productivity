rm(list=ls())

#[jogsmith@ucsc.edu](mailto:jogsmith@ucsc.edu)

#load packages and set wd
library(tidyverse)
library(googlesheets4)
library(janitor)


################################################################################
# Format time as 24-hour HHMM
# e.g., 830 -> "0830", 1346 -> "1346", 13:46 -> "1346"

format_time_24 <- function(x) {
  
  # If Google Sheets reads the column as a date/time
  if (inherits(x, "POSIXt")) {
    return(ifelse(
      is.na(x),
      NA_character_,
      format(x, "%H%M")
    ))
  }
  
  # If Google Sheets reads the column as a time/difftime
  if (inherits(x, "difftime")) {
    secs <- as.numeric(x, units = "secs")
    
    return(ifelse(
      is.na(secs),
      NA_character_,
      sprintf(
        "%02d%02d",
        floor(secs / 3600) %% 24,
        floor((secs %% 3600) / 60)
      )
    ))
  }
  
  # If Google Sheets reads the column as numeric
  if (is.numeric(x)) {
    
    return(ifelse(
      is.na(x),
      NA_character_,
      
      # Google Sheets can store time as a fraction of a day
      ifelse(
        x >= 0 & x < 1,
        sprintf(
          "%02d%02d",
          floor(x * 24),
          floor((((x * 24) - floor(x * 24)) * 60) + 1e-8)
        ),
        
        # Otherwise assume HHMM was entered directly
        sprintf("%04d", as.integer(round(x)))
      )
    ))
  }
  
  # If already character
  x <- str_trim(as.character(x))
  
  case_when(
    is.na(x) | x == "" ~ NA_character_,
    
    # Convert "13:46" or "13:46:00" to "1346"
    str_detect(x, "^\\d{1,2}:\\d{2}(:\\d{2})?$") ~
      str_replace_all(str_sub(x, 1, 5), ":", ""),
    
    # Convert 830 to 0830, while leaving 1346 as 1346
    str_detect(x, "^\\d{1,4}$") ~
      str_pad(x, width = 4, side = "left", pad = "0"),
    
    TRUE ~ x
  )
}


sheet_id <- "129OyWLuE64lXgXDu7VDU1EHLWpjZ0C5bverqpEhO9S8"

# gid of the specific tab y

tab_gid <- 380545784

tab_name <- sheet_properties(sheet_id) %>%
  filter(id == tab_gid) %>%
  pull(name)

tab_name

################################################################################
#Inspect

header <- read_sheet(
  sheet_id,
  sheet = tab_name,
  range = "A5:AB5",
  col_names = FALSE
) %>%
  unlist(use.names = FALSE) %>%
  janitor::make_clean_names()

header

################################################################################
#read dta starting at row 7

fronds <- read_sheet(
  sheet_id,
  sheet = tab_name,
  range = "A7:AB",
  col_names = header
) %>%
  clean_names()


################################################################################
# Format time as 24-hour HHMM

fronds <- fronds %>%
  mutate(
    time = format_time_24(time)
  )


################################################################################
#drop shorcut

fronds <- fronds %>%
  select(-contains("windows"))

################################################################################
#Convert to long format

fronds_long <- fronds %>%
  
  # Make F1-F12 character because some contain "*"
  
  mutate(
    across(
      matches("^f([1-9]|1[0-2])$"),
      as.character
    )
  ) %>%
  
  # Convert F1-F12 from columns to rows
  
  # Empty/NA cells do not generate rows
  
  pivot_longer(
    cols = matches("^f([1-9]|1[0-2])$"),
    names_to = "frond",
    values_to = "frond_length",
    values_drop_na = TRUE
  ) %>%
  
  # Remove anything else representing an empty cell
  
  filter(
    !is.na(frond_length),
    str_trim(frond_length) != "",
    str_trim(frond_length) != "NULL"
  ) %>%
  
  # Create variables
  
  mutate(
    damaged_stipe = if_else(
      str_detect(frond_length, fixed("*")),
      "yes",
      "no"
    ),
    
    frond_length = parse_number(frond_length),
    
    frond = parse_number(frond)
  )

################################################################################
#Inspect

glimpse(fronds_long)

#View(fronds_long)

################################################################################
#write back to original

sheet_write(
  fronds_long,
  ss = sheet_id,
  sheet = "fronds_pull_down_long"
)
