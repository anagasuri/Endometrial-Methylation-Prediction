############################################################
# Create Table 1: Cohort metadata / sample characteristics
#
# Table sections:
#   1. Sample size
#   2. Original menstrual cycle phase annotations
#   3. Menstrual cycle phase grouping used for binary analysis
#   4. Grouped endometriosis stage among cases
#   5. Institute for analysis
#   6. Batch
#
# Uses columns available in:
# SH-Data Annotation_07.07.20_1.csv
############################################################

# Install if needed:
# install.packages(c("data.table", "dplyr", "gridExtra"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

# ==========================================================
# Input files
# ==========================================================

ANNOT_CSV <- "SH-Data Annotation_07.07.20_1.csv"
KEEP_IDS_TXT <- "R01_Study_IDs_to_include_220622.txt"

# ==========================================================
# Output directory
# ==========================================================

OUT_DIR <- "/mnt/efs/home/ubuntu/Table1_demographics"

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

# ==========================================================
# Load annotation data
# ==========================================================

annot <- fread(
  ANNOT_CSV,
  data.table = FALSE
)

# Remove accidental whitespace from column names
colnames(annot) <- trimws(colnames(annot))

# Load sample IDs to retain
keep_ids <- readLines(
  KEEP_IDS_TXT,
  warn = FALSE
)

keep_ids <- trimws(keep_ids)

# Remove blank lines and possible header
keep_ids <- keep_ids[
  keep_ids != "" &
    keep_ids != "Study.ID" &
    keep_ids != "Study ID"
]

# Standardize Study ID formatting
annot$`Study ID` <- trimws(
  as.character(annot$`Study ID`)
)

# Retain only samples included in the analysis
annot_sub <- annot %>%
  filter(`Study ID` %in% keep_ids)

# ==========================================================
# Clean variables and create analysis groupings
# ==========================================================

table1_df <- annot_sub %>%
  mutate(
    # ------------------------------------------------------
    # Endometriosis case-control status
    # ------------------------------------------------------
    Endometriosis_raw = trimws(
      as.character(`Endometriosis (Yes/No)`)
    ),
    
    Endometriosis = factor(
      Endometriosis_raw,
      levels = c("No", "Yes"),
      labels = c("Control", "Endometriosis")
    ),
    
    # ------------------------------------------------------
    # Original menstrual cycle phase annotation
    # ------------------------------------------------------
    Cycle_phase_raw = toupper(
      trimws(
        as.character(`Cycle phase for Analysis`)
      )
    ),
    
    Cycle_phase_original = case_when(
      Cycle_phase_raw == "PE" ~ "PE",
      Cycle_phase_raw == "ESE" ~ "ESE",
      Cycle_phase_raw == "MSE" ~ "MSE",
      Cycle_phase_raw == "LSE" ~ "LSE",
      Cycle_phase_raw == "SE" ~ "SE",
      Cycle_phase_raw == "MENSTRUAL" ~ "Menstrual",
      TRUE ~ NA_character_
    ),
    
    Cycle_phase_original = factor(
      Cycle_phase_original,
      levels = c(
        "PE",
        "ESE",
        "MSE",
        "LSE",
        "SE",
        "Menstrual"
      )
    ),
    
    # ------------------------------------------------------
    # Binary menstrual cycle phase grouping used in modeling
    #
    # PE                    -> Proliferative
    # ESE, MSE, LSE, and SE -> Secretory
    # Menstrual             -> Excluded
    # ------------------------------------------------------
    Cycle_phase_binary = case_when(
      Cycle_phase_raw == "PE" ~ "Proliferative",
      
      Cycle_phase_raw %in% c(
        "ESE",
        "MSE",
        "LSE",
        "SE"
      ) ~ "Secretory",
      
      Cycle_phase_raw == "MENSTRUAL" ~
        "Excluded (menstrual)",
      
      TRUE ~ NA_character_
    ),
    
    Cycle_phase_binary = factor(
      Cycle_phase_binary,
      levels = c(
        "Proliferative",
        "Secretory",
        "Excluded (menstrual)"
      )
    ),
    
    # ------------------------------------------------------
    # Grouped endometriosis stage
    # ------------------------------------------------------
    Stage_group_raw = toupper(
      trimws(
        as.character(
          `Endometriosis stage grouped (I-II), (III-IV)`
        )
      )
    ),
    
    Endometriosis_stage_grouped = case_when(
      # Controls do not have an endometriosis stage
      Endometriosis_raw == "No" ~ NA_character_,
      
      Stage_group_raw %in% c(
        "I-II",
        "I/II",
        "STAGE I-II",
        "STAGE I/II"
      ) ~ "Stage I/II",
      
      Stage_group_raw %in% c(
        "III-IV",
        "III/IV",
        "STAGE III-IV",
        "STAGE III/IV"
      ) ~ "Stage III/IV",
      
      Stage_group_raw %in% c(
        "UNKNOWN",
        "UNK",
        "NA",
        "N/A",
        ""
      ) ~ "Unknown",
      
      # Any case with a missing grouped stage is unknown
      Endometriosis_raw == "Yes" &
        is.na(Stage_group_raw) ~ "Unknown",
      
      TRUE ~ NA_character_
    ),
    
    Endometriosis_stage_grouped = factor(
      Endometriosis_stage_grouped,
      levels = c(
        "Stage I/II",
        "Stage III/IV",
        "Unknown"
      )
    ),
    
    # ------------------------------------------------------
    # Institute and batch
    # ------------------------------------------------------
    Institute = factor(
      trimws(
        as.character(`Institute for Analysis`)
      )
    ),
    
    Batch = factor(
      trimws(
        as.character(Batch)
      )
    )
  ) %>%
  select(
    Endometriosis,
    Cycle_phase_original,
    Cycle_phase_binary,
    Endometriosis_stage_grouped,
    Institute,
    Batch
  )

# ==========================================================
# Verify expected total sample count
# ==========================================================

if (nrow(table1_df) != 984) {
  warning(
    paste0(
      "Expected 984 samples, but found ",
      nrow(table1_df),
      ". Check the Study ID filtering."
    )
  )
}

# ==========================================================
# Helper function: format n (%)
# ==========================================================

n_pct <- function(x, level) {
  
  # Denominator is all nonmissing values in the supplied group
  total <- sum(!is.na(x))
  
  # Number belonging to the requested level
  n <- sum(
    as.character(x) == level,
    na.rm = TRUE
  )
  
  if (total == 0) {
    return("0 (NA%)")
  }
  
  sprintf(
    "%d (%.1f%%)",
    n,
    100 * n / total
  )
}

# ==========================================================
# Build Table 1
# ==========================================================

rows <- list()

# ----------------------------------------------------------
# Add sample-size row
# ----------------------------------------------------------

add_sample_size <- function() {
  
  rows[[length(rows) + 1]] <<- data.frame(
    Characteristic = "Sample size",
    
    Overall = as.character(
      nrow(table1_df)
    ),
    
    Controls = as.character(
      sum(
        table1_df$Endometriosis == "Control",
        na.rm = TRUE
      )
    ),
    
    Endometriosis = as.character(
      sum(
        table1_df$Endometriosis == "Endometriosis",
        na.rm = TRUE
      )
    ),
    
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------
# Add categorical section
# ----------------------------------------------------------

add_categorical <- function(
    label,
    var,
    case_only = FALSE
) {
  
  # Add section heading
  rows[[length(rows) + 1]] <<- data.frame(
    Characteristic = label,
    Overall = "",
    Controls = "",
    Endometriosis = "",
    stringsAsFactors = FALSE
  )
  
  x <- table1_df[[var]]
  
  # Preserve specified factor-level order
  if (is.factor(x)) {
    levels_var <- levels(x)
  } else {
    levels_var <- sort(
      unique(
        as.character(
          x[!is.na(x)]
        )
      )
    )
  }
  
  # Remove factor levels that do not appear in the data
  levels_var <- levels_var[
    levels_var %in% as.character(x)
  ]
  
  for (lvl in levels_var) {
    
    control_values <- x[
      table1_df$Endometriosis == "Control"
    ]
    
    case_values <- x[
      table1_df$Endometriosis == "Endometriosis"
    ]
    
    rows[[length(rows) + 1]] <<- data.frame(
      Characteristic = paste0(
        "  ",
        lvl
      ),
      
      Overall = n_pct(
        x,
        lvl
      ),
      
      Controls = if (case_only) {
        "—"
      } else {
        n_pct(
          control_values,
          lvl
        )
      },
      
      Endometriosis = n_pct(
        case_values,
        lvl
      ),
      
      stringsAsFactors = FALSE
    )
  }
}

# ----------------------------------------------------------
# Assemble sections
# ----------------------------------------------------------

add_sample_size()

add_categorical(
  label = "Menstrual cycle phase (original annotation)",
  var = "Cycle_phase_original"
)

add_categorical(
  label = "Menstrual cycle phase (binary analysis)",
  var = "Cycle_phase_binary"
)

add_categorical(
  label = "Endometriosis stage (cases only)",
  var = "Endometriosis_stage_grouped",
  case_only = TRUE
)

add_categorical(
  label = "Institute for analysis",
  var = "Institute"
)

add_categorical(
  label = "Batch",
  var = "Batch"
)

# Combine all table rows
table1_summary <- bind_rows(rows)

# ==========================================================
# Create publishable table theme
# ==========================================================

table_theme <- ttheme_default(
  base_size = 10,
  
  core = list(
    fg_params = list(
      hjust = 0,
      x = 0.03
    ),
    
    bg_params = list(
      fill = "white",
      col = "black",
      lwd = 0.8
    )
  ),
  
  colhead = list(
    fg_params = list(
      fontface = "bold",
      hjust = 0.5
    ),
    
    bg_params = list(
      fill = "grey90",
      col = "black",
      lwd = 1.2
    )
  )
)

# ==========================================================
# Create table grob
# ==========================================================

table_grob <- tableGrob(
  table1_summary,
  rows = NULL,
  theme = table_theme
)

# ==========================================================
# Save publishable PNG
# ==========================================================

png(
  filename = file.path(
    OUT_DIR,
    "Table1_cohort_metadata_publishable.png"
  ),
  
  width = 2400,
  height = 2800,
  res = 300
)

grid.newpage()
grid.draw(table_grob)

dev.off()

# ==========================================================
# Save summary table as CSV
# ==========================================================

write.csv(
  table1_summary,
  
  file = file.path(
    OUT_DIR,
    "Table1_cohort_metadata_summary.csv"
  ),
  
  row.names = FALSE
)

# ==========================================================
# Save cleaned input data as CSV
# ==========================================================

write.csv(
  table1_df,
  
  file = file.path(
    OUT_DIR,
    "Table1_cohort_metadata_input_data.csv"
  ),
  
  row.names = FALSE
)

# ==========================================================
# Print validation checks
# ==========================================================

cat("\n========================================\n")
cat("TOTAL SAMPLES INCLUDED\n")
cat("========================================\n")

print(
  nrow(table1_df)
)

cat("\n========================================\n")
cat("CASE-CONTROL COUNTS\n")
cat("========================================\n")

print(
  table(
    table1_df$Endometriosis,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("ORIGINAL CYCLE PHASE COUNTS\n")
cat("========================================\n")

print(
  table(
    table1_df$Cycle_phase_original,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("BINARY CYCLE PHASE COUNTS\n")
cat("========================================\n")

print(
  table(
    table1_df$Cycle_phase_binary,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("GROUPED ENDOMETRIOSIS STAGE COUNTS\n")
cat("========================================\n")

print(
  table(
    table1_df$Endometriosis_stage_grouped,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("INSTITUTE COUNTS\n")
cat("========================================\n")

print(
  table(
    table1_df$Institute,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("BATCH COUNTS\n")
cat("========================================\n")

print(
  table(
    table1_df$Batch,
    useNA = "ifany"
  )
)

# ==========================================================
# Check expected key counts
# ==========================================================

expected_counts <- c(
  Total = 984,
  Controls = 347,
  Endometriosis = 637,
  Proliferative = 473,
  Secretory = 461,
  Menstrual_excluded = 50,
  Stage_I_II = 344,
  Stage_III_IV = 286,
  Stage_unknown = 7
)

observed_counts <- c(
  Total = nrow(table1_df),
  
  Controls = sum(
    table1_df$Endometriosis == "Control",
    na.rm = TRUE
  ),
  
  Endometriosis = sum(
    table1_df$Endometriosis == "Endometriosis",
    na.rm = TRUE
  ),
  
  Proliferative = sum(
    table1_df$Cycle_phase_binary == "Proliferative",
    na.rm = TRUE
  ),
  
  Secretory = sum(
    table1_df$Cycle_phase_binary == "Secretory",
    na.rm = TRUE
  ),
  
  Menstrual_excluded = sum(
    table1_df$Cycle_phase_binary == "Excluded (menstrual)",
    na.rm = TRUE
  ),
  
  Stage_I_II = sum(
    table1_df$Endometriosis_stage_grouped == "Stage I/II",
    na.rm = TRUE
  ),
  
  Stage_III_IV = sum(
    table1_df$Endometriosis_stage_grouped == "Stage III/IV",
    na.rm = TRUE
  ),
  
  Stage_unknown = sum(
    table1_df$Endometriosis_stage_grouped == "Unknown",
    na.rm = TRUE
  )
)

count_check <- data.frame(
  Category = names(expected_counts),
  Expected = as.integer(expected_counts),
  Observed = as.integer(observed_counts),
  Match = expected_counts == observed_counts,
  row.names = NULL
)

cat("\n========================================\n")
cat("EXPECTED VERSUS OBSERVED COUNTS\n")
cat("========================================\n")

print(count_check)

if (!all(count_check$Match)) {
  warning(
    paste(
      "One or more observed counts do not match",
      "the expected manuscript counts.",
      "Review the printed count-check table."
    )
  )
}

# ==========================================================
# Print output information
# ==========================================================

cat("\n========================================\n")
cat("FILES WRITTEN TO\n")
cat("========================================\n")

print(OUT_DIR)

print(
  list.files(
    OUT_DIR,
    full.names = TRUE
  )
)

cat("\n========================================\n")
cat("PREVIEW OF TABLE 1\n")
cat("========================================\n")

print(
  table1_summary,
  row.names = FALSE
)