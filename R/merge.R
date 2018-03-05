# ------------------------------------------------------------------
# Merging the metabolite data across ages  
# ------------------------------------------------------------------

## Run this script as part of setup.R

stopifnot(exists("datafile_metabs"))

# Add in a unique identifier
datafile_metabs <- datafile_metabs %>%
	mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
	dplyr::select(u_ID, everything())

# ------------------------------------------------------------------
# Alteration of metabolite names 
# ------------------------------------------------------------------

## NB: This step just changes the metabolite names in the data.frame so they are 
##     more readable in tables and graphs

# Name metabolites with appropriate labels in the dataset
fine_mets <- c("DAG", "PC", "ApoA1", "ApoB", "FALen", "UnsatDeg", "DHA", "LA", "CLA", "FAw3", "FAw6", "PUFA", "MUFA", "SFA", "Glc", "Lac", "Pyr", "Cit", "Ala", "Gln", "His", "Ile", "Leu", "Val", "Phe", "Tyr", "Ace", "AcAce", "bOHBut", "Crea", "Alb", "Gp")
# Name the lipoproteins that don't have a size suffix
diff_lipos <- c("VLDLD", "LDLD", "HDLD", "VLDLC", "LDLC", "HDLC", "HDL2C", "HDL3C", "VLDLTG", "LDLTG", "HDLTG")
new_diff_lipos <- paste(
  str_extract(diff_lipos, "V?[HIL]DL"),
  str_sub(str_extract(diff_lipos, "DL.+"), 3, 10),
  sep = "-"
  )
# Put together a vector of how the wrongly named non-lipo-protein metabolites should be named - in order of the metabolites in the column names of the df
new_non_lipos <- c("Serum-C", "Remnant-C", "Est-C", "Free-C", "Serum-TG", "DAG/TG", "Tot-PG", "TG/PG", "Tot-Cho", "ApoB/ApoA1", "Tot-FA", "DHA/FA", "LA/FA", "CLA/FA", "FAw3/FA", "FAw6/FA", "PUFA/FA", "MUFA/FA", "SFA/FA")
df_list <- list()
# F7 = age 7
# TF3 = age 15
# TF4 = age 17
suffixes <- c("_F7", "_TF3", "_TF4")
# Make a loop to extract data at all 3 ages and change the metabolite names
for (i in suffixes) {
  df <- dplyr::select(datafile_metabs, u_ID, matches(paste(i, "$", sep = ""), ignore.case = FALSE))
  colnames(df) <- str_replace(colnames(df), i, "")
  lipos <- colnames(df)[matches("[HL]DL", vars = colnames(df))] %>%
    .[!(. %in% diff_lipos)]
  sizes <- str_c(c("X+L", "L", "M" ,"S", "X+S"), collapse = "|")
  new_lipos <- paste(
    str_to_lower(str_extract(lipos, sizes)), #First part
    str_extract(lipos, "V?[HL]DL"), #Second part
    str_sub(str_extract(lipos, "DL.+"), 3, 10), #Third part
    sep = "-"
    )
  IDLs <- colnames(df)[matches("IDL", vars = colnames(df))]
  new_IDLs <- paste(
    str_extract(IDLs, "IDL"), #Second part
    str_sub(str_extract(IDLs, "DL.+"), 3, 10), #Third part
    sep = "-"
    ) 
  colnames(df)[colnames(df) %in% diff_lipos] <- new_diff_lipos
  colnames(df)[colnames(df) %in% lipos] <- new_lipos
  colnames(df)[colnames(df) %in% IDLs] <- new_IDLs
  non_lipos <- colnames(df)[!(colnames(df) %in% new_lipos | colnames(df) %in% new_diff_lipos | colnames(df) %in% new_IDLs)] %>%
    .[-1] %>%
    .[!(. %in% fine_mets)]

  colnames(df)[colnames(df) %in% non_lipos] <- new_non_lipos

  df_list[[i]] <- df
}

# ------------------------------------------------------------------
# Merging 3 age groups
# ------------------------------------------------------------------

# Remove observations with missing values
df_main <- df_list[[1]] %>%
  .[complete.cases(.), ]
# Add in age
df_main[["age"]] <- 7

# df_tf3_gluc <- df_list[[2]] %>%
#   select(u_ID, glucose, insulin)

# df_main <- left_join(df_main, df_tf3_gluc)

df_tf3 <- df_list[[2]] %>%
  .[complete.cases(.), ] %>%
  .[!(.[[1]] %in% df_main[[1]]), ]

df_tf3[["age"]] <- 15

df_main <- rbind(df_main, df_tf3)

df_tf4 <- df_list[[3]] %>%
  .[complete.cases(.), ] %>%
  .[!(.[[1]] %in% df_main[[1]]), ]

# Add in glucose and insulin to age 17 as missing - was only measured at age 15
# df_tf4[["glucose"]] <- NA
# df_tf4[["insulin"]] <- NA
df_tf4[["age"]] <- 17
df_main <- rbind(df_main, df_tf4)
# Test to see if an individual is represented more than once in the dataset
stopifnot(length(unique(df_main[[1]])) == nrow(df_main))

HDL_sub <- c("HDL-2C", "HDL-3C")
colnames(df_main)[colnames(df_main) %in% HDL_sub] <- c("HDL2-C", "HDL3-C")

# Extract metabolite names and metabolites that aren't ratios
mnames <- colnames(df_main)[-c(1, ncol(df_main))]

nr_mnames <- mnames[-grep("_P|/", mnames)]

df_main_not_transformed <- df_main
# rank normalise the metabolites
for (i in mnames) {
  df_main[[i]] <- as.numeric(df_main[[i]])
  df_main[[i]] <- rntransform(formula = df_main[[i]], data = df_main, family = gaussian)
}


print("Metabolite names changed and data at 3 ages merged")



