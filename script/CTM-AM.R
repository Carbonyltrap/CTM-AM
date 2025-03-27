library(shiny)
library(DT)
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(readxl)
library(writexl)
library(stringr)
library(shinyWidgets)
library(bslib)
library(openxlsx)


options(shiny.maxRequestSize = 500 * 1024^2)

ui <- navbarPage(
  title = "CTM-AM",
  
  # ---- Library Tab ----
  tabPanel("Library",
           sidebarLayout(
             sidebarPanel(
               selectInput("carbonyl", "Select Carbonyl Compound", 
                           choices = c("MGO", "Glyoxal", "3-Deoxyglucosone"),
                           selected = "MGO"),
               checkboxGroupInput("modifications", "Include Modifications",
                                  choices = c("-H2" = "H2", "-H2O" = "H2O"),
                                  selected = c("H2", "H2O")),
               sliderInput("max_n", "Max n for Combination", min = 1, max = 10, value = 3),
               sliderInput("ppm", "PPM Error", min = 1, max = 100, value = 50),
               actionButton("generate", "Generate Library"),
               downloadButton("download_csv", "Download CSV")
             ),
             mainPanel(
               h4("Generated Combination Library"),
               DT::dataTableOutput("lib_table")
             )
           )
  ),
  
  # ---- Preprocessing Tab ----
  tabPanel("Preprocessing",
           sidebarLayout(
             sidebarPanel(
               fileInput("ms_file", "Upload MS Data"),
               actionButton("preprocess", "Run Preprocessing"),
               hr(),
               sliderInput("mz_tolerance", "m/z Tolerance (ppm)", min = 1, max = 20, value = 5),
               sliderInput("rt_tolerance", "Retention Time Tolerance (min)", min = 0.1, max = 2, value = 0.5, step = 0.1),
               checkboxInput("enable_classification", "Enable S/P Classification", value = TRUE),
               actionButton("group_classify", "Run Grouping and Classification")
             ),
             mainPanel(
               uiOutput("file_selector_ui"),
               hr(),
               DT::dataTableOutput("pre_table")
             )
           )
  ),
  
  # ---- Mining Tab ----
  tabPanel("Mining",
           sidebarLayout(
             sidebarPanel(
               fileInput("database_file_input", "Upload Carbonyl Combination Results File (CSV)"),
               sliderInput("rt_range", "Retention Time Range (min)", min = 0, max = 20, value = c(1, 19), step = 0.5),
               actionButton("run_mining1", "Run Mining - MS1 Level"),
               hr(),
               fileInput("mgf_file", "Select MGF File"),
               numericInput("noise_threshold", "Noise Threshold", value = 1000),
               sliderInput("ppm_tolerance", "PPM Tolerance", min = 1, max = 50, value = 5, step = 1),
               
               checkboxGroupInput("similarity_metrics", "Similarity Metrics",
                                  choices = c("Jaccard Similarity" = "jaccard",
                                              "Weighted Cosine Similarity" = "cosine"),
                                  selected = c("jaccard", "cosine")),
               
               actionButton("run_mining2", "Run Mining - MS2 Level")
             ),
             mainPanel(
               uiOutput("ms2_file_selector_ui"),
               DT::dataTableOutput("mining_table")
             )
           )
  ),
  
  # ---- Annotation Tab ----
  # Annotation Tab
  tabPanel("Annotation",
           sidebarLayout(
             sidebarPanel(
               h4("MGF Matching"),
               fileInput("mgf_file_annotation", "Upload MGF File (.mgf)", accept = c(".mgf")),
               radioButtons("match_target", "Select Target for MGF Matching:",
                            choices = c("Candidate Compounds" = "S", "Candidate Products" = "P"), selected = "S"),
               actionButton("run_mgf_match", "Match MGF"),
               downloadButton("download_matched_mgf", "Download Matched MGF"),
               
               hr(),
               h4("Annotation Processing"),
               
               fileInput("structure_csv", "Upload Structure Identifications (CSV)", accept = ".csv"),
               fileInput("canopus_csv", "Upload Canopus Summary (CSV)", accept = ".csv"),
               fileInput("formula_csv", "Upload Formula Identifications (CSV)", accept = ".csv"),
               
               checkboxInput("enable_structure_id", "Structure Identification", value = TRUE),
               checkboxInput("enable_confidence_filter", "Molecular Formula Matching", value = TRUE),
               actionButton("run_annotation", "Submit Annotation Processing"),
               hr(),
               uiOutput("annot_file_selector_ui")
             ),
             mainPanel(
               h4("Annotation Results"),
               uiOutput("annot_file_selector_ui"), 
               DT::dataTableOutput("annot_table")
             )
           )
  )
  
)

server <- function(input, output, session) {
  
  # ---- Library Tab Logic ----
  observeEvent(input$generate, {
    S <- 0
    compound_mass <- switch(input$carbonyl,
                            "MGO" = 72.0210,
                            "Glyoxal" = 58.0055,
                            "3-Deoxyglucosone" = 162.0528)
    H2_mass <- 2.0157
    H2O_mass <- 18.0106
    
    include_H2 <- "H2" %in% input$modifications
    include_H2O <- "H2O" %in% input$modifications
    
    calculate_ppm_range <- function(value, ppm) {
      error <- value * ppm / 1e6
      c(value - error, value + error)
    }
    
    results <- data.frame(
      n = integer(),
      value_type = character(),
      value = double(),
      lower_limit = double(),
      upper_limit = double(),
      stringsAsFactors = FALSE
    )
    
    generate_combinations <- function(n, current_n, i_count, j_count) {
      if (current_n == n) {
        value <- S + n * compound_mass - i_count * H2_mass - j_count * H2O_mass
        value_type <- paste0("S+", n, "*", input$carbonyl, 
                             if (include_H2 && i_count > 0) paste0("-", i_count, "*H2") else "",
                             if (include_H2O && j_count > 0) paste0("-", j_count, "*H2O") else "")
        range <- calculate_ppm_range(value, input$ppm)
        results <<- rbind(results, data.frame(
          n = n,
          value_type = value_type,
          value = round(value, 4),
          lower_limit = round(range[1], 4),
          upper_limit = round(range[2], 4),
          stringsAsFactors = FALSE
        ))
      } else {
        if (include_H2) generate_combinations(n, current_n + 1, i_count + 1, j_count)
        if (include_H2O) generate_combinations(n, current_n + 1, i_count, j_count + 1)
        generate_combinations(n, current_n + 1, i_count, j_count)
      }
    }
    
    for (n in 1:input$max_n) {
      generate_combinations(n, 0, 0, 0)
    }
    
    results <- results[!duplicated(results$value), ]
    
    output$lib_table <- DT::renderDataTable({
      DT::datatable(results,
                    options = list(pageLength = 10, autoWidth = TRUE),
                    class = 'cell-border stripe hover compact',
                    rownames = FALSE)
    })
    
    output$download_csv <- downloadHandler(
      filename = function() {
        paste0("Library_", input$carbonyl, "_n", input$max_n, "_ppm", input$ppm, ".csv")
      },
      content = function(file) {
        write.csv(results, file, row.names = FALSE)
      }
    )
  })
  
  # ---- Preprocessing Tab Logic ----
  observeEvent(input$preprocess, {
    req(input$ms_file)
    
    withProgress(message = 'Preprocessing uploaded file...', value = 0, {
      df <- read.csv(input$ms_file$datapath)
      process_data <- function(df) {
        col_names <- colnames(df)
        base_cols <- c("row.ID", "row.m.z", "row.retention.time")
        e_cols <- grep("Ext_E", col_names, value = TRUE)
        if (length(e_cols) == 0) return(list())
        e_prefixes <- unique(gsub("_(MGO_)?24h_ddms30.mzML.Peak.area", "", e_cols))
        processed_data_list <- list()
        for (prefix in e_prefixes) {
          e_group_cols <- grep(paste0("^", prefix, "_"), col_names, value = TRUE)
          if (length(e_group_cols) > 1) {
            mgo_col <- grep("MGO", e_group_cols, value = TRUE)
            non_mgo_cols <- setdiff(e_group_cols, mgo_col)
            sorted_group_cols <- c(non_mgo_cols, mgo_col)
            combined_df <- df[, c(base_cols, sorted_group_cols)]
            row_has_nonzero <- rowSums(combined_df[, -c(1:3)] != 0) > 0
            filtered_df <- combined_df[row_has_nonzero, ]
            if (nrow(filtered_df) == 0) next
            sorted_df <- filtered_df[order(filtered_df$`row.m.z`), ]
            processed_data_list[[prefix]] <- sorted_df
          }
        }
        return(processed_data_list)
      }
      processed_data_list <- process_data(df)
      output_folder <- "output_data/2_1_Data_Ext_MZminedata_separate_SPname"
      dir_create(output_folder)
      for (name in names(processed_data_list)) {
        final_data <- processed_data_list[[name]]
        standardized_name <- sub("Ext2", "Ext", name)
        output_path <- file.path(output_folder, paste0(standardized_name, "_processed_Ext_Mzminedata_quant.csv"))
        write.csv(final_data, output_path, row.names = FALSE)
      }
      incProgress(1)
    })
  })
  
  observeEvent(input$group_classify, {
    withProgress(message = 'Processing and grouping data...', value = 0, {
      input_folder <- "output_data/2_1_Data_Ext_MZminedata_separate_SPname"
      output_folder <- "output_data/2_2_Data_Ext_MZminedata_separate_SPname"
      dir_create(output_folder)
      calculate_ppm <- function(observed, theoretical) {
        abs(observed - theoretical) / theoretical * 1e6
      }
      process_grouped_data <- function(df, mz_tol, rt_tol, classify) {
        mgo_cols <- grep("_MGO_", colnames(df), value = TRUE)
        non_mgo_cols <- setdiff(grep("Ext_E", colnames(df), value = TRUE), mgo_cols)
        df <- df %>%
          mutate(row.m.z = as.numeric(row.m.z), row.retention.time = as.numeric(row.retention.time)) %>%
          arrange(row.m.z, row.retention.time) %>%
          mutate(
            prev_mz = lag(row.m.z),
            prev_rt = lag(row.retention.time),
            mz_diff = calculate_ppm(row.m.z, prev_mz),
            rt_diff = abs(row.retention.time - prev_rt),
            Peakgroup_id = cumsum(is.na(prev_mz) | mz_diff > mz_tol | rt_diff > rt_tol)
          ) %>%
          select(-prev_mz, -prev_rt, -mz_diff, -rt_diff)
        if (classify) {
          group_summary <- df %>%
            group_by(Peakgroup_id) %>%
            summarise(
              Total_NonMGO = sum(across(all_of(non_mgo_cols)), na.rm = TRUE),
              Total_MGO = sum(across(all_of(mgo_cols)), na.rm = TRUE),
              .groups = "drop"
            ) %>%
            mutate(
              Label = case_when(
                (Total_NonMGO == 0 & Total_MGO > 0) ~ "P",
                ((Total_NonMGO - Total_MGO) > (0.05 * Total_NonMGO) & Total_NonMGO > 0) ~ "S",
                TRUE ~ ""
              )
            )
          df <- df %>%
            left_join(group_summary, by = "Peakgroup_id") %>%
            select(-Total_NonMGO, -Total_MGO)
        } else {
          df$Label <- ""
        }
        return(df)
      }
      csv_files <- dir_ls(input_folder, regexp = "\\.csv$")
      for (file_path in csv_files) {
        df <- read_csv(file_path, col_types = cols(row.m.z = col_double(), row.retention.time = col_double()))
        final_data <- process_grouped_data(df, input$mz_tolerance, input$rt_tolerance, input$enable_classification)
        write_csv(final_data, path(output_folder, basename(file_path)))
        incProgress(1 / length(csv_files))
      }
      processed_files <- dir_ls(output_folder, regexp = "\\.csv$")
      output$file_selector_ui <- renderUI({
        tagList(
          tags$label("Select File to View Result", style = "font-weight:bold;"),
          selectInput("selected_file", NULL, choices = basename(processed_files), width = '100%')
        )
      })
    })
  })
  
  observeEvent(input$selected_file, {
    req(input$selected_file)
    file_path <- file.path("output_data/2_2_Data_Ext_MZminedata_separate_SPname", input$selected_file)
    df <- read_csv(file_path)
    output$pre_table <- DT::renderDataTable({
      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE),
                    class = 'cell-border stripe hover compact', rownames = FALSE)
    })
  })
  
  # ---- Mining - MS1 Level ----
  observeEvent(input$run_mining1, {
    req(input$database_file_input)
    withProgress(message = "Running Mining - MS1 Level...", value = 0, {
      data_folder <- "./output_data/2_2_Data_Ext_MZminedata_separate_SPname/"
      output_folder <- "./output_data/3_1_Filter_Ext_MGO_MS1_SPmatch_mzRT/"
      dir_create(output_folder)
      
      database_path <- input$database_file_input$datapath
      database <- read.csv(database_path)
      database$lower_limit <- as.numeric(database$lower_limit)
      database$upper_limit <- as.numeric(database$upper_limit)
      
      data_files <- list.files(data_folder, pattern = "\\.csv$", full.names = TRUE)
      total_files <- length(data_files)
      progress_per_file <- 1 / total_files
      
      for (data_file in data_files) {
        tryCatch({
          E_number <- gsub(".*_E(\\d+)_.*", "\\1", data_file)
          
          S_peak_area_col1 <- paste0("X20240730_Ext_E", E_number, "_24h_ddms30.mzML.Peak.area")
          S_peak_area_col2 <- paste0("X20240730_Ext_E", E_number, "_MGO_24h_ddms30.mzML.Peak.area")
          
          data <- read.csv(data_file)
          if (nrow(data) == 0) next
          
          if (!(S_peak_area_col1 %in% names(data)) || !(S_peak_area_col2 %in% names(data))) next
          
          S_data <- data[data$Label == "S", ]
          P_data <- data[data$Label == "P", ]
          
          if (nrow(S_data) == 0 | nrow(P_data) == 0) next
          
          filtered_pairs <- data.frame(
            Peakgroup_id = integer(),
            S_row.ID = integer(),
            S_value = numeric(),
            S_retention_time = numeric(),
            S_Peak.area1 = numeric(),
            S_Peak.area2 = numeric(),
            P_row.ID = integer(),
            P_value = numeric(),
            P_retention_time = numeric(),
            P_Peak.area1 = numeric(),
            P_Peak.area2 = numeric(),
            n = integer(),
            value_type = character(),
            stringsAsFactors = FALSE
          )
          
          for (i in seq_along(S_data$`row.m.z`)) {
            s_val <- S_data$`row.m.z`[i]
            s_row_ID <- S_data$row.ID[i]
            s_peakgroup_id <- S_data$Peakgroup_id[i]
            
            for (j in seq_along(P_data$`row.m.z`)) {
              p_val <- P_data$`row.m.z`[j]
              p_row_ID <- P_data$row.ID[j]
              
              if (p_val > s_val) {
                difference <- p_val - s_val
                within_range <- database$lower_limit <= difference & difference <= database$upper_limit
                if (any(within_range, na.rm = TRUE)) {
                  matched_row <- which(within_range)[1]
                  n_val <- database$n[matched_row]
                  value_type_val <- as.character(database$value_type[matched_row])
                  
                  S_info <- S_data[i, c("row.retention.time", S_peak_area_col1, S_peak_area_col2)]
                  P_info <- P_data[j, c("row.retention.time", S_peak_area_col1, S_peak_area_col2)]
                  
                  rt_diff <- abs(S_info[1, 1] - P_info[1, 1])
                  if (rt_diff < 6) {
                    filtered_pairs <- rbind(filtered_pairs, data.frame(
                      Peakgroup_id = s_peakgroup_id,
                      S_row.ID = s_row_ID,
                      S_value = s_val,
                      S_retention_time = S_info[1, 1],
                      S_Peak.area1 = S_info[1, 2],
                      S_Peak.area2 = S_info[1, 3],
                      P_row.ID = p_row_ID,
                      P_value = p_val,
                      P_retention_time = P_info[1, 1],
                      P_Peak.area1 = P_info[1, 2],
                      P_Peak.area2 = P_info[1, 3],
                      n = n_val,
                      value_type = value_type_val,
                      stringsAsFactors = FALSE
                    ))
                  }
                }
              }
            }
          }
          
          if (nrow(filtered_pairs) > 0) {
            filtered_pairs <- filtered_pairs %>%
              filter(S_retention_time >= input$rt_range[1], S_retention_time <= input$rt_range[2],
                     P_retention_time >= input$rt_range[1], P_retention_time <= input$rt_range[2],
                     S_Peak.area1 != 0, P_Peak.area1 == 0,
                     !(n >= 3 & P_retention_time > S_retention_time)) %>%
              mutate(
                match_label = case_when(
                  value_type == "S+1*MGO-0*H2-1*H2O" & 
                    P_retention_time < {
                      subset_times <- P_retention_time[value_type %in% c("S+1*MGO-0*H2-0*H2O", "S+1*MGO-1*H2-0*H2O")]
                      if (length(subset_times) > 0) min(subset_times, na.rm = TRUE) else Inf
                    } ~ "High confidence",
                  value_type == "S+1*MGO-0*H2-0*H2O" & P_retention_time < S_retention_time ~ "High confidence",
                  value_type == "S+1*MGO-0*H2-1*H2O" & P_retention_time < S_retention_time ~ "High confidence",
                  value_type == "S+1*MGO-1*H2-0*H2O" & S_retention_time < P_retention_time ~ "High confidence",
                  TRUE ~ "Mid confidence"
                )
              ) %>% ungroup()
            
            wb <- createWorkbook()
            addWorksheet(wb, "Sheet1")
            writeData(wb, "Sheet1", filtered_pairs)
            output_file <- file.path(output_folder, paste0("SP_E", E_number, "_X20240730_Ext_combined_mgo.xlsx"))
            saveWorkbook(wb, output_file, overwrite = TRUE)
          }
          
        }, error = function(e) {
          message("Error processing file: ", data_file, " - ", e$message)
        })
        incProgress(progress_per_file)
      }
      
      result_files <- list.files(path = output_folder, pattern = "\\.xlsx$", full.names = FALSE)
      output$result_file_selector_ui <- renderUI({
        tagList(
          tags$label("Select Result File to View", style = "font-weight:bold;"),
          selectInput("selected_result_file", NULL, choices = result_files, width = '100%')
        )
      })
    })
  })
  
  observeEvent(input$selected_result_file, {
    req(input$selected_result_file)
    file_path <- file.path("./output_data/3_1_Filter_Ext_MGO_MS1_SPmatch_mzRT/", input$selected_result_file)
    tryCatch({
      df <- read_excel(file_path, sheet = 1)
      output$mining_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
      })
    }, error = function(e) {
      showNotification(paste("Error reading file:", input$selected_result_file), type = "error")
    })
  })
  
  
  # ---- Mining - MS2 Level ----
  observeEvent(input$run_mining2, {
    req(input$mgf_file)
    
    withProgress(message = "Running Mining - MS2 Level...", value = 0, {
      library(openxlsx)
      
      folder_path <- "./output_data/3_1_Filter_Ext_MGO_MS1_SPmatch_mzRT"
      mgf_path <- input$mgf_file$datapath
      output_folder <- "./output_data/3_2_1_Filter_Ext_MGO_MS2_SPmatch_mzRT"
      filtered_output_folder <- "./output_data/3_2_2_Filter_Ext_MGO_MS2_SPmatch_mzRT_Filtered_Data"
      dir_create(output_folder)
      dir_create(filtered_output_folder)
      
      noise_threshold <- input$noise_threshold
      ppm_tolerance <- input$ppm_tolerance
      
      read_mgf_file <- function(file_path) {
        spectra <- list()
        mgf_lines <- readLines(file_path)
        current_spectrum <- NULL
        current_feature_id <- NULL
        
        for (line in mgf_lines) {
          if (grepl("BEGIN IONS", line, ignore.case = TRUE)) {
            current_spectrum <- list("m/z" = numeric(), "intensity" = numeric())
            current_feature_id <- NULL
          } else if (grepl("END IONS", line, ignore.case = TRUE)) {
            if (!is.null(current_feature_id)) {
              valid_indices <- which(current_spectrum[["intensity"]] >= noise_threshold)
              current_spectrum[["m/z"]] <- current_spectrum[["m/z"]][valid_indices]
              current_spectrum[["intensity"]] <- current_spectrum[["intensity"]][valid_indices]
              spectra[[current_feature_id]] <- current_spectrum
            }
          } else if (!is.null(current_spectrum)) {
            if (grepl("^FEATURE_ID=", line, ignore.case = TRUE)) {
              current_feature_id <- sub("^FEATURE_ID=", "", line)
            } else {
              parts <- strsplit(line, "\\s+")[[1]]
              if (length(parts) >= 2) {
                current_spectrum[["m/z"]] <- c(current_spectrum[["m/z"]], as.numeric(parts[1]))
                current_spectrum[["intensity"]] <- c(current_spectrum[["intensity"]], as.numeric(parts[2]))
              }
            }
          }
        }
        return(spectra)
      }
      
      spectra_list <- read_mgf_file(mgf_path)
      incProgress(0.1, detail = paste0("Parsed MGF - ", length(spectra_list), " spectra"))
      
      match_peaks <- function(s_mzs, s_intensities, p_mzs, p_intensities, ppm_tolerance, noise_threshold) {
        if (length(s_mzs) == 0 || length(p_mzs) == 0) {
          return(list(mz_s = NA, mz_p = NA, s_intensity = NA, p_intensity = NA))
        }
        
        valid_s_indices <- which(s_intensities >= noise_threshold)
        valid_p_indices <- which(p_intensities >= noise_threshold)
        s_mzs <- s_mzs[valid_s_indices]
        s_intensities <- s_intensities[valid_s_indices]
        p_mzs <- p_mzs[valid_p_indices]
        p_intensities <- p_intensities[valid_p_indices]
        
        matched_mz_s <- matched_mz_p <- matched_s_intensity <- matched_p_intensity <- c()
        
        for (mz in s_mzs) {
          ppm_error_limit <- (ppm_tolerance / 1e6) * mz
          p_match_idx <- which(abs(p_mzs - mz) <= ppm_error_limit)
          
          if (length(p_match_idx) > 0) {
            best_p_idx <- p_match_idx[which.min(abs(p_mzs[p_match_idx] - mz))]
            matched_mz_s <- c(matched_mz_s, mz)
            matched_mz_p <- c(matched_mz_p, p_mzs[best_p_idx])
            matched_s_intensity <- c(matched_s_intensity, s_intensities[which(s_mzs == mz)][1])
            matched_p_intensity <- c(matched_p_intensity, p_intensities[best_p_idx])
          }
        }
        
        if (length(matched_mz_s) == 0) {
          return(list(mz_s = NA, mz_p = NA, s_intensity = NA, p_intensity = NA))
        }
        
        return(list(
          mz_s = paste(sprintf("%.4f", matched_mz_s), collapse = ", "),
          mz_p = paste(sprintf("%.4f", matched_mz_p), collapse = ", "),
          s_intensity = paste(sprintf("%.2f", matched_s_intensity), collapse = ", "),
          p_intensity = paste(sprintf("%.2f", matched_p_intensity), collapse = ", ")
        ))
      }
      
      compute_jaccard_similarity <- function(s_mzs, p_mzs, tolerance) {
        if (length(s_mzs) == 0 || length(p_mzs) == 0) return(0)
        matched_peaks <- unique(s_mzs[abs(outer(s_mzs, p_mzs, "-")) <= tolerance])
        matched_count <- length(matched_peaks)
        total_unique_peaks <- length(unique(c(s_mzs, p_mzs)))
        return(ifelse(total_unique_peaks > 0, matched_count / total_unique_peaks, 0))
      }
      
      calculate_weighted_cosine_similarity <- function(mz_values, intensity_s, intensity_p) {
        if (length(mz_values) == 0 || length(intensity_s) == 0 || length(intensity_p) == 0 ||
            length(mz_values) != length(intensity_s) || length(mz_values) != length(intensity_p)) {
          return(NA)
        }
        mz_sqrt <- sqrt(mz_values)
        intensity_s_sqrt <- sqrt(intensity_s)
        intensity_p_sqrt <- sqrt(intensity_p)
        numerator <- sum(mz_sqrt * intensity_s_sqrt * mz_sqrt * intensity_p_sqrt)
        denominator_s <- sqrt(sum((mz_sqrt * intensity_s_sqrt)^2))
        denominator_p <- sqrt(sum((mz_sqrt * intensity_p_sqrt)^2))
        similarity <- numerator / (denominator_s * denominator_p)
        return(ifelse(is.na(similarity) | is.infinite(similarity), NA, similarity))
      }
      
      excel_files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
      total_files <- length(excel_files)
      
      for (i in seq_along(excel_files)) {
        excel_file_path <- excel_files[i]
        excel_data <- read_excel(excel_file_path, sheet = 1)
        
        s_ids <- as.character(excel_data$S_row.ID)
        p_ids <- as.character(excel_data$P_row.ID)
        
        matched_mz_s_list <- matched_mz_p_list <- matched_intensity_s_list <- matched_intensity_p_list <- vector("character", length(s_ids))
        jaccard_similarities <- numeric(length(s_ids))
        
        for (j in seq_along(s_ids)) {
          s_spectrum <- spectra_list[[s_ids[j]]]
          p_spectrum <- spectra_list[[p_ids[j]]]
          
          if (!is.null(s_spectrum) && !is.null(p_spectrum)) {
            matched_result <- match_peaks(s_spectrum[["m/z"]], s_spectrum[["intensity"]],
                                          p_spectrum[["m/z"]], p_spectrum[["intensity"]],
                                          ppm_tolerance, noise_threshold)
            
            matched_mz_s_list[j] <- matched_result$mz_s
            matched_mz_p_list[j] <- matched_result$mz_p
            matched_intensity_s_list[j] <- matched_result$s_intensity
            matched_intensity_p_list[j] <- matched_result$p_intensity
            
            jaccard_similarities[j] <- compute_jaccard_similarity(
              s_spectrum[["m/z"]], p_spectrum[["m/z"]],
              ppm_tolerance * mean(c(s_spectrum[["m/z"]], p_spectrum[["m/z"]])) / 1e6
            )
          }
        }
        
        excel_data$Matched_mz_S <- matched_mz_s_list
        excel_data$Matched_mz_P <- matched_mz_p_list
        excel_data$Matched_mz_Intensity_S <- matched_intensity_s_list
        excel_data$Matched_mz_Intensity_P <- matched_intensity_p_list
        excel_data$Jaccard_Similarity <- jaccard_similarities
        
        excel_data$Weighted_Cosine_Similarity <- sapply(seq_len(nrow(excel_data)), function(k) {
          mz_values_s <- as.numeric(unlist(strsplit(excel_data$Matched_mz_S[k], ", ")))
          mz_values_p <- as.numeric(unlist(strsplit(excel_data$Matched_mz_P[k], ", ")))
          intensity_s <- as.numeric(unlist(strsplit(excel_data$Matched_mz_Intensity_S[k], ", ")))
          intensity_p <- as.numeric(unlist(strsplit(excel_data$Matched_mz_Intensity_P[k], ", ")))
          if (any(is.na(mz_values_s)) || any(is.na(mz_values_p)) ||
              any(is.na(intensity_s)) || any(is.na(intensity_p))) {
            return(NA)
          }
          calculate_weighted_cosine_similarity((mz_values_s + mz_values_p) / 2, intensity_s, intensity_p)
        })
        
        write_xlsx(excel_data, file.path(output_folder, basename(excel_file_path)))
        
        filtered_data <- excel_data %>%
          filter(!is.na(Matched_mz_S) & Matched_mz_S != "" &
                   sapply(strsplit(Matched_mz_S, ", "), length) >= 3 &
                   !is.na(Weighted_Cosine_Similarity) & Weighted_Cosine_Similarity >= 0.7) %>%
          mutate(n = as.numeric(n)) %>%
          filter(!(n >= 3 & (P_retention_time >= S_retention_time | (S_retention_time - P_retention_time) < 2))) %>%
          group_by(P_value) %>%
          slice_max(order_by = Weighted_Cosine_Similarity, with_ties = FALSE) %>%  
          ungroup() %>%
          arrange(Peakgroup_id, P_value, P_retention_time)
        
        
        write_xlsx(filtered_data, file.path(filtered_output_folder, basename(excel_file_path)))
        incProgress(1 / total_files, detail = paste("Processed", i, "of", total_files))
      }
      
      ms2_result_files <- list.files(
        path = filtered_output_folder,
        pattern = "\\.xlsx$",
        full.names = FALSE
      )
      
      output$ms2_file_selector_ui <- renderUI({
        tagList(
          tags$label("Select MS2 Result File to View", style = "font-weight:bold;"),
          selectInput("selected_ms2_file", NULL, choices = ms2_result_files, width = '100%')
        )
      })
    })
  })
  

  observeEvent(input$selected_ms2_file, {
    req(input$selected_ms2_file)
    file_path <- file.path("./output_data/3_2_2_Filter_Ext_MGO_MS2_SPmatch_mzRT_Filtered_Data", input$selected_ms2_file)
    
    tryCatch({
      df <- read_excel(file_path, sheet = 1)
      output$mining_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
      })
    }, error = function(e) {
      showNotification(paste("Error reading MS2 file:", input$selected_ms2_file), type = "error")
    })
  })
  
  # ---- Annotation ----
  
  # ---- 1. MGF Matching ----
  observeEvent(input$run_mgf_match, {
    withProgress(message = "Matching MGF...", value = 0, {
      target <- input$match_target
      excel_folder_path <- "./output_data/3_2_2_Filter_Ext_MGO_MS2_SPmatch_mzRT_Filtered_Data"
      req(input$mgf_file_annotation) 
      mgf_path <- input$mgf_file_annotation$datapath
      
      
      if (target == "S") {
        id_column <- "S_row.ID"
        matched_mgf_path <- "./output_data/matched_output.mgf"
        unmatched_mgf_path <- "./output_data/unmatched_output.mgf"
      } else {
        id_column <- "P_row.ID"
        matched_mgf_path <- "./output_data/P_matched_output.mgf"
        unmatched_mgf_path <- "./output_data/P_unmatched_output.mgf"
      }
      
      row_id_list <- unique(unlist(lapply(
        list.files(excel_folder_path, pattern = "\\.xlsx$", full.names = TRUE),
        function(file) {
          tryCatch({
            excel_data <- readxl::read_excel(file, sheet = 1, col_types = "text")
            if (id_column %in% colnames(excel_data)) {
              return(as.character(excel_data[[id_column]]))
            }
          }, error = function(e) {
            return(NULL)
          })
        }
      )))
      
      mgf_lines <- readLines(mgf_path)
      matched_mgf_content <- unmatched_mgf_content <- c()
      matched_count <- unmatched_count <- 0
      capture_spectrum <- FALSE
      current_spectrum <- c()
      
      for (line in mgf_lines) {
        if (grepl("BEGIN IONS", line)) {
          capture_spectrum <- TRUE
          current_spectrum <- c(line)
        } else if (grepl("END IONS", line)) {
          current_spectrum <- c(current_spectrum, line)
          feature_id_line <- grep("FEATURE_ID=", current_spectrum, value = TRUE)
          if (length(feature_id_line) > 0) {
            feature_id <- trimws(sub("FEATURE_ID=", "", feature_id_line))
            if (feature_id %in% row_id_list) {
              matched_mgf_content <- c(matched_mgf_content, current_spectrum)
              matched_count <- matched_count + 1
            } else {
              unmatched_mgf_content <- c(unmatched_mgf_content, current_spectrum)
              unmatched_count <- unmatched_count + 1
            }
          } else {
            unmatched_mgf_content <- c(unmatched_mgf_content, current_spectrum)
            unmatched_count <- unmatched_count + 1
          }
          capture_spectrum <- FALSE
        } else if (capture_spectrum) {
          current_spectrum <- c(current_spectrum, line)
        }
      }
      
      writeLines(matched_mgf_content, matched_mgf_path)
      writeLines(unmatched_mgf_content, unmatched_mgf_path)
      
      matched_msg <- if (target == "S") {
        paste0("Candidate Compounds Matching Complete:\nMatched spectra: ", matched_count,
               " Unmatched spectra: ", unmatched_count)
      } else {
        paste0("Candidate Products Matching Complete:\nMatched spectra: ", matched_count,
               " Unmatched spectra: ", unmatched_count)
      }
      
      showNotification(matched_msg, type = "message", duration = 10)
      
      output$download_matched_mgf <- downloadHandler(
        filename = function() if (target == "S") "matched_output.mgf" else "P_matched_output.mgf",
        content = function(file) file.copy(matched_mgf_path, file)
      )
      
      output$download_unmatched_mgf <- downloadHandler(
        filename = function() if (target == "S") "unmatched_output.mgf" else "P_unmatched_output.mgf",
        content = function(file) file.copy(unmatched_mgf_path, file)
      )
      
      incProgress(1)
    })
  })
  
  # ---- 2. Annotation Processing ----
  observeEvent(input$run_annotation, {
    withProgress(message = "Running Annotation Processing...", value = 0, {

      req(input$structure_csv)
      req(input$canopus_csv)
      req(input$formula_csv)
      
      csv_file <- input$structure_csv$datapath
      target_file <- input$canopus_csv$datapath
      formula_file <- input$formula_csv$datapath
      
      folder_in <- "./output_data/3_2_2_Filter_Ext_MGO_MS2_SPmatch_mzRT_Filtered_Data"
      folder_full <- "./output_data/4_Annotation_Ext_structure_library2"
      folder_filtered <- "./output_data/4_Annotation_Ext_structure_library3"
      
      dir.create(folder_full, recursive = TRUE, showWarnings = FALSE)
      dir.create(folder_filtered, recursive = TRUE, showWarnings = FALSE)
      
      data1 <- read.csv(csv_file) %>%
        dplyr::select(mappingFeatureId, InChI, name, smiles, ConfidenceScoreExact, ConfidenceScoreApproximate) %>%
        dplyr::mutate(mappingFeatureId = as.character(mappingFeatureId))
      
      target_data <- readr::read_csv(target_file) %>%
        dplyr::mutate(mappingFeatureId = as.character(mappingFeatureId)) %>%
        dplyr::select(mappingFeatureId, precursorFormula, molecularFormula, adduct, `NPC#class`, `NPC#superclass`)
      
      formula_data <- readr::read_csv(formula_file) %>%
        dplyr::mutate(mappingFeatureId = as.character(mappingFeatureId)) %>%
        dplyr::select(mappingFeatureId, molecularFormula) %>%
        dplyr::rename(P_sirius_molecularFormula = molecularFormula)
      
      formula_lookup <- list("MGO" = "C3H4O2", "H2O" = "H2O", "H2" = "H2")
      
      parse_formula <- function(formula) {
        elements <- stringr::str_extract_all(formula, "[A-Z][a-z]?\\d*")[[1]]
        parsed <- setNames(as.numeric(gsub("[^0-9]", "", elements)), gsub("\\d", "", elements))
        parsed[is.na(parsed)] <- 1
        return(parsed)
      }
      
      extract_carbon_count <- function(formula_vector) {
        sapply(formula_vector, function(formula) {
          if (is.na(formula)) return(NA)
          match <- stringr::str_match(formula, "C(\\d*)")
          if (is.na(match[2]) || match[2] == "") return(1)
          return(as.numeric(match[2]))
        })
      }
      
      adjust_formula <- function(molecularFormula, value_type) {
        base_formula <- parse_formula(molecularFormula)
        modifications <- stringr::str_extract_all(value_type, "[+-]?\\d*\\*?[A-Za-z0-9]+")[[1]]
        for (mod in modifications) {
          sign <- ifelse(stringr::str_detect(mod, "^-"), -1, 1)
          mod <- stringr::str_replace(mod, "^[+-]", "")
          count <- as.numeric(stringr::str_extract(mod, "^\\d+"))
          if (is.na(count)) count <- 1
          molecule <- stringr::str_replace(mod, "^\\d+\\*", "")
          formula <- formula_lookup[[molecule]]
          if (!is.null(formula)) {
            parsed_mod <- parse_formula(formula)
            for (element in names(parsed_mod)) {
              base_formula[element] <- base_formula[element] + sign * count * parsed_mod[element]
            }
          }
        }
        updated_formula <- paste0(names(base_formula), ifelse(base_formula > 1, base_formula, ""), collapse = "")
        return(updated_formula)
      }
      
      files <- list.files(folder_in, pattern = "\\.xlsx$", full.names = TRUE)
      
      for (i in seq_along(files)) {
        incProgress(1 / length(files), detail = paste("Processing", basename(files[i])))
        df <- readxl::read_xlsx(files[i]) %>% dplyr::mutate(S_row.ID = as.character(S_row.ID))
        
        annotated <- df
        
        if (input$enable_structure_id) {
          annotated <- annotated %>%
            dplyr::left_join(data1, by = c("S_row.ID" = "mappingFeatureId")) %>%
            dplyr::mutate(name = if_else(is.na(name), paste0("unnamed_", S_row.ID), name),
                          ConfidenceScoreExact = dplyr::coalesce(ConfidenceScoreExact, 0),
                          ConfidenceScoreApproximate = dplyr::coalesce(ConfidenceScoreApproximate, 0)) %>%
            dplyr::left_join(target_data, by = c("S_row.ID" = "mappingFeatureId")) %>%
            dplyr::mutate(precursorFormula = dplyr::coalesce(precursorFormula, "unknown"),
                          molecularFormula = dplyr::coalesce(molecularFormula, "unknown"),
                          adduct = dplyr::coalesce(adduct, "unknown"),
                          `NPC#class` = dplyr::coalesce(`NPC#class`, "unknown"),
                          `NPC#superclass` = dplyr::coalesce(`NPC#superclass`, "unknown")) %>%
            dplyr::mutate(P_molecularFormula = ifelse(!is.na(molecularFormula) & !is.na(value_type),
                                                      mapply(adjust_formula, molecularFormula, value_type),
                                                      molecularFormula)) %>%
            dplyr::mutate(P_row.ID = as.character(P_row.ID)) %>%
            dplyr::left_join(formula_data, by = c("P_row.ID" = "mappingFeatureId")) %>%
            dplyr::mutate(P_sirius_molecularFormula = dplyr::coalesce(P_sirius_molecularFormula, "unknown"))
        }
        
        if (input$enable_confidence_filter) {
          annotated <- annotated %>%
            dplyr::mutate(P_label = dplyr::case_when(
              P_molecularFormula == P_sirius_molecularFormula ~ "high confidence",
              extract_carbon_count(P_molecularFormula) == extract_carbon_count(P_sirius_molecularFormula) ~ "mid confidence",
              TRUE ~ "low confidence"
            ))
        } else {
          annotated$P_label <- "not classified"
        }
        
        writexl::write_xlsx(annotated, file.path(folder_full, paste0("SP_", basename(files[i]))))
        filtered <- annotated %>% dplyr::filter(P_label != "low confidence")
        writexl::write_xlsx(filtered, file.path(folder_filtered, paste0("SP_filtered_", basename(files[i]))))
      }
      
      full_files <- list.files(folder_full, pattern = "\\.xlsx$", full.names = FALSE)
      filtered_files <- list.files(folder_filtered, pattern = "\\.xlsx$", full.names = FALSE)
      
      output$annot_file_selector_ui <- renderUI({
        files <- filtered_files 
        tagList(
          tags$label("Select Annotation Result File:", style = "font-weight:bold;"),
          selectInput("selected_annot_file", NULL, choices = files, width = '100%')
        )
      })
    })
  })
  
  observeEvent(input$selected_annot_file, {
    req(input$selected_annot_file)
    folder <- "./output_data/4_Annotation_Ext_structure_library3"
    file_path <- file.path(folder, input$selected_annot_file)
    tryCatch({
      df <- readxl::read_excel(file_path, sheet = 1)
      output$annot_table <- DT::renderDataTable({
        DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
      })
    }, error = function(e) {
      showNotification(paste("Error reading file:", input$selected_annot_file), type = "error")
    })
  })
  
}

shinyApp(ui, server)
