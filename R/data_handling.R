import_QS_files <- function(paths,
                            Ch1_is_mutation = TRUE,
                            annotations = NULL,
                            sample_annotations = NULL,
                            merge_wells = "no",
                            include_qs_concentrations = FALSE) {

  # Split input into files and directories
  file_paths <- paths[utils::file_test("-f", paths)]
  dir_paths <- paths[utils::file_test("-d", paths)]

  # Get files
  load_files_df <- suppressWarnings(readr::read_csv(file_paths, id = "FilePath", show_col_types = FALSE))

  # Get .csv files from directories
  dir_files <- list.files(dir_paths, pattern = ".csv", full.names = TRUE)
  load_dirs_df <- suppressWarnings(readr::read_csv(dir_files, id = "FilePath", show_col_types = FALSE))

  # TODO: Ch1 and Ch2 in wide format. include concentrations

  # Bind data

  df <- dplyr::bind_rows(load_files_df, load_dirs_df)

  ch1_df <- df %>%
    dplyr::filter(
      grepl("Ch1", .data$TargetType)
    ) %>%
    dplyr::mutate(
      Ch1TargetType = stringr::str_remove(.data$TargetType, pattern = "Ch1")
    ) %>%
    dplyr::select(-c("TargetType"))

  ch2_df <- df %>%
    dplyr::filter(
      grepl("Ch2", .data$TargetType)
    ) %>%
    dplyr::mutate(
      Ch2TargetType = stringr::str_remove(.data$TargetType, pattern = "Ch2")
    ) %>%
    dplyr::select("FilePath", "Well", "ExptType", "Experiment", "Sample", "Ch2TargetType")

  df <- dplyr::full_join(ch1_df, ch2_df, by = c("FilePath", "Well", "ExptType", "Experiment", "Sample"))

  df <- df %>%
    dplyr::mutate(
      FileName = basename(.data$FilePath)
    ) %>%
    dplyr::select(
      "FileName",
      "Well", "Sample", "Ch1TargetType", "Ch2TargetType", "Target",
      "Ch1+Ch2-", "Ch1-Ch2+", "Ch1-Ch2-", "Ch1+Ch2+",
      "AcceptedDroplets", "FilePath", "MergedWells"
    ) %>%
    dplyr::rename(
      MutantOnlyDroplets = ifelse(Ch1_is_mutation, "Ch1+Ch2-", "Ch1-Ch2+"),
      WildtypeOnlyDroplets = ifelse(Ch1_is_mutation, "Ch1-Ch2+", "Ch1+Ch2-"),
      DoubleNegativeDroplets = "Ch1-Ch2-",
      DoublePositiveDroplets = "Ch1+Ch2+",
      TotalDroplets = .data$AcceptedDroplets
    ) %>%
    dplyr::mutate(
      NumberOfMergedWells = stringr::str_count(ifelse(is.na(.data$MergedWells), "", .data$MergedWells), pattern = ",") + 1
    )

  if (merge_wells == "Yes") {
    df <- df %>%
      dplyr::filter(
        !grepl("M", .data$Well)
      ) %>%
      dplyr::group_by(.data$Sample) %>%
      dplyr::summarise(
        FileName = paste0(unique(.data$FileName), collapse = ","),
        Target = paste0(unique(.data$Target), collapse = ","),
        Ch1TargetType = paste0(unique(.data$Ch1TargetType), collapse = ","),
        Ch2TargetType = paste0(unique(.data$Ch2TargetType), collapse = ","),
        WildtypeOnlyDroplets = sum(.data$WildtypeOnlyDroplets),
        MutantOnlyDroplets = sum(.data$MutantOnlyDroplets),
        DoubleNegativeDroplets = sum(.data$DoubleNegativeDroplets),
        DoublePositiveDroplets = sum(.data$DoublePositiveDroplets),
        TotalDroplets = sum(.data$TotalDroplets),
        NumberOfMergedWells = dplyr::n(),
        MergedWells = paste0(c("(", paste0(c(.data$Well), collapse = ","), ")"), collapse = ""),
        .groups = "drop"
      )
  } else if (merge_wells == "QS") {
    merged_samples <- df %>%
      dplyr::filter(grepl("M", .data$Well)) %>%
      dplyr::pull(.data$Sample) %>%
      unique()

    df <- df %>%
      dplyr::filter(grepl("M", .data$Well) | !.data$Sample %in% merged_samples)
  }

  if (!is.null(annotations)) {
    df <- df %>% dplyr::bind_cols(annotations)
  }

  if (!is.null(sample_annotations)) {
    df <- df %>% dplyr::left_join(sample_annotations, by = "Sample")
  }

  return(df)
}



paths <- list.files("~/Downloads/Tenna example/data/", full.names = TRUE)

import_QS_files(paths)

annotations <- data.frame(
  project = "Horse",
  analysis_time = Sys.time()
)

sample_annotations <- data.frame(
  Sample = c("Plasma_A", "Plasma_B"),
  my_annotation = c("My favorite", "Just another sample")
)

annotated_import <- import_QS_files(paths = paths, annotations = annotations, sample_annotations = sample_annotations)

# View(annotated_import)

annotated_import <- import_QS_files(paths = paths, annotations = annotations, sample_annotations = sample_annotations, merge_wells = "Yes")
