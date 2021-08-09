test_that("columns", {
  path_to_my_test_sample <- system.file(
    "extdata/test_data/",
    c(
      "patient_1.csv",
      "patient_2.csv"
    ),
    package = "castle"
  )

  single_import_patient <- import_QS_files(path_to_my_test_sample[1])
  multi_import_patient <- import_QS_files(path_to_my_test_sample)
  merge_wells_yes_import_patient <- import_QS_files(path_to_my_test_sample, merge_wells = "yes")
  merge_wells_no_import_patient <- import_QS_files(path_to_my_test_sample, merge_wells = "no")
  merge_files_import_patient <- import_QS_files(path_to_my_test_sample, merge_wells = "yes", merge_files = TRUE)

  expected_columns <- c(
    "FileName",
    "Well", "Sample",
    "Ch1TargetType", "Ch2TargetType", "Target",
    "MutantOnlyDroplets", "WildtypeOnlyDroplets", "DoubleNegativeDroplets", "DoublePositiveDroplets", "TotalDroplets",
    "MergedWells", "NumberOfMergedWells"
  )

  single_observed_columns <- colnames(single_import_patient)
  multi_observed_columns <- colnames(multi_import_patient)
  merge_wells_yes_observed_columns <- colnames(merge_wells_yes_import_patient)
  merge_wells_no_observed_columns <- colnames(merge_wells_no_import_patient)
  merge_files_observed_columns <- colnames(merge_files_import_patient)


  expect_true(setequal(expected_columns, single_observed_columns))
  expect_true(setequal(expected_columns, multi_observed_columns))
  expect_true(setequal(expected_columns, merge_wells_yes_observed_columns))
  expect_true(setequal(expected_columns, merge_wells_no_observed_columns))
  expect_true(setequal(expected_columns, merge_files_observed_columns))
})

test_that("Wrong function calls", {
  path_to_my_test_samples <- system.file(
    "extdata/test_data/", c("patient_1.csv"),
    package = "castle"
  )

  expect_error(
    import_QS_files("do/not/exist/1"),
    regexp = "do/not/exist/1.*do not exist"
  )
  expect_error(
    import_QS_files(c("do/not/exist/1", "do/not/exist/2")),
    regexp = "do/not/exist/1.*do/not/exist/2.*do not exist"
  )

  expect_error(
    import_QS_files(c("")),
    regexp = "''.*do not exist"
  )

  expect_error(
    import_QS_files(path_to_my_test_samples, merge_wells = "unusable_string"),
    regexp = "merge_wells should be 'yes', 'no', 'qs' or 'none'"
  )
})

test_that("annotations", {
  path_to_my_test_samples <- system.file(
    "extdata/test_data/", c("patient_1.csv", "patient_2.csv"),
    package = "castle"
  )

  ann_import_df <-
    import_QS_files(
      path_to_my_test_samples,
      annotations = data.frame(my_column = "test")
    )
  expect_true(
    "my_column" %in% colnames(ann_import_df)
  )

  ann_import_df <-
    import_QS_files(
      path_to_my_test_samples,
      annotations = list(my_column = "test")
    )
  expect_true(
    "my_column" %in% colnames(ann_import_df)
  )
})

test_that("sample annotations", {
  path_to_my_test_samples <- system.file(
    "extdata/test_data/", c("patient_1.csv", "patient_2.csv"),
    package = "castle"
  )

  # Import samples with annotations
  sample_ann_df <-
    data.frame(
      Sample = c("Patient1 PlasmaA", "Patient2 PlasmaA"),
      my_sample_column = c("P1", "P2")
    )
  sammple_ann_import_df <- import_QS_files(path_to_my_test_samples,
    sample_annotations = sample_ann_df
  )

  expect_true(
    "my_sample_column" %in% colnames(sammple_ann_import_df)
  )
  expect_true(
    all(sammple_ann_import_df %>%
      dplyr::filter(.data$Sample == "Patient1 PlasmaA") %>%
      dplyr::pull(.data$my_sample_column) == "P1") &
      all(sammple_ann_import_df %>%
        dplyr::filter(.data$Sample == "Patient2 PlasmaA") %>%
        dplyr::pull(.data$my_sample_column) == "P2")
  )
  expect_error(
    import_QS_files(path_to_my_test_samples,
      sample_annotations = data.frame(a = "a")
    ),
    regexp = "'sample_annotations' does not include a column 'Sample'."
  )
})

test_that("single import", {
  path_to_my_test_samples <- system.file(
    "extdata/test_data/",
    c(
      "patient_1.csv"
    ),
    package = "castle"
  )

  single_import_patient <- import_QS_files(path_to_my_test_samples)
  merge_wells_yes_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "yes")
  merge_wells_no_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "no")
  merge_wells_QS_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "qs")

  expect_snapshot(single_import_patient %>% data.frame())
  expect_snapshot(merge_wells_yes_import_patient %>% data.frame())
  expect_snapshot(merge_wells_no_import_patient %>% data.frame())
  expect_snapshot(merge_wells_QS_import_patient %>% data.frame())
})

test_that("multi import", {
  path_to_my_test_samples <- system.file(
    "extdata/test_data/",
    c(
      "patient_1.csv",
      "patient_2.csv",
      "patient_2_extra.csv"
    ),
    package = "castle"
  )

  multi_import_patient <- import_QS_files(path_to_my_test_samples)
  merge_wells_yes_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "yes")
  merge_wells_no_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "no")
  merge_wells_and_files_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "yes", merge_files = TRUE)
  merge_wells_QS_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = "qs")

  expect_equal(
    nrow(multi_import_patient), 10
  )
  expect_equal(
    nrow(merge_wells_yes_import_patient), 4
  )
  expect_equal(
    nrow(merge_wells_no_import_patient), 7
  )
  expect_equal(
    nrow(merge_wells_and_files_import_patient), 3
  )
  expect_equal(
    nrow(merge_wells_QS_import_patient), 4
  )


  expect_snapshot(multi_import_patient %>% data.frame())
  expect_snapshot(merge_wells_yes_import_patient %>% data.frame())
  expect_snapshot(merge_wells_no_import_patient %>% data.frame())
  expect_snapshot(merge_wells_and_files_import_patient %>% data.frame())
  expect_snapshot(merge_wells_QS_import_patient %>% data.frame())
})
