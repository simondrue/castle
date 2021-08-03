test_that("columns", {
  path_to_my_test_sample <- system.file(
    "extdata/test_samples",
    c(
      "patient_1_plasma_KRAS_G12D.csv",
      "patient_2_plasma_KRAS_G12D.csv"
    ),
    package = "castle"
  )

  single_import_patient <- import_QS_files(path_to_my_test_sample[1])
  multi_import_patient <- import_QS_files(path_to_my_test_sample)
  merge_wells_import_patient <- import_QS_files(path_to_my_test_sample, merge_wells = TRUE)
  merge_files_import_patient <- import_QS_files(path_to_my_test_sample, merge_wells = TRUE, merge_files = TRUE)

  expected_columns <- c(
    "FileName",
    "Well", "Sample",
    "Ch1TargetType", "Ch2TargetType", "Target",
    "MutantOnlyDroplets", "WildtypeOnlyDroplets", "DoubleNegativeDroplets", "DoublePositiveDroplets", "TotalDroplets",
    "MergedWells", "NumberOfMergedWells"
  )

  single_observed_columns = colnames(single_import_patient)
  multi_observed_columns = colnames(multi_import_patient)
  merge_wells_observed_columns = colnames(merge_wells_import_patient)
  merge_files_observed_columns = colnames(merge_files_import_patient)


  expect_true(setequal(expected_columns,single_observed_columns))
  expect_true(setequal(expected_columns,multi_observed_columns))
  expect_true(setequal(expected_columns,merge_wells_observed_columns))
  expect_true(setequal(expected_columns,merge_files_observed_columns))
})

test_that("snapshot", {
  path_to_my_test_samples <- system.file(
    "tests/testthat/test_data/",
    c(
      "patient_1.csv",
      "patient_2.csv"
    ),
    package = "castle"
  )

  single_import_patient <- import_QS_files(path_to_my_test_samples[1])
  multi_import_patient <- import_QS_files(path_to_my_test_samples)
  merge_wells_import_patient <- import_QS_files(path_to_my_test_samples, merge_wells = TRUE)

  expect_snapshot(single_import_patient %>% data.frame())
  expect_snapshot(multi_import_patient %>% data.frame())
  expect_snapshot(merge_wells_import_patient %>% data.frame())
})

