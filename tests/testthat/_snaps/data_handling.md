# snapshot - single import

    Code
      single_import_patient %>% data.frame()
    Output
             FileName           Sample Well Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv Patient1 PlasmaA  A01       Unknown       Unknown KRAS G12D mut
      3 patient_1.csv Patient1 PlasmaA  B01       Unknown       Unknown KRAS G12D mut
      4 patient_1.csv          Control  C01       Unknown       Unknown KRAS G12D mut
        DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1                  70449                17647                 20
      2                  11147                 2771                  1
      3                  11871                 3002                  2
      4                  11871                 3002                  2
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                      2         88118                   2   (A01,B01)
      2                      0         13919                   1        <NA>
      3                      0         14875                   1        <NA>
      4                      0         14875                   1        <NA>

---

    Code
      merge_wells_import_patient %>% data.frame()
    Output
             FileName           Sample Well Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv          Control  C01       Unknown       Unknown KRAS G12D mut
        DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1                  23018                 5773                  3
      2                  11871                 3002                  2
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                      0         28794                   2   (A01,B01)
      2                      0         14875                   1        <NA>

---

    Code
      merge_wells_QS_import_patient %>% data.frame()
    Output
             FileName           Sample Well Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv          Control  C01       Unknown       Unknown KRAS G12D mut
        DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1                  70449                17647                 20
      2                  11871                 3002                  2
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                      2         88118                   2   (A01,B01)
      2                      0         14875                   1        <NA>

# snapshot - multi import

    Code
      multi_import_patient %>% data.frame()
    Output
                    FileName           Sample Well Ch1TargetType Ch2TargetType
      1        patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown
      2        patient_1.csv Patient1 PlasmaA  A01       Unknown       Unknown
      3        patient_1.csv Patient1 PlasmaA  B01       Unknown       Unknown
      4        patient_1.csv          Control  C01       Unknown       Unknown
      5        patient_2.csv Patient2 PlasmaA  M01       Unknown       Unknown
      6        patient_2.csv Patient2 PlasmaA  A01       Unknown       Unknown
      7        patient_2.csv Patient2 PlasmaA  B01       Unknown       Unknown
      8  patient_2_extra.csv Patient2 PlasmaA  M01       Unknown       Unknown
      9  patient_2_extra.csv Patient2 PlasmaA  A01       Unknown       Unknown
      10 patient_2_extra.csv Patient2 PlasmaA  B01       Unknown       Unknown
                Target DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1  KRAS G12D mut                  70449                17647                 20
      2  KRAS G12D mut                  11147                 2771                  1
      3  KRAS G12D mut                  11871                 3002                  2
      4  KRAS G12D mut                  11871                 3002                  2
      5  KRAS G12D mut                  70449                17647                 20
      6  KRAS G12D mut                  11147                 2771                  1
      7  KRAS G12D mut                  11871                 3002                  2
      8  KRAS G12D mut                  70449                17647                 20
      9  KRAS G12D mut                  11147                 2771                  1
      10 KRAS G12D mut                  11871                 3002                  2
         DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                       2         88118                   2   (A01,B01)
      2                       0         13919                   1        <NA>
      3                       0         14875                   1        <NA>
      4                       0         14875                   1        <NA>
      5                       2         88118                   2   (A01,B01)
      6                       0         13919                   1        <NA>
      7                       0         14875                   1        <NA>
      8                       2         88118                   2   (A01,B01)
      9                       0         13919                   1        <NA>
      10                      0         14875                   1        <NA>

---

    Code
      merge_wells_import_patient %>% data.frame()
    Output
                   FileName           Sample Well Ch1TargetType Ch2TargetType
      1       patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown
      2       patient_2.csv Patient2 PlasmaA  M02       Unknown       Unknown
      3 patient_2_extra.csv Patient2 PlasmaA  M03       Unknown       Unknown
      4       patient_1.csv          Control  C01       Unknown       Unknown
               Target DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1 KRAS G12D mut                  23018                 5773                  3
      2 KRAS G12D mut                  23018                 5773                  3
      3 KRAS G12D mut                  23018                 5773                  3
      4 KRAS G12D mut                  11871                 3002                  2
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                      0         28794                   2   (A01,B01)
      2                      0         28794                   2   (A01,B01)
      3                      0         28794                   2   (A01,B01)
      4                      0         14875                   1        <NA>

---

    Code
      merge_wells_and_files_import_patient %>% data.frame()
    Output
                                 FileName           Sample Well Ch1TargetType
      1                     patient_1.csv Patient1 PlasmaA  M01       Unknown
      2 patient_2.csv,patient_2_extra.csv Patient2 PlasmaA  M02       Unknown
      3                     patient_1.csv          Control  C01       Unknown
        Ch2TargetType        Target DoubleNegativeDroplets WildtypeOnlyDroplets
      1       Unknown KRAS G12D mut                  23018                 5773
      2       Unknown KRAS G12D mut                  46036                11546
      3       Unknown KRAS G12D mut                  11871                 3002
        MutantOnlyDroplets DoublePositiveDroplets TotalDroplets NumberOfMergedWells
      1                  3                      0         28794                   2
      2                  6                      0         57588                   4
      3                  2                      0         14875                   1
              MergedWells
      1         (A01,B01)
      2 (A01,B01,A01,B01)
      3              <NA>

---

    Code
      merge_wells_QS_import_patient %>% data.frame()
    Output
                   FileName           Sample Well Ch1TargetType Ch2TargetType
      1       patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown
      2       patient_1.csv          Control  C01       Unknown       Unknown
      3       patient_2.csv Patient2 PlasmaA  M01       Unknown       Unknown
      4 patient_2_extra.csv Patient2 PlasmaA  M01       Unknown       Unknown
               Target DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1 KRAS G12D mut                  70449                17647                 20
      2 KRAS G12D mut                  11871                 3002                  2
      3 KRAS G12D mut                  70449                17647                 20
      4 KRAS G12D mut                  70449                17647                 20
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                      2         88118                   2   (A01,B01)
      2                      0         14875                   1        <NA>
      3                      2         88118                   2   (A01,B01)
      4                      2         88118                   2   (A01,B01)

