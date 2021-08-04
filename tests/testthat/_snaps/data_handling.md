# snapshot

    Code
      single_import_patient %>% data.frame()
    Output
             FileName           Sample Well Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv Patient1 PlasmaA  A01       Unknown       Unknown KRAS G12D mut
      3 patient_1.csv Patient1 PlasmaA  B01       Unknown       Unknown KRAS G12D mut
        DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1                  70449                17647                 20
      2                  11147                 2771                  1
      3                  11871                 3002                  2
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells
      1                      2         88118                   6
      2                      0         13919                   1
      3                      0         14875                   1
                      MergedWells
      1 (A01,B01,C01,D01,E01,F01)
      2                      <NA>
      3                      <NA>

---

    Code
      multi_import_patient %>% data.frame()
    Output
             FileName           Sample Well Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv Patient1 PlasmaA  A01       Unknown       Unknown KRAS G12D mut
      3 patient_1.csv Patient1 PlasmaA  B01       Unknown       Unknown KRAS G12D mut
      4 patient_2.csv Patient2 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      5 patient_2.csv Patient2 PlasmaA  A01       Unknown       Unknown KRAS G12D mut
      6 patient_2.csv Patient2 PlasmaA  B01       Unknown       Unknown KRAS G12D mut
        DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1                  70449                17647                 20
      2                  11147                 2771                  1
      3                  11871                 3002                  2
      4                  70449                17647                 20
      5                  11147                 2771                  1
      6                  11871                 3002                  2
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells
      1                      2         88118                   6
      2                      0         13919                   1
      3                      0         14875                   1
      4                      2         88118                   6
      5                      0         13919                   1
      6                      0         14875                   1
                      MergedWells
      1 (A01,B01,C01,D01,E01,F01)
      2                      <NA>
      3                      <NA>
      4 (A01,B01,C01,D01,E01,F01)
      5                      <NA>
      6                      <NA>

---

    Code
      merge_wells_import_patient %>% data.frame()
    Output
             FileName           Sample Well Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv Patient1 PlasmaA  M01       Unknown       Unknown KRAS G12D mut
      2 patient_2.csv Patient2 PlasmaA  M02       Unknown       Unknown KRAS G12D mut
        DoubleNegativeDroplets WildtypeOnlyDroplets MutantOnlyDroplets
      1                  23018                 5773                  3
      2                  23018                 5773                  3
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells
      1                      0         28794                   2   (A01,B01)
      2                      0         28794                   2   (A01,B01)

