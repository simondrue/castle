# snapshot

    Code
      single_import_patient %>% data.frame()
    Output
             FileName Well           Sample Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv  M01 Patient1 PlasmaA       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv  A01 Patient1 PlasmaA       Unknown       Unknown KRAS G12D mut
      3 patient_1.csv  B01 Patient1 PlasmaA       Unknown       Unknown KRAS G12D mut
        MutantOnlyDroplets WildtypeOnlyDroplets DoubleNegativeDroplets
      1                 20                17647                  70449
      2                  1                 2771                  11147
      3                  2                 3002                  11871
        DoublePositiveDroplets TotalDroplets               MergedWells
      1                      2         88118 (A01,B01,C01,D01,E01,F01)
      2                      0         13919                      <NA>
      3                      0         14875                      <NA>
        NumberOfMergedWells
      1                   6
      2                   1
      3                   1

---

    Code
      multi_import_patient %>% data.frame()
    Output
             FileName Well           Sample Ch1TargetType Ch2TargetType        Target
      1 patient_1.csv  M01 Patient1 PlasmaA       Unknown       Unknown KRAS G12D mut
      2 patient_1.csv  A01 Patient1 PlasmaA       Unknown       Unknown KRAS G12D mut
      3 patient_1.csv  B01 Patient1 PlasmaA       Unknown       Unknown KRAS G12D mut
      4 patient_2.csv  M01 Patient2 PlasmaA       Unknown       Unknown KRAS G12D mut
      5 patient_2.csv  A01 Patient2 PlasmaA       Unknown       Unknown KRAS G12D mut
      6 patient_2.csv  B01 Patient2 PlasmaA       Unknown       Unknown KRAS G12D mut
        MutantOnlyDroplets WildtypeOnlyDroplets DoubleNegativeDroplets
      1                 20                17647                  70449
      2                  1                 2771                  11147
      3                  2                 3002                  11871
      4                 20                17647                  70449
      5                  1                 2771                  11147
      6                  2                 3002                  11871
        DoublePositiveDroplets TotalDroplets               MergedWells
      1                      2         88118 (A01,B01,C01,D01,E01,F01)
      2                      0         13919                      <NA>
      3                      0         14875                      <NA>
      4                      2         88118 (A01,B01,C01,D01,E01,F01)
      5                      0         13919                      <NA>
      6                      0         14875                      <NA>
        NumberOfMergedWells
      1                   6
      2                   1
      3                   1
      4                   6
      5                   1
      6                   1

---

    Code
      merge_wells_import_patient %>% data.frame()
    Output
                  Sample      FileName        Target Ch1TargetType Ch2TargetType
      1 Patient1 PlasmaA patient_1.csv KRAS G12D mut       Unknown       Unknown
      2 Patient2 PlasmaA patient_2.csv KRAS G12D mut       Unknown       Unknown
        WildtypeOnlyDroplets MutantOnlyDroplets DoubleNegativeDroplets
      1                 5773                  3                  23018
      2                 5773                  3                  23018
        DoublePositiveDroplets TotalDroplets NumberOfMergedWells MergedWells Well
      1                      0         28794                   2   (A01,B01)  M01
      2                      0         28794                   2   (A01,B01)  M02

