# training model - snapshot

    Code
      train_simple_ddpcr_model(training_samples = training_samples)
    Output
      $a_est
      [1] 5.125098e-07
      
      $b_est
      [1] 8.088402e-07
      
      $c_est
      [1] 0.1132381
      
      $l_est_vec
      [1] 0.3036824 0.5045560 0.6623755 0.0000000
      

# simulation - training - test

    Code
      multiple_test_res
    Output
        N_d_neg N_d_pos N_WT_only N_M_only mutant_molecules_per_droplet
      1    7521     644      4889      946                 0.1006261466
      2    8466     152      5252      130                 0.0004842163
        wildtype_molecules_per_droplet     p_val test_statistic mutation_detected
      1                      0.5028811 0.0000000   3172.7099905              TRUE
      2                      0.4877604 0.6735607      0.1774652             FALSE
        allele_frequency total_mutant_molecules total_wildtype_molecules
      1     0.1667356131            1408.766052                 7040.335
      2     0.0009917494               6.779028                 6828.645
        mutant_molecules_per_droplet_CI_lower mutant_molecules_per_droplet_CI_upper
      1                            0.09480199                           0.106645448
      2                            0.00000000                           0.002857597
        allele_frequency_CI_lower allele_frequency_CI_upper
      1                 0.1586158               0.174964410
      2                 0.0000000               0.005824486
        total_mutant_molecules_CI_lower total_mutant_molecules_CI_upper
      1                        1327.228                      1493.03627
      2                           0.000                        40.00636
        wildtype_molecules_per_droplet_CI_lower
      1                               0.4895578
      2                               0.4746937
        wildtype_molecules_per_droplet_CI_upper
      1                               0.5164573
      2                               0.5010782

---

    Code
      no_CIs_res
    Output
        N_d_neg N_d_pos N_WT_only N_M_only mutant_molecules_per_droplet
      1    7521     644      4889      946                 0.1006261466
      2    8466     152      5252      130                 0.0004842163
        wildtype_molecules_per_droplet     p_val test_statistic mutation_detected
      1                      0.5028811 0.0000000   3172.7099905              TRUE
      2                      0.4877604 0.6735607      0.1774652             FALSE
        allele_frequency total_mutant_molecules total_wildtype_molecules
      1     0.1667356131            1408.766052                 7040.335
      2     0.0009917494               6.779028                 6828.645

