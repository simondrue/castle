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
        N_d_neg N_d_pos N_WT_only N_M_only        r_est     l_est pval_r_leq_0
      1    7521     644      4889      946 0.1006261466 0.5028811    0.0000000
      2    8466     152      5252      130 0.0004842163 0.4877604    0.6735607
        r_LLR_test_stat is_tumor_positive total_tumor_molecules_expected
      1    3172.7099905              TRUE                    1408.766052
      2       0.1774652             FALSE                       6.779028
        total_droplets r_CI_lower  r_CI_upper total_tumor_molecules_CI_lower
      1          14000 0.09480199 0.106645448                       1327.228
      2          14000 0.00000000 0.002857597                          0.000
        total_tumor_molecules_CI_upper l_CI_lower l_CI_upper
      1                     1493.03627  0.4895578  0.5164573
      2                       40.00636  0.4746937  0.5010782

---

    Code
      no_CIs_res
    Output
        N_d_neg N_d_pos N_WT_only N_M_only        r_est     l_est pval_r_leq_0
      1    7521     644      4889      946 0.1006261466 0.5028811    0.0000000
      2    8466     152      5252      130 0.0004842163 0.4877604    0.6735607
        r_LLR_test_stat is_tumor_positive total_tumor_molecules_expected
      1    3172.7099905              TRUE                    1408.766052
      2       0.1774652             FALSE                       6.779028
        total_droplets
      1          14000
      2          14000

