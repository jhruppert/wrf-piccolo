## Notes on simulations and testing

First WRF tests of 2-domain setup (5km, 1km):
- crashes of d02:
    - crashing after about 30 min with no indication why
    - a test with GEFS ICs/BCs ran with wdamp turned on and epssm increased from 0.1 to 0.2
    - therefore, now testing downscaled run to isolate the key factor:
        - wrf_fine: wdamp on, epssm unchanged (0.1)
            - outcome: CRASHED EARLY
        - wrf_fine_epssm: wdamp off, epssm increased to 0.2
            - outcome: OUT OF TIME
        - wrf_fine_both: wdamp on, epssm increased to 0.2
            - outcome: OUT OF TIME
        - Conclusion: epssm is the key!! Will now go with 0.2 for d02, which = that of d01
- d01:
    - first run completed in 04:27 over 40 nodes = 23k core hours = 115k for 5 ens members
    - will try 33 nodes next
- WRF for d02 completed in 
over 33 nodes = xxk core hours = xxxk for 5 ens members
