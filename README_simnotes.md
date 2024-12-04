## Notes on simulations and testing

First WRF tests of 2-domain setup (5km, 1km):
- Problems/testing:
    - d02 crashing:
        - crashing after about 30 min with no indication why
        - a test with GEFS ICs/BCs ran with wdamp turned on and epssm increased from 0.1 to 0.2
        - therefore, now testing downscaled run to isolate the key factor:
            - wrf_fine: wdamp on, epssm unchanged (0.1)
                - outcome: CRASHED EARLY
            - wrf_fine_epssm: wdamp off, epssm increased to 0.2
                - outcome: OUT OF TIME
            - wrf_fine_both: wdamp on, epssm increased to 0.2
                - outcome: OUT OF TIME
            - ***Conclusion:*** epssm is the key!! Will now go with 0.2 for d02, which = that of d01
- tests of run time:
    - d01:
        - first run completed in 04:27 over 40 nodes = 23k core hours = 115k for 5 ens members
        - completed in 04:50 with 33 nodes = 18k core hours = 86k for 5 members
    - d02:
        - completed in 12+10=22 with 33 nodes = 93k core hours = 465k for 5 members
    - total cost for one 4-day case = 551k core hours
