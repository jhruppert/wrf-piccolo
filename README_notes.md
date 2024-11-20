## Notes on simulations and testing

First test of outer domain from 2-domain setup (5km, 1km) completed in 04:27 over 40 nodes.

## wrf-ndown steps

Assumes doing a two-domain setup with one domain, ndown, then the second.

1. Run WPS for both domains
    - produces: met_em files for d01, d02
2. Run REAL for coarse domain in wrf_coarse
    - produces: wrfinput_d01, wrflowinput_d01, wrfbdy_d01
3. Run WRF for coarse domain in wrf_coarse
    - change: max_dom = 1
4. Run REAL in wrf_fine for both domains to get new ICs for d02 init time
5. Run NDOWN in wrf_fine:
    - in prep:
        - rename wrfinput_d02 to wrfndi_d02
        - rename wrflowinput_d02 to wrflowinput_d01
        - link wrf_coarse/wrfout_* to wrf_fine
        - add namelist settings:
            - io_form_auxinput2=2
            - set interval_seconds to history output frequency of d01
            - max_dom=2
    - produces: wrfinput_d02, wrfbdy_d02
6. Run WRF on fine domain
    - rename wrfinput_d02 to wrfinput_d01 and wrfbdy_d02 to wrfbdy_d01
    - adapt the namelist to shift d02 settings to d01 position
    - namelist settings:
        - max_dom=1
        - in &bdy_control:
            - have_bcs_moist=.true.
            - have_bcs_scalar=.true.
        
