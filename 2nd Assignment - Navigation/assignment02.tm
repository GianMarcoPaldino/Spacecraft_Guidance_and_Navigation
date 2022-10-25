KPL/MK

Description:
   Meta-kernel for assignment 01
   
Kernel content:

 - LSK kernel: 
    naif0012.tls -> Standard leap second kernel

 - SPK kernels:
    de425s.bsp -> Planetary ephemeris
    stations.bsp -> List the stations

 - PCK kernels: 
    pck00010.tpc -> planetary constant kernel (TEXT kernel)
    earth_latest_high_prec.bpc -> high precision Earth planetary constants
                                 (nutation, polar motion). Binary kernel.
                                 Updated daily.

 - FK kernels:
    earth_fixed.tf -> Defines an alias for a EARTH_FIXED frame, which can 
                      be called for transformations.
                      EARTH_FIXED can map either IAU_EARTH or ITRF93.
    stations.tf -> Frame kernels defining stations ID and their corresponding 
                   topocentric frames.


 Notes:
   This material was prepared to support the course 'Satellite Guidance
   and Navigation', AY 2021/2022.


 NB: this kernel was generated on Windows PC, file paths and line-ending 
     shall be changed on MacOS and linux.

\begindata

    PATH_VALUES = ( 'kernels' )
    PATH_SYMBOLS = ( 'KERNELS' )
    KERNELS_TO_LOAD = (
                       '$KERNELS/naif0012.tls',
                       '$KERNELS/pck00010.tpc',
                       '$KERNELS/earth_latest_high_prec.bpc',
					   '$KERNELS/earth_fixed.tf',
					   '$KERNELS/de425s.bsp',
                       '$KERNELS/gm_de431.tpc',
					   '$KERNELS/stations.tf',
					   '$KERNELS/stations.bsp',
					   '$KERNELS/stations_ex3.tf',
					   '$KERNELS/stations_ex3.bsp'
                      )

\begintext

KERNELS: 
(downloaded from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/)