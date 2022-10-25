KPL/FK
 
   FILE: stations.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2021-11-23T15:01:55
   PINPOINT DEFINITIONS FILE: stations_latlon.def
   PINPOINT PCK FILE:         pck00010.tpc
   PINPOINT SPK FILE:         stations.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'MILANO'
   NAIF_BODY_CODE                      += 399101
 
   NAIF_BODY_NAME                      += 'WELLINGTON'
   NAIF_BODY_CODE                      += 399102
 
   NAIF_BODY_NAME                      += 'LA-SILLA'
   NAIF_BODY_CODE                      += 399103
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame MILANO_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame MILANO_TOPO is centered at the
      site MILANO, which has Cartesian coordinates
 
         X (km):                  0.4421005637378E+04
         Y (km):                  0.7124529090551E+03
         Z (km):                  0.4526578207835E+04
 
      and planetodetic coordinates
 
         Longitude (deg):         9.1546100000000
         Latitude  (deg):        45.5012200000000
         Altitude   (km):         0.2000000000295E-01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_MILANO_TOPO                   =  1399101
   FRAME_1399101_NAME                  =  'MILANO_TOPO'
   FRAME_1399101_CLASS                 =  4
   FRAME_1399101_CLASS_ID              =  1399101
   FRAME_1399101_CENTER                =  399101
 
   OBJECT_399101_FRAME                 =  'MILANO_TOPO'
 
   TKFRAME_1399101_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399101_SPEC                =  'ANGLES'
   TKFRAME_1399101_UNITS               =  'DEGREES'
   TKFRAME_1399101_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399101_ANGLES              =  (   -9.1546100000000,
                                             -44.4987800000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame WELLINGTON_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame WELLINGTON_TOPO is centered at the
      site WELLINGTON, which has Cartesian coordinates
 
         X (km):                 -0.4779894399227E+04
         Y (km):                  0.4377829599012E+03
         Z (km):                 -0.4186283224013E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       174.7669700000000
         Latitude  (deg):       -41.2843700000000
         Altitude   (km):         0.1170000000005E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_WELLINGTON_TOPO               =  1399102
   FRAME_1399102_NAME                  =  'WELLINGTON_TOPO'
   FRAME_1399102_CLASS                 =  4
   FRAME_1399102_CLASS_ID              =  1399102
   FRAME_1399102_CENTER                =  399102
 
   OBJECT_399102_FRAME                 =  'WELLINGTON_TOPO'
 
   TKFRAME_1399102_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399102_SPEC                =  'ANGLES'
   TKFRAME_1399102_UNITS               =  'DEGREES'
   TKFRAME_1399102_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399102_ANGLES              =  ( -174.7669700000000,
                                            -131.2843700000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame LA-SILLA_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame LA-SILLA_TOPO is centered at the
      site LA-SILLA, which has Cartesian coordinates
 
         X (km):                  0.1838367489900E+04
         Y (km):                 -0.5258770040539E+04
         Z (km):                 -0.3100360107992E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       -70.7313300000000
         Latitude  (deg):       -29.2611700000000
         Altitude   (km):         0.2399999999999E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_LA-SILLA_TOPO                 =  1399103
   FRAME_1399103_NAME                  =  'LA-SILLA_TOPO'
   FRAME_1399103_CLASS                 =  4
   FRAME_1399103_CLASS_ID              =  1399103
   FRAME_1399103_CENTER                =  399103
 
   OBJECT_399103_FRAME                 =  'LA-SILLA_TOPO'
 
   TKFRAME_1399103_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399103_SPEC                =  'ANGLES'
   TKFRAME_1399103_UNITS               =  'DEGREES'
   TKFRAME_1399103_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399103_ANGLES              =  ( -289.2686700000000,
                                            -119.2611700000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file stations_latlon.def
--------------------------------------------------------------------------------
 
begindata
 
   SITES         = ( 'MILANO',
                     'WELLINGTON',
                     'LA-SILLA' )
 
   MILANO_CENTER = 399
   MILANO_FRAME  = 'EARTH_FIXED'
   MILANO_IDCODE = 399101
   MILANO_LATLON = ( 45.50122, 9.15461, 0.020 )
   MILANO_UP     = 'Z'
   MILANO_NORTH  = 'X'
 
   WELLINGTON_CENTER = 399
   WELLINGTON_FRAME  = 'EARTH_FIXED'
   WELLINGTON_IDCODE = 399102
   WELLINGTON_LATLON = ( -41.28437, 174.76697, 0.117 )
   WELLINGTON_UP     = 'Z'
   WELLINGTON_NORTH  = 'X'
 
   LA-SILLA_CENTER = 399
   LA-SILLA_FRAME  = 'EARTH_FIXED'
   LA-SILLA_IDCODE = 399103
   LA-SILLA_LATLON = ( -29.26117, -70.73133, 2.400 )
   LA-SILLA_UP     = 'Z'
   LA-SILLA_NORTH  = 'X'
 
begintext
 
begintext
 
[End of definitions file]
 
