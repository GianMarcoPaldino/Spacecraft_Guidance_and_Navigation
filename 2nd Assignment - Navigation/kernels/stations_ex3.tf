KPL/FK
 
   FILE: .\stations_latlon_ex3.fk
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2021-11-25T18:14:55
   PINPOINT DEFINITIONS FILE: .\stations_latlon_ex3.def
   PINPOINT PCK FILE:         .\pck00010.tpc
   PINPOINT SPK FILE:         .\stations_latlon_ex3.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'NEW_NORCIA'
   NAIF_BODY_CODE                      += 399104
 
   NAIF_BODY_NAME                      += 'MALARGUE'
   NAIF_BODY_CODE                      += 399105
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame NEW_NORCIA_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame NEW_NORCIA_TOPO is centered at the
      site NEW_NORCIA, which has Cartesian coordinates
 
         X (km):                 -0.2414063804502E+04
         Y (km):                  0.4907870242295E+04
         Z (km):                 -0.3270605163873E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       116.1914678000000
         Latitude  (deg):       -31.0482254000000
         Altitude   (km):         0.2519999999996E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_NEW_NORCIA_TOPO               =  1399104
   FRAME_1399104_NAME                  =  'NEW_NORCIA_TOPO'
   FRAME_1399104_CLASS                 =  4
   FRAME_1399104_CLASS_ID              =  1399104
   FRAME_1399104_CENTER                =  399104
 
   OBJECT_399104_FRAME                 =  'NEW_NORCIA_TOPO'
 
   TKFRAME_1399104_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399104_SPEC                =  'ANGLES'
   TKFRAME_1399104_UNITS               =  'DEGREES'
   TKFRAME_1399104_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399104_ANGLES              =  ( -116.1914678000000,
                                            -121.0482254000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame MALARGUE_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame MALARGUE_TOPO is centered at the
      site MALARGUE, which has Cartesian coordinates
 
         X (km):                  0.1823334827791E+04
         Y (km):                 -0.4850439313812E+04
         Z (km):                 -0.3708962261707E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       -69.3981934000000
         Latitude  (deg):       -35.7760086000000
         Altitude   (km):         0.1550000000001E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_MALARGUE_TOPO                 =  1399105
   FRAME_1399105_NAME                  =  'MALARGUE_TOPO'
   FRAME_1399105_CLASS                 =  4
   FRAME_1399105_CLASS_ID              =  1399105
   FRAME_1399105_CENTER                =  399105
 
   OBJECT_399105_FRAME                 =  'MALARGUE_TOPO'
 
   TKFRAME_1399105_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399105_SPEC                =  'ANGLES'
   TKFRAME_1399105_UNITS               =  'DEGREES'
   TKFRAME_1399105_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399105_ANGLES              =  ( -290.6018066000000,
                                            -125.7760086000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file .\stations_latlon_ex3.def
--------------------------------------------------------------------------------
 
begindata
 
   SITES         = ( 'NEW_NORCIA',
                     'MALARGUE' )
 
   NEW_NORCIA_CENTER = 399
   NEW_NORCIA_FRAME  = 'EARTH_FIXED'
   NEW_NORCIA_IDCODE = 399104
   NEW_NORCIA_LATLON = (-31.0482254, 116.1914678, 0.252)
   NEW_NORCIA_UP     = 'Z'
   NEW_NORCIA_NORTH  = 'X'
 
   MALARGUE_CENTER = 399
   MALARGUE_FRAME  = 'EARTH_FIXED'
   MALARGUE_IDCODE = 399105
   MALARGUE_LATLON = (-35.7760086, -69.3981934, 1.55)
   MALARGUE_UP     = 'Z'
   MALARGUE_NORTH  = 'X'
 
begintext
 
begintext
 
[End of definitions file]
 
