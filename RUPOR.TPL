#TMC_RT_H

#define L_waveg @ ( 800.00)
#define W_file  @ r7
#define W_freq  @ (  10.00)
#define W_eps @ ( 4. )
#define W_waveg @ (  400.00)
#define W_input @ (   36.00)
#define W_delta @ (   1.0)
#define W_r @ (   40.00)
#define W_type @ MAGNETIC

#define a1 @ ( W_input )
#define l_rupor @ ( 140 )
#define h_rupor @ ( 64 )

#STEP

#PARM
 ANGLE_UNIT radian;
 FREQ_UNIT GHz;
 LONG_UNIT mm;
 TIME_UNIT ns;
 DELTA W_delta;
 TIME 0.; 10.;
 X_MIN  0.0;
 X_MAX L_waveg;
 Y_MIN  0.0;
 Y_MAX W_waveg;
 FREQ  W_freq;
#END_PARM

#TOPOLOGY

 BLOCK 7;
  INPUT_X  W_type; -a1/2; a1/2; 0; 0.05;
 END_B

 BLOCK 5;
  POLYGON_STAT W_type;
   L   0;   a1/2;
   L   2*a1;   a1/2;
   L   2*a1+(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*((l_rupor)/4);
   L   2*a1+2*(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*(2*(l_rupor)/4);
   L   2*a1+3*(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*(3*(l_rupor)/4);
   L   2*a1+l_rupor;   a1/2+h_rupor;
  END_B
 BLOCK 6;
  POLYGON_STAT W_type;
   L   0;   -a1/2;
   L   2*a1;   -a1/2;
   L   2*a1+(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*((l_rupor)/4);
   L   2*a1+2*(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*(2*(l_rupor)/4);
   L   2*a1+3*(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*(3*(l_rupor)/4);
   L   2*a1+l_rupor;   -a1/2-h_rupor;
  END_B

 BLOCK 55;
!  RECT_STAT [eps1]1.; 0; 100; 0; 100;
  FILE eps1; 0; 100; 0; 100;
 END_B

 BLOCK 56;
!  RECT_STAT [eps1]1.; 0; 100; 0; 100;
  FILE eps1; 0; 100; 0; 100;
 END_B





 BLOCK 8;
  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0.0; W_waveg-2*W_delta;
 END_B
 BLOCK 9;
  RECT_STAT ABSORBER; 1*W_delta; L_waveg-2*W_delta; W_delta; 5*W_delta;
 END_B
 BLOCK 10;
  RECT_STAT ABSORBER; 1*W_delta; L_waveg; W_waveg-5*W_delta; W_waveg;
 END_B
 BLOCK 11;
  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; W_waveg/2-a1/2-4*W_delta;
 END_B
 BLOCK 12;
  RECT_STAT ABSORBER; 0.0; 5*W_delta; W_waveg/2+a1/2+4*W_delta; W_waveg;
 END_B

#END_TOPOLOGY

#LINK_LIST
 T  5; 0.0; W_waveg/2;
 T  6; 0.0; W_waveg/2-W_delta;
 T  7; 0.0; W_waveg/2;
 T  8; 0.0; 0.0;
 T  9; 0.0; 0.0;
 T 10; 0.0; 0.0;
 T 11; 0.0; 0.0;
 T 12; 0.0; 0.0;
 T 55; 100; 100;
 T 56; 150; 100;
#END_LINK

#OUTPUT
 FILE W_file;
 FIELDS;
 TOPOLOGY;
 FIELD_DISTRIBUTION_M  W_file; W_freq; 8; 9;
#END_OUTPUT

#END_STEP


#EOF


