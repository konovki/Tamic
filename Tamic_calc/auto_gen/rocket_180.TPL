#TMC_RT_H
#define L_waveg @ ( 80000.00)
#define W_file  @ rocket_180
#define W_freq  @ (  0.5)
#define W_eps @ ( 4. )
#define W_time @ ( 500. )
#define W_waveg @ (  50000.00)
#define W_input @ (  W_waveg-1000.00)
#define W_delta @ (   20.0)
#define W_r @ (   40.00)
#define W_type @ MAGNETIC
#define a1 @ ( W_input )
#define a11 @ ( 20 )
#define l_rupor @ ( 10 )
#define h_rupor @ ( 0 )

#STEP

#PARM
 ANGLE_UNIT radian;
 FREQ_UNIT GHz;
 LONG_UNIT mm;
 TIME_UNIT ns;
 DELTA W_delta;
 TIME 0.; W_time;
 X_MIN  0.0;
 X_MAX L_waveg;
 Y_MIN  0.0;
 Y_MAX W_waveg;
 FREQ  W_freq;
#END_PARM

#TOPOLOGY

 BLOCK 7;
  INPUT_X  W_type; -a1/2; a1/2; 0; W_time;
 END_B

 BLOCK 5;
  POLYGON_STAT W_type;
   L   0;   a1/2;
   L   2*a11;   a1/2;
   L   2*a11+(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*((l_rupor)/4);
   L   2*a11+2*(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*(2*(l_rupor)/4);
   L   2*a11+3*(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*(3*(l_rupor)/4);
   L   2*a11+l_rupor;   a1/2+h_rupor;
  END_B
 BLOCK 6;
  POLYGON_STAT W_type;
   L   0;   -a1/2;
   L   2*a11;   -a1/2;
   L   2*a11+(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*((l_rupor)/4);
   L   2*a11+2*(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*(2*(l_rupor)/4);
   L   2*a11+3*(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*(3*(l_rupor)/4);
   L   2*a11+l_rupor;   -a1/2-h_rupor;
  END_B

 BLOCK 55;
  FILE W_file; 0; 10000; 0; 10000;
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
 T 55; L_waveg/2; 30000.0;
#END_LINK

#OUTPUT
 FILE W_file;
 FIELDS;
 TOPOLOGY;
 FIELD_DISTRIBUTION_M  W_file; W_freq; 490; 495;
#END_OUTPUT

#END_STEP

#EOF