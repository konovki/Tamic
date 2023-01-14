#TMC_RT_H

!#define L_input @ (   1.40)
#define W_input @ (   1.00)
#define W_delta @ (   0.02)
#define W_freq @ (  0.6 )

#define a1 @ ( W_input )
#define a2 @ ( 0.2 )

#define h1 @ ( 0.880 )
#define h2 @ ( 6.807 )
#define h3 @ ( 1.493 )
#define h4 @ ( 2.100 )

#define d1  @ ( 2.800 )
#define d2  @ ( 3.715 )
#define d3  @ ( 2.900 )

#define d_a @ ( (d1+d2)/2 )
#define h_a @ ( h1*(d_a-d1)/(d2-d1) )

#define d_a1 @ ( d_a + 6*(W_delta) )
#define h_a1 @ ( h1*(d_a1-d1)/(d2-d1) )

#define l_a @ ( 0.354 )
#define x_a @ ( h1*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )
#define y_a @ ( ((d2-d1)/2)*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )

#define l_r @ ( 0.470 )
#define x_r @ ( h1*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )
#define y_r @ ( ((d2-d1)/2)*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )
#define  d_r @ ( 0.177 )
#define xd_r @ ( ((d2-d1)/2)*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )
#define yd_r @ ( h1*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )
#define x_shift_fakel @ 1.0

#define x_shift @ ( 7.000 )



#define W_type1 @ MAGNETIC
#define W_type  @ METAL


#define W_wavegPlus @ ( 0 )
#define W_waveg @ (  h1 + h2 + h3 + h4 + 10 + W_wavegPlus )
#define L_waveg @ (  d2 + 10 + x_shift )


#STEP

#PARM
 ANGLE_UNIT radian;
 FREQ_UNIT GHz;
 LONG_UNIT m;
 TIME_UNIT ns;
 DELTA W_delta;
 TIME 0.; 500.;
 X_MIN  0.0;
 X_MAX L_waveg;
 Y_MIN  0.0;
 Y_MAX W_waveg;
 FREQ  W_freq;
#END_PARM

#TOPOLOGY

 BLOCK 1;
  INPUT_X  W_type1 ; h_a-2*W_delta; h_a; 0; 720; 1; 0.;
 END_B

 BLOCK 101;
  INPUT_X  W_type1 ; h_a1-2*W_delta; h_a1; 0; 720; 1; 3.141592653589;
 END_B


 BLOCK 4;
  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0.0; W_waveg-2*W_delta;
 END_B
 BLOCK 5;
  RECT_STAT ABSORBER; 1*W_delta; L_waveg-2*W_delta; W_delta; 5*W_delta;
 END_B
 BLOCK 6;
  RECT_STAT ABSORBER; 1*W_delta; L_waveg; W_waveg-5*W_delta; W_waveg;
 END_B
 BLOCK 7;
  RECT_STAT ABSORBER; 0.0; 4*W_delta; 0.0; W_waveg;
 END_B

 BLOCK 8;
  POLYGON_STAT W_type1;
   L   -(d1)/2;           0.0;
   L   -(d2)/2;            h1;
   L   -(d2)/2;         h1+h2;
   L   -(d3)/2;      h1+h2+h3;
   L       0.0;   h1+h2+h3+h4;
   L    (d3)/2;      h1+h2+h3;
   L    (d2)/2;         h1+h2;
   L    (d2)/2;            h1;
   L    (d1)/2;           0.0;
   L   -(d1)/2;           0.0;
 END_B

 BLOCK 9;
  POLYGON_STAT W_type1;
   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;
   L  (L_waveg)/2+x_shift/2-d_a/2-x_a;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a;
   L (L_waveg)/2+x_shift/2-d_a1/2-x_a; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a-W_delta;
   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;
 END_B
 BLOCK 10;
  POLYGON_STAT W_type1;
   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;
   L  (L_waveg)/2+x_shift/2-d_a/2-x_r;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r;
   L (L_waveg)/2+x_shift/2-d_a1/2-x_r; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r-W_delta;
   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;
 END_B


 !BLOCK 1033;
 ! FILE fakel ; -50; 50; 0.; 100.;
! END_B


#END_TOPOLOGY

#LINK_LIST
 T  1; (L_waveg)/2+x_shift/2-d_a/2; (W_waveg)/2-(h1+h2+h3+h4)/2 + (W_wavegPlus)/2;
 T  101; (L_waveg)/2+x_shift/2-d_a1/2; (W_waveg)/2-(h1+h2+h3+h4)/2 + (W_wavegPlus)/2;
 T  4; 0.0; 0.0;
 T  5; 0.0; 0.0;
 T  6; 0.0; 0.0;
 T  7; 0.0; 0.0;
 T  8; (L_waveg)/2+x_shift/2; (W_waveg)/2-(h1+h2+h3+h4)/2 + (W_wavegPlus)/2;
 T  9; 0.0; 0.0 + (W_wavegPlus)/2;
 T 10; -xd_r; yd_r + (W_wavegPlus)/2;
 !T 1033; (L_waveg)/2+x_shift/2; (W_waveg)/2-(h1+h2+h3+h4)/2 ;
#END_LINK

#OUTPUT
 FILE fregat_bez_f;
 FIELDS;
 TOPOLOGY;
 FIELD_DISTRIBUTION_M fregat_bez_f; W_freq; 490; 495;
#END_OUTPUT

#END_STEP


#EOF


