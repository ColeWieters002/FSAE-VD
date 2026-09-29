%% VEHICLE PARAMETERS
% FSAE Vehicle Dynamics
%
% STANDARD UNITS:
% Geometry       = in
% Force          = lbf
% Moment         = ft*lbf
% Spring Rate    = lbf/in
% Roll Stiffness = ft*lbf/rad
% Mass           = slug
% Inertia        = slug*ft^2
% Velocity       = ft/s
% Acceleration   = ft/s^2
% Angles         = rad internally
% Display Angles = deg

%% FILE PATHS

ThisFile = mfilename("fullpath");
ThisDir = fileparts(ThisFile);


%% CONSTANTS

Gravity = 32.174; % ft/s^2

RAD2DEG = 180 / pi;
DEG2RAD = pi / 180;

MPH2FPS = 1.46667;
FPS2MPH = 1 / MPH2FPS;


%% CONVERSIONS

IN2FT = 1 / 12;
LBM2SLUG = 1 / Gravity;

%% MASS / INERTIA

% Original measured values converted from kg to lbm
Vehicle_lbm = 388.0131; % lbm
Driver_lbm  = 167.5511;  % lbm

TotalMass_lbm = Vehicle_lbm + Driver_lbm;

UnsprungMass_lbm = 92.5940;  % lbm
SprungMass_lbm = TotalMass_lbm - UnsprungMass_lbm;

% Convert lbm to slugs for dynamics
Vehicle_slug = Vehicle_lbm * LBM2SLUG;
Driver_slug = Driver_lbm * LBM2SLUG;
TotalMass_slug = TotalMass_lbm * LBM2SLUG;

UnsprungMass_slug = UnsprungMass_lbm * LBM2SLUG;
SprungMass_slug = SprungMass_lbm * LBM2SLUG;

% Front static weight distribution
WeightDist = 0.48;

% Inertia
RollInertia_RA = 15.0; % slug*ft^2,
YawInertia = 67.8557; % slug*ft^2

% SprungPitchInertia

% Calculate Static Weight
VehicleWeight_lbf = Vehicle_slug * Gravity;
DriverWeight_lbf = Driver_slug * Gravity;
TotalWeight_lbf = TotalMass_slug * Gravity;
FrontWeight_lbf = TotalWeight_lbf * WeightDist;
RearWeight_lbf  = TotalWeight_lbf * (1 - WeightDist);
StaticFZ_FL = FrontWeight_lbf / 2;
StaticFZ_FR = FrontWeight_lbf / 2;
StaticFZ_RL = RearWeight_lbf / 2;
StaticFZ_RR = RearWeight_lbf / 2;

%% CG

CG_in = 11.25;         % in
CG_ft = CG_in * IN2FT; % ft


%% WHEEL SPACING

Wheelbase_in = 60.5;                % in

FTrackwidth_in = 47.0;              % in
RTrackwidth_in = 46.0;              % in

% Feet versions used in rigid-body dynamics
Wheelbase_ft = Wheelbase_in * IN2FT;

FTrackwidth_ft = FTrackwidth_in * IN2FT;
RTrackwidth_ft = RTrackwidth_in * IN2FT;


%% CG LONGITUDINAL POSITION

% WeightDist = fraction of static weight on FRONT axle
%
% Front weight fraction = b / L
% Rear weight fraction  = a / L

b_in = WeightDist * Wheelbase_in;
a_in = Wheelbase_in - b_in;

b_ft = b_in * IN2FT;
a_ft = a_in * IN2FT;


%% ACTUATION

FrontSpringRate = 820;       % lbf/in
RearSpringRate  = 640;       % lbf/in

FrontRollStiffness = 18875;  % ft*lbf/rad
RearRollStiffness  = 22285;  % ft*lbf/rad

FrontHeaveStiffness = 400;   % lbf/in
RearHeaveStiffness  = 375;   % lbf/in

TotalRollStiffness = FrontRollStiffness + RearRollStiffness;

FrontRollStiffnessDistribution = FrontRollStiffness / TotalRollStiffness;

RearRollStiffnessDistribution = RearRollStiffness / TotalRollStiffness;

%% Roll Damper Kinematics

FrontRollMotionRatio = 0.77;
RearRollMotionRatio  = 0.1213;

FrontRollDamperLever_in = FrontRollMotionRatio * FTrackwidth_in / 2;
RearRollDamperLever_in  = RearRollMotionRatio  * RTrackwidth_in / 2;


%% CAMBER

CamberBounds = [0, -2, -3];   % deg

% Corresponding suspension travel
CamberTravel_in = [-1, 0, 1]; % in

FrontStaticCamber_deg = 0; %deg
RearStaticCamber_deg = 0; % deg

FrontCamberGainRoll_deg_deg = -0.57; %deg/deg 
RearCamberGainRoll_deg_deg = -0.49; %deg/deg 

FrontCamberGainHeave_deg_in = -0; %deg/in 
RearCamberGainHeave_deg_in = -0; %deg/in


%% ROLL CENTER

FrontRollCenter_in = -0.02; % in
RearRollCenter_in  = 2.474; % in

FrontRollCenter_ft = FrontRollCenter_in * IN2FT;
RearRollCenter_ft  = RearRollCenter_in * IN2FT;

RollAxisHeight_ft = (FrontRollCenter_ft*b_ft + RearRollCenter_ft*a_ft) / Wheelbase_ft;
RollMomentArm_ft = CG_ft - RollAxisHeight_ft;


%% STEERING

Ackerman = 0.0;

FrontStaticToe_deg = 0.0;
RearStaticToe_deg  = 0.0;

FrontToeGainRoll_deg_deg = 0; %deg/deg
RearToeGainRoll_deg_deg = 0; %deg/deg

FrontToeGainHeave_deg_in = 0; %deg/in
RearToeGainHeave_deg_in = 0; %deg/in


%% CASTER

FrontCaster_deg = 8.23;
RearCaster_deg  = 2.81;


%% KPI

FrontKPI_deg = 0.0;
RearKPI_deg  = 14.9;


%% AERODYNAMICS

AirDensity = 0.002377; % slug/ft^3

AeroArea_ft2 = 10.7639; % ft^2

CL = -3.75;
CD = 1.4;

AeroBalance = 0.40; % Front


%% TIRES

TirePressure_psi = 12.0; % psi

RollResistance = 0.4;

TireModel = fullfile(ThisDir,"Tires","Hoosier_16x75_10_R20.tir");


























%% MULTIMATIC DAMPER DATA
% Source data:
%   Velocity = mm/s
%   Force    = N
%
% Converted database:
%   Velocity = in/s
%   Force    = lbf
%
% Matrix format:
%   ROW    = damper velocity
%   COLUMN = adjuster position 0-11

N2LBF     = 0.224808943;
MMPS2INPS = 1 / 25.4;

%% Velocity Breakpoints

DamperVelocity_mm_s = [10 25 50 100 150 200 250 300 350 400 450 500];

DamperVelocity_in_s = DamperVelocity_mm_s * MMPS2INPS;


%% ================================================================
%  HEAVE DAMPER
%  ================================================================

% ---------------- VC01 ----------------

Heave.VC01.Compression_N = [ ...
17 17 17 18 19 20 23 28 41 68 114 179;
23 24 27 30 35 43 56 78 116 173 244 328;
44 50 58 69 84 106 137 179 235 304 385 475;
117 135 157 185 219 261 311 370 437 512 594 683;
215 244 278 318 364 416 475 540 612 689 772 859;
325 362 405 452 505 564 628 697 771 849 932 1019;
440 484 532 585 643 705 772 844 919 999 1082 1168;
556 605 658 715 776 841 910 983 1060 1140 1223 1310;
673 726 782 842 906 973 1044 1118 1195 1276 1359 1445;
790 846 904 967 1032 1101 1173 1248 1326 1407 1491 1577;
903 964 1025 1089 1156 1227 1300 1376 1454 1535 1619 1705;
1003 1073 1143 1210 1279 1350 1424 1500 1579 1661 1744 1830];

Heave.VC01.Extension_N = [ ...
17 17 18 18 19 21 23 28 43 85 159 249;
25 27 29 33 38 46 59 78 108 163 235 319;
50 57 66 78 94 115 143 179 224 279 342 416;
137 156 179 206 238 276 320 369 425 486 552 623;
246 275 308 345 386 433 483 538 598 662 729 800;
363 399 439 482 529 580 634 693 754 819 887 958;
482 523 567 614 665 719 776 836 899 965 1034 1104;
599 644 691 742 795 851 910 971 1035 1102 1171 1242;
715 762 812 865 920 978 1038 1100 1165 1232 1301 1372;
828 878 930 984 1041 1100 1161 1224 1289 1357 1426 1497;
939 991 1044 1100 1158 1218 1279 1343 1409 1477 1546 1617;
1048 1101 1156 1213 1271 1332 1395 1459 1525 1593 1662 1733];

% ---------------- VC02 ----------------

Heave.VC02.Compression_N = [ ...
13 16 17 18 20 22 26 31 40 54 75 101;
17 19 22 28 38 50 59 72 88 108 131 156;
22 27 35 48 66 87 111 135 157 179 203 228;
39 48 63 82 105 131 158 187 217 248 280 311;
60 73 89 110 134 161 190 220 250 282 314 347;
82 98 116 137 162 188 217 246 277 309 341 374;
107 123 142 164 189 215 243 272 303 334 366 399;
132 149 169 192 216 242 270 299 329 360 391 424;
158 176 197 219 244 270 297 326 355 386 417 449;
185 204 225 248 272 298 325 354 383 413 444 476;
213 233 254 277 301 327 354 382 411 441 472 503;
242 262 284 307 332 357 384 412 441 471 501 532];

Heave.VC02.Extension_N = [ ...
11 14 16 16 18 21 24 29 35 44 59 82;
16 17 18 21 29 42 55 71 91 112 137 165;
19 21 25 34 50 75 104 135 164 196 230 265;
30 37 47 63 87 118 154 193 234 276 320 363;
46 56 71 92 119 151 187 227 269 313 357 402;
64 78 96 120 148 180 216 255 297 340 384 430;
84 101 121 146 175 208 244 282 323 365 409 453;
105 123 146 172 201 234 270 308 348 390 433 477;
126 146 170 197 227 260 295 333 373 414 456 500;
147 169 193 221 252 285 320 358 397 438 480 523;
169 191 217 245 276 309 345 382 421 461 503 545;
191 214 240 269 300 333 369 406 444 485 526 568];

% ---------------- VC03 ----------------

Heave.VC03.Compression_N = [ ...
16 16 17 17 18 19 21 25 34 54 97 160;
19 21 23 26 30 37 48 64 89 126 177 239;
31 36 43 53 66 85 110 142 180 227 281 342;
70 85 104 127 157 192 233 279 330 386 445 507;
121 145 173 207 245 289 337 389 445 504 567 631;
177 208 243 282 326 375 427 482 541 603 666 733;
236 271 310 354 402 453 508 566 626 689 754 821;
296 333 376 422 472 526 582 642 703 767 833 900;
356 396 440 488 540 595 653 713 775 840 906 974;
417 459 504 553 606 661 720 781 844 909 976 1044;
478 521 568 618 671 727 785 846 910 975 1042 1110;
540 584 632 682 736 792 850 911 974 1040 1106 1175];

Heave.VC03.Extension_N = [ ...
16 16 16 17 18 19 23 29 41 63 103 166;
18 20 21 24 29 38 53 76 112 159 218 286;
27 31 38 48 64 86 119 162 215 277 346 422;
58 71 90 116 149 191 241 299 364 435 510 590;
101 123 152 189 233 284 343 407 476 550 628 709;
151 180 216 259 309 366 429 497 569 644 723 805;
204 238 279 326 380 440 505 575 649 725 805 887;
258 297 341 391 447 508 575 646 720 798 878 960;
313 355 402 454 512 574 641 712 786 864 944 1027;
367 412 461 515 574 637 704 775 849 927 1007 1089;
421 468 519 575 635 698 766 837 911 987 1067 1148;
475 524 576 633 694 758 826 897 970 1047 1125 1206];


%% ================================================================
%  ROLL DAMPER
%  ================================================================

% ---------------- VC01 ----------------

Roll.VC01.Compression_N = [ ...
12 13 15 15 16 16 16 17 18 21 33 61;
16 17 17 18 18 19 21 24 31 42 64 95;
20 21 23 25 27 31 37 46 60 79 106 138;
35 39 44 50 58 68 81 98 120 145 175 208;
57 64 73 83 95 110 129 150 175 202 233 267;
84 94 106 120 136 154 175 199 225 255 286 320;
115 127 141 158 176 197 220 246 273 303 335 369;
147 162 178 196 217 239 264 290 319 349 382 415;
180 197 215 235 257 281 306 334 363 394 426 460;
214 232 252 273 296 321 347 376 405 436 469 503;
249 268 289 311 335 361 388 416 446 478 511 545;
283 304 325 349 374 400 427 456 487 518 551 585];

% Roll VC01 compression and extension are identical in source data.
Roll.VC01.Extension_N = Roll.VC01.Compression_N;

% ---------------- VC02 ----------------

Roll.VC02.Compression_N = [ ...
13 15 16 16 17 19 22 27 33 42 53 64;
16 17 19 21 25 31 39 47 56 66 77 89;
20 23 27 32 39 47 57 67 79 90 103 115;
31 36 43 51 60 70 80 92 104 116 129 142;
42 49 57 66 76 86 98 109 121 134 147 160;
53 61 70 79 89 100 112 124 136 149 161 175;
64 72 81 91 102 113 125 137 149 161 174 187;
74 83 92 103 113 125 136 148 161 173 186 199;
84 93 103 114 125 136 148 160 172 185 198 211;
95 104 114 124 135 147 159 171 183 196 209 222;
105 115 124 135 146 157 169 181 194 207 219 232;
116 125 135 146 157 168 180 192 204 217 230 243];

Roll.VC02.Extension_N = Roll.VC02.Compression_N;

% ---------------- VC03 ----------------

Roll.VC03.Compression_N = [ ...
11 13 15 15 16 16 17 19 22 30 44 64;
16 16 17 18 19 21 25 31 42 56 75 96;
19 21 22 25 28 34 43 55 70 89 109 132;
31 35 40 46 54 65 78 95 113 134 156 179;
48 54 62 71 82 96 111 128 147 168 190 214;
68 76 86 97 110 125 141 159 178 199 221 244;
89 99 110 123 137 152 169 188 207 228 250 273;
111 122 134 148 163 179 197 215 235 256 278 300;
133 145 158 173 188 205 223 242 262 283 304 327;
155 168 182 197 213 231 249 268 288 309 330 352;
178 191 206 222 238 255 274 293 313 334 356 378;
201 215 230 246 262 280 299 318 338 359 381 403];

Roll.VC03.Extension_N = Roll.VC03.Compression_N;


%% ================================================================
%  CORNER DAMPER
%  ================================================================

% ---------------- VC01 ----------------

Corner.VC01.Compression_N = [ ...
17 17 17 18 18 19 20 23 27 35 56 86;
24 25 27 30 33 39 47 59 74 96 125 160;
48 54 61 69 80 93 110 131 156 185 219 256;
124 138 154 172 192 215 241 270 301 335 372 410;
216 235 257 280 306 333 363 395 429 465 502 542;
312 335 360 387 415 445 477 511 546 583 621 661;
408 434 461 490 520 552 585 620 656 694 732 772;
503 531 560 590 622 655 689 725 761 799 838 878;
597 626 657 688 721 755 790 826 863 901 940 980;
690 721 752 784 818 853 888 925 962 1000 1040 1080;
783 814 846 879 913 948 984 1021 1059 1097 1137 1177;
874 906 939 973 1007 1043 1079 1116 1154 1193 1233 1273];

Corner.VC01.Extension_N = [ ...
17 17 18 18 19 20 21 24 29 38 57 83;
25 27 29 32 36 43 51 62 78 98 123 153;
52 58 66 75 86 99 115 134 157 182 211 243;
132 146 161 178 197 218 241 267 294 323 354 387;
223 241 261 282 305 330 356 384 413 443 475 508;
314 336 358 382 407 433 461 490 520 551 584 617;
403 427 451 476 503 531 559 589 620 652 684 718;
490 515 540 567 594 623 652 683 714 746 779 812;
574 600 626 654 682 711 741 772 803 836 869 902;
656 682 710 738 767 796 827 858 889 922 955 989;
736 763 791 819 849 879 909 941 973 1005 1038 1072;
813 841 870 899 928 959 990 1021 1053 1086 1119 1153];

% ---------------- VC02 ----------------

Corner.VC02.Compression_N = [ ...
17 19 23 30 37 47 61 79 103 129 159 190;
24 32 45 62 83 106 128 150 174 200 227 256;
40 55 74 97 123 151 180 210 241 272 304 336;
73 93 117 142 170 199 230 261 293 325 358 391;
102 124 149 176 204 234 264 296 328 360 393 426;
128 152 178 205 234 263 294 325 357 390 423 456;
154 178 204 232 261 291 321 353 385 417 450 483;
179 204 230 258 287 317 348 379 411 443 476 509;
204 229 256 284 313 343 374 405 437 469 501 534;
229 255 282 310 339 369 400 431 463 495 527 560;
254 281 308 336 365 395 426 457 489 521 553 586;
281 307 335 363 392 422 453 484 516 548 580 613];

Corner.VC02.Extension_N = [ ...
17 19 22 29 36 45 59 77 98 122 148 174;
23 30 41 56 74 94 114 133 155 177 201 226;
38 51 67 86 108 131 156 181 207 233 260 287;
66 83 103 124 148 172 197 223 250 277 305 333;
89 108 129 152 175 200 226 252 279 306 333 361;
110 130 151 174 198 223 249 275 302 329 356 384;
128 149 171 194 218 243 269 295 322 349 376 404;
145 166 189 212 236 261 287 313 340 367 394 422;
161 183 205 229 253 278 304 330 357 384 411 438;
177 199 221 245 269 294 320 346 373 400 427 454;
192 214 237 261 285 310 336 362 388 415 442 470;
207 229 252 276 300 325 351 377 404 430 457 485];

% ---------------- VC03 ----------------

Corner.VC03.Compression_N = [ ...
17 17 17 18 18 19 20 23 27 35 49 68;
24 25 27 30 33 39 45 54 65 81 99 121;
47 52 57 65 73 83 96 110 127 146 167 190;
113 123 135 148 162 178 195 214 234 255 277 301;
186 200 215 231 248 266 285 306 327 349 372 396;
261 277 294 312 330 349 370 391 413 435 459 483;
335 352 370 389 409 429 450 471 494 517 540 565;
408 426 445 465 485 506 527 549 572 595 619 643;
480 499 518 538 559 580 602 624 647 671 695 719;
551 571 591 611 632 654 676 698 722 745 769 793;
622 642 662 683 705 726 749 771 795 818 842 867;
693 713 734 755 776 798 821 844 867 891 915 939];

Corner.VC03.Extension_N = [ ...
17 17 18 18 19 20 21 24 29 37 49 65;
25 27 29 32 36 41 48 56 67 81 97 116;
50 55 61 68 77 87 98 111 126 143 161 180;
117 127 138 150 163 178 193 209 226 244 264 284;
189 201 215 229 244 260 277 294 312 331 351 371;
258 272 287 303 319 336 353 371 390 409 429 449;
325 340 356 372 389 407 424 443 462 481 501 522;
389 405 422 438 456 474 492 511 530 550 570 590;
451 468 485 502 520 538 557 576 595 615 635 655;
512 529 546 564 582 600 619 638 658 678 698 718;
571 588 606 624 642 660 679 699 718 738 758 779;
628 646 664 682 700 719 738 758 777 797 818 838];


%% ================================================================
%  CONVERT ALL FORCE TABLES TO LBF
%  ================================================================

ValveCodes = {'VC01','VC02','VC03'};

for i = 1:length(ValveCodes)

    VC = ValveCodes{i};

    Heave.(VC).Compression_lbf = ...
        Heave.(VC).Compression_N * N2LBF;

    Heave.(VC).Extension_lbf = ...
        Heave.(VC).Extension_N * N2LBF;

    Roll.(VC).Compression_lbf = ...
        Roll.(VC).Compression_N * N2LBF;

    Roll.(VC).Extension_lbf = ...
        Roll.(VC).Extension_N * N2LBF;

    Corner.(VC).Compression_lbf = ...
        Corner.(VC).Compression_N * N2LBF;

    Corner.(VC).Extension_lbf = ...
        Corner.(VC).Extension_N * N2LBF;

end

clear i VC ValveCodes