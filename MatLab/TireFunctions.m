function Tire = TireFunctions()
%% TIRE FUNCTIONS
% Pacejka MF5.2 Tire Model
%
% EXTERNAL UNITS:
% Slip Angle     = deg
% Slip Ratio     = decimal
% Vertical Load  = lbf
% Camber         = deg
% Lateral Force  = lbf
% Long. Force    = lbf
% Align. Moment  = ft*lbf
%
% INTERNAL MAGIC FORMULA UNITS:
% Force          = N
% Moment         = N*m
% Angle          = rad

Tire.FY = @FY;
Tire.FX = @FX;
Tire.MZ = @MZ;

end


%% SCALING FACTORS

function L = ScalingFactors()

L.LFZO = 1.0;

L.LCZ = 1.0;

L.LCX = 1.0;
L.LMUX = 1.0;
L.LEX = 1.0;
L.LKX = 1.0;
L.LHX = 1.0;
L.LVX = 1.0;
L.LGAX = 1.0;

L.LCY = 1.0;
L.LMUY = 0.85;
L.LEY = 1.0;
L.LKY = 1.0;
L.LHY = 1.0;
L.LVY = 1.0;
L.LGAY = 1.0;

L.LGAZ = 1.0;
L.LTR = 1.0;
L.LRES = 1.0;
L.LVMX = 1.0;
L.LMY = 1.0;

L.LXAL = 1.0;
L.LYKA = 1.0;
L.LVYKA = 1.0;
L.LS = 1.0;

end


%% LATERAL FORCE

function FY_lbf = FY(X, p, fz0_lbf)

L = ScalingFactors();

% Conversions
LBF2N = 4.4482216153;
N2LBF = 1 / LBF2N;

% Inputs
SA = X(1);       % deg
FZ_lbf = X(2);   % lbf
IA = X(3);       % deg

% Convert to MF5.2 internal units
FZ = FZ_lbf * LBF2N;
fz0 = fz0_lbf * LBF2N;

ALPHA = deg2rad(SA);
GAMMA = deg2rad(IA);

% Parameters
PCY1 = p(1);
PDY1 = p(2);
PDY2 = p(3);
PDY3 = p(4);
PEY1 = p(5);
PEY2 = p(6);
PEY3 = p(7);
PEY4 = p(8);
PKY1 = p(9);
PKY2 = p(10);
PKY3 = p(11);
PHY1 = p(12);
PHY2 = p(13);
PHY3 = p(14);
PVY1 = p(15);
PVY2 = p(16);
PVY3 = p(17);
PVY4 = p(18);

% Magic Formula
FZ0 = fz0 * L.LFZO;
DFZ = (FZ - FZ0) / FZ0;
GAMY = GAMMA * L.LGAY;

SHY = (PHY1 + PHY2 * DFZ) * L.LHY + PHY3 * GAMY;
ALPHAY = ALPHA + SHY;
SVY = FZ * ((PVY1 + PVY2 * DFZ) * L.LVY + (PVY3 + PVY4 * DFZ) * GAMY) * L.LMUY;

CY = PCY1 * L.LCY;
MUY = (PDY1 + PDY2 * DFZ) * (1.0 - PDY3 * GAMY^2) * L.LMUY;
DY = MUY * FZ;

KY = PKY1 * fz0 * sin(2.0 * atan(FZ / (PKY2 * fz0 * L.LFZO))) * (1.0 - PKY3 * abs(GAMY)) * L.LFZO * L.LKY;
BY = KY / (CY * DY + 1e-9);
EY = (PEY1 + PEY2 * DFZ) * (1.0 - (PEY3 + PEY4 * GAMY) * sign(ALPHAY)) * L.LEY;

FY_N = DY * sin(CY * atan(BY * ALPHAY - EY * (BY * ALPHAY - atan(BY * ALPHAY)))) + SVY;

% Output in lbf
FY_lbf = FY_N * N2LBF;

end


%% ALIGNING MOMENT

function MZ_ftlb = MZ(X, q, p, fz0_lbf, r0_ft)

L = ScalingFactors();

% Conversions
LBF2N = 4.4482216153;
FT2M = 0.3048;
NM2FTLB = 0.7375621493;

% Inputs
SA = X(1);       % deg
FZ_lbf = X(2);   % lbf
IA = X(3);       % deg

% Convert to MF5.2 internal units
FZ = FZ_lbf * LBF2N;
fz0 = fz0_lbf * LBF2N;
r0 = r0_ft * FT2M;

ALPHA = deg2rad(SA);
GAMMA = deg2rad(IA);

% Lateral Force Parameters
PCY1 = p(1);
PDY1 = p(2);
PDY2 = p(3);
PDY3 = p(4);
PEY1 = p(5);
PEY2 = p(6);
PEY3 = p(7);
PEY4 = p(8);
PKY1 = p(9);
PKY2 = p(10);
PKY3 = p(11);
PHY1 = p(12);
PHY2 = p(13);
PHY3 = p(14);
PVY1 = p(15);
PVY2 = p(16);
PVY3 = p(17);
PVY4 = p(18);

% Aligning Moment Parameters
QBZ1 = q(1);
QBZ2 = q(2);
QBZ3 = q(3);
QBZ4 = q(4);
QBZ5 = q(5);
QBZ9 = q(6);
QBZ10 = q(7);
QCZ1 = q(8);
QDZ1 = q(9);
QDZ2 = q(10);
QDZ3 = q(11);
QDZ4 = q(12);
QDZ6 = q(13);
QDZ7 = q(14);
QDZ8 = q(15);
QDZ9 = q(16);
QEZ1 = q(17);
QEZ2 = q(18);
QEZ3 = q(19);
QEZ4 = q(20);
QEZ5 = q(21);
QHZ1 = q(22);
QHZ2 = q(23);
QHZ3 = q(24);
QHZ4 = q(25);

% Lateral Force Calculation
FZ0 = fz0 * L.LFZO;
DFZ = (FZ - FZ0) / FZ0;

GAMY = GAMMA * L.LGAY;
GAMZ = GAMMA * L.LGAZ;

SHY = (PHY1 + PHY2 * DFZ) * L.LHY + PHY3 * GAMY;
SVY = FZ * ((PVY1 + PVY2 * DFZ) * L.LVY + (PVY3 + PVY4 * DFZ) * GAMY) * L.LMUY;
ALPHAY = ALPHA + SHY;

KY = PKY1 * fz0 * sin(2.0 * atan(FZ / (PKY2 * fz0 * L.LFZO))) * (1.0 - PKY3 * abs(GAMY)) * L.LFZO * L.LKY;
CY = PCY1 * L.LCY;
MUY = (PDY1 + PDY2 * DFZ) * (1.0 - PDY3 * GAMY^2) * L.LMUY;
DY = MUY * FZ;
BY = KY / (CY * DY + 1e-9);
EY = (PEY1 + PEY2 * DFZ) * (1.0 - (PEY3 + PEY4 * GAMY) * sign(ALPHAY)) * L.LEY;

FY0 = DY * sin(CY * atan(BY * ALPHAY - EY * (BY * ALPHAY - atan(BY * ALPHAY)))) + SVY;

% Pneumatic Trail
SHT = QHZ1 + QHZ2 * DFZ + (QHZ3 + QHZ4 * DFZ) * GAMZ;
ALPHAT = ALPHA + SHT;

SHF = SHY + SVY / (KY + 1e-9);
ALPHAR = ALPHA + SHF;

BT = (QBZ1 + QBZ2 * DFZ + QBZ3 * DFZ^2) * (1.0 + QBZ4 * GAMZ + QBZ5 * abs(GAMZ)) * L.LKY / L.LMUY;
CT = QCZ1;
DT = FZ * (QDZ1 + QDZ2 * DFZ) * (1.0 + QDZ3 * GAMZ + QDZ4 * GAMZ^2) * (r0 / fz0) * L.LTR;
ET = (QEZ1 + QEZ2 * DFZ + QEZ3 * DFZ^2) * (1.0 + (QEZ4 + QEZ5 * GAMZ) * (2.0 / pi) * atan(BT * CT * ALPHAT));

% Pacejka MF5.2 limiter
ET = min(ET, 1.0);

BR = QBZ9 * L.LKY / L.LMUY + QBZ10 * BY * CY;
DR = FZ * ((QDZ6 + QDZ7 * DFZ) * L.LRES + (QDZ8 + QDZ9 * DFZ) * GAMZ) * r0 * L.LMUY;

TRAIL = DT * cos(CT * atan(BT * ALPHAT - ET * (BT * ALPHAT - atan(BT * ALPHAT)))) * cos(ALPHA);
MZR = DR * cos(atan(BR * ALPHAR)) * cos(ALPHA);

MZ_Nm = -TRAIL * FY0 + MZR;

% Output in ft*lbf
MZ_ftlb = MZ_Nm * NM2FTLB;

end


%% LONGITUDINAL FORCE

function FX_lbf = FX(X, p, fz0_lbf)

L = ScalingFactors();

% Conversions
LBF2N = 4.4482216153;
N2LBF = 1 / LBF2N;

% Inputs
SR = X(1);
FZ_lbf = X(2);
IA = X(3);

% Convert to MF5.2 internal units
FZ = FZ_lbf * LBF2N;
fz0 = fz0_lbf * LBF2N;

KAPPA = SR;
GAMMA = deg2rad(IA);

% Parameters
PCX1 = p(1);
PDX1 = p(2);
PDX2 = p(3);
PDX3 = p(4);
PEX1 = p(5);
PEX2 = p(6);
PEX3 = p(7);
PEX4 = p(8);
PKX1 = p(9);
PKX2 = p(10);
PKX3 = p(11);
PHX1 = p(12);
PHX2 = p(13);
PVX1 = p(14);
PVX2 = p(15);

% Magic Formula
FZ0 = fz0 * L.LFZO;
DFZ = (FZ - FZ0) / FZ0;
GAMX = GAMMA * L.LGAX;

SHX = (PHX1 + PHX2 * DFZ) * L.LHX;
KAPPAX = KAPPA + SHX;

SVX = FZ * (PVX1 + PVX2 * DFZ) * L.LVX * L.LMUX;

CX = PCX1 * L.LCX;
MUX = (PDX1 + PDX2 * DFZ) * (1.0 - PDX3 * GAMX^2) * L.LMUX;
DX = MUX * FZ;

KX = FZ * (PKX1 + PKX2 * DFZ) * exp(PKX3 * DFZ) * L.LKX;
BX = KX / (CX * DX + 1e-9);
EX = (PEX1 + PEX2 * DFZ + PEX3 * DFZ^2) * (1.0 - PEX4 * sign(KAPPAX)) * L.LEX;

FX_N = DX * sin(CX * atan(BX * KAPPAX - EX * (BX * KAPPAX - atan(BX * KAPPAX)))) + SVX;

% Output in lbf
FX_lbf = FX_N * N2LBF;

end