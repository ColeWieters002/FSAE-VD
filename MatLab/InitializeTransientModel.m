%% INITIALIZE TRANSIENT VEHICLE MODEL

clear;
clc;

%% VEHICLE PARAMETERS

VehicleParameters;
[RS3GPS_Speed, RS3LateralAcc, RS3YawRate, RS3SteeringSignal] = RS3Input("RS3 Data/RS3 Dobbins Skidpad2.xlsx");




%% LOAD TIRE DATA

TireData = LoadTire(TireModel);

%% TIRE GENERAL PARAMETERS

Tire.FNOMIN_lbf = TireData.FNOMIN_lbf;
Tire.R0_ft = TireData.R0_ft;

%% TIRE SCALING FACTORS

Tire.LFZO = 1;
Tire.LCZ  = 1;

Tire.LCX  = 1;
Tire.LMUX = 1;
Tire.LEX  = 1;
Tire.LKX  = 1;
Tire.LHX  = 1;
Tire.LVX  = 1;
Tire.LGAX = 1;

Tire.LCY  = 1;
Tire.LMUY = 0.85;
Tire.LEY  = 1;
Tire.LKY  = 1;
Tire.LHY  = 1;
Tire.LVY  = 1;
Tire.LGAY = 1;

Tire.LGAZ = 1;
Tire.LTR  = 1;
Tire.LRES = 1;
Tire.LVMX = 1;
Tire.LMY  = 1;

Tire.LXAL  = 1;
Tire.LYKA  = 1;
Tire.LVYKA = 1;
Tire.LS    = 1;


%% LATERAL FORCE COEFFICIENTS

p = TireData.FY_Params;

Tire.PCY1 = p(1);

Tire.PDY1 = p(2);
Tire.PDY2 = p(3);
Tire.PDY3 = p(4);

Tire.PEY1 = p(5);
Tire.PEY2 = p(6);
Tire.PEY3 = p(7);
Tire.PEY4 = p(8);

Tire.PKY1 = p(9);
Tire.PKY2 = p(10);
Tire.PKY3 = p(11);

Tire.PHY1 = p(12);
Tire.PHY2 = p(13);
Tire.PHY3 = p(14);

Tire.PVY1 = p(15);
Tire.PVY2 = p(16);
Tire.PVY3 = p(17);
Tire.PVY4 = p(18);


%% LONGITUDINAL FORCE COEFFICIENTS

p = TireData.FX_Params;

Tire.PCX1 = p(1);

Tire.PDX1 = p(2);
Tire.PDX2 = p(3);
Tire.PDX3 = p(4);

Tire.PEX1 = p(5);
Tire.PEX2 = p(6);
Tire.PEX3 = p(7);
Tire.PEX4 = p(8);

Tire.PKX1 = p(9);
Tire.PKX2 = p(10);
Tire.PKX3 = p(11);

Tire.PHX1 = p(12);
Tire.PHX2 = p(13);

Tire.PVX1 = p(14);
Tire.PVX2 = p(15);


%% ALIGNING MOMENT COEFFICIENTS

p = TireData.MZ_Params;

Tire.QBZ1  = p(1);
Tire.QBZ2  = p(2);
Tire.QBZ3  = p(3);
Tire.QBZ4  = p(4);
Tire.QBZ5  = p(5);
Tire.QBZ9  = p(6);
Tire.QBZ10 = p(7);

Tire.QCZ1 = p(8);

Tire.QDZ1 = p(9);
Tire.QDZ2 = p(10);
Tire.QDZ3 = p(11);
Tire.QDZ4 = p(12);
Tire.QDZ6 = p(13);
Tire.QDZ7 = p(14);
Tire.QDZ8 = p(15);
Tire.QDZ9 = p(16);

Tire.QEZ1 = p(17);
Tire.QEZ2 = p(18);
Tire.QEZ3 = p(19);
Tire.QEZ4 = p(20);
Tire.QEZ5 = p(21);

Tire.QHZ1 = p(22);
Tire.QHZ2 = p(23);
Tire.QHZ3 = p(24);
Tire.QHZ4 = p(25);

clear p

%% BICYCLE MODEL STATIC TIRE LOADS

FZF_lbf = FrontWeight_lbf / 2;
FZR_lbf = RearWeight_lbf / 2;


%% INITIAL CAMBER

FrontCamber_deg = 0;
RearCamber_deg = 0;