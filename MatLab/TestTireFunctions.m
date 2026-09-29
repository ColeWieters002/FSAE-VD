clear;
clc;

%% LOAD VEHICLE

VehicleParameters;


%% LOAD TIRE

Tire = TireFunctions();
TireData = LoadTire(TireModel);


%% DISPLAY TIRE DATA

fprintf("\n====================================\n");
fprintf("TIRE DATA\n");
fprintf("====================================\n");

fprintf("FNOMIN = %.4f lbf\n", TireData.FNOMIN_lbf);
fprintf("R0     = %.4f ft\n", TireData.R0_ft);


%% TEST CONDITION

SA = 0.0;       % deg
FZ = 150.0;     % lbf
IA = 0.0;      % deg

X = [SA, FZ, IA];


%% CALCULATE FORCES

FY = Tire.FY(X, TireData.FY_Params, TireData.FNOMIN_lbf);

MZ = Tire.MZ(X, TireData.MZ_Params, TireData.FY_Params, TireData.FNOMIN_lbf, TireData.R0_ft);


%% RESULTS

fprintf("\n====================================\n");
fprintf("TIRE TEST\n");
fprintf("====================================\n");

fprintf("Slip Angle    = %+.2f deg\n", SA);
fprintf("Vertical Load = %.2f lbf\n", FZ);
fprintf("Camber        = %+.2f deg\n", IA);

fprintf("\nFY = %+.6f lbf\n", FY);
fprintf("MZ = %+.6f ft*lbf\n", MZ);