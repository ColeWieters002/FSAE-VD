clear;
clc;
close all;

%% LOAD

VehicleParameters;

Tire = TireFunctions();
TireData = LoadTire(TireModel);


%% TEST CONDITIONS

SlipAngle = -15:0.25:15;     % deg
SlipRatio = -0.20:0.005:0.20;
FZ = 150;                    % lbf
IA = -2;                     % deg

FY = zeros(size(SlipAngle));
MZ = zeros(size(SlipAngle));
FX = zeros(size(SlipRatio));


%% FY AND MZ SWEEP

for i = 1:length(SlipAngle)

    X = [SlipAngle(i), FZ, IA];

    FY(i) = Tire.FY(X, TireData.FY_Params, TireData.FNOMIN_lbf);
    MZ(i) = Tire.MZ(X, TireData.MZ_Params, TireData.FY_Params, TireData.FNOMIN_lbf, TireData.R0_ft);

end


%% FX SWEEP

for i = 1:length(SlipRatio)

    X = [SlipRatio(i), FZ, IA];

    FX(i) = Tire.FX(X, TireData.FX_Params, TireData.FNOMIN_lbf);

end


%% FY GRAPH

figure;

plot(SlipAngle, FY, 'LineWidth', 2);
grid on;

xlabel('Slip Angle (deg)');
ylabel('Lateral Force (lbf)');
title('Slip Angle vs Lateral Force');


%% MZ GRAPH

figure;

plot(SlipAngle, MZ, 'LineWidth', 2);
grid on;

xlabel('Slip Angle (deg)');
ylabel('Aligning Moment (ft*lbf)');
title('Slip Angle vs Aligning Moment');


%% FX GRAPH

figure;

plot(SlipRatio, FX, 'LineWidth', 2);
grid on;

xlabel('Slip Ratio');
ylabel('Longitudinal Force (lbf)');
title('Slip Ratio vs Longitudinal Force');


%% FY RESULTS

[MaxFY, MaxFYIndex] = max(FY);
[MinFY, MinFYIndex] = min(FY);


%% FX RESULTS

[MaxFX, MaxFXIndex] = max(FX);
[MinFX, MinFXIndex] = min(FX);


%% PRINT RESULTS

fprintf("\n====================================\n");
fprintf("TIRE SWEEP RESULTS\n");
fprintf("====================================\n");

fprintf("\nLATERAL FORCE\n");
fprintf("Maximum FY = %+.3f lbf at %+.2f deg\n", MaxFY, SlipAngle(MaxFYIndex));
fprintf("Minimum FY = %+.3f lbf at %+.2f deg\n", MinFY, SlipAngle(MinFYIndex));

fprintf("\nLONGITUDINAL FORCE\n");
fprintf("Maximum FX = %+.3f lbf at %+.3f slip ratio\n", MaxFX, SlipRatio(MaxFXIndex));
fprintf("Minimum FX = %+.3f lbf at %+.3f slip ratio\n", MinFX, SlipRatio(MinFXIndex));