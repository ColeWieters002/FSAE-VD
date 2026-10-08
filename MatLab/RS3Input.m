function [RS3GPS_Speed, RS3LateralAcc, RS3YawRate, RS3SteeringSignal] = RS3Input(path)

    data = readmatrix(path, ...
        "Sheet", "2.csv", ...
        "Range", "A18:Q2157");

    RS3Time_s = data(:,1);
    RS3Time_s = RS3Time_s - RS3Time_s(1);

    GPS_kph  = data(:,2);
    LatAcc_g = data(:,15);
    Yaw_deg  = data(:,16);

    Steering_deg = ...
        (data(:,17) - 2410) * ...
        (130 / (4390 - 2410)) / 4.8;

    RS3GPS_Speed       = timeseries(GPS_kph, RS3Time_s);
    RS3LateralAcc      = timeseries(LatAcc_g, RS3Time_s);
    RS3YawRate         = timeseries(Yaw_deg, RS3Time_s);
    RS3SteeringSignal  = timeseries(Steering_deg, RS3Time_s);

end