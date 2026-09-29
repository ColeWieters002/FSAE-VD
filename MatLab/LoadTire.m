function TireData = LoadTire(tir_path)
%% LOAD TIRE
% Reads a Pacejka .tir file and loads the coefficients needed by
% TireFunctions.m
%
% OUTPUT UNITS:
% FNOMIN = lbf
% R0     = ft

%% CONVERSIONS

LBF2N = 4.4482216153;
FT2M = 0.3048;


%% PARAMETER NAMES

FY_Names = ["PCY1","PDY1","PDY2","PDY3","PEY1","PEY2","PEY3","PEY4", ...
            "PKY1","PKY2","PKY3","PHY1","PHY2","PHY3", ...
            "PVY1","PVY2","PVY3","PVY4"];

MZ_Names = ["QBZ1","QBZ2","QBZ3","QBZ4","QBZ5","QBZ9","QBZ10","QCZ1", ...
            "QDZ1","QDZ2","QDZ3","QDZ4","QDZ6","QDZ7","QDZ8","QDZ9", ...
            "QEZ1","QEZ2","QEZ3","QEZ4","QEZ5","QHZ1","QHZ2","QHZ3","QHZ4"];

FX_Names = ["PCX1","PDX1","PDX2","PDX3","PEX1","PEX2","PEX3","PEX4", ...
            "PKX1","PKX2","PKX3","PHX1","PHX2","PVX1","PVX2"];


%% PARSE .TIR FILE

Sections = ParseTIR(tir_path);


%% LATERAL COEFFICIENTS

TireData.FY_Params = zeros(1, length(FY_Names));

for i = 1:length(FY_Names)
    TireData.FY_Params(i) = Sections.LATERAL_COEFFICIENTS.(FY_Names(i));
end


%% ALIGNING MOMENT COEFFICIENTS

TireData.MZ_Params = zeros(1, length(MZ_Names));

for i = 1:length(MZ_Names)
    TireData.MZ_Params(i) = Sections.ALIGNING_COEFFICIENTS.(MZ_Names(i));
end


%% LONGITUDINAL COEFFICIENTS

TireData.FX_Params = zeros(1, length(FX_Names));

for i = 1:length(FX_Names)
    TireData.FX_Params(i) = Sections.LONGITUDINAL_COEFFICIENTS.(FX_Names(i));
end


%% NOMINAL LOAD

FNOMIN_N = Sections.VERTICAL.FNOMIN;

TireData.FNOMIN_lbf = FNOMIN_N / LBF2N;


%% UNLOADED RADIUS

R0_m = Sections.DIMENSION.UNLOADED_RADIUS;

TireData.R0_ft = R0_m / FT2M;


%% FILE INFORMATION

TireData.Path = tir_path;

end


%% TIR PARSER

function Sections = ParseTIR(path)

fid = fopen(path, 'r');

if fid == -1
    error('Could not open tire file: %s', path);
end

Sections = struct();
CurrentSection = '';

while ~feof(fid)

    Line = strtrim(fgetl(fid));

    % Remove comments
    CommentLocation = strfind(Line, '!');

    if ~isempty(CommentLocation)
        Line = strtrim(Line(1:CommentLocation(1)-1));
    end

    % Skip empty lines
    if isempty(Line)
        continue
    end

    % Detect section
    if Line(1) == '[' && Line(end) == ']'

        CurrentSection = upper(strtrim(Line(2:end-1)));
        CurrentSection = matlab.lang.makeValidName(CurrentSection);

        if ~isfield(Sections, CurrentSection)
            Sections.(CurrentSection) = struct();
        end

        continue
    end

    % Skip anything outside a section
    if isempty(CurrentSection)
        continue
    end

    % Find equals sign
    EqualLocation = strfind(Line, '=');

    if isempty(EqualLocation)
        continue
    end

    Key = strtrim(Line(1:EqualLocation(1)-1));
    ValueText = strtrim(Line(EqualLocation(1)+1:end));

    Key = matlab.lang.makeValidName(Key);

    % Remove quotes
    ValueText = strrep(ValueText, '''', '');
    ValueText = strrep(ValueText, '"', '');

    % Convert to number if possible
    ValueNumber = str2double(ValueText);

    if ~isnan(ValueNumber)
        Value = ValueNumber;
    else
        Value = ValueText;
    end

    Sections.(CurrentSection).(Key) = Value;

end

fclose(fid);

end
