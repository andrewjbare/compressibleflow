% Perfect gas isentropic flow relations
% Functionality based on Compressible Aerodynamics Calculator
% https://devenport.aoe.vt.edu/aoe3114/calc.html
% Andrew Bare

% --Usage--
% isentropic(gamma, inputName, inputValues[, outputName])

% --Arguments--
% gamma: Constant ratio of specific heats
% inputName: One of:
%   mach - Mach number
%   TT0 - Total temperature ratio
%   PP0 - Total pressure ratio
%   RR0 - Total density ratio
%   AAsup - Area ratio (supersonic)
%   AAsub - Area ratio (subsonic)
%   mu - Mach angle (degrees)
%   nu - Prandtl-Meyer angle (degrees)
% inputValues: Value of inputName, may be scalar, row or column vector
% outputName (optional): Name of output; see below for output names

% --Output--
% If outputName is not provided, out is a table of values in the order:
% Value (outputName)
% Mach number (Mach)
% Mach angle (mu)
% Prandtl-Meyer angle (PM)
% Stagnation (total) pressure ratio p/p_0 (PP0)
% Stagnation (total) density ratio rho/rho_0 (RR0)
% Stagnation (total) temperature ratio T/T_0 (TT0)
% Sonic (critical) pressure ratio P/P* (PPs)
% Sonic (critical) density ratio rho/rho* (RRs)
% Sonic (critical) temperature ratio T/T* (TTs)
% Area ratio A/A* (AAs)
% If outputName is provided, return only a vertical array of the
% value(s) of outputName in order of inputValues (in parentheses above)

function out = isentropic(gamma, inputName, inputValues, varargin)
    % Possible outputs from function
    possibleOutputs = ["mach"; "TT0"; "PP0"; "RR0"; "AAsup";...
        "AAsub"; "mu"; "nu"];

    if gamma <= 1
        error('gamma must be greater than 1');
    end

    [inputHeight, inputWidth] = size(inputValues);
    if inputHeight == 1 && inputWidth ~= 1
        inputValues = inputValues';
    end
    if inputHeight ~= 1 && inputWidth ~= 1
        error('Matrix of inputValue not accepted.');
    end

    if any(inputValues <= 0)
        error('InputValue must only contain positive values');
    end

    % Set up table of values to return
    outputs = possibleOutputs;
    outputs(outputs == inputName) = [];
    if length(outputs) == length(possibleOutputs)
        error('inputName invalid; only specify one of valid names.');
    end
    switch inputName
        case "mach"
            outputTable = tableFromMach(gamma, inputValues);
        case "TT0"
            if any(inputValues >= 1)
                error('TT0 must be between 0 and 1');
            end
            M = TT0_M(gamma,inputValues);
            outputTable = tableFromMach(gamma, M);
        case "PP0"
            if any(inputValues >= 1)
                error('PP0 must be between 0 and 1');
            end
            M = PP0_M(gamma,inputValues);
            outputTable = tableFromMach(gamma, M);
        case "RR0"
            if any(inputValues >= 1)
                error('RR0 must be between 0 and 1');
            end
            M = RR0_M(gamma,inputValues);
            outputTable = tableFromMach(gamma, M);
        case "AAsup"
            if any(inputValues <= 1)
                error('AAsup must be greater than 1');
            end
            [~,M] = AAs_M(gamma,inputValues);
            outputTable = tableFromMach(gamma, M);
        case "AAsub"
            if any(inputValues <= 1)
                error('AAsub must be greater than 1');
            end
            M = AAs_M(gamma,inputValues);
            outputTable = tableFromMach(gamma, M);
        case "mu"
            if any(inputValues <= 1)
                error('Mach angle mu must be between 0 and 90');
            end
            M = mu_M(inputValues);
            outputTable = tableFromMach(gamma, M);
        case "nu"
            if any(inputValues <= 1)
                error(['Prandtl-Meyer angle nu must be between 0 and' ...
                    '130.454076']);
            end
            M = PM_M(inputValues);
            outputTable = tableFromMach(gamma, M);
        otherwise
            error('InputName invalid.');
    end

    nVarargs = length(varargin);
    if nVarargs == 0
        out = outputTable;
    else % nVarargs == 1
        outputName = varargin{1};
        out = outputTable.(outputName);
    end

    function outputTable = tableFromMach(gamma, M)
        outputTable = table(M, ...
            M_mu(M), ...
            M_PM(gamma, M), ...
            M_PP0(gamma, M), ...
            M_RR0(gamma, M), ...
            M_TT0(gamma, M), ...
            M_PPs(gamma, M), ...
            M_RRs(gamma, M), ...
            M_TTs(gamma, M), ...
            M_AAs(gamma, M));
        outputTable.Properties.VariableNames = ...
            ["Mach","mu","PM","PP0","RR0","TT0","PPs","RRs","TTs","AAs"];
    end

    % For any inputName besides mach, get mach number and then use the mach
    % number formulas to compute the rest of the values.
    
    % Parameter Mach number
    function out = M_mu(M)
        out = asind(1./M);
    end

    function PM = M_PM(gamma, M)
    PM = sqrt((gamma+1)./(gamma-1)).*...
        atand(sqrt(((gamma-1)./(gamma+1)).*(M.^2-1)))...
        -atand(sqrt(M.^2-1));
    end

    function TT0 = M_TT0(gamma, M)
        TT0 = 1./(1+(M.^2).*((gamma-1)./2));
    end
    
    function PP0 = M_PP0(gamma, M)
        PP0 = (1+((gamma-1)./2).*M.^2).^-(gamma./(gamma-1));
    end
    
    function RR0 = M_RR0(gamma, M)
        RR0 = (1+((gamma-1)./2).*M.^2).^-(1./(gamma-1));
    end

    function PPs = M_PPs(gamma, M)
        PPs = M_PP0(gamma, M)./((2./(gamma+1))).^(gamma./(gamma-1));
    end

    function RRs = M_RRs(gamma, M)
        RRs = M_RR0(gamma, M)./((2./(gamma+1))).^(1./(gamma-1));
    end

    function TTs = M_TTs(gamma, M)
        TTs = M_TT0(gamma, M)./(2./(gamma+1));
    end

    function AAs = M_AAs(gamma, M)
        AAs = sqrt((1./(M.^2)).*((2./(gamma+1)) ...
            .*(1+((M.^2).*(gamma-1))./2)).^((gamma+1)/(gamma-1)));
    end

    % inputName to Mach Number Relations

    function M = TT0_M(gamma,TT0)
        M = (sqrt(2).*sqrt(1-TT0))./sqrt(TT0.*(gamma-1));
    end

    function M = PP0_M(gamma, PP0)
        M = (sqrt(2).*sqrt(PP0.^(1/gamma)-PP0))./sqrt(PP0.*(gamma-1));
    end

    function M = RR0_M(gamma, RR0)
        M = (sqrt(2).*sqrt(RR0.^(1-gamma)-1))./sqrt(gamma-1);
    end

    % Use Newton's method
    function [M_sub, M_sup] = AAs_M(gamma, AAs)
        P = 2./(gamma+1);
        Q = 1-P;
        E = 1./Q;
        % sub sup
        R = [AAs.^2 AAs.^(2.*Q./P)];
        a = [P.^(1./Q) Q.^(1./P)];
        r = (R-1)./(2.*a);
        Xn = 1./((1+r)+sqrt(r.*(r+2)));
        for n = 1:7
            Xi = Xn;
            Xn(:,1) = P.*(Xi(:,1)-1)./(1-R(:,1).*(P+Q.*Xi(:,1)).^(-P./Q));
            Xn(:,2) = Q.*(Xi(:,2)-1)./(1-R(:,2).*(Q+P.*Xi(:,2)).^(-Q./P));
        end
        M_sub = sqrt(Xn(:,1));
        M_sup = 1./sqrt(Xn(:,2));
    end

    function M = mu_M(mu)
        M = 1./sind(mu);
    end

    % Use Hall (1975) method
    function M = PM_M(nu)
        nu = deg2rad(nu);
        A = 1.3604;
        B = 0.0962;
        C = -0.5127;
        D = -0.6722;
        E = -0.3278;
        nu_0 = 0.5*pi*(sqrt(6)-1);
        y = (nu./nu_0).^(2/3);
        M = (1+A.*y+B.*y.^2+C.*y.^3)./(1+D.*y+E.*y.^2);
    end

end