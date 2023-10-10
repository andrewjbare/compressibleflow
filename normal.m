function out = normal(gamma, inputName, inputValues, varargin)
% Perfect gas normal shock relations
% Functionality based on Compressible Aerodynamics Calculator
% https://devenport.aoe.vt.edu/aoe3114/calc.html
% Andrew Bare
%
% --Usage--
% normal(gamma, inputName, inputValues[, outputName])
%
% --Arguments--
% gamma: Constant ratio of specific heats
% inputName: One of:
%   M1 - Upstream Mach number
%   M2 - Downstream Mach number
%   P2P1 - Ratio of downstream to upstream static pressure
%   R2R1 - Ratio of downstream to upstream density
%   T2T1 - Ratio of downstream to upstream temperature
%   P02P01 - Ratio of downstream to upstream stagnation (total) pressure
%   P1P02 - Ration of upstream static pressure to downstream
%       stagnation (total) pressure
% inputValues: Values of inputName, may be scalar, row or column vector
% outputName (optional): Name of output; see below for output names
%
% --Output--
% If outputName is not provided, output is a table of values in the order:
%   Value (outputName)
%   Upstream Mach number (M1)
%   Downstream Mach number (M2)
%   Ratio of downstream to upstream stagnation (total) pressure (P02P01)
%   Ratio of upstream static pressure to downstream
%       stagnation (total) pressure (P1P02)
%   Ratio of downstream to upstream static pressure (P2P1)
%   Ratio of downstream to upstream density (R2R1)
%   Ratio of downstream to upstream temperature (T2T1)
% If outputName is provided, return only a vertical array of the
% value(s) of outputName in order of inputValues (in parentheses above)

    % --- Input validation & conditioning ---
    if ~isnumeric(inputValues) || ~isnumeric(gamma)
        error('gamma & inputValues must be numeric')
    end
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
        error('InputValue must only contain positive values.');
    end
    possibleInputs = ["M1"; "M2"; "P02P01"; "P1P02";
        "P2P1"; "R2R1"; "T2T1"];
    if isempty(intersect(possibleInputs,inputName))
        error('inputName invalid.');
    end

    % --- Get table of output values & validate inputValues ---
    switch inputName
        case "M1"
            if any(inputValues <= 1)
                error('M1 must be greater than 1');
            end
            outputTable = tableFromMachUpstream(gamma, inputValues);
        case "M2"
            if any(inputValues <= 0.37796447) || any(inputValues >= 1)
                error('M2 must be between 0.37796447 and 1');
            end
            M = M1_M2(gamma,inputValues);
            outputTable = tableFromMachUpstream(gamma, M);
        case "P2P1"
            if any(inputValues <= 1)
                error('P2P1 must be greater than 1');
            end
            M = M1_P2P1(gamma,inputValues);
            outputTable = tableFromMachUpstream(gamma, M);
        case "R2R1"
            if any(inputValues <= 1) || any(inputValues >= 6)
                error('R2R1 must be between 1 and 6');
            end
            M = M1_R2R1(gamma,inputValues);
            outputTable = tableFromMachUpstream(gamma, M);
        case "T2T1"
            if any(inputValues <= 1)
                error('T2T1 must be greater than 1');
            end
            M = M1_T2T1(gamma,inputValues);
            outputTable = tableFromMachUpstream(gamma, M);
        case "P02P01"
            if any(inputValues <= 0) || any(inputValues >= 1)
                error('P02P01 must be between 0 and 1');
            end
            M = M1_P02P01(gamma,inputValues);
            outputTable = tableFromMachUpstream(gamma, M);
        case "P1P02"
            if any(inputValues <= 0) || any(inputValues >= 1)
                error('P1P02 must be between 0 and 1');
            end
            M = M1_P1P02(gamma,inputValues);
            outputTable = tableFromMachUpstream(gamma, M);
        otherwise
            error('InputName invalid.');
    end

    % Return whole table if outputName not provided. Otherwise return
    % single column as array
    nVarargs = length(varargin);
    if nVarargs == 0
        out = outputTable;
    else % nVarargs == 1
        outputName = varargin{1};
        out = outputTable.(outputName);
    end

    function outputTable = tableFromMachUpstream(gamma, M)
        outputTable = table(M, ...
            M2_M1(gamma, M), ...
            P2P1_M1(gamma, M), ...
            R2R1_M1(gamma, M), ...
            T2T1_M1(gamma, M), ...
            P02P01_M1(gamma, M), ...
            P1P02_M1(gamma, M));
        outputTable.Properties.VariableNames = ...
            ["M1","M2","P2P2","R2R1","T2T1","P02P01","P1P02"];
    end

    % For any inputName besides mach, get Mach number and then use the Mach
    % number formulas to compute the rest of the values.

    % M2 from M1
    function M2 = M2_M1(gamma, M1)
        M2 = sqrt((1+(((gamma-1)./2)).*M1.^2)./ ...
            (gamma.*M1.^2-(gamma-1)./2));
    end

    % p2/p1
    function P2P1 = P2P1_M1(gamma, M1)
        P2P1 = (2.*gamma.*M1.^2-gamma+1)./(gamma+1);
    end

    % rho2/rho1 and u1/u2
    function R2R1 = R2R1_M1(gamma, M1)
        R2R1 = ((gamma+1).*(M1.^2))./(2+(gamma-1).*(M1.^2));
    end
    
    % T2/T1 (and h2/h1)
    function T2T1 = T2T1_M1(gamma, M1)
        T2T1 = (1+(M1.^2-1).*(2.*gamma)./(gamma+1)).*((2+(gamma-1).*...
            (M1.^2))./((gamma+1).*(M1.^2)));
    end
    
    % P1/P02
    function P1P02 = P1P02_M1(gamma, M1)
        P1P02 = 1./(((((gamma+1).*(M1.^2))./2).^(gamma./(gamma-1))).*...
            ((gamma+1)./(2.*gamma.*(M1.^2)-(gamma-1))).^(1/(gamma-1)));
    end
    
    % P02/P01
    function P02P01 = P02P01_M1(gamma, M1)
        P02P01 = ((((gamma+1).*(M1.^2))./ ...
            ((gamma-1).*(M1.^2)+2)).^(gamma./(gamma-1))).*...
            ((gamma+1)./(2.*gamma.*(M1.^2)-(gamma-1))).^(1/(gamma-1));
    end

    % inputName to upstream Mach number Relations

    function M1 = M1_M2(gamma, M2)
        M1 = sqrt(gamma.*M2.^2-M2.^2+2)./sqrt(2.*gamma.*M2.^2-gamma+1);
    end

    function M1 = M1_P2P1(gamma, P2P1)
        M1 = sqrt(gamma.*P2P1+gamma+P2P1-1)./(sqrt(2).*sqrt(gamma));
    end

    function M1 = M1_R2R1(gamma, R2R1)
        M1 = (sqrt(2).*sqrt(R2R1))./sqrt(-gamma.*R2R1+gamma+R2R1+1);
    end

    % 0.0001 < T2T1 < 1000
    function M1 = M1_T2T1(gamma, T2T1)
        M1 = zeros(length(T2T1),1);
        for i = 1:length(T2T1)
            fun = @(M) T2T1_M1(gamma, M) - T2T1(i);
            M1(i) = fzero(fun, [1 72]);
        end
        
    end
    
    % 0.0001 < P02/P01 < 1
    function M1 = M1_P02P01(gamma, P02P01)
        M1 = zeros(length(P02P01),1);
        for i = 1:length(P02P01)
            fun = @(M) P02P01_M1(gamma, M) - P02P01(i);
            M1(i) = fzero(fun,[1 21]);
        end
        
    end

    % 0.0001 < P1/P02 < 0.52828178
    function M1 = M1_P1P02(gamma, P02P1)
        M1 = zeros(length(P02P1),1);
        for i = 1:length(P02P1)
            fun = @(M) P1P02_M1(gamma, M) - P02P1(i);
            M1(i) = fzero(fun, [1 89]);
        end
    end
end