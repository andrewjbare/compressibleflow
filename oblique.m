function out = oblique(gamma, M1, inputName, inputValues, varargin)
% Perfect gas oblique shock relations
% Functionality based on Compressible Aerodynamics Calculator
% https://devenport.aoe.vt.edu/aoe3114/calc.html
% Andrew Bare
%
% --Usage--
% oblique(gamma, M1, inputName, inputValues[, outputName])
%
% --Arguments--
% gamma: Constant ratio of specific heats
% M1: Upstream Mach number
% inputName: One of
%   Weak shock turn angle (weak), degrees
%   Strong shock turn angle (strong), degrees
%   Wave angle (wave), degrees
%   Normal component of upstream Mach number (M1n)
% inputValues: Values of inputName, may be scalar, row or column vector
% outputName (optional): Name of output; see below for output names
%
% --Output--
% If outputName is not provided, output is a table of values in the order:
%   Value (outputName)
%   Downstream Mach number (M2), degrees
%   Flow turning angle (theta), degrees
%   Shock wave angle (beta), degrees
%   Ratio of downstream to upstream static pressure (P2P1)
%   Ratio of downstream to upstream density (R2R1)
%   Ratio of downstream to upstream temperature (T2T1)
%   Ratio of downstream to upstream stagnation (total) pressure (P02P01)
%   Normal component of upstream Mach number (M1n)
%   Normal component of downstream Mach number (M2n)
% If outputName is provided, return only a vertical array of the
% value(s) of outputName in order of inputValues (in parentheses above)

    % --- Input validation & conditioning ---
    if ~isnumeric(inputValues) || ~isnumeric(gamma)
        error('gamma & inputValues must be numeric');
    end
    if gamma <= 1
        error('gamma must be greater than 1');
    end
    if isrow(inputValues)
        inputValues = inputValues';
    end
    if isrow(M1)
        M1 = M1';
    end
    if ~isvector(inputValues) && ~isscalar(inputValues)
        error('Matrix of inputValues not accepted.');
    end
    if ~isscalar(inputValues) && ~isscalar(M1)
        error('Either M1 or inputValues must be scalar.');
    end
    if any(inputValues <= 0)
        error('InputValue must only contain positive values.');
    end
    if any(M1 <= 0)
        error('M1 must only contain positive values.');
    end
    possibleInputs = ["weak", "strong", "wave", "M1n"];
    if isempty(intersect(possibleInputs,inputName))
        error('inputName invalid.');
    end

    % --- Get table of output values & validate inputValues ---
    switch inputName
        case "weak"
            if any(inputValues <= 0)
                error('Turn angle must be greater than 0');
            end
            [~,beta] = beta_thetaM(gamma, inputValues, M1);
            M1n = M1n_beta(M1, beta);
            outputTable = tableFromUpstreamNormalMach(gamma, M1, M1n);
        case "strong"
            if any(inputValues <= 0)
                error('Turn angle must be greater than 0');
            end
            beta = beta_thetaM(gamma, inputValues, M1);
            M1n = M1n_beta(M1, beta);
            outputTable = tableFromUpstreamNormalMach(gamma, M1, M1n);
        case "wave"
            if any(inputValues <= mu_M1(M1))
                error('Wave angle must be greater than Mach angle');
            end
            M1n = M1n_beta(M1, inputValues);
            outputTable = tableFromUpstreamNormalMach(gamma, M1, M1n);
        case "M1n"
            if any(inputValues >= M1)
                error('M1n must be between 0 and M1');
            end
            outputTable = tableFromUpstreamNormalMach(gamma, M1, ...
                inputValues);
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

    function outputTable = tableFromUpstreamNormalMach(gamma, M1, M1n)
        if isscalar(M1)
            M1 = repmat(M1,length(M1n),1);
        end
        if isscalar(M1n)
            M1n = repmat(M1n,length(M1),1);
        end

        outputTable = table(M1, ...
            M2_M2n(beta_M1n(M1, M1n), ...
                thetabetaM(gamma, beta_M1n(M1,M1n), M1), ...
                M2n_M1n(gamma, M1n)), ...
            thetabetaM(gamma, beta_M1n(M1,M1n), M1), ...
            beta_M1n(M1,M1n), ...
            P2P1_M1n(gamma, M1n), ...
            R2R1_M1n(gamma, M1n), ...
            T2T1_M1n(gamma, M1n), ...
            P02P01_M1n(gamma, M1n), ...
            M1n, ...
            M2n_M1n(gamma, M1n));
        outputTable.Properties.VariableNames = ...
            ["M1","M2","theta","beta","P2P1","R2R1", ...
                "T2T1","P02P01_M1","M1n","M2n"];
    end

    function mu = mu_M1(M1)
        mu = asind(1./M1);
    end

    % For any inputName besides M1n, get M1n & compute values from there

    % M2n from M1n
    function M2n = M2n_M1n(gamma, M1n)
        M2n = sqrt((1+(((gamma-1)./2)).*M1n.^2)./ ...
            (gamma.*M1n.^2-(gamma-1)./2));
    end

    function M2 = M2_M2n(beta, theta, M2n)
        M2 = M2n./(sind(beta-theta));
    end

    function theta = thetabetaM(gamma, beta, M1)
        theta = atand(2.*cotd(beta).*(M1.^2.*sind(beta).^2-1)./ ...
            (2+(gamma+cosd(2.*beta)).*M1.^2));
    end

    function beta = beta_M1n(M1, M1n)
        beta = asind(M1n./M1);
    end


    % p2/p1
    function P2P1 = P2P1_M1n(gamma, M1n)
        P2P1 = (2.*gamma.*M1n.^2-gamma+1)./(gamma+1);
    end

    % rho2/rho1 and u1/u2
    function R2R1 = R2R1_M1n(gamma, M1n)
        R2R1 = ((gamma+1).*(M1n.^2))./(2+(gamma-1).*(M1n.^2));
    end
    
    % T2/T1 (and h2/h1)
    function T2T1 = T2T1_M1n(gamma, M1n)
        T2T1 = (1+(M1n.^2-1).*(2.*gamma)./(gamma+1)).*((2+(gamma-1).*...
            (M1n.^2))./((gamma+1).*(M1n.^2)));
    end

    function P02P01 = P02P01_M1n(gamma, M1n)
        P02P01 = ((((gamma+1).*(M1n.^2))./ ...
            ((gamma-1).*(M1n.^2)+2)).^(gamma./(gamma-1))).*...
            ((gamma+1)./(2.*gamma.*(M1n.^2)-(gamma-1))).^(1/(gamma-1));
    end

    % --- Functions for finding M1n ---
    function M1n = M1n_beta(M1, beta)
        M1n = M1.*sind(beta);
    end

    function [beta_strong, beta_weak] = beta_thetaM(gamma, theta, M)
        b = -(M.^2+2)./(M.^2)-gamma.*(sind(theta).^2);
        c = (2.*(M.^2)+1)./M.^4+ ...
            (((gamma+1).^2)./4+(gamma-1)./M.^2).*sind(theta).^2;
        d = -(cosd(theta).^2)./(M.^4);
        phi = acos(((9/2).*b.*c-b.^3-(27/2).*d)./ ...
            ((b.^2-3.*c).^(3/2)));
        beta_strong = asind(sqrt(-b./3+(2/3).* ...
            (cos((phi+0*pi)/3)).*(sqrt(b.^2-3*c))));
        beta_weak = asind(sqrt(-b./3+(2/3).* ...
            (cos((phi+4*pi)/3)).*(sqrt(b.^2-3*c))));
end

    
end