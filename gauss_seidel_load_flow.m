%% ================================================================
%              GAUSS-SEIDEL LOAD FLOW ANALYSIS
%                   Power Systems Assignment
%
%  Generalized MATLAB program for load-flow analysis using the
%  Gauss-Seidel iterative method.
%
%  Bus Types:
%       1 -> Slack Bus
%       2 -> PQ Bus
%       3 -> PV Bus
%
%  The program:
%       1. Forms Ybus
%       2. Performs Gauss-Seidel iterations
%       3. Handles PQ and PV buses
%       4. Checks convergence
%       5. Calculates bus voltages and powers
%       6. Calculates line power flows
%       7. Calculates transmission losses
%       8. Checks power balance
%       9. Plots convergence
%
%  Test system: 4-bus example from lecture material
% ================================================================

clc;
clear;
close all;

%% ================================================================
% 1. SYSTEM DATA
% ================================================================

% ------------------------------------------------
% BUS DATA
%
% Columns:
% [Bus Type P Q Vspec Angle Qmin Qmax Vmin Vmax]
%
% Type:
%   1 = Slack
%   2 = PQ
%   3 = PV
% ------------------------------------------------

busData = [
    1   1    0       0       1.05    0     -Inf    Inf    0.90   1.10;
    2   3   -0.45    0       1.00    0     -Inf    Inf    0.90   1.10;
    3   2   -0.51   -0.25     0       0     -Inf    Inf    0.90   1.10;
    4   2   -0.60   -0.30     0       0     -Inf    Inf    0.90   1.10
];

% ------------------------------------------------
% TRANSMISSION LINE DATA
%
% Columns:
% [FromBus ToBus R X B]
%
% B = total line charging susceptance.
% B = 0 means line charging is neglected.
% ------------------------------------------------

lineData = [
    1   2   0.08   0.20   0;
    1   4   0.05   0.10   0;
    2   3   0.04   0.12   0;
    3   4   0.04   0.14   0
];

%% ================================================================
% 2. GAUSS-SEIDEL PARAMETERS
% ================================================================

% Acceleration factor used in lecture example
alpha = 1.2;

% Convergence tolerance
tol = 1e-5;

% Maximum number of iterations
maxIter = 100;

%% ================================================================
% 3. SYSTEM INFORMATION
% ================================================================

NB = size(busData,1);
NL = size(lineData,1);

% Identify bus types
slackBus = find(busData(:,2) == 1);
PQbuses  = find(busData(:,2) == 2);
PVbuses  = find(busData(:,2) == 3);

% Check slack bus
if isempty(slackBus)
    error('No slack bus has been specified.');
end

if length(slackBus) > 1
    error('More than one slack bus has been specified.');
end

slackBus = slackBus(1);

%% ================================================================
% 4. DISPLAY SYSTEM INFORMATION
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('              GAUSS-SEIDEL LOAD FLOW ANALYSIS\n');
fprintf('=================================================================\n');

fprintf('\nSYSTEM INFORMATION\n');
fprintf('-----------------------------------------------------------------\n');

fprintf('Number of buses              : %d\n', NB);
fprintf('Number of transmission lines : %d\n', NL);
fprintf('Slack bus                    : %d\n', slackBus);

fprintf('PQ buses                     : ');
fprintf('%d ', PQbuses);
fprintf('\n');

fprintf('PV buses                     : ');

if isempty(PVbuses)
    fprintf('None');
else
    fprintf('%d ', PVbuses);
end

fprintf('\n');

fprintf('Acceleration factor (alpha) : %.2f\n', alpha);
fprintf('Convergence tolerance       : %.1e\n', tol);
fprintf('Maximum iterations          : %d\n', maxIter);

%% ================================================================
% 5. FORM YBUS
% ================================================================

Ybus = buildYbus(busData, lineData);

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                         Y-BUS MATRIX\n');
fprintf('=================================================================\n');

for i = 1:NB

    for k = 1:NB

        fprintf('%12.6f %+.6fi   ', ...
            real(Ybus(i,k)), imag(Ybus(i,k)));

    end

    fprintf('\n');

end

%% ================================================================
% 6. INITIALIZE BUS VOLTAGES
% ================================================================

V = ones(NB,1);

for i = 1:NB

    busType = busData(i,2);

    if busType == 1

        % Slack bus
        V(i) = busData(i,5) * ...
            exp(1j*deg2rad(busData(i,6)));

    elseif busType == 3

        % PV bus
        V(i) = busData(i,5) * ...
            exp(1j*deg2rad(busData(i,6)));

    else

        % PQ bus - flat start
        V(i) = 1 + 1j*0;

    end

end

%% ================================================================
% 7. GAUSS-SEIDEL ITERATION
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                    GAUSS-SEIDEL ITERATION\n');
fprintf('=================================================================\n');

fprintf('%10s %25s\n', ...
    'Iteration','Maximum Voltage Change');

fprintf('-----------------------------------------------------------------\n');

% Store convergence history
convergenceHistory = zeros(maxIter,1);

% Current bus types
currentType = busData(:,2);

converged = false;

for iter = 1:maxIter

    % Store old voltages
    Vold = V;

    %% ------------------------------------------------------------
    % Update each non-slack bus
    % ------------------------------------------------------------

    for i = 1:NB

        % Slack bus voltage remains fixed
        if i == slackBus
            continue;
        end

        %% ========================================================
        % PQ BUS
        % ========================================================

        if currentType(i) == 2

            P = busData(i,3);
            Q = busData(i,4);

            % Calculate sum(Yik*Vk), k ~= i
            sumYV = 0;

            for k = 1:NB

                if k ~= i
                    sumYV = sumYV + Ybus(i,k)*V(k);
                end

            end

            % Gauss-Seidel voltage equation
            Vcalc = (1/Ybus(i,i)) * ...
                ((P - 1j*Q)/conj(V(i)) - sumYV);

            % Apply acceleration factor
            Vnew = V(i) + alpha*(Vcalc - V(i));

            % Voltage magnitude limits
            Vmag = abs(Vnew);

            Vmin = busData(i,9);
            Vmax = busData(i,10);

            if Vmag < Vmin

                Vnew = Vmin * exp(1j*angle(Vnew));

            elseif Vmag > Vmax

                Vnew = Vmax * exp(1j*angle(Vnew));

            end

            V(i) = Vnew;

        %% ========================================================
        % PV BUS
        % ========================================================

        elseif currentType(i) == 3

            P = busData(i,3);

            % ----------------------------------------------------
            % Calculate reactive power
            % ----------------------------------------------------

            I = 0;

            for k = 1:NB
                I = I + Ybus(i,k)*V(k);
            end

            S = V(i)*conj(I);

            Qcalc = imag(S);

            % Reactive power limits
            Qmin = busData(i,7);
            Qmax = busData(i,8);

            % ----------------------------------------------------
            % Check reactive power limits
            % ----------------------------------------------------

            if Qcalc < Qmin

                Q = Qmin;
                currentType(i) = 2;

            elseif Qcalc > Qmax

                Q = Qmax;
                currentType(i) = 2;

            else

                Q = Qcalc;
                currentType(i) = 3;

            end

            % ----------------------------------------------------
            % Calculate voltage
            % ----------------------------------------------------

            sumYV = 0;

            for k = 1:NB

                if k ~= i
                    sumYV = sumYV + Ybus(i,k)*V(k);
                end

            end

            Vcalc = (1/Ybus(i,i)) * ...
                ((P - 1j*Q)/conj(V(i)) - sumYV);

            % ----------------------------------------------------
            % Maintain specified voltage magnitude for PV bus
            % ----------------------------------------------------

            if currentType(i) == 3

                Vcalc = busData(i,5) * ...
                    exp(1j*angle(Vcalc));

            end

            % Apply acceleration
            Vnew = V(i) + alpha*(Vcalc - V(i));

            % Enforce specified magnitude
            if currentType(i) == 3

                Vnew = busData(i,5) * ...
                    exp(1j*angle(Vnew));

            end

            V(i) = Vnew;

        end

    end

    %% ------------------------------------------------------------
    % Calculate maximum voltage change
    % ------------------------------------------------------------

    deltaV = abs(V - Vold);

    maxDeltaV = max(deltaV);

    convergenceHistory(iter) = maxDeltaV;

    fprintf('%10d %25.10e\n', ...
        iter, maxDeltaV);

    %% ------------------------------------------------------------
    % Convergence check
    % ------------------------------------------------------------

    if maxDeltaV < tol

        converged = true;
        break;

    end

end

% Remove unused convergence-history entries
convergenceHistory = convergenceHistory(1:iter);

fprintf('-----------------------------------------------------------------\n');

%% ================================================================
% 8. CONVERGENCE STATUS
% ================================================================

if converged

    fprintf('\n');
    fprintf('RESULT: LOAD FLOW CONVERGED SUCCESSFULLY\n');
    fprintf('Number of iterations = %d\n', iter);
    fprintf('Final maximum voltage change = %.6e pu\n', maxDeltaV);

else

    fprintf('\n');
    fprintf('RESULT: LOAD FLOW DID NOT CONVERGE\n');
    fprintf('Maximum iterations reached = %d\n', maxIter);

end

%% ================================================================
% 9. CALCULATE BUS POWERS
% ================================================================

Sbus = zeros(NB,1);

for i = 1:NB

    I = 0;

    for k = 1:NB
        I = I + Ybus(i,k)*V(k);
    end

    Sbus(i) = V(i)*conj(I);

end

Pbus = real(Sbus);
Qbus = imag(Sbus);

%% ================================================================
% 10. BUS RESULTS
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                         BUS RESULTS\n');
fprintf('=================================================================\n');

fprintf('%6s %10s %14s %16s %14s %14s\n', ...
    'Bus','Type','|V| (pu)','Angle (deg)','P (pu)','Q (pu)');

fprintf('-----------------------------------------------------------------\n');

for i = 1:NB

    if busData(i,2) == 1

        typeName = 'Slack';

    elseif busData(i,2) == 2

        typeName = 'PQ';

    else

        typeName = 'PV';

    end

    fprintf('%6d %10s %14.6f %16.6f %14.6f %14.6f\n', ...
        busData(i,1), ...
        typeName, ...
        abs(V(i)), ...
        rad2deg(angle(V(i))), ...
        Pbus(i), ...
        Qbus(i));

end

%% ================================================================
% 11. FINAL COMPLEX VOLTAGES
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                    FINAL COMPLEX VOLTAGES\n');
fprintf('=================================================================\n');

fprintf('%6s %18s %18s\n', ...
    'Bus','Real(V)','Imag(V)');

fprintf('-----------------------------------------------------------------\n');

for i = 1:NB

    fprintf('%6d %18.8f %18.8f\n', ...
        busData(i,1), ...
        real(V(i)), ...
        imag(V(i)));

end

%% ================================================================
% 12. SLACK BUS POWER
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                       SLACK BUS POWER\n');
fprintf('=================================================================\n');

fprintf('Slack Bus = %d\n', slackBus);

fprintf('P_slack   = %.8f pu\n', ...
    Pbus(slackBus));

fprintf('Q_slack   = %.8f pu\n', ...
    Qbus(slackBus));

%% ================================================================
% 13. LINE POWER FLOWS
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                         LINE FLOWS\n');
fprintf('=================================================================\n');

fprintf('%6s %6s %14s %14s %14s %14s\n', ...
    'From','To','Pij','Qij','Pji','Qji');

fprintf('-----------------------------------------------------------------\n');

lineResults = zeros(NL,8);

totalLoss = 0;

for n = 1:NL

    from = lineData(n,1);
    to   = lineData(n,2);

    R = lineData(n,3);
    X = lineData(n,4);
    B = lineData(n,5);

    % Series impedance
    Z = R + 1j*X;

    % Series admittance
    y = 1/Z;

    % Bus voltages
    Vi = V(from);
    Vk = V(to);

    % Current from i to k
    Iik = (Vi - Vk)*y + Vi*(1j*B/2);

    % Current from k to i
    Iki = (Vk - Vi)*y + Vk*(1j*B/2);

    % Complex power flows
    Sik = Vi*conj(Iik);
    Ski = Vk*conj(Iki);

    % Line loss
    Sloss = Sik + Ski;

    totalLoss = totalLoss + Sloss;

    % Store results
    lineResults(n,:) = [
        from
        to
        real(Sik)
        imag(Sik)
        real(Ski)
        imag(Ski)
        real(Sloss)
        imag(Sloss)
    ]';

    fprintf('%6d %6d %14.6f %14.6f %14.6f %14.6f\n', ...
        from, to, ...
        real(Sik), imag(Sik), ...
        real(Ski), imag(Ski));

end

%% ================================================================
% 14. TRANSMISSION LOSSES
% ================================================================

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                    TRANSMISSION LOSSES\n');
fprintf('=================================================================\n');

fprintf('%6s %6s %16s %18s\n', ...
    'From','To','P Loss (pu)','Q Loss (pu)');

fprintf('-----------------------------------------------------------------\n');

for n = 1:NL

    fprintf('%6d %6d %16.8f %18.8f\n', ...
        lineResults(n,1), ...
        lineResults(n,2), ...
        lineResults(n,7), ...
        lineResults(n,8));

end

fprintf('-----------------------------------------------------------------\n');

fprintf('Total Real Power Loss     = %.8f pu\n', ...
    real(totalLoss));

fprintf('Total Reactive Power Loss = %.8f pu\n', ...
    imag(totalLoss));

%% ================================================================
% 15. POWER BALANCE
% ================================================================

totalGeneration = sum(Pbus(Pbus > 0));

totalLoad = -sum(Pbus(Pbus < 0));

realLoss = real(totalLoss);

powerBalanceError = ...
    abs((totalGeneration - totalLoad) - realLoss);

fprintf('\n');
fprintf('=================================================================\n');
fprintf('                        POWER BALANCE\n');
fprintf('=================================================================\n');

fprintf('Total Generation = %.8f pu\n', totalGeneration);
fprintf('Total Load       = %.8f pu\n', totalLoad);
fprintf('Total Loss       = %.8f pu\n', realLoss);

fprintf('Generation - Load = %.8f pu\n', ...
    totalGeneration - totalLoad);

fprintf('\nPower balance error = %.3e pu\n', ...
    powerBalanceError);

%% ================================================================
% 16. CONVERGENCE PLOT
% ================================================================

figure;

semilogy(1:iter, convergenceHistory, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5);

grid on;

xlabel('Iteration Number');

ylabel('Maximum Voltage Change |\DeltaV| (pu)');

title('Gauss-Seidel Convergence');

%% ================================================================
% 17. FINAL SUMMARY
% ================================================================

fprintf('\n');

fprintf('                         KEY RESULTS\n');

if converged
    fprintf('Load-flow status          : CONVERGED\n');
else
    fprintf('Load-flow status          : NOT CONVERGED\n');
end

fprintf('Number of iterations      : %d\n', iter);

fprintf('Final voltage deviation   : %.6e pu\n', maxDeltaV);

fprintf('Acceleration factor       : %.2f\n', alpha);

fprintf('Convergence tolerance     : %.1e\n', tol);

fprintf('\nSlack bus generation:\n');

fprintf('    P = %.6f pu\n', ...
    Pbus(slackBus));

fprintf('    Q = %.6f pu\n', ...
    Qbus(slackBus));

fprintf('\nTransmission losses:\n');

fprintf('    P_loss = %.6f pu\n', ...
    real(totalLoss));

fprintf('    Q_loss = %.6f pu\n', ...
    imag(totalLoss));

fprintf('\nPower balance error       : %.3e pu\n', ...
    powerBalanceError);

fprintf('=================================================================\n');
fprintf('                         PROGRAM END\n');
fprintf('=================================================================\n');


%% ================================================================
%                         LOCAL FUNCTION
% ================================================================

function Ybus = buildYbus(busData,lineData)

    % Number of buses
    NB = size(busData,1);

    % Initialize Ybus
    Ybus = zeros(NB,NB);

    % Number of transmission lines
    NL = size(lineData,1);

    % --------------------------------------------------------------
    % Process every transmission line
    % --------------------------------------------------------------

    for n = 1:NL

        from = lineData(n,1);
        to   = lineData(n,2);

        R = lineData(n,3);
        X = lineData(n,4);
        B = lineData(n,5);

        % Series impedance
        Z = R + 1j*X;

        % Series admittance
        y = 1/Z;

        % ----------------------------------------------------------
        % Diagonal elements
        % ----------------------------------------------------------

        Ybus(from,from) = ...
            Ybus(from,from) + y + 1j*B/2;

        Ybus(to,to) = ...
            Ybus(to,to) + y + 1j*B/2;

        % ----------------------------------------------------------
        % Off-diagonal elements
        % ----------------------------------------------------------

        Ybus(from,to) = ...
            Ybus(from,to) - y;

        Ybus(to,from) = ...
            Ybus(to,from) - y;

    end

end