function deltaV = compute_overshoot(mean_den, GmU, mD, Ee, Vtarget)

% Computes the necessary overshoot voltage (deltaV) added to Vtarget
% to account for synaptic conductance effects at a distal dendritic site.

% Inputs:
% mean_den  - synaptic density (per unit length)
% GmU       - membrane conductance in nS/um^2
% mD        - dendrite diameter in um
% Ee        - synaptic reversal potential (e.g., 40 mV)
% Vtarget   - desired depolarization at distal site (e.g., 8 mV)

% Output:
% deltaV    - required voltage overshoot (in mV)

    Gmpd = GmU * pi * mD; % passive membrane conductance

    % Try a range of candidate overshoots
    delta_V_candidates = 0:0.01:5; % in mV

    for dV = delta_V_candidates
        Vdist_try  = (Vtarget + dV) * 1 / mean_den;     % adjusted input voltage per synapse
        Idist_try  = Vdist_try * Gmpd;                  % estimated current
        gsynT_try  = Idist_try / Ee;                    % total synaptic conductance
        Vsyn_try   = (Idist_try / (1 / mean_den)) / (Gmpd + gsynT_try / (1 / mean_den)); % actual synaptic voltage

        if abs(Vsyn_try - Vtarget) ==0 % within 0.01 mV of target
            deltaV = dV;
            return
        end
    end

    % If no suitable deltaV found, return NaN
    warning('No overshoot found within tested range.');
    deltaV = NaN;
end
