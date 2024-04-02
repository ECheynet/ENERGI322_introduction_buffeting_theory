function [S_response,Hmeca] = getResponse_FD(meanU,Su,rho,A,Cd,f,M,K,C,modelAAF)
% This function computes the frequency-domain response of a Single Degree 
% of Freedom (SDOF) system subjected to wind loading. 
% It incorporates the Aerodynamic Admittance Function (AAF) 
% to adjust the input wind load Power Spectral Density (PSD) and 
% calculates the mechanical admittance to derive the system's response spectrum.
% 
% #### Syntax
% 
% ```matlab
% [S_response, Hmeca] = getResponse_FD(meanU, Su, rho, A, Cd, f, M, K, C, modelAAF)
% ```
% 
% #### Inputs
% 
% - `meanU`: Mean wind velocity.
% - `Su`: Power Spectral Density (PSD) of the along-wind velocity component.
% - `rho`: Air density.
% - `A`: Cross-sectional area of the structure subject to the wind load.
% - `Cd`: Drag coefficient.
% - `f`: Frequency vector for the analysis.
% - `M`: Mass of the SDOF system.
% - `K`: Stiffness of the SDOF system.
% - `C`: Damping coefficient of the SDOF system.
% - `modelAAF`: String specifying the AAF model to be applied. Currently, 'Liepman' is implemented.
% 
% #### Outputs
% 
% - `S_response`: PSD of the SDOF system displacement response.
% - `Hmeca`: Mechanical admittance of the system.
% 
% Author: E Cheynet - UiB - last modified 02-04-2024


if strcmpi(modelAAF,'Liepman')
    B = sqrt(A);
    fr = f*B/meanU;
    % Define the square of Aerodynamic Admittance Function (AAF), phi
    phi2 = 1./(1 + 2*pi^2*fr); % Adjust according to the correct formula if needed
else
    warning('modelAAF not implemented yet'),
    B = sqrt(A);
    fr = f*B./meanU;
    phi2 = ones(size(fr));
end



SFb = (rho.*A.*Cd.*meanU)^2.*Su.*phi2;
omega_i = 1i.*2*pi*f;
Ctot = C + rho.*A*Cd.*meanU;
Hmeca = 1./(K + Ctot.*omega_i+M.*(omega_i).^2);
S_response = SFb.*abs(Hmeca.^2);


end