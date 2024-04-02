function [newU] = applyAAF_TD(oldU,f,meanU,B,model)
% % This function applies the Aerodynamic Admittance Function (AAF) in the time domain to a given wind velocity time series. It is designed to modify the input wind velocity data based on the specified AAF model, particularly useful for wind engineering analyses where aerodynamic effects on structures are considered.
% 
% #### Syntax
% 
% ```matlab
% newU = applyAAF_TD(oldU, f, meanU, B, model)
% ```
% #### Inputs
% 
% - `oldU`: The original wind velocity time series vector.
% - `f`: Frequency vector associated with the wind velocity data.
% - `meanU`: Mean wind velocity.
% - `B`: Characteristic length (e.g., building width or height impacting the flow).
% - `model`: String specifying the AAF model to be applied. Currently, only 'Liepman' is implemented.
% 
% #### Outputs
% 
% - `newU`: The modified wind velocity time series after applying the AAF.
% 
% #### Description of Main Steps
% 
% 1. **Model Selection**: Based on the input `model` string, the function calculates the Aerodynamic Admittance Function squared (\(\phi^2\)) for the 'Liepman' model. If the model is not recognized, it defaults to a unitary admittance function (no modification).
% 2. **Spectrum Modification**: The function calculates the Fourier Transform of the input velocity, adjusts its magnitude spectrum according to the AAF, and then converts it back to the time domain.
% 
% Author: E Cheynet - UiB - last modified 02-04-2024



if strcmpi(model,'Liepman')
    fr = f*B/meanU;
    % Define the square of Aerodynamic Admittance Function (AAF), phi
    phi2 = 1./(1 + 2*pi^2*fr); % Adjust according to the correct formula if needed
else
    warning('models not implemented yet')
    fr = f*B./meanU;
    phi2 = ones(size(fr));
end

% Extend phi to negative frequencies to match the FFT output dimensions
phi_full = [phi2, fliplr(phi2(2:end-1))];

U = fft(oldU(:));
S_u = abs(U(:)).^2; % Velocity spectrum
S_f = S_u(:) .* phi_full(:); % Load spectrum
% Convert the load spectrum back to time domain
f_t = ifft(sqrt(S_f(:)).*exp(1j*angle(U(:))));
newU = real(f_t);

end