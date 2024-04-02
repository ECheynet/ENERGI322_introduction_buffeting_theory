function [rx,drx,ddrx] = getResponse1D_ODE45(t,u_old,meanU,A,Cd,rho,method,modelAAF,M,K,C)
% This function calculates the dynamic response of a Single Degree of Freedom (SDOF) system under wind loading, using the `ode45` solver for differential equations in MATLAB. It supports both linear and non-linear methods for aerodynamic force calculation and applies the Aerodynamic Admittance Function (AAF) to the wind velocity input data.
% #### Syntax
% 
% ```matlab
% [rx, drx, ddrx] = getResponse1D_ODE45(t, u_old, meanU, A, Cd, rho, method, modelAAF, M, K, C)
% ```
% #### Inputs
% 
% - `t`: Time vector (s).
% - `u_old`: Original wind velocity vector (m/s) before applying AAF.
% - `meanU`: Mean wind velocity (m/s).
% - `A`: Cross-sectional area of the structure (m²).
% - `Cd`: Drag coefficient.
% - `rho`: Air density (kg/m³).
% - `method`: String specifying the method for aerodynamic force calculation ('linear' or 'non-linear').
% - `modelAAF`: Model name for the Aerodynamic Admittance Function (AAF).
% - `M`: Mass of the SDOF system (kg).
% - `K`: Stiffness of the SDOF system (N/m).
% - `C`: Damping coefficient of the SDOF system (N·s/m).
% 
% #### Outputs
% 
% - `rx`: Displacement response of the SDOF system (m).
% - `drx`: Velocity response of the SDOF system (m/s).
% - `ddrx`: Acceleration response of the SDOF system (m/s²), calculated numerically.
% 
% #### Description of Main Steps
% 
% 1. **AAF Application**: The function applies the specified AAF model to the original wind velocity data to account for the frequency-dependent aerodynamic effects.
% 2. **Preallocation Check**: Adjusts the damping coefficient for the linear method by accounting for aerodynamic damping.
% 3. **Solving Differential Equation**: Utilizes `ode45` to solve the mass-damping-stiffness system differential equations, based on the interpolated wind velocity and the specified method for aerodynamic force calculation.
% 4. **Interpolation**: Interpolates the `ode45` solution back to the original time vector for displacement and velocity.
% 5. **Acceleration Calculation**: Calculates acceleration by numerically differentiating the velocity response.
% 
% #### Helper Functions
% 
% - `massDampingStiffnessSystem`: Defines the differential equations for the mass-damping-stiffness system given the current state, system properties, and interpolated wind velocity.
% - `getFaero`: Calculates the aerodynamic force based on the current relative wind velocity, system velocity, and specified method (linear or non-linear).
% 
% Author: E Cheynet - UiB - last modified 02-04-2024


fs = 1/median(diff(t));
N=numel(t);
% Frequency vector
f0 = 1/t(end);
f = 0:f0:(fs/2+f0);
B = sqrt(A);
[u] = applyAAF_TD(u_old,f,meanU,B,modelAAF);

% Preallocation
if strcmpi(method,'linear')
    C = C + rho*A*Cd*meanU;
end

% Initial conditions
tSpan = [t(1), t(end)]; % Define the time span for the simulation
Y0 = [0; 0]; % Initial conditions: [Initial displacement, Initial velocity]

% Define options for ode45, if necessary
options = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Assuming 't' is your time vector and 'u' is your wind velocity vector
u_interp = @(t_query) interp1(t, u, t_query, 'linear', 'extrap');
% Solve the system
[T,Y] = ode45(@(t,Y) massDampingStiffnessSystem(t,Y,M,K,C,u_interp,meanU,A,Cd,rho,method), tSpan, Y0, options);

D_interp = interp1(T, Y(:,1), t, 'linear');
V_interp = interp1(T, Y(:,2), t, 'linear');

% Extract displacement and velocity from the solution
rx = D_interp;
drx = V_interp;

% If you need acceleration, differentiate velocity or recalculate it from forces
ddrx = diff(drx)./diff(T); % This is a simple numerical derivative, consider using a more accurate method

function dYdt = massDampingStiffnessSystem(t, Y, M, K, C, u_interp, meanU, A, Cd, rho, method)
    D = Y(1); % Displacement
    V = Y(2); % Velocity
    u_current = u_interp(t); % Current wind velocity from interpolated function
    F_aero = getFaero(u_current, meanU, V, A, Cd, rho, method); % Compute aerodynamic force
    
    % Calculate acceleration
    A = M\(F_aero - C*V - K*D); 
    
    % Derivatives of displacement and velocity
    dYdt = [V; A];
end


function [Faero]= getFaero(u_prime,meanU,drx,A,Cd,rho,method)


COEFF = 1/2*rho*A*Cd;
if strcmpi(method,'linear')
    % Vrel2 =(meanU.^2 + 2*meanU.*u_prime - 2*meanU.*drx);
    Vrel2 =(meanU.^2 + 2*meanU.*u_prime);
elseif strcmpi(method,'non-linear')
    Vrel2 = (meanU + u_prime - drx).^2;
else
    error('unknown method');
end
Faero = COEFF.*Vrel2;



end


end