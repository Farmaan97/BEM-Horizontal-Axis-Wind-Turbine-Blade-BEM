%% Input Variable

Uinf  = input('Enter the freestream velocity (m/s): \n');     % freestream veloctiy (m/s)
R     = 6.25;           % Radius of Rotor (m)                 % Input 1
N     = input('Enter the Rotor Speed (rpm): \n');             % Rotor Speed (rpm)
N_rad = 2*pi*N/60;      % Rotor Speed (rad/s)

%% Blade Segment

r = round(xlsread('C:\Users\ajmal\Desktop\Turbine\New Turbine Simulation\BEM Result.xlsx','BEM Result','C2:C18'),2); %segments (m) %Input 2
c = round(xlsread('C:\Users\ajmal\Desktop\Turbine\New Turbine Simulation\BEM Result.xlsx','BEM Result','D2:D18'),2); %chord length of the airfoil at segement (m) %Input 3
theta_t = round(xlsread('C:\Users\ajmal\Desktop\Turbine\New Turbine Simulation\BEM Result.xlsx','BEM Result','E2:E18'),2); %twist angle (deg) %Input 4

%% Initialization Vairable

a     = zeros(length(r),1); % axial induction factor inital assumption
at    = zeros(length(r),1); % angular induction factor
T     = zeros(length(r),1); %Thrust
Q     = zeros(length(r),1); %Torque
P     = zeros(length(r),1); %Power
B     = 3;                  % no of blades
rho   = 1.225;              % density of air @ 25 degCel

%% Setup Details

n = 400;    % no of iterations
e = 0.0001;  % tolerance limit

%% Problem Solution
% For Loop for cycling through each blade element
for j = 1:length(r)
    % For Loop for no of iterations
    for i= 1:n
        %% Disk Velocity and wake Velocity
        Vd = (1-a(j)).*Uinf;
        Vw = (1-2*a(j)).*Uinf;
        
        %% relative velocity and angle of attack determination
        t       = (Uinf/(N_rad.*r(j)))*((1-a(j))./(1+at(j)));
        phi     = atan(t);          % rad
        phi_d   = 180/pi()*phi;     % deg
        alpha_d = phi_d - theta_t(j);  % deg
        alpha   = alpha_d*pi()/180; % rad

        %% Coefficient of Lift and Drag,  and Relative Velocity
        Cl = -0.0091*(alpha_d).^2 + 0.2695*(alpha_d) - 0.1205;
        Cd = -0.0001*(alpha_d).^2 - 0.0051*alpha_d + 0.0527;
        Vr = sqrt(((Uinf.*(1 - a(j))).^2) + ((N_rad*r(j)*(1 + at(j))).^2)); % relative velocity
        
        %% Normal Force Coefficient and Tangential Force Coefficient
        Cn = Cl*cos(phi) + Cd*sin(phi);     % normal force coefficient
        Ct = Cl*sin(phi) - Cd*cos(phi);     % tangential force coefficient
        
        %% Pradntl Tip Loss Factor
        f = (B*(R-r(j)))./(2.*r(j)*sin(phi));
        F = (2/pi)*(acos(exp(-f)));
        
        %% Solidity Constant
        sigmar = (B*c(j))/(2*pi.*r(j));
        Za     = sigmar*Cn + 4*F*(sin(phi)).^2;
        Zt     = N_rad.*r(j)*Cn + 8*F*pi.*r(j)*(sin(phi)).^2;
        
        %% Axial and Tangential Induction Factor re-calculation
        anew   = (sigmar*Cn)/Za;
        atnew  = (Uinf*Ct)/Zt;
        
        %fprintf('a%d = %.4f\n',i,anew)
        %fprintf('at%d = %.4f\n',i,atnew)
        %% Test Tolerance Limit for Loop Break
        if abs(anew-a(j))<e & abs(atnew-at(j))<e
            break;
        end
        %% Reassigning Induction factor for next iteration of loop    
        a(j)  =  anew;
        at(j) = atnew;
    end
    dT = @(r) 2*F*rho*(Uinf^2).*a(j).*(1-a(j)).*2*pi().*r;
    T(j) = integral(dT,0,r(j)); %thrust
    dQ = @(r) 2*rho*Uinf*N_rad*F.*at(j).*(1-a(j)).*2*pi().*(r.^3);
    Q(j) = integral(dQ,0,r(j)); %Torque
    P(j) = N_rad*Q(j);  %Power
end
T(isnan(T))=0;
Q(isnan(Q))=0;
P(isnan(P))=0;
P_Comp = mean(P); %Power Result