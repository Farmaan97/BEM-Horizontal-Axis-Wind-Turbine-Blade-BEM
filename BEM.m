%% Input Variable

Uinf = 11; %freestream veloctiy (m/s)
R = 6.25; %Radius of Rotor (m) %Input 1
gamma = 90; %Rotor Speed (rpm)
gamma_rad = 2*pi()*gamma/60; %Rotor Speed (rad/s)

%% Blade Segment

r = 1.56; %segments (m) %Input 2
c = 0.63; %chord length of the airfoil at segement (m) %Input 3
theta_t = 14.22; %twist angle (deg) %Input 4

%% Initialization Vairable

a = 0; %axial induction factor inital assumption
at = 0; %angular induction factor
B = 3; %no of blades
rho = 1.225; %density of air @ 25 degCel

%% Setup Details

n = 400 %no of iterations
e = 1e-03 %tolerance limit

%% Problem Solution

for i= 1:n
Vd = (1-a)*Uinf
Vw = (1-2*a)*Uinf

t = (Uinf/(gamma_rad*r))*((1-a)/(1+at));
phi = atan(t); %rad
phi_d = 180/pi()*phi; %deg
alpha_d = phi_d - theta_t %deg
alpha = alpha_d*pi()/180 %rad

Cl = -0.0091*(alpha_d)^2 + 0.2695*(alpha_d) - 0.1205
Cd = -0.0001*(alpha_d)^2 - 0.0051*alpha_d + 0.0527
Vr = sqrt(((Uinf*(1-a))^2)+((gamma_rad*r*(1+at))^2)); %relative velocity

Cn = Cl*cos(phi)+Cd*sin(phi); %normal force coefficient
Ct = Cl*sin(phi)-Cd*cos(phi); %tangential force coefficient

f = (B*(R-r))/(2*r*sin(phi));
F = (2/pi())*(acos(exp(-f)))

sigmar = (B*c)/(2*pi()*r);
Za = sigmar*Cn+4*F*(sin(phi))^2;
Zt = gamma_rad*r*Cn+8*F*pi()*r*(sin(phi))^2;
anew = (sigmar*Cn)/Za;
atnew = (Uinf*Ct)/Zt;
fprintf('a%d = %.4f\n',i,anew)
fprintf('at%d = %.4f\n',i,atnew)
    if abs(anew-a)<e && abs(atnew-at)<e
        break;
    end
a = anew;
at = atnew;
end

fprintf('\n')
fprintf('F = %.4f \n',F)
a = anew;
at = atnew;
fprintf('converged a = %.4f \n', a)
fprintf('converged at = %.4f \n', at)

%resulting values
dT = @(r) 2*F*rho*(Uinf^2)*a*(1-a)*2*pi().*r;
T = integral(dT,0,r) %thrust
dQ = @(r) 2*F*at*(1-a)*rho*Uinf*gamma_rad*(r.^3)*2*pi();
Q = integral(dQ,0,r)
P = gamma_rad*Q

