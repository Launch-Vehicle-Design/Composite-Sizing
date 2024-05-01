clear;
clc;

%% Current Sizing
% Casing - 20 layers
% Stage 2 tanks - X layers

%% Graphite-Polymer Composite Properties
% https://www.rockwestcomposites.com/14019-d-group
E1 = 165E9; % Pa
E2 = 8.56E9; % Pa
E3 = 8.56E9; % Pa
nu23 = 0.458;
nu13 = 0.326;
nu12 = 0.326;
G23 = 3.2E9; % Pa
G13 = 4.39E9; % Pa
G12 =  4.39E9; % Pa
alpha1 = -0.01800E-6; % 1/C
alpha2 = 24.3E-6; % 1/C
alpha3 = 24.3E-6; % 1/C
beta1 = 146E-6; % /%M
beta2 = 4770E-6; % /%M
beta3 = 4770E-6; % /%M

% Stengths
Xt = 3190E6; % Fiber tensile strength
Xc = 1440E6; % Fiber compression strength
Yt = 82E6; % Transverse tensile strength
Yc = 164E6; % Transverse compression strength
SS = 141E6; % in-plane shear

%% Laminate Geometry
thick = 0.00011684; % m
h = thick*ones(1,8); % Individual layer thickness

thetavector = 1:90;

for thetaloop = 1:90
% thetan = [-thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop]; % Layer fiber angle
thetan = [-thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop -thetaloop thetaloop]; % Layer fiber angle
deltaT = 0; % Layer Temperature Change
deltaM = 0; % Layer Moisture Change %

thetan = deg2rad(thetan);

H = sum(h); % Total laminate thickness
N = length(h); % Number of layers
Layer = zeros(2*N,1);
n = 1;
while n <= N
    Layer(2*n-1) = n;
    Layer(2*n) = n;
    n = n + 1;
end
z = zeros(2*N,1);
n = 2;
z(1) = -H/2;
z(2) = z(1)+h(1);
while n <= N
    z(2*n-1) = z(2*(n-1));
    z(2*n) = z(2*n-1)+h(n);
    n = n + 1;
end

%% Input Errors
if length(h) ~= length(thetan)
    error('Invalid Laminate Geometry')
end
n = 1;
while n <= N
    if h(n) <= 0
        error('Cannot Have Negative Layer Thickness')
    end
    if thetan(n) > 90*pi/180
        error('Fiber Angle Cannot Be Greater Than 90 Degrees')
    elseif thetan(n) < -90*pi/180
        error('Fiber Angle Cannot Be Less Than 90 Degrees')
    end
    n = n + 1;
end

%% Strains and Curvatures
epsx0 = 0;
epsy0 = 0;
gamxy0 = 0;
Kx0 = 0; % 1/m
Ky0 = 0; % 1/m
Kxy0 = 0; % 1/m
epsx = @(x,y,z) epsx0 + z*Kx0;
epsy = @(x,y,z) epsy0 + z*Ky0;
gamxy = @(x,y,z) gamxy0 + z*Kxy0;

eps0 = [epsx0;epsy0;gamxy0;Kx0;Ky0;Kxy0];

%% Laminate Loads

%% Stage 1 MaxQ
% R = 0.3;
% Ixx = pi*R^3*H;
% p = 9.5e6; % Internal Pressure (Pa)
% 
% vsigmax = -45000/(2*pi*R*H) + abs(-270000)*R/Ixx + p*R/(2*H); % (N/m) axial
% vsigmay = p*R/H; % (N/m) hoop
% vtauxy = (5200 + 1250)/(2*pi*R*H); % (N/m) shear

%% Stage 2 Tanks MaxQ
R = 0.3;
Ixx = pi*R^3*H;
p = 2.5e6; % Internal Pressure (Pa)

% Al liner (0.5mm thick) (Only carries pressure)
tAl = 0.0005; % m

vsigmax = -7000/(2*pi*R*(H)) + abs(1000)*R/Ixx + p*R/(2*(H)); % (N/m) axial
vsigmay = p*R/(H); % (N/m) hoop
vtauxy = (700 + 300)/(2*pi*R*(H)); % (N/m) shear

Nx = vsigmax*H; % N/m
Ny = vsigmay*H; % N/m
Nxy = vtauxy*H; % N/m
Mx = 0; % Nm/m
My = 0; % Nm/m
Mxy = 0; % Nm/m

Loads = [Nx;Ny;Nxy;Mx;My;Mxy];

if norm(eps0) > 0
    if norm(Loads) > 0
        error('Can Only Input Strains/Curvatures OR Loads')
    end
end

%% Compliance & Stiffness Matrix
S = [1/E1     -nu12/E1   -nu13/E1   0       0      0
    -nu12/E1  1/E2       -nu23/E2   0       0      0
    -nu13/E1  -nu23/E2   1/E3       0       0      0
    0         0          0          1/G23   0      0
    0         0          0          0       1/G13  0
    0         0          0          0       0      1/G12];

C = inv(S);

Q = [C(1,1)-C(1,3)^2/C(3,3)        C(1,2)-C(1,3)*C(2,3)/C(3,3)    0         0 0 0
    C(1,2)-C(1,3)*C(2,3)/C(3,3)    C(2,2)-C(2,3)^2/C(3,3)         0         0 0 0
    0                              0                              0         0 0 0
    0                              0                              0         0 0 0
    0                              0                              0         0 0 0
    0                              0                              0         0 0 C(6,6)];

%% ABD Matrix (and abd)
n = 1;
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

    NxThat = 0;
    NyThat = 0;
    NxyThat = 0;
    MxThat = 0;
    MyThat = 0;
    MxyThat = 0;
    NxMhat = 0;
    NyMhat = 0;
    NxyMhat = 0;
    MxMhat = 0;
    MyMhat = 0;
    MxyMhat = 0;

while n<=N
    theta = thetan(n);
    Qbar11 = Q(1,1)*cos(theta)^4+2*(Q(1,2)+2*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(2,2)*sin(theta)^4;
    Qbar12 = (Q(1,1)+Q(2,2)-4*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(1,2)*(sin(theta)^4+cos(theta)^4);
    Qbar16 = (Q(1,1)-Q(1,2)-2*Q(6,6))*sin(theta)*cos(theta)^3+(Q(1,2)-Q(2,2)+2*Q(6,6))*sin(theta)^3*cos(theta);
    Qbar22 = Q(1,1)*sin(theta)^4+2*(Q(1,2)+2*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(2,2)*cos(theta)^4;
    Qbar26 = (Q(1,1)-Q(1,2)-2*Q(6,6))*sin(theta)^3*cos(theta)+(Q(1,2)-Q(2,2)+2*Q(6,6))*sin(theta)*cos(theta)^3;
    Qbar66 = (Q(1,1)+Q(2,2)-2*Q(1,2)-2*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(6,6)*(sin(theta)^4+cos(theta)^4);

    Qbar11n(n) = Qbar11;
    Qbar12n(n) = Qbar12;
    Qbar16n(n) = Qbar16;
    Qbar22n(n) = Qbar22;
    Qbar26n(n) = Qbar26;
    Qbar66n(n) = Qbar66;

    alphax(n) = alpha1*cos(theta)^2+alpha2*sin(theta)^2;
    alphay(n) = alpha1*sin(theta)^2+alpha2*cos(theta)^2;
    alphaxy(n) = 2*(alpha1-alpha2)*cos(theta)*sin(theta);

    betax(n) = beta1*cos(theta)^2+beta2*sin(theta)^2;
    betay(n) = beta1*sin(theta)^2+beta2*cos(theta)^2;
    betaxy(n) = 2*(beta1-beta2)*cos(theta)*sin(theta);

    NxThat = NxThat + (Qbar11n(n)*alphax(n)+Qbar12n(n)*alphay(n)+Qbar16n(n)*alphaxy(n))*(z(2*n)-z(2*n-1));
    NyThat = NyThat + (Qbar12n(n)*alphax(n)+Qbar22n(n)*alphay(n)+Qbar26n(n)*alphaxy(n))*(z(2*n)-z(2*n-1));
    NxyThat = NxyThat + (Qbar16n(n)*alphax(n)+Qbar26n(n)*alphay(n)+Qbar66n(n)*alphaxy(n))*(z(2*n)-z(2*n-1));
    MxThat = MxThat + 1/2*(Qbar11n(n)*alphax(n)+Qbar12n(n)*alphay(n)+Qbar16n(n)*alphaxy(n))*(z(2*n)^2-z(2*n-1)^2);
    MyThat = MyThat + 1/2*(Qbar12n(n)*alphax(n)+Qbar22n(n)*alphay(n)+Qbar26n(n)*alphaxy(n))*(z(2*n)^2-z(2*n-1)^2);
    MxyThat = MxyThat + 1/2*(Qbar16n(n)*alphax(n)+Qbar26n(n)*alphay(n)+Qbar66n(n)*alphaxy(n))*(z(2*n)^2-z(2*n-1)^2);

    NxMhat = NxMhat + (Qbar11n(n)*betax(n)+Qbar12n(n)*betay(n)+Qbar16n(n)*betaxy(n))*(z(2*n)-z(2*n-1));
    NyMhat = NyMhat + (Qbar12n(n)*betax(n)+Qbar22n(n)*betay(n)+Qbar26n(n)*betaxy(n))*(z(2*n)-z(2*n-1));
    NxyMhat = NxyMhat + (Qbar16n(n)*betax(n)+Qbar26n(n)*betay(n)+Qbar66n(n)*betaxy(n))*(z(2*n)-z(2*n-1));
    MxMhat = MxMhat + 1/2*(Qbar11n(n)*betax(n)+Qbar12n(n)*betay(n)+Qbar16n(n)*betaxy(n))*(z(2*n)^2-z(2*n-1)^2);
    MyMhat = MyMhat + 1/2*(Qbar12n(n)*betax(n)+Qbar22n(n)*betay(n)+Qbar26n(n)*betaxy(n))*(z(2*n)^2-z(2*n-1)^2);
    MxyMhat = MxyMhat + 1/2*(Qbar16n(n)*betax(n)+Qbar26n(n)*betay(n)+Qbar66n(n)*betaxy(n))*(z(2*n)^2-z(2*n-1)^2);
    
    A(1,1) = A(1,1) + Qbar11n(n)*(z(2*n)-z(2*n-1));
    A(2,2) = A(2,2) + Qbar22n(n)*(z(2*n)-z(2*n-1));
    A(3,3) = A(3,3) + Qbar66n(n)*(z(2*n)-z(2*n-1));
    A(1,2) = A(1,2) + Qbar12n(n)*(z(2*n)-z(2*n-1));
    A(2,1) = A(1,2);
    A(1,3) = A(1,3) + Qbar16n(n)*(z(2*n)-z(2*n-1));
    A(3,1) = A(1,3);
    A(2,3) = A(2,3) + Qbar26n(n)*(z(2*n)-z(2*n-1));
    A(3,2) = A(2,3);

    B(1,1) = B(1,1) + (1/2)*Qbar11n(n)*(z(2*n)^2-z(2*n-1)^2);
    B(2,2) = B(2,2) + (1/2)*Qbar22n(n)*(z(2*n)^2-z(2*n-1)^2);
    B(3,3) = B(3,3) + (1/2)*Qbar66n(n)*(z(2*n)^2-z(2*n-1)^2);
    B(1,2) = B(1,2) + (1/2)*Qbar12n(n)*(z(2*n)^2-z(2*n-1)^2);
    B(2,1) = B(1,2);
    B(1,3) = B(1,3) + (1/2)*Qbar16n(n)*(z(2*n)^2-z(2*n-1)^2);
    B(3,1) = B(1,3);
    B(2,3) = B(2,3) + (1/2)*Qbar26n(n)*(z(2*n)^2-z(2*n-1)^2);
    B(3,2) = B(2,3);

    D(1,1) = D(1,1) + (1/3)*Qbar11n(n)*(z(2*n)^3-z(2*n-1)^3);
    D(2,2) = D(2,2) + (1/3)*Qbar22n(n)*(z(2*n)^3-z(2*n-1)^3);
    D(3,3) = D(3,3) + (1/3)*Qbar66n(n)*(z(2*n)^3-z(2*n-1)^3);
    D(1,2) = D(1,2) + (1/3)*Qbar12n(n)*(z(2*n)^3-z(2*n-1)^3);
    D(2,1) = D(1,2);
    D(1,3) = D(1,3) + (1/3)*Qbar16n(n)*(z(2*n)^3-z(2*n-1)^3);
    D(3,1) = D(1,3);
    D(2,3) = D(2,3) + (1/3)*Qbar26n(n)*(z(2*n)^3-z(2*n-1)^3);
    D(3,2) = D(2,3);

    n = n + 1;
end

% %% Check if Laminate is Balanced and/or Symmetric
% if isequal(h,flip(h)) == true
%     if isequal(thetan,flip(thetan)) == true
%         B = zeros(3,3);
%         % MxThat = 0;
%         % MyThat = 0;
%         % MxyThat = 0;
%         disp('Composite Laminate is Symmetric')
%     else
%         disp('Composite Laminate is NOT Symmetric')
%     end
% end
% 
% balance = thetan*180/pi;
% bal = 0;
% if length(thetan) > 1
%     if rem(length(thetan),2) == 0
%         while length(balance) > 0
%             angle = balance(1);
%             angle_index = find(balance == angle,1);
%             if find(balance == -angle,1) > 0
%                 balance_index = find(balance == -angle,1);
%                 balance(balance_index) = [];
%                 balance(angle_index) = [];
%             else
%                 balance = [];
%                 disp('Composite Laminate is NOT Balanced')
%                 bal = 1;
%             end
%         end
%     end
%     if bal == 0
%         disp('Composite Laminate is Balanced')
%         A(1,3) = 0;
%         A(2,3) = 0;
%         A(3,1) = 0;
%         A(3,2) = 0;
%         % NxyThat = 0;
%     end
% end

%% Display ABD and abd
format shortG
ABD = [A B
       B D];
A;
B;
D;
abd = inv(ABD);
a = [abd(1,1) abd(1,2) abd(1,3)
    abd(2,1) abd(2,2) abd(2,3)
    abd(3,1) abd(3,2) abd(3,3)];
b = [abd(1,4) abd(1,5) abd(1,6)
    abd(2,4) abd(2,5) abd(2,6)
    abd(3,4) abd(3,5) abd(3,6)];
d = [abd(4,4) abd(4,5) abd(4,6)
    abd(5,4) abd(5,5) abd(5,6)
    abd(6,4) abd(6,5) abd(6,6)];

%% Effective Engineering Properties
Exbar = (A(1,1)*A(2,2)-A(1,2)^2)/(A(2,2)*H);
Eybar = (A(1,1)*A(2,2)-A(1,2)^2)/(A(1,1)*H);
Gxybar = A(3,3)/H;
nuxybar = A(1,2)/A(2,2);
nuyxbar = A(1,2)/A(1,1);
alphaxbar = (A(2,2)*NxThat - A(1,2)*NyThat)/(A(1,1)*A(2,2)-A(1,2)^2);
alphaybar = (A(1,1)*NyThat - A(1,2)*NxThat)/(A(1,1)*A(2,2)-A(1,2)^2);
alphaxybar = a(1,3)*NxThat + a(2,3)*NyThat + a(3,3)*NxyThat;

format shortG
table(Exbar,Eybar,Gxybar,nuxybar,nuyxbar);
% disp('Effective Engineering Properties');
% disp('**For Balanced-Symmetric Laminate Only!**');

%% Laminate Loads Output
NxThat;
NyThat;
NxyThat;
MxThat;
MyThat;
MxyThat;

NxT = NxThat*deltaT;
NyT = NyThat*deltaT;
NxyT = NxyThat*deltaT;
MxT = MxThat*deltaT;
MyT = MyThat*deltaT;
MxyT = MxyThat*deltaT;

NxMhat;
NyMhat;
NxyMhat;
MxMhat;
MyMhat;
MxyMhat;

NxM = NxMhat*deltaM;
NyM = NyMhat*deltaM;
NxyM = NxyMhat*deltaM;
MxM = MxMhat*deltaM;
MyM = MyMhat*deltaM;
MxyM = MxyMhat*deltaM;

LoadsT = [NxT; NyT; NxyT; MxT; MyT; MxyT];
LoadsM = [NxM; NyM; NxyM; MxM; MyM; MxyM];

Loads = Loads + LoadsT + LoadsM;

if norm(eps0) > 0
    Loads = ABD*eps0;
    % eps0 = abd*Loads;
elseif norm(Loads) > 0
    eps0 = abd*Loads;
    % Loads = ABD*eps0;
end

epsx0 = eps0(1);
epsy0 = eps0(2);
gamxy0 = eps0(3);
Kx0 = eps0(4); % 1/m
Ky0 = eps0(5); % 1/m
Kxy0 = eps0(6); % 1/m

%% Stress and Strain Calculations
n = 1;
epsxy_top = zeros(3,N);
epsxy_bot = zeros(3,N);
sigxy_top = zeros(3,N);
sigxy_bot = zeros(3,N);

eps12_top = zeros(3,N);
eps12_bot = zeros(3,N);
sig12_top = zeros(3,N);
sig12_bot = zeros(3,N);

epsx = zeros(2*N,1);
epsy = zeros(2*N,1);
gamxy = zeros(2*N,1);
sigx = zeros(2*N,1);
sigy = zeros(2*N,1);
tauxy = zeros(2*N,1);

eps1 = zeros(2*N,1);
eps2 = zeros(2*N,1);
gam12 = zeros(2*N,1);
sig1 = zeros(2*N,1);
sig2 = zeros(2*N,1);
tau12 = zeros(2*N,1);

while n<=N
    theta = thetan(n);
    Qbar11 = Q(1,1)*cos(theta)^4+2*(Q(1,2)+2*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(2,2)*sin(theta)^4;
    Qbar12 = (Q(1,1)+Q(2,2)-4*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(1,2)*(sin(theta)^4+cos(theta)^4);
    Qbar16 = (Q(1,1)-Q(1,2)-2*Q(6,6))*sin(theta)*cos(theta)^3+(Q(1,2)-Q(2,2)+2*Q(6,6))*sin(theta)^3*cos(theta);
    Qbar22 = Q(1,1)*sin(theta)^4+2*(Q(1,2)+2*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(2,2)*cos(theta)^4;
    Qbar26 = (Q(1,1)-Q(1,2)-2*Q(6,6))*sin(theta)^3*cos(theta)+(Q(1,2)-Q(2,2)+2*Q(6,6))*sin(theta)*cos(theta)^3;
    Qbar66 = (Q(1,1)+Q(2,2)-2*Q(1,2)-2*Q(6,6))*sin(theta)^2*cos(theta)^2+Q(6,6)*(sin(theta)^4+cos(theta)^4);

    Qbar11n(n) = Qbar11;
    Qbar12n(n) = Qbar12;
    Qbar16n(n) = Qbar16;
    Qbar22n(n) = Qbar22;
    Qbar26n(n) = Qbar26;
    Qbar66n(n) = Qbar66;


    Qbar = [Qbar11 Qbar12 Qbar16
            Qbar12 Qbar22 Qbar26
            Qbar16 Qbar26 Qbar66];
    
    epsTM = [-alphax(n)*deltaT - betax(n)*deltaM
            -alphay(n)*deltaT - betay(n)*deltaM
            -alphaxy(n)*deltaT - betaxy(n)*deltaM];

    epsxy_top(:,n) = [epsx0+z(2*n-1)*Kx0
                    epsy0+z(2*n-1)*Ky0
                    gamxy0+z(2*n-1)*Kxy0];
    epsxy_bot(:,n) = [epsx0+z(2*n)*Kx0
                    epsy0+z(2*n)*Ky0
                    gamxy0+z(2*n)*Kxy0];

    sigxy_top(:,n) = Qbar*(epsxy_top(:,n) + epsTM);
    sigxy_bot(:,n) = Qbar*(epsxy_bot(:,n) + epsTM);

    epsx(2*n-1) = epsxy_top(1,n);
    epsx(2*n) = epsxy_bot(1,n);
    epsy(2*n-1) = epsxy_top(2,n);
    epsy(2*n) = epsxy_bot(2,n);
    gamxy(2*n-1) = epsxy_top(3,n);
    gamxy(2*n) = epsxy_bot(3,n);

    sigx(2*n-1) = sigxy_top(1,n);
    sigx(2*n) = sigxy_bot(1,n);
    sigy(2*n-1) = sigxy_top(2,n);
    sigy(2*n) = sigxy_bot(2,n);
    tauxy(2*n-1) = sigxy_top(3,n);
    tauxy(2*n) = sigxy_bot(3,n);
    
    T = [cos(theta)^2           sin(theta)^2            2*cos(theta)*sin(theta)
        sin(theta)^2            cos(theta)^2            -2*cos(theta)*sin(theta)
        -cos(theta)*sin(theta)  cos(theta)*sin(theta)   cos(theta)^2-sin(theta)^2];
      
    sig12_top(:,n) = T*sigxy_top(:,n);
    sig12_bot(:,n) = T*sigxy_bot(:,n);

    sig1(2*n-1) = sig12_top(1,n);
    sig1(2*n) = sig12_bot(1,n);
    sig2(2*n-1) = sig12_top(2,n);
    sig2(2*n) = sig12_bot(2,n);
    tau12(2*n-1) = sig12_top(3,n);
    tau12(2*n) = sig12_bot(3,n);
    
    S_red = [1/E1     -nu12/E1   0
        -nu12/E1  1/E2       0
        0       0            1/G12];

    eps12_top(:,n) = S_red*sig12_top(:,n)+[alpha1*deltaT+beta1*deltaM;alpha2*deltaT+beta2*deltaM;0];
    eps12_bot(:,n) =S_red*sig12_bot(:,n)+[alpha1*deltaT+beta1*deltaM;alpha2*deltaT+beta2*deltaM;0];


    eps1(2*n-1) = eps12_top(1,n);
    eps1(2*n) = eps12_bot(1,n);
    eps2(2*n-1) = eps12_top(2,n);
    eps2(2*n) = eps12_bot(2,n);
    gam12(2*n-1) = eps12_top(3,n);
    gam12(2*n) = eps12_bot(3,n);

    n = n + 1;

end

%% Output Tables
thetan = rad2deg(thetan);
for n = 1:N
    thetatable(n) = thetan(n);
end
thetatable = thetatable';
format shortG
table(Layer,z,epsx,epsy,gamxy);
table(Layer,z,sigx,sigy,tauxy);
table(Layer,z,eps1,eps2,gam12);
table(Layer,z,sig1,sig2,tau12);

%% Maximum Stress Failure
for n = 1:N
    if sig1(2*n-1) > 0
        failsig1(2*n-1) = abs(sig1(2*n-1)/Xt);
    else
        failsig1(2*n-1) = abs(sig1(2*n-1)/Xc);
    end
    
    if sig1(2*n) > 0
        failsig1(2*n) = abs(sig1(2*n)/Xt);
    else
        failsig1(2*n) = abs(sig1(2*n)/Xc);
    end
    
    if sig2(2*n-1) > 0
        failsig2(2*n-1) = abs(sig2(2*n-1)/Yt);
    else
        failsig2(2*n-1) = abs(sig2(2*n-1)/Yc);
    end
    
    if sig2(2*n) > 0
        failsig2(2*n) = abs(sig2(2*n)/Yt);
    else
        failsig2(2*n) = abs(sig2(2*n)/Yc);
    end

    failtau12(2*n-1) = abs(tau12(2*n-1)/SS);
    failtau12(2*n) = abs(tau12(2*n)/SS);
end
format shortG
Layer = Layer(1:2:end);
Maxsig1(thetaloop) = max(failsig1)';
Maxsig2(thetaloop) = max(failsig2)';
Maxtau12(thetaloop) = max(failtau12)';
% table(table(thetatable,failsig1,failsig2,failtau12),'VariableNames',{'Maximum Stress Failure'});
% Nfailsig1 = 1./failsig1./1E6;
% Nfailsig2 = 1./failsig2./1E6;
% Nfailtau12 = 1./failtau12./1E6;
% Nfailsig1 = Nfailsig1;
% Nfailsig2 = Nfailsig2;
% Nfailtau12 = Nfailtau12;
end
% table(table(Layer,thetatable,Nfailsig1,Nfailsig2,Nfailtau12),'VariableNames',{'Maximum Stress Failure - Load to Failure (MN)'});

figure(1)
hold on
semilogy(thetavector,Maxsig1)
semilogy(thetavector,Maxsig2)
semilogy(thetavector,Maxtau12)
yline(1)
legend('sigma1','sigma2','tau12')
ylabel('Failure Criteria')
xlabel('Fiber Angle (Deg)')
% title('Maximum Stress Failure')
hold off

disp(H*1000)

% %% Tsai-Wu Failure
% F1 = (1/Xt - 1/Xc);
% F2 = (1/Yt - 1/Yc);
% F11 = 1/(Xt*Xc);
% F22 = 1/(Yt*Yc);
% F66 = (1/SS)^2;
% 
% Nload = 1;
% tsaiwutens = zeros(N,1);
% while max(tsaiwutens) < 1
% for n = 1:N
%     F11sig1tens(n,1) = F11*(Nload*sig1(2*n-1))^2;
%     F22sig2tens(n,1) = F22*(Nload*sig2(2*n-1))^2;
%     F66tau12tens(n,1) = F66*(Nload*tau12(2*n-1))^2;
%     F1sig1tens(n,1) = F1*Nload*sig1(2*n-1);
%     F2sig2tens(n,1) = F2*Nload*sig2(2*n-1);
%     sqrtF11F22tens(n,1) = sqrt(F11*F22)*Nload*sig1(2*n-1)*Nload*sig2(2*n-1);
%     tsaiwutens(n,1) = F11sig1tens(n,1) + F22sig2tens(n,1) + F66tau12tens(n,1) + F1sig1tens(n,1) + F2sig2tens(n,1) + sqrtF11F22tens(n,1);
% end
% Nload = Nload+1;
% end
% Ntens = Nload*ones(N,1)./1E6;
% 
% Nload = 1;
% tsaiwucomp = zeros(N,1);
% while max(tsaiwucomp) < 1
% for n = 1:N
%     F11sig1comp(n,1) = F11*(Nload*sig1(2*n-1))^2;
%     F22sig2comp(n,1) = F22*(Nload*sig2(2*n-1))^2;
%     F66tau12comp(n,1) = F66*(Nload*tau12(2*n-1))^2;
%     F1sig1comp(n,1) = F1*Nload*sig1(2*n-1);
%     F2sig2comp(n,1) = F2*Nload*sig2(2*n-1);
%     sqrtF11F22comp(n,1) = sqrt(F11*F22)*Nload*sig1(2*n-1)*Nload*sig2(2*n-1);
%     tsaiwucomp(n,1) = F11sig1comp(n,1) + F22sig2comp(n,1) + F66tau12comp(n,1) + F1sig1comp(n,1) + F2sig2comp(n,1) + sqrtF11F22comp(n,1);
% end
% Nload = Nload-1;
% end
% Ncomp = Nload*ones(N,1)./1E6;
% format shortG
% tsaiwufauluretens = table(Layer,thetatable,Ntens,F11sig1tens,F22sig2tens,F66tau12tens,F1sig1tens,F2sig2tens,sqrtF11F22tens,tsaiwutens);
% tsaiwufaulurecomp = table(Layer,thetatable,Ncomp,F11sig1comp,F22sig2comp,F66tau12comp,F1sig1comp,F2sig2comp,sqrtF11F22comp,tsaiwucomp);