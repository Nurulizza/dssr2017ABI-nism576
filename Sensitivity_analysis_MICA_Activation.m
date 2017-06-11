function [t,y,prm] = Sensitivity_analysis_MICA_Activation()

%%Initial conditions(uM) 
MICA_0 = 10;
DAP10_0 = 10.97;
D_0 = 10^-4.1816;
D1_0 = 10^-3.3614;
D2_0 = 10^-4.9460;
D3_0 = 0.0016;
C1_0 = 10^-1.9727;
C2_0 = 10^-2.4014;
C3_0 = 10^-3.5698; 
C4_0 = 10^-4.8143;
Vav1_0 = 0.19;
SFK_0 = 237.79;
CD45_0 = 39.18;
X_0 = [MICA_0,DAP10_0,D_0,D1_0,D2_0,D3_0,C1_0,C2_0,C3_0,C4_0,Vav1_0,SFK_0,CD45_0];
    
%%Define parameters
prm(1) = 10^-3.2984;
prm(2) = 10^-1.7843;
prm(3) = 10^-1.8213;
prm(4) = 10^-3.2350;
prm(5) = 10^-3.3972;
prm(6) = 10^-0.0537;
prm(7) = 10^-2.4588;
prm(8) = 10^-4.1380;
prm(9) = 10^-2.9776;
prm(10) = 10^-1.7006;
prm(11) = 10^-2.2110;
prm(12) = 10^-1.9272;
prm(13) = 10^-3.3338;
prm(14) = 10^-1.9775;
prm(15) = 10^-4.2713;
prm(16) = 10^-2.9600;

tspan = [0 300];
[T,Y] = ode15s(@ODE_MICA,tspan,X_0,[],prm);

prm(1) = 10^-3.2984;
prm(2) = 10^-1.7843;
prm(3) = 10^-1.8213;
prm(4) = 10^-3.2350;
prm(5) = 10^-3.3972;
prm(6) = 10^-0.0537;
prm(7) = 10^-2.4588;
prm(8) = 10^-4.1380;
prm(9) = 10^-2.9776;
prm(10) = 10^-1.7006;
prm(11) = 10^-2.2110;
prm(12) = 10^-1.9272;
prm(13) = 10^-3.3338;
prm(14) = 10^-1.9775;
prm(15) = 10^-4.2713;
prm(16) = 10^-2.9600;

tspan = [0 300];
[t,y] = ode15s(@(t,y)ODE_MICA(t,y,prm),tspan+T(end),Y(end,:));
T = [T;t(2:end)];
Y = [Y;y(2:end,:)];
      
figure();
plot(T,Y(:,6));
xlabel('time');
ylabel('y');
end

function dXdT = ODE_MICA(~,X,prm)

%%Define parameters
k_r1 = prm(1);
k_f1 = prm(2);
k_f2 = prm(3);
k_r2 = prm(4);
k_f5 = prm(5);
k_f3 = prm(6);
k_f4 = prm(7);
k_r4 = prm(8);
k_f6 = prm(9);
k_r6 = prm(10);
k_f7 = prm(11);
k_r7 = prm(12);
k_f10 = prm(13);
k_f8 = prm(14);
k_f9 = prm(15);
k_r9 = prm(16);

%%Variables
MICA = X(1);
DAP10 = X(2);
D = X(3);
D1 = X(4);
D2 = X(5);
D3 = X(6);
C1 = X(7);
C2 = X(8);
C3 = X(9); 
C4 = X(10);
Vav1 = X(11);
SFK = X(12);
CD45 = X(13);

%%The differential equations
dXdT = zeros(13,1);
% dXdT(MICA) = -J1;
dXdT(1) = k_r1*D - k_f1*MICA*DAP10;
% dXdT(DAP10) = -J1;
dXdT(2) = k_r1*D - k_f1*MICA*DAP10;
% dXdT(D) = J1 - J2 + J5;
dXdT(3) = -k_r1*D + k_f1*MICA*DAP10 - k_f2*D*SFK + k_r2*C1 + k_f5*C2;
% dXdT(D1) = J3 - J4 - J6;
dXdT(4) = k_f3*C1 - k_f4*D1*CD45 + k_r4*C2 - k_f6*Vav1*D1 + k_r6*D2; 
% dXdT(D2) = J6 - J7 + J10;
dXdT(5) = k_f6*Vav1*D1 - k_r6*D2 - k_f7*D2*SFK + k_r7*C3 + k_f10*C4;
% dXdT(D3) = J8 - J9;
dXdT(6) = k_f8*C3 - k_f9*D3*CD45 + k_r9*C4;
% dXdT(C1) = J2 - J3;
dXdT(7) = k_f2*D*SFK - k_r2*C1 - k_f3*C1;
% dXdT(C2) = J4 - J5;
dXdT(8) = k_f4*D1*CD45 - k_r4*C2 - k_f5*C2;
% dXdT(C3) = J7 - J8;
dXdT(9) = k_f7*D2*SFK - k_r7*C3 - k_f8*C3;
% dXdT(C4) = J9 -J10;
dXdT(10) = k_f9*D3*CD45 - k_r9*C4 - k_f10*C4;
% dXdT(Vav1) = -J6;
dXdT(11) = -k_f6*Vav1*D1 + k_r6*D2;
% dXdT(SFK) = -J2 + J3 -J7 +J8;
dXdT(12) = -k_f2*D*SFK + k_r2*C1 + k_f3*C1 - k_f7*D2*SFK + k_r7*C3 + k_f8*C3;
% dXdT(CD45) = -J4 + J5 - J9 + J10;
dXdT(13) = -k_f4*D1*CD45 + k_r4*C2 + k_f5*C2 - k_f9*D3*CD45 + k_r9*C4 + k_f10*C4;

end