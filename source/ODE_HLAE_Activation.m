function dXdT = ODE_HLAE_Activation(~,X,prm)
%%Define parameters
k_r1 = 10^prm(1);
k_f1 = 10^prm(2);
k_f2 = 10^prm(3);
k_r2 = 10^prm(4);
k_f5 = 10^prm(5);
k_f3 = 10^prm(6);
k_f4 = 10^prm(7);
k_r4 = 10^prm(8);
k_f6 = 10^prm(9);
k_r6 = 10^prm(10);
k_f7 = 10^prm(11);
k_r7 = 10^prm(12);
k_f10 = 10^prm(13);
k_f8 = 10^prm(14);
k_f9 = 10^prm(15);
k_r9 = 10^prm(16);

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