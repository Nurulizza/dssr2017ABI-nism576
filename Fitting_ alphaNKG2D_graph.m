function [prm_opt,fval,exitflag,output] = Fitting_alphaNKG2D_graph()

%%assessment point 
exp =  [0 0.0064 0.0193 0.0257 0.0578 0.0358 0.0257;
        0 0.1477 0.5906 0.9508 1.0143 0.5007 0.2889;
        0 0.1925 0.6741 1.6435 1.9003 1.1106 0.6805]; %%virtual experiment
time = [0;30;90;180;300;450;600];

tspan = [min(time),max(time)];
nvars = 23;

%Assigning weights
% w1 = 2.0;
% w2 = 1.5;
% w3 = 1.5;

%%Initial conditions(uM) 
MICA_01 = 10;
MICA_02 = 100;
MICA_03 = 1000;
DAP10_0 = 10.97;
D3_0 = 0.0016;
Vav1_0 = 0.19;
SFK_0 = 237.79;
CD45_0 = 39.18;
SHP1_0 = 3.47;

LB = [-8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  -8];
UB = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

objectivesToPlot = [1,2,3]; % Set the objective you wanted to plot
PlotFcn = @(options,state,flag)gaplotpareto(options,state,flag,objectivesToPlot);
options = gaoptimset;
options = gaoptimset(options, 'OutputFcn',@gaoutputfcn);
options = gaoptimset(options,'Display','iter');
options = gaoptimset(options,'Display','final');
options = gaoptimset(options,'PlotFcn',PlotFcn);
[prm_opt,fval,exitflag,output] = gamultiobj(@minimize,nvars,[],[],[],[],LB,UB,[],options);

%%Solve it with ode solver
function fitnessFcn = minimize(prm)
    
D_0 = 10^prm(18);
D1_0 = 10^prm(19);
D2_0 = 10^prm(20);
C1_0 = 10^prm(21);
C2_0 = 10^prm(22);
C3_0 = 10^prm(23); 
C4_0 = 10^prm(24);

X_01 = [MICA_01,DAP10_0,D_0,D1_0,D2_0,D3_0,C1_0,C2_0,C3_0,C4_0,Vav1_0,SFK_0,CD45_0,SHP1_0]; 
sol1 = ode15s(@ODE_HLAE_Activation,tspan,X_01,[],prm(1:17));
y1 = deval(sol1, time);
fitnessFcn(1) = sum((y1(6,:) - exp(1,:)).^2);
fitnessFcn(1) = w1*fitnessFcn(1);

X_02 = [MICA_02,DAP10_0,D_0,D1_0,D2_0,D3_0,C1_0,C2_0,C3_0,C4_0,Vav1_0,SFK_0,CD45_0,SHP1_0];
sol2 = ode15s(@ODE_HLAE_Activation,tspan,X_02,[],prm(1:17));
y2 = deval(sol2, time);
fitnessFcn(2) = sum((y2(6,:) - exp(2,:)).^2);
fitnessFcn(2) = w2*fitnessFcn(2);

X_03 = [MICA_03,DAP10_0,D_0,D1_0,D2_0,D3_0,C1_0,C2_0,C3_0,C4_0,Vav1_0,SFK_0,CD45_0,SHP1_0];
sol3 = ode15s(@ODE_HLAE_Activation,tspan,X_03,[],prm(1:17));
y3 = deval(sol3, time);
fitnessFcn(3) = sum((y3(6,:) - exp(3,:)).^2);
fitnessFcn(3) = w3*fitnessFcn(3);
end

%Function for GA output function (optional)
function [state, options,optchanged] = gaoutputfcn(options,state,flag)
optchanged = false;

switch flag
case 'init'
disp('Starting the algorithm');
case {'iter','interrupt'}
% valueFitt = state.Best;
disp('Iterating ...')
fname=[pwd,'\',num2str(state.Generation),'.mat'];
save(fname,'state')
case 'done'
disp('Performing final task');
fname=[pwd,'\',num2str(state.Generation),'.mat'];
save(fname,'state')
end 
end

% Plot ODE with optimal parameter
% [T,Y] = ode15s(@ODE_HLAE_Activation, tspan, X_0,[],prm_opt);
% filename='GA_HLAE_ACTIVATION.xlsx';%%Save data in Excel file
% A=table(T,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6));
% writetable(A,filename,'sheet',1)
% figure
% plot(time, exp)
% hold on
% plot(T,Y(:,6), 'r--', 'linewidth', 2)%%Column 6 for D3
% legend ('Exp Data','prm_opt');
% xlabel('Time');
% ylabel('dxdt');
end

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
k_f20 = 10^prm(17); 

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
dXdT(6) = k_f8*C3 - k_f9*D3*CD45 + k_r9*C4 - k_f20*SHP1*D3;
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