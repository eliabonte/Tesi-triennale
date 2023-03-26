%% parametri modello quarter car
clc
clear
close

%dati modello quarter car
mu = 70; %%unsprung mass
ms = 450; %%sprung mass
k = 27000; %%suspension stiffness
kt = 300000; %%tire stiffness
c_min = 800; %%minimum damping
k0 = 38000; %%saturation slope
f_max = 3000; %%maximum saturation level

%definisco il sistema dinamico lineare, caso con uMR = 0
a_32 = -k/ms;
a_34 = -c_min/ms;
a_41 = kt/mu;
a_42 = -(k/ms)-((kt+k)/mu);
a_44 = -(c_min*(ms+mu))/(mu*ms);
A = [0 0 1 0; 0 0 0 1; 0 a_32  0 a_34; a_41 a_42 0 a_44];
B = [0; 0; 0; -kt/mu];
C = [0 -k/ms 0 -c_min/ms];
D = 0;

%valori per generare la funzione di trasferimento Gr
er = 0.7;
v = 25; % m/s
lc = 20; %risoluzione spaziale
wr = 2*pi*v/lc;
% gaussian noise -> Gr -> zr

%Gaussian noise
lung_strada = 2500;
t_tot = lung_strada/v; %considero percorso lungo 2.5km
dt = 0.001; %tempo di campionamento
t = 0:dt:t_tot;
n = randn(length(t),1);
N = [t',n];

uMR = 0;
out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
%road surface
output_zr = getElement(out.yout, 'zr');
zr = output_zr.Values.Data(:);

strada = 0:(dt/t_tot)*lung_strada:lung_strada;
figure(1)
plot(strada,zr,'LineWidth',0.05)
title('Profilo stradale')
xlabel('[m]')
ylabel('[m]')
xlim([100 200])
ylim([-0.025 0.02])


out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
output_zs2punti_l = getElement(out.yout, 'zs_2punti_l');
zs_2punti_l = output_zs2punti_l.Values.Data(:);
figure;
plot(t,zs_2punti_nl,'LineWidth',0.05)
hold on
plot(t,zs_2punti_l,'LineWidth',0.05)
title('output zs_2punti')
xlabel('[m]')
ylabel('[m/s^2]')


%% Performance migliore -> caso passivo (uMR = 0)

n = randn(length(t),1);
N = [t',n];
uMR = 0;
out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
zs_2punti_nl = output_zs2punti_nl.Values.Data(:);

Jmin = rms(zs_2punti_nl)


%% indice di performance J in funzione di valori di uMR = [0, fmax]

uMR_values = 0:50:f_max;

J = zeros(1,length(uMR_values));

for i = 1 : length(uMR_values)
    
    uMR = uMR_values(i);
    
    out = sim("QuarterCarModel_nonLin_vs_Lin.slx"); %simulazione del modello quarter car

    %estraggo l'output zs_2punti
    output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
    zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
    %calcolo indice di performance
    J(i) = rms(zs_2punti_nl);
end

figure(2)
plot(uMR_values,J)
title('Indice J per valori differenti di u_M_R')
ylabel('J [m/s^2]')
xlabel('u_M_R [N]');
grid on;

%% Plot grafico forza sospensione - velocità di stroke
clc
n = randn(length(t),1);
N = [t',n];

uMR_val_princ = [0 f_max/4 f_max/2 f_max*(3/4) f_max];

fd = zeros(length(t),5);

for i = 1 : length(uMR_val_princ)
    
    uMR = uMR_val_princ(i);
    
    out = sim("QuarterCarModel_nonLin_vs_Lin.slx"); %simulazione del modello quarter car
    
    %estraggo l'output (zs_punto-zu_punto)
    outputZpunto = getElement(out.yout, 'zTilde_punto');
    Zpunto_notsort = outputZpunto.Values.Data(:);
    Zpunto = sort(Zpunto_notsort);
    %estraggo l'output fd
    outputfd = getElement(out.yout, 'fd');
    fd_notsort(:,i) = outputfd.Values.Data(:);
    fd(:,i) = sort(fd_notsort(:,i));
    
    hold on
    figure(3)
    if i == 1
        plot(Zpunto,fd(:,i), '--r','LineWidth',2)
    elseif i == 5
        plot(Zpunto,fd(:,i), '--k','LineWidth',2)
    else
        plot(Zpunto,fd(:,i))
    end    
    title('Relazione forza sospensione - velocità')
    xlabel('zs(punto) - zu(punto) [m/s]')
    ylabel('f_d [N]')
    xlim([-1.5 1.5])
    ylim([-4000 4000])
    
end

legend('u_M_R = 0', 'u_M_R = fmax/4', 'u_M_R = fmax/2', 'u_M_R = (fmax*3/4)', 'u_M_R = fmax')
grid on
hold off

%% FUNZIONE DI TRASFERIMENTO MODELLO LINEARE E NON LINEARE
clc

fs = 1/dt;

uMR = 0;
out = sim("QuarterCarModel_nonLin_vs_Lin.slx"); %simulazione del modello quarter car

%%Funzione di trasferimento sistema lineare
figure(4)
sys = ss(A,B,C,D);
G = tf(sys);
bode(G);

%%Confronto Funzioni di trasferimento sistema non lineare vs lineare con
%%funzione tfestimate()

%estraggo l'output zs_2punti del sistema non lineare
zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
y = zs2punti_nl.Values.Data(:);
%estraggo l'output zs_2punti del sistema lineare
zs2punti_l = getElement(out.yout, 'zs_2punti_l');
yl = zs2punti_l.Values.Data(:);
%estraggo l'input zr
zr = getElement(out.yout, 'zr');
x = zr.Values.Data(:);

[H_nl,f_nl] = tfestimate(x,y,[],[],[],fs); %stima della funzione di trasferimento sist NL
[H_l,f_l] = tfestimate(x,yl,[],[],[],fs); %stima della funzione di trasferimento sist L

figure(5)
subplot(2,1,1)
semilogx(2*pi*f_nl,mag2db(abs(H_nl)))
hold on
semilogx(2*pi*f_l,mag2db(abs(H_l)))
title('Diagramma di Bode FdT - Modello Non Lineare vs Modello Lineare')
legend('Modello Non Lineare', 'Modello Lineare')
ylabel('|H|[dB]')
xlabel('Frequenza [rad/s]');
xlim([1 10^3])
ylim([0 60])
grid
hold off
subplot(2,1,2)
semilogx(2*pi*f_nl,rad2deg(angle(H_nl)))
hold on
semilogx(2*pi*f_l,rad2deg(angle(H_l)))
title('Diagramma di Bode FdT - Modello Non Lineare vs Modello Lineare')
legend('Non Linear Model', 'Linear Model')
ylabel('Fase [deg]')
xlabel('Frequenza [rad/s]');
xlim([1 10^3])
ylim([-200 200])
yticks([-180 -90 0 90 180])
grid
hold off

%% Ottimizzazione Bayesiana per il valore uMR
%%clc
n = randn(length(t),1);
N = [t',n];
assignin('base','N',N);

I_max = 20;
uMR_opt = optimizableVariable('uMR',[0,f_max]);
bo_uMR = bayesopt(@myfun,uMR_opt,...
                  'MaxObjectiveEvaluations', I_max, ...
                  'AcquisitionFunctionName', 'expected-improvement', ...
                  'explorationratio', 0.5, ...
                  'NumSeedPoints', 4,...
                  'PlotFcn',[]);

%% Ottimizzazione Bayesiana per legge di controllo LINEARE
%vettore variabili di stato x = [xs delta xs_punto delta_punto]
clc

k_molt = 1000;

simulazioni = length(k_molt);

iter = zeros(simulazioni,1);
Jopt = zeros(simulazioni,1);
k1_t = zeros(simulazioni,1); k2_t = zeros(simulazioni,1); k3_t = zeros(simulazioni,1); k4_t = zeros(simulazioni,1); 
best_J = 2;
bound_value = zeros(simulazioni,1);

%n = randn(length(t),1);
%N = [t',n];
uMR = 0;
out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
Jmin_k = rms(zs_2punti_nl);

for i = 1:simulazioni
 
    f_max = 3000;
    uMR_values = 0:150:f_max;

    k_range = calcolorange_k(uMR_values,f_max,k_molt(i));

    k1 = optimizableVariable('k1',[-k_range(1),k_range(1)]);
    k2 = optimizableVariable('k2',[-k_range(2),k_range(2)]);
    k3 = optimizableVariable('k3',[-k_range(3),k_range(3)]);
    k4 = optimizableVariable('k4',[-k_range(4),k_range(4)]);

    I_max = 150;
    bo_lin = bayesopt(@fun_lin,[k1,k2,k3,k4],...
                     'MaxObjectiveEvaluations', I_max, ...
                     'AcquisitionFunctionName', 'expected-improvement', ...
                     'explorationratio', 0.5, ...
                     'NumSeedPoints', 4,...
                     'PlotFcn',[]);
                 
    %iterazioni = 1:bo_lin.NumObjectiveEvaluations;
    %J_visited = bo_lin.ObjectiveTrace;
    %J_min_visited = bo_lin.ObjectiveMinimumTrace;
    %J_min_est = bo_lin.EstimatedObjectiveMinimumTrace;

    %plot_fcnJ(iterazioni,J_visited, J_min_visited, J_min_est, 7);
    
    %{
    uMR = 0;
    out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
    output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
    zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
    %}
 
    
    iter(i) = bo_lin.NumObjectiveEvaluations;
    Jopt(i) = bo_lin.MinEstimatedObjective;
    k1_t(i) = bo_lin.XAtMinEstimatedObjective.k1;
    k2_t(i) = bo_lin.XAtMinEstimatedObjective.k2;
    k3_t(i) = bo_lin.XAtMinEstimatedObjective.k3;
    k4_t(i) = bo_lin.XAtMinEstimatedObjective.k4; 
    bound_value(i) = max(k_range);
    
    %salvo il set di parametri migliore
    if bo_lin.MinEstimatedObjective < best_J
        best_J = bo_lin.MinEstimatedObjective;
        best_Ktilde = [k1_t(i) k2_t(i) k3_t(i) k4_t(i)];
        best_x = n;
    end
end

%k_tilde = [bo_lin.XAtMinEstimatedObjective.k1 bo_lin.XAtMinEstimatedObjective.k2 bo_lin.XAtMinEstimatedObjective.k3 bo_lin.XAtMinEstimatedObjective.k4];

Jmin_k_t = Jmin_k + zeros(simulazioni,1);
table_results = table(iter, k1_t, k2_t, k3_t, k4_t, bound_value, Jopt, Jmin_k_t);


%% grafico della FdT sospensione controllata da legge lineare(set di parametri migliori) vs FdT sospensione passiva ottima

fs = 1/dt;
N = [t',n];
k_tilde = best_Ktilde;
out = sim("controllo_leggelineare.slx");
zr = getElement(out.yout, 'zr');
x = zr.Values.Data(:);
xs2punti = getElement(out.yout, 'xs_2punti');
y = xs2punti.Values.Data(:);
xs2punti_lineare = getElement(out.yout, 'xs_2punti_lineare');
y_lineare = xs2punti_lineare.Values.Data(:);
[H_leggelin,f_leggelin] = tfestimate(x,y,6000,[],[],fs); %stima della funzione di trasferimento, modello con legge di controllo lineare
[H_lineare,f_lineare] = tfestimate(x,y_lineare,[],[],[],fs); %stima della funzione di trasferimento, sistema dinamico lineare

d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.1,'DesignMethod','butter');

H_filt_lin = filtfilt(d1,mag2db(abs(H_leggelin)));

figure(7)
semilogx(2*pi*f_leggelin,H_filt_lin,'-r','LineWidth', 2.5)
hold on
semilogx(2*pi*f_lineare,mag2db(abs(H_lineare)),'-g','LineWidth', 1.5)
title('Bode Diagram Transfer Function - Legge di controllo lineare vs Configurazione passiva')
legend('Modello con legge di controllo lineare', 'Modello con configurazione passiva')
ylabel('Amplitude [dB]')
xlabel('Frequency [rad/s]');
xlim([1 10^3])
ylim([0 60])
grid
hold off

%{
figure(9)
semilogx(2*pi*f_leggelin,rad2deg(angle(H_leggelin)))
hold on
semilogx(2*pi*f_lineare,rad2deg(angle(H_lineare)))
title('Bode Diagram Transfer Function - Non Linear Model vs Linear Model')
legend('Non Linear Model', 'Linear Model')
ylabel('Phase [deg]')
xlabel('Frequency [rad/s]');
xlim([1 10^3])
ylim([-200 200])
yticks([-180 -90 0 90 180])
grid
hold off
%}

%% Ottimizzazione Bayesiana per legge di controllo (SOLO) QUADRATICA
clc
p_molt = 10000;

simulazioni = length(p_molt);

iter = zeros(simulazioni,1);
Jopt_quad = zeros(simulazioni,1);
uMR_op = zeros(simulazioni,1);
p11_t = zeros(simulazioni,1); p13_t = zeros(simulazioni,1); p14_t = zeros(simulazioni,1); 
p33_t = zeros(simulazioni,1); p34_t = zeros(simulazioni,1); p44_t = zeros(simulazioni,1);
bound_value_p = zeros(simulazioni,1);

%%Jmin = zeros(simulazioni,1);
best_J = 3;

N = [t',n];
uMR = 0;
out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
Jmin_p = rms(zs_2punti_nl);

f_max = 3000;
uMR_values = 0:150:f_max;

for i = 1:simulazioni
    
    Px_range = calcolorange_Pnn(uMR_values,f_max, p_molt(i));
    
    p11 = optimizableVariable('p11',[-Px_range(1),Px_range(1)]);
    p13 = optimizableVariable('p13',[-Px_range(2),Px_range(2)]);
    p14 = optimizableVariable('p14',[-Px_range(3),Px_range(3)]);
    p33 = optimizableVariable('p33',[-Px_range(4),Px_range(4)]);
    p34 = optimizableVariable('p34',[-Px_range(5),Px_range(5)]);
    p44 = optimizableVariable('p44',[-Px_range(6),Px_range(6)]);
    
    
    Px = [p11 p13 p14 p33 p34 p44];
    
    
    I_max = 150;
    bo_quad = bayesopt(@fun_quad,Px,...
                      'MaxObjectiveEvaluations', I_max, ...
                      'AcquisitionFunctionName', 'expected-improvement', ...
                      'explorationratio', 0.5, ...
                      'NumSeedPoints', 4, ...   
                      'PlotFcn', []);
    %{
    iterazioni = 1:bo_quad.NumObjectiveEvaluations;
    J_visited = bo_quad.ObjectiveTrace;
    J_min_visited = bo_quad.ObjectiveMinimumTrace;
    J_min_est = bo_quad.EstimatedObjectiveMinimumTrace;

    plot_fcnJ(iterazioni,J_visited, J_min_visited, J_min_est, 8);

    
    uMR = 0;
    out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
    output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
    zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
    Jmin(i) = rms(zs_2punti_nl);
    %}
                      
    iter(i) = bo_quad.NumObjectiveEvaluations;
    Jopt_quad(i) = bo_quad.MinEstimatedObjective;
    p11_t(i) = bo_quad.XAtMinEstimatedObjective.p11;
    p13_t(i) = bo_quad.XAtMinEstimatedObjective.p13;
    p14_t(i) = bo_quad.XAtMinEstimatedObjective.p14;
    p33_t(i) = bo_quad.XAtMinEstimatedObjective.p33;
    p34_t(i) = bo_quad.XAtMinEstimatedObjective.p34; 
    p44_t(i) = bo_quad.XAtMinEstimatedObjective.p34;
    
    bound_value_p(i) = max(Px_range);
    
    
    if bo_quad.MinEstimatedObjective < best_J
        best_J = bo_quad.MinEstimatedObjective;
        best_elementsP = [p11_t(i) 0 p13_t(i) p14_t(i); 0 0 0 0; p13_t(i) 0 p33_t(i) p34_t(i); p14_t(i) 0 p34_t(i) p44_t(i)];
        best_x = n;
    end
    
end

Jmin_p_t = Jmin_p + zeros(simulazioni,1);

table_results_pxdiag = table(iter, p11_t, p13_t, p14_t, p33_t, p34_t, p44_t, bound_value_p, Jopt_quad, Jmin_p_t);

%% grafico della FdT sospensione controllata da legge quadratica(set di parametri migliori) vs FdT sospensione passiva ottima

fs = 1/dt;
N = [t',n];
Px = best_elementsP;

out = sim("controllo_leggequadratica.slx");
zr = getElement(out.yout, 'zr');
x = zr.Values.Data(:);

xs2punti = getElement(out.yout, 'xs_2punti');
y_quad = xs2punti.Values.Data(:);

xs2punti_lineare = getElement(out.yout, 'xs_2punti_lineare');
y_lineare = xs2punti_lineare.Values.Data(:);

[H_leggequad,f_leggequad] = tfestimate(x,y_quad,6000,[],[],fs); %stima della funzione di trasferimento, legge quadratica
[H_lineare,f_lineare] = tfestimate(x,y_lineare,[],[],[],fs); %stima della funzione di trasferimento, sistema dinamico lineare

H_filt_quad = filtfilt(d1,mag2db(abs(H_leggequad)));

figure(8)
semilogx(2*pi*f_leggequad,H_filt_quad,'-b','LineWidth', 2.5)
hold on
semilogx(2*pi*f_lineare,mag2db(abs(H_lineare)),'-g','LineWidth', 1.5)
title('Bode Diagram Transfer Function - Legge di controllo quadratica vs Configurazione passiva')
legend('Modello con legge di controllo quadratica', 'Modello con configurazione passiva')
ylabel('Amplitude [dB]')
xlabel('Frequency [rad/s]');
xlim([1 10^3])
ylim([0 60])
grid
hold off

%%
figure(9)
semilogx(2*pi*f_leggequad,H_filt_quad,'-b','LineWidth', 2.5)
hold on
semilogx(2*pi*f_leggelin,H_filt_lin,'-r','LineWidth', 2.5)
hold on
semilogx(2*pi*f_lineare,mag2db(abs(H_lineare)),'-g','LineWidth', 1.5)
title('Legge di controllo quadratica vs Legge di controllo lineare vs Configurazione passiva')
legend('Modello con legge di controllo quadratica','Modello con legge di controllo lineare', 'Modello con configurazione passiva')
ylabel('Amplitude [dB]')
xlabel('Frequency [rad/s]');
xlim([1 10^3])
ylim([0 60])
grid
hold off


%% Ottimizzazione Bayesiana per legge di controllo LINEARE + QUADRATICA
clc
p_molt = 100000;

simulazioni = length(p_molt);

iter = zeros(simulazioni,1);
Jopt_linquad = zeros(simulazioni,1);
uMR_op = zeros(simulazioni,1);
k1_t = zeros(simulazioni,1); k3_t = zeros(simulazioni,1); k4_t = zeros(simulazioni,1); 
p13_t = zeros(simulazioni,1); p14_t = zeros(simulazioni,1); p34_t = zeros(simulazioni,1);
bound_value_p = zeros(simulazioni,1);
bound_value_k = zeros(simulazioni,1);

%%Jmin = zeros(simulazioni,1);
best_J = 3;

N = [t',n];
uMR = 0;
out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
Jmin_p = rms(zs_2punti_nl);

f_max = 3000;
uMR_values = 0:150:f_max;

for i = 1:simulazioni
    
    k_range = calcolorange_k(uMR_values,f_max,p_molt(i));
    Px_range = calcolorange_Pnn(uMR_values,f_max, p_molt(i));
    
    k1 = optimizableVariable('k1',[-k_range(1),k_range(1)]);
    k3 = optimizableVariable('k3',[-k_range(3),k_range(3)]);
    k4 = optimizableVariable('k4',[-k_range(4),k_range(4)]);
    
    p13 = optimizableVariable('p13',[-Px_range(2),Px_range(2)]);
    p14 = optimizableVariable('p14',[-Px_range(3),Px_range(3)]);
    p34 = optimizableVariable('p34',[-Px_range(5),Px_range(5)]);
    
    
    
    param = [p13 p14 p34 k1 k3 k4];
    
    I_max = 150;
    bo_lin_quad = bayesopt(@fun_lin_e_quad,param,...
                      'MaxObjectiveEvaluations', I_max, ...
                      'AcquisitionFunctionName', 'expected-improvement', ...
                      'explorationratio', 0.5, ...
                      'NumSeedPoints', 4, ...   
                      'PlotFcn', []);
    %{
    iterazioni = 1:bo_quad.NumObjectiveEvaluations;
    J_visited = bo_quad.ObjectiveTrace;
    J_min_visited = bo_quad.ObjectiveMinimumTrace;
    J_min_est = bo_quad.EstimatedObjectiveMinimumTrace;

    plot_fcnJ(iterazioni,J_visited, J_min_visited, J_min_est, 8);

    
    uMR = 0;
    out = sim("QuarterCarModel_nonLin_vs_Lin.slx");
    output_zs2punti_nl = getElement(out.yout, 'zs_2punti_nl');
    zs_2punti_nl = output_zs2punti_nl.Values.Data(:);
    Jmin(i) = rms(zs_2punti_nl);
    %}
                      
    iter(i) = bo_lin_quad.NumObjectiveEvaluations;
    Jopt_linquad(i) = bo_lin_quad.MinEstimatedObjective;
    p13_t(i) = bo_lin_quad.XAtMinEstimatedObjective.p13;
    p14_t(i) = bo_lin_quad.XAtMinEstimatedObjective.p14; 
    p34_t(i) = bo_lin_quad.XAtMinEstimatedObjective.p34;
    k1_t(i) = bo_lin_quad.XAtMinEstimatedObjective.k1;
    k3_t(i) = bo_lin_quad.XAtMinEstimatedObjective.k3;
    k4_t(i) = bo_lin_quad.XAtMinEstimatedObjective.k4; 
    
    bound_value_p(i) = max(Px_range);
    bound_value_k(i) = max(k_range);
    
    
    if bo_lin_quad.MinEstimatedObjective < best_J
        best_J = bo_lin_quad.MinEstimatedObjective;
        best_elementsP = [0 0 p13_t(i) p14_t(i); 0 0 0 0; p13_t(i) 0 0 p34_t(i); p14_t(i) 0 p34_t(i) 0];
        best_Ktilde = [k1_t(i) 0 k3_t(i) k4_t(i)];
        best_x = n;
    end
    
end

Jmin_p_t = Jmin_p + zeros(simulazioni,1);

table_results_fin = table(iter, k1_t, k3_t, k4_t, p13_t, p14_t, p34_t, bound_value_p, bound_value_k,Jopt_linquad, Jmin_p_t);


%% grafico della FdT sospensione controllata da legge lineare+legge quadratica

fs = 1/dt;
N = [t',n];
Px = best_elementsP;
k_tilde = best_Ktilde;

out = sim("controllo_leggelin_e_quad.slx");
zr = getElement(out.yout, 'zr');
x = zr.Values.Data(:);

xs2punti = getElement(out.yout, 'xs_2punti');
y = xs2punti.Values.Data(:);

xs2punti_lineare = getElement(out.yout, 'xs_2punti_lineare');
y_lineare = xs2punti_lineare.Values.Data(:);


[H_leggefin,f_leggefin] = tfestimate(x,y,6000,[],[],fs); %stima della funzione di trasferimento, legge quadratica
[H_lineare,f_lineare] = tfestimate(x,y_lineare,[],[],[],fs); %stima della funzione di trasferimento, sistema dinamico lineare

H_filt_fin = filtfilt(d1,mag2db(abs(H_leggefin)));


figure(8)
semilogx(2*pi*f_leggefin,H_filt_fin,'-','color',[0.9290 0.6940 0.1250],'LineWidth', 2.5)
hold on
semilogx(2*pi*f_leggequad,H_filt_quad,'-b','LineWidth', 2.5)
hold on
semilogx(2*pi*f_leggelin,H_filt_lin,'-r','LineWidth', 2.5)
hold on
semilogx(2*pi*f_lineare,mag2db(abs(H_lineare)),'-','color',[0.4660 0.6740 0.1880],'LineWidth', 1.5)
title('Bode Diagram Transfer Function - Legge di controllo quadratica vs Configurazione passiva')
legend('Modello con legge di controllo lin+quad','Modello con legge di controllo quadratica','Modello con legge di controllo lineare', 'Modello con configurazione passiva')
ylabel('Amplitude [dB]')
xlabel('Frequency [rad/s]');
xlim([1 10^3])
ylim([0 60])
grid
hold off