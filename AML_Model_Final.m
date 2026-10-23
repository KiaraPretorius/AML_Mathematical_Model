close all

% White backgrounds for figures
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')


% Parameter values:

% Proliferation rates (h^-1)
p.r_H = 0.02888113;
p.r_S = 0.00163093;
p.r_R = 0.00062166;
p.r_B = 0.00309028;

% Natural cell death rates (h^-1)
p.gamma_H = 0.02888113;
p.gamma_S = 0.00031803;
p.gamma_R = 0.00002611;
p.gamma_B = 0.00309028;

% Cytarabine-induced cell killing rate (h^-1)
p.beta_H = 0.01444;
p.beta_S = 0.05020;
p.beta_R = 0.00170;
p.beta_B = 0.00677;

% Unitless parameter values
p.alpha = 2.21176e-9;
p.K = 1;
p.H0 = 0.05;
p.phi = 0.2497;


% Initial conditions:
% H, S, R, B
y0 = [0.05; 0.10; 0.01; 0.84];
tEnd = 17520;


% Chemotherapy settings:
chemoON = false;
secondChemoON = false;
chemoStart1 = 309.9656;
chemoEnd1 = chemoStart1 + 168;
chemoStart2 = 5280.0388;
chemoEnd2 = chemoStart2 + 168;

opts = odeset('RelTol',1e-9,'AbsTol',1e-12, ...
              'NonNegative',1:4,'MaxStep',1);


% Solve ODEs:
if ~chemoON
    [t,y] = ode45(@(t,y) AML_ODE(t,y,p,0),[0 tEnd],y0,opts);
else
    % Before chemotherapy
    [t1,y1] = ode45(@(t,y) AML_ODE(t,y,p,0), ...
                    [0 chemoStart1],y0,opts);

    % First chemotherapy cycle
    [t2,y2] = ode45(@(t,y) AML_ODE(t,y,p,1), ...
                    [chemoStart1 chemoEnd1],y1(end,:)',opts);

    if secondChemoON
        % Between cycles
        [t3,y3] = ode45(@(t,y) AML_ODE(t,y,p,0), ...
                        [chemoEnd1 chemoStart2],y2(end,:)',opts);
        % Second chemotherapy cycle
        [t4,y4] = ode45(@(t,y) AML_ODE(t,y,p,1), ...
                        [chemoStart2 chemoEnd2],y3(end,:)',opts);
        % After second cycle
        [t5,y5] = ode45(@(t,y) AML_ODE(t,y,p,0), ...
                        [chemoEnd2 tEnd],y4(end,:)',opts);

        t = [t1; t2(2:end); t3(2:end); t4(2:end); t5(2:end)];
        y = [y1; y2(2:end,:); y3(2:end,:); y4(2:end,:); y5(2:end,:)];
    else
        % After first cycle
        [t3,y3] = ode45(@(t,y) AML_ODE(t,y,p,0), ...
                        [chemoEnd1 tEnd],y2(end,:)',opts);
        t = [t1; t2(2:end); t3(2:end)];
        y = [y1; y2(2:end,:); y3(2:end,:)];
    end
end

H = y(:,1);
S = y(:,2);
R = y(:,3);
B = y(:,4);
N = H + S + R + B;


% Find when total blast percentage crosses 20%:
blasts = 100*(H+S+R)./N;
threshold = 20;

crossings = find(blasts(1:end-1) < threshold & blasts(2:end) >= threshold);

for i = 1:length(crossings)
    k = crossings(i);

    blastTime = interp1(blasts(k:k+1),t(k:k+1),threshold);

    fprintf('Crossing %d: %.4f hours\n',i,blastTime)
end


% Figures 1-4:

% Figure 1: Four populations
figure(1)
h1 = plot(t,H,'b-','LineWidth',1.5);
hold on
h2 = plot(t,S,'m-','LineWidth',1.5);
h3 = plot(t,R,'r-','LineWidth',1.5);
h4 = plot(t,B,'-','Color',[0.08 0.61 0.40],'LineWidth',1.5);
xlabel('Time (hours)')
ylabel('Normalised population')
title('Bone Marrow Population Dynamics in Acute Myeloid Leukaemia')
ylim([0 1])
grid on
if chemoON
    h5 = xregion(chemoStart1,chemoEnd1,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
    if secondChemoON
        xregion(chemoStart2,chemoEnd2,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5, ...
            'HandleVisibility','off');
    end
    legend([h1 h2 h3 h4 h5],{'H','S','R','B','C'},'Location','eastoutside')
else
    legend([h1 h2 h3 h4],{'H','S','R','B'},'Location','eastoutside')
end

% Figure 2: Total marrow population
figure(2)
plot(t,N,'k-','LineWidth',1.5)
hold on
yline(p.K,'r--','K = 1','LineWidth',1.5)
xlabel('Time (hours)')
ylabel('Total normalised population')
title('Total Normalised Bone Marrow Population')
ylim([0 1.1])
grid on

% Figure 3: Percentage of total bone marrow
figure(3)
plot(t,100*H./max(N,eps),'b-','LineWidth',1.5)
hold on
plot(t,100*S./max(N,eps),'m-','LineWidth',1.5)
plot(t,100*R./max(N,eps),'r-','LineWidth',1.5)
plot(t,100*B./max(N,eps),'-','Color',[0.08 0.61 0.40],'LineWidth',1.5)
xlabel('Time (hours)')
ylabel('Percentage of modelled BM cells (%)')
title('Relative Bone Marrow Cell Composition')
legend('H','S','R','B','Location','eastoutside')
ylim([0 100])
grid on

% Figure 4: Cytarabine schedule
figure(4)
if ~chemoON
    tC = [0 tEnd];
    Cvalues = [0 0];
elseif secondChemoON
    tC = [0 chemoStart1 chemoStart1 chemoEnd1 chemoEnd1 ...
          chemoStart2 chemoStart2 chemoEnd2 chemoEnd2 tEnd];
    Cvalues = [0 0 1 1 0 0 1 1 0 0];
else
    tC = [0 chemoStart1 chemoStart1 chemoEnd1 chemoEnd1 tEnd];
    Cvalues = [0 0 1 1 0 0];
end
plot(tC,Cvalues,'k-','LineWidth',1.5)
xlabel('Time (hours)')
ylabel('Chemotherapy C(t)')
title('Cytarabine Administration Schedule')
xlim([0 tEnd])
ylim([-0.1 1.1])
grid on


% Local ODE function:
function dydt = AML_ODE(~,y,p,C)
H = y(1);
S = y(2);
R = y(3);
B = y(4);

L = S + R;
N = H + S + R + B;
G = max(0,1-L/p.K);
V = max(0,p.K-N);

% Net leukaemic expansion responsible for healthy-cell displacement:
D = max(0,p.r_S*S*G+p.r_R*R*G-p.gamma_S*S-p.gamma_R*R);

% Healthy-cell displacement and recovery:
if H+B > 0
    D_H = D*H/(H+B);
    D_B = D*B/(H+B);
    Q_H = p.r_H*p.phi*max(p.H0*G-H,0)*B/(H+B);
    recovery_B = p.r_B*V*B/(H+B);
else
    D_H = 0;
    D_B = 0;
    Q_H = 0;
    recovery_B = 0;
end


% Differential equations:

dHdt = p.r_H*H-p.gamma_H*H+Q_H-p.beta_H*C*H-D_H;

dSdt = p.r_S*S*G-p.alpha*S*G-p.gamma_S*S-p.beta_S*C*S;

dRdt = p.r_R*R*G+p.alpha*S*G-p.gamma_R*R-p.beta_R*C*R;

dBdt = p.r_B*B+recovery_B-p.gamma_B*B-p.beta_B*C*B-D_B-Q_H;

dydt = [dHdt; dSdt; dRdt; dBdt];
end
