clear

%parameters
param=csvread('parametersNew.csv');

%number of runs
S = 100;
%periods
T = 876;
parfor s=1:S
    disp(s)
    feature('jit','off')
    feature('accel', 'off')
    [Y_real,defGDP,debGDP,labforce1,labforce2,Lev_a,Lev_ab,Lev_al,Lev_af,Lev_ak,defaults,defaults_k,dividendsB,average_interest_rate,average_interest_rate_k,Hdemand,CPI_b,CPI_l,dprob,lock_duration,detect_new,inf_detected,distance,inf_share2,sus_share2,rec_share2,healthcare,serious_share,inf_cum,inf_new,sus_share,inf_share,rec_share,dead_share,gdp_deflator, Investment,EXP, consumption, Un,stock_bonds,wages_t, desired_consumption, totK,I,totE,bankruptcy_rate,et] = model_policy(s,T, param);
    %this is the data which will be saved for multiple runs
    Y(:,s)  = Y_real;
    P(:,s)  = gdp_deflator;
    Ig(:,s) = I;
    Inv(:,s)  = Investment;
    C(:,s)  = consumption;
    DC(:,s) = desired_consumption;
    U(:,s) = 1-Un;
    W(:,s) = wages_t;
    elapsed(s) = et;
    G(:,s)=EXP;
    B(:,s)=stock_bonds;
    Infected(:,s)=inf_cum;
    Infected_new(:,s)=inf_new;
    Dead(:,s)=dead_share;
    Serious(:,s)=serious_share;
    Distance(:,s)=distance;
    Detected(:,s)=inf_detected;
    Detected_new(:,s)=detect_new;
    Lock_duration(s)=lock_duration;
    Detectprob(:,s)=dprob;
    bankrupt(:,s)=bankruptcy_rate;
    P_l(:,s)=CPI_l;
    P_b(:,s)=CPI_b;
    HealthD(:,s)=Hdemand;
    rate(:,s)=average_interest_rate;
    rate_k(:,s)=average_interest_rate_k;
    bankrupt2(:,s)=defaults+defaults_k;
    Lev_tot(:,s)=Lev_a;
    Lev_f(:,s)=Lev_af;
    Lev_k(:,s)=Lev_ak;
    Lev_b(:,s)=Lev_ab;
    Lev_l(:,s)=Lev_al;
    lab1(:,s)=labforce1;
    lab2(:,s)=labforce2;
    DefGDP(:,s)=defGDP;
    DebGDP(:,s)=debGDP;
end