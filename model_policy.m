
function  [Y_real,defGDP,debGDP,labforce1,labforce2,Lev_a,Lev_ab,Lev_al,Lev_af,Lev_ak,defaults,defaults_k,dividendsB,average_interest_rate,average_interest_rate_k,Hdemand,CPI_b,CPI_l,dprob,lock_duration,detect_new,inf_detected,distance,inf_share2,sus_share2,rec_share2,healthcare,serious_share,inf_cum,inf_new,sus_share,inf_share,rec_share,dead_share,gdp_deflator, Investment,EXP, consumption, Un,stock_bonds,wages_t, DC, totK,I,totE,bankruptcy_rate,et] = model_policy(seed, T, par)

tic


%Parameters
Bk=1;                               %no. of banks
F=250; %200;%                         %no. of C-firms
F_l=150;                               % no. of "luxury" firms
F_b=100;                                 % no. of basic firms
W=2500; % 3000;%                        %no. of workers
N = 50;% 50;%                          %no. of capital producing firms
r_f=0.01;%   0.015;% 0.005          %general refinancing rate
r_f=r_f/3;                         % dividing by 12 to get a rough weekly rate
r_d=r_f/2;                          % deposit rate (markdown on risk-free rate)
%r_d=0;

%Parameters
z_c=par(1);                              %no. of aplications in consumption good market
z_k = par(2);                         %no. of aplications in capiatal good market
z_e = par(3);                            %number of job applications  

%%main seed
rng(seed)

%%create shocks et al
shock_pk = rand(T,N);
shock_p  = rand(T,F);
prob_k = rand(T,F);
permutations_consumption_b = NaN(T,W+N+F,z_c+1); %might happen that two are the same...
permutations_consumption_l = NaN(T,W+N+F,z_c+1); %need a set of permutations for both l and b goods; size needs to be z_c+1; one element will be replaced by largest firm from last period below
permutations_un = NaN(T,W,z_e);
permutations_capital = NaN(T,F,z_k+1);

for tt = 1:T
    for ii = 1:W+N+F
        permutations_consumption_b(tt,ii,:) = randperm(F_b,z_c+1);
        permutations_consumption_l(tt,ii,:) = F_b+randperm(F_l,z_c+1);
    end
    for ii = 1:W
        permutations_un(tt,ii,:) = randperm(F+N,z_e);
    end
    for ii = 1:F
        permutations_capital(tt,ii,:) = randperm(N,z_k+1);
    end
    
end

largest_b=randi([1 F_b],1,W+F+N);   %create initial "favourite" firm for each agent
largest_l=randi([(F_b+1) F],1,W+F+N);
largest_k=randi([1 N],1,F);

seeds_unemployed = randi(2^30,T);
seeds_capital = randi(2^30,T);
seeds_consumption = randi(2^30,T);
seeds_surplusk = randi(2^30,T,N);
seeds_surplus = randi(2^30,T,F);
seeds_infection2=randi(2^30,T*4);
seeds_hdem2=randi(2^30,T*4);
shocks_healthcare2=rand(T*4,W);
seeds_susceptible=randi(2^30,T*4);

%%initialise household network & SIR
age=zeros(1,W);
rage=rand(W,1);         %randomise agents' age --> age structure needs to be calibrated properly!
%age structure (15% young, 65% middle-aged, 20% old --> can be changed!)
age(rage<0.15)=1;       %young
age(rage>0.15)=2;       %middle-aged
age(rage>0.8)=3;        %old 
serious=zeros(1,W);     % if agent is infected, do they develop serious symptoms (==1)? set randomly according to age below
epidemic=0;                         %switch epidemic off (0) or on (1)
social_distancing=0;
dead_assets=0;
c_shock_l=ones(1,W);
c_shock_b=ones(1,W);
shock_l=2/3;
shock_b=1.2;
distanced=zeros(1,W);

if epidemic==1
shocks_postlock=rand(T,F);
shocks_distancing=rand(T*4,W);
seeds_connections=randi(2^30,T*4);
seeds_consampling=randi(2^30,T*4);
seeds_consampling2=randi(2^30,T*4);
seeds_hdem=randi(2^30,T*4);
shocks_healthcare=rand(T*4,W);
shocks_death=rand(T*4,W);
shocks_detection=rand(T*4,W);    

max_dur=randi([4,6],1,W);              % maximum duration of the disease
for i=1:length(serious)
    if age(i)==1
       if rand<0.01     %probabilities of getting serious symptoms --> parameter!
          serious(i)=1;
       end
    end
    if age(i)==2
       if rand<0.025     % parameter!
          serious(i)=1;
       end
    end
    if age(i)==3
       if rand<0.2      %parameter!
          serious(i)=1;
       end
    end
end

patient0=randperm(W,5);  %choose a few agents (5 in this case) to be the initial cases of the disease

con_mat=zeros(W,W);     %adjacency matrix for the network
ncons=((W*(W-1)/2));   %maximum number of possible connections
%parameters!
perm_cons=round(ncons/750); %number of permanent connections 
m_randcons=ncons/2800;      %mean of new random connections per period
sd_randcons=ncons/25000;     %sd of above
counter=0;
ordering=randperm(W,W);
for i=1:W                       %This loop is to make sure each agent has at least 1 permanent connection to another (could be changed)
   agent=ordering(i);
   if sum(con_mat(agent,:))==0   %check whether agent has already one permanent connection
      con=agent;                %choose agent to connect with
      while con==agent          %in case you accidentally choose yourself
      con=randperm(W,1);
      end
      con_mat(agent,con)=3;     % permanent connections are marked with a "3"
      con_mat(con,agent)=3;
      counter=counter+1;
   end
end

while counter<perm_cons      %make additional permanent connections until desired number is reached
con=randperm(W,2);           %choose 2 agents
    if con_mat(con(1),con(2))==0    %make sure connection does not already exist!
        con_mat(con(1),con(2))=3;
        con_mat(con(2),con(1))=3;
        counter=counter+1;
    end
end

con_mat=triu(con_mat,1);        %need the upper triangular part of the matrix -->
[row, col]=find(con_mat==3);    %--> in order to get a list of the permanent connections without double counting
fix_cons=[row col];
con_mat=(con_mat+con_mat') - eye(size(con_mat,1)).*diag(con_mat);       %copy the upper triangular part to the lower

end


%parameters
active=ones(1,W);       %which agents are active?
active(age==3)=0;       %Note that here (unless there is an epidemic) old and inactive codincide! i.e. there are no younger inactive! Could be changed if desired
labforce=sum(active);   %size of labour force
sickpay=zeros(1,W);     %to indicate which agent receives sickpay
dead=zeros(1,W);        %to record dead agents
sus_share=zeros(1,T*4);   %share of susceptible
inf_share=zeros(1,T*4);   %share of infected
inf_detected=zeros(1,T*4);   %share of infected & detected
rec_share=zeros(1,T*4);   %share of recovered
dead_share=zeros(1,T*4);  %share dead
inf_cum=zeros(1,T*4);     %cumulative cases
inf_new=zeros(1,T*4);     %new cases
distance=zeros(1,T*4);
ddetect=zeros(1,T*4);
detect_new=(zeros(1,T*4));
av_inf_new=zeros(1,T);
serious_share=zeros(1,T*4); %share of infections with serious symptoms
inf_share2=zeros(1,T*4);   %share of infected (normal disease)
sus_share2=zeros(1,T*4);   %share of susceptible (normal disease)
rec_share2=zeros(1,T*4);   %share of recovered (normal disease)

lockdown_exp=0;         %dummy for lockdown
CIGS=0;
labour_transfer=zeros(1,F);
labour_transfer_k=zeros(1,N);
liquidity_transfers=0;
liquidity_transfer=zeros(1,F);
liquidity_transfer_k=zeros(1,N);
credit_guarantees=0;
loans_g=0;
deb_g=zeros(1,F);
deb_g_k=zeros(1,N);
equity_support=0;
equity_injection=zeros(1,F);
equity_injection_k=zeros(1,N);
gov_ownership=zeros(1,F);
gov_ownership_k=zeros(1,N);
div_g=zeros(1,F);
div_g_k=zeros(1,N);
income_support=0;
%parameter
lf_share=round(F_l/3);
deactivate_f=F_b+1:F_b+lf_share;
active_f=ones(1,F);     %1=firm is active, 0=firm is inactive (lockdown)
active_k=ones(1,N);
homeoffice=zeros(1,F);
%parameters
constraint_k=ones(1,N);         %degree to which K-production is constrained during lockdown
constraint_f=ones(1,F);         %degree to which L-firms that implement smart working reduce production
lockdown=0;
lock_start=0;
reduce_cons=1;          
lock_duration=0;
lock_ended=0;
lock_min_duration=3;
lock_end_threshold=1;
post_lock_adjustment=1/3;
dis_cost_lockdown=-2;
dis_cost_init=2;
distancing_persistence=0.7;
reduce_cons_lockdown=0.25;
transmission_rate=0.185;
distancing_effect=0.4;
lockdown_threshold=5;
lock_constraint=0.9;

%"normal" disease
infprob2=0.0012;
susprob=0.1;
max_dur2=4*ones(1,W);
duration2=zeros(1,W);
susceptible2=ones(1,W);
infected2=zeros(1,W);
recovered2=zeros(1,W);

health_demanded=zeros(1,W);
health_demanded2=zeros(1,W);
health_supplied=zeros(1,W);
health_supplied2=zeros(1,W);
Hdemand=zeros(1,T*4);

h_1=0.55;
h_2=0.1;
h_3=0.06;

susceptible=ones(1,W);  %susceptible=1; not susceptible=0
infected=zeros(1,W);    
detected=zeros(1,W);    % has agent's infection been detected?
duration=zeros(1,W);    % if agent is infected how many periods will they remain ill?
recovered=zeros(1,W);   % has agent recovered?

distancing=zeros(1,W);
dprob=ones(1,T*4);
if epidemic==1
dis_index=-3*ones(1,W);
dis_prob=zeros(1,W);
intensity=[1 1.5 1];
dis_cost=2;
dis_threshold=1;

%parameters!
detectprob=0.01; % probability of being detected if positive (reflects number of tests conducted etc.)
detect_adjustment=0.0003;
dieprob=0.02;  % what is probability of dying (scaled by age below; conditional on serious==1)
dprob=detectprob*ones(1,T*4);
end

%parameters
xi = par(4);                          %memory parameter human wealth
chi = par(5);                            %fraction of wealth devoted to consumption
q_adj = par(6);                            %quantity adjustment parameter
p_adj = par(7);                            %price adjustment parameter    
mu =  par(8);                           %bank's gross mark-up
eta = par(9);                         %capital depreciation
Iprob=par(10);                         %probability of investing
phi =  par(11);                        %bank's leverage parameter
theta=par(12);                         %rate of debt reimbursment
delta = par(13);                        %memory parameter in the capital utilization rate
alpha=2/9;                              %labour productivity
k = 1/9;                            %capital productivity
div =par(16);                            %share of dividends
div_B=par(16);                                %bank dividends
barX=par(17);                          %desired capital utilization
inventory_depreciation = par(18);%0.3;           %rate at which capital firms' inventories depreciate
b1 = par(19);   
b2 = par(20);                       %Parameters for risk evaluation by banks
b_k1 = par(21); 
b_k2 = par(22); 

%parameters
interest_rate = par(23);
subsidy = par(24);                              %unemployment subsidy --> now roughly in line with net replacement rate for Italy acc. to OECD
tax_rate = par(25);                     %tax rate --> needs to be changed in line with unemployment & inactive subsidy & healthcare spending

%parameters
wage_update_up = par(26);
wage_update_down = par(27);
u_target = par(28);

%Government
Gov = 0;
%parameters
healthcare=zeros(1,T);                     %healthcare spending time series
healthcare(:)=0.04*labforce*alpha;         %healthcare spending is 4% of full employment output (needs to be calibrated but note that this shouldn't be calibrated to total healthcare spending over GDP!)

%phillips curve
wb=1/3;                                % initial wage rate


%else government will use unemployment subsidies as line below.
%parameters
bond_interest_rate = interest_rate;   %% ricordare di aggiungere anche questi con liquidità
unemployment_subsidy_init = subsidy;
inactive_sub=1.2;                     %subsidy received by inactives (as fraction of unemployment benefit) --> in line with pensions for Italy

G=zeros(1,T);                       %government expenditures on subsidies
Healthcare=zeros(1,T);              %nominal healthcare spending
TA=zeros(1,T);                      %government income
GB=zeros(1,T);                      %governament budget  GB = TA - subsidies - EXP - Healthcare - interest
GB(1)=-(1000+1/3*2800+5*250+5*50);
EXP = Gov*ones (1,T);       % spesa pubblica EROGABILE, update erogata in fondo
EXP(1,1)=0;               % inizializzazione, vedi sotto
pub_exp_c= zeros(1,T);    % fu il valore totale EROGATO alle consumption firms 
%parameter
tax_rate_f = par(29);        %taxes on dividends --> needs to be calibrated in line with gov. spending on subsidies & healthcare
tax_rate_k = par(29);
tax_rate_b = par(29);
exp_c=zeros(F,T);         % valore EROGATO individualmente per updating liquidity
quota_exp_c=zeros(F,T);   % quota singola erogabile impresa per totale
quota_exp_c(:,1)=1/F;     % time 1: occhio se cambi iniziando con t=2
public_dem_c=zeros(F,T);  % domanda pubblica per le consumption
quota_health=zeros(F+N,T); %share of healthcare spending going to each firm
quota_health(:,1)=1/(F+N);
RIC=ones(1,F);
RIC_k=ones(1,N);
bonds = zeros(1,T);
bonds(1)=-GB(1);             
stock_bonds=zeros(1,T);
stock_bonds(1)=bonds(1);
average_interest_rate=zeros(1,T);
average_interest_rate_k=zeros(1,T);
%%Bank's Parameters
b = [b1;b2];
b_k = [b_k1;b_k2];
                           
%%Initial conditions
%capital firm
Leff_k = zeros(1,N);
Y_k =zeros(1,N);
Y_prev_k=3/3*ones(1,N);
Y_kd=Y_prev_k;
P_k=3*ones(1,N);
A_k =5*ones(1,N);
A_control=A_k;
A_init_k=A_k;
liquidity_k = A_k;
De_k=1/3*ones(1,N);
deb_k=zeros(1,N);
price_k=zeros(1,(T+1));
price_k(1:2)=mean(P_k);                %capital price index
Q_k=zeros(1,N);                     %sales 
Ftot_k=zeros(1,N);
interest_r_k=zeros(1,N);
initialA_k=zeros(1,T);
%firms
value_investments=zeros(1,F);
investment=zeros(1,F);              %pshysical capital acquired in the period
K=15*ones(1,F);
A=5*ones(1,F)+ K*price_k(1); 
A_init=A;
liquidity=A-K*price_k(1);           %firm liquid resources
capital_value = K*price_k(1);
P=3*ones(1,F);                        %prices                  
Y_prev=5/3*ones(1,F);                     %past production
Yd=Y_prev;
Q=zeros(1,F);                       %sales   
Leff=zeros(1,F);                    %current employees
De=1/3*ones(1,F);                       %expected demand
deb=zeros(1,F);                     %firm debts
barK=K;
barYK=Y_prev/k;
x = barYK./barK;
interest_r=zeros(1,F);
%households
w=zeros(1,W);                       %earned wages
PA=1/3*ones(1,W+F+N);                 %household personal assets (saving stock)
Oc=zeros(1,W);                      %employment status: Oc(i)=j --> worker j employed by firm i; if i=0, j is unemployed
totE=zeros(1,T);                    %bank equity time series
E=1000*ones(1,Bk);                     %bank's equity
loans=sum(deb);                            %total loans
totE(1:2)=E;
%macro time series
consumption=zeros(1,T);             %total cunsumption                     
price=zeros(1,(T+1));                %consumer price index time series
price(1:2)=P(1);
Un=zeros(1,T);                      %unemployment
wages_t = NaN(1,T);                 %wages    
dividends=zeros(1,T);               %total dividends
Ftot=zeros(1,F);                    %borrowings
defaults=zeros(1,T);                %number of bankruptcies
profitsB=zeros(T,Bk);               %bank profit time series
dividendsB=zeros(T,Bk);
unsatisfiedDemand = zeros(1,T);
totK = zeros(1,T);
Investment = NaN(1,T);
defaults_k=zeros(1,T); 
bankruptcy_rate=zeros(1,T);
CPI_l=zeros(1,T);
CPI_l(1)=mean(P);                         %CPI for luxury goods
CPI_b=zeros(1,T);
CPI_b(1)=mean(P);                         %CPI for basic goods
DC=zeros(1,T);
Lev_a=zeros(1,T);
Lev_af=zeros(1,T);
Lev_ab=zeros(1,T);
Lev_al=zeros(1,T);
Lev_ak=zeros(1,T);
labforce1=zeros(1,T);
labforce2=zeros(1,T);
defGDP=zeros(1,T);
debGDP=zeros(1,T);


Y_nominal_k = NaN(1,T);
Y_nominal_c = NaN(1,T);
Y_nominal_tot= NaN(1,T);
Y_real =NaN(1,T);
Y_real(1)=0;
gdp_deflator = NaN(1,T);
I = NaN(1,T);
dividends_income_k = zeros(1,N);
dividends_income   = zeros(1,F);
permanent_income = 1/3*ones(F+W+N,1)./price(1);
permanent_income(W+1:W+F)=0.1/price(1);

money(1:2) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)+totE(1);

%%%HEREAFTER THE MODEL STARTS RUNNING FOR T PERIODS

 for t=2:T 

    if lockdown==1 && t-lock_start>=lock_min_duration && mean(detect_new(tt-1:tt))<=lock_end_threshold
        lockdown_exp=0;
        lockdown=0;                 %activate all firms
        sumY=sum(Yd(F_b+1:F_b+F_l));
        divisor=F_l-sum(active_f==0);
        meanY=sumY/divisor;
        Y_prev(active_f==0)=meanY;
        Yd(active_f==0)=(Y_prev(active_f==0));
        active_f(:)=1;
        active_k(:)=1;
        constraint_k(:)=1;
        lock_duration=tt-lock_start2;
        lock_ended=1;
    end
    
    if lock_ended==1
       for i=1:F_l
           if homeoffice(F_b+i)==1 && shocks_postlock(t,i)<post_lock_adjustment
              homeoffice(F_b+i)=0;
              constraint_f(F_b+i)=1;
           end
       end
    end
    
    
     %define shares by which gov expenditure (discretionary and health) is
     %distributed between firms
     if lockdown==0
     quota_exp_c(:,t)= (RIC./sum(RIC));  % quota di mercato calcolata come share entrate
     quota_health(1:F,t)=(RIC./(sum(RIC)+sum(RIC_k)));
     quota_health(F+1:F+N,t)=(RIC_k./(sum(RIC)+sum(RIC_k)));
     else
     rev=RIC;
     revk=RIC_k;
     rev(active_f==0)=0;
     revk(active_k==0)=0;
     quota_exp_c(:,t)= (rev./sum(rev));  % quota di mercato calcolata come share entrate
     quota_health(1:F,t)=(rev./(sum(rev)+sum(revk)));
     quota_health(F+1:F+N,t)=(revk./(sum(rev)+sum(revk)));
     end

    if t==648 && epidemic ==1    %infect the first individuals
    infected(patient0)=1;
    inf_cum(t*4)=inf_cum(t*4)+length(patient0);
    susceptible(patient0)=0;
    end
    Oc(dead==1)=-1;             %make sure that dead are not counted as employed or unemployed
    
    totK(t)=sum(K);
    
    %firms define expected demands and bid prices
    %note that given search and matching in the demand, the demand
    %perceived by the firms is overestimating actual demand. 
    stock = Y_prev-Yd;
    
    for i=1:F
        if active_f(i)==0   %if firm is in lockdown, produce nothing and keep price unchanged
            De(i)=0;
        else
        if stock(i)<=0 && P(i)>=price(t)
            De(i) =  Y_prev(i) + (-stock(i))*q_adj;
        elseif stock(i)<=0 && P(i)<price(t)
            De(i)=Y_prev(i);
            P(i)=P(i)*(1+shock_p(t,i)*p_adj);
            
        elseif stock(i)>0 &&P(i)>price(t)
            De(i)=Y_prev(i); 
            P(i)=P(i)*(1-shock_p(t,i)*p_adj);
           
        elseif stock(i)>0 && P(i)<=price(t)           
             De(i) = Y_prev(i) - stock(i)*q_adj;
        end    
        
        if De(i)<alpha                                                  %in order to hire at least one worker
            De(i)=alpha;    
        end
        end
    end
    
    
    %CAPITAL PRODUCTION DECISION
    %%capital good is durable, inventories depreciates
    inv_dep=inventory_depreciation*Y_k;
    inv_dep=inv_dep.*P_k;
    inventory_k=(1-inventory_depreciation)*Y_k;
    invent_start=inventory_k.*P_k;
    
    %for K-firms everything is as above for C-firms!
    stock_k = Y_prev_k-Y_kd; %virtual variation of inventories
    
    for i=1:N
        if active_k(i)==0
             De_k(i)=0;
        else
        if stock_k(i)<=0 && P_k(i)>=price_k(t)       
             De_k(i) = (Y_prev_k(i) + (-stock_k(i))*q_adj)-inventory_k(i);          
        elseif stock_k(i)<=0 && P_k(i)<price_k(t)
            De_k(i)=Y_prev_k(i);
            P_k(i)=P_k(i)*(1+shock_pk(t,i)*p_adj);      
        elseif stock_k(i)>0 && P_k(i)>price_k(t)
            De_k(i)=Y_prev_k(i); 
            P_k(i)=P_k(i)*(1-shock_pk(t,i)*p_adj);  
        elseif stock_k(i)>0 && P_k(i)<=price_k(t)
             De_k(i) = (Y_prev_k(i) - (stock_k(i))*q_adj)-inventory_k(i);
        end    
        
        if De_k(i)<alpha                                                  %in order to hire at least one worker
            De_k(i)=alpha;    
        end
        end
    end    
    
   
    
    %% investments
    %check who is allowed to invest and formulate investment demand
    prob=prob_k(t,:);
    K_dem=zeros(1,F);
    K_des=zeros(1,F);
    K_des(prob<Iprob) = barYK(prob<Iprob)/barX;
    %only used capital depreciates
    depreciation = barX*K_des*eta*1/(Iprob);
    K_dem(prob<Iprob)= K_des(prob<Iprob) - K(prob<Iprob)+depreciation(prob<Iprob);

    K_dem(K_dem<0)=0;
    
    %labour requirement (consumption good)
    Ld=min(ceil(De./alpha), ceil(K.*k/alpha));    %%demand for labour is the minimum between the actual production needs and the employable workers given the capital                                                 %labour demand approximated by excess
    wages=wb*Ld;   
   
    %the consumption good firms looks at the prices on the capital market
    %have an idea of the investment costs
    p_k=price_k(t);
    
   
    %labour requirement (capital good)
    Ld_k=ceil(De_k./alpha);
    wages_k = wb*Ld_k;

    
    if CIGS==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    labour_transfer(:)=0;
    labour_transfer_k(:)=0;
       for i=1:F
           if Leff(i)>Ld(i)
              labour_transfer(i)=(Leff(i)-Ld(i))*wb;
              liquidity(i)=liquidity(i)+labour_transfer(i);
              A(i)=A(i)+labour_transfer(i);
           end
       end
       for j=1:N
           if Leff_k(j)>Ld_k(j)
              labour_transfer_k(j)=(Leff_k(j)-Ld_k(j))*wb;
              liquidity_k(j)=liquidity_k(j)+labour_transfer_k(j);
              A_k(j)=A_k(j)+labour_transfer_k(j);
           end
       end
    else
        labour_transfer(:)=0;
        labour_transfer_k(:)=0;
    end



%% CREDIT MARKET OPENS
   
    Ftot(:)=0;
    Ftot_k(:)=0;
   
        
    %%THE FIRMS ASK FOR LIQUIDITY TO BUY THE CAPITAL AT A PRICE EQUAL TO
    %%THE ACTUAL AVERAGE CAPITAL PRICE
    
    %compute financial gap
    B=wages+K_dem*p_k-liquidity;                                        %financial gap          
    B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
    B(B<0)=0;
    B_k(B_k<0)=0;
    
    %compute leverage and use it to define bankruptcy probability and
    %interest rates
    lev=(deb+B)./(A+B+deb);                                             %leverage
    lev_k=(deb_k+B_k)./(A_k+B_k+deb_k);
    
    loan_applications=find(B>0);                                        %only firms with positive credit demand will apply for an additional loan
    loan_applications_k=find(B_k>0); 
    %evaluate bankruptcy probability and expected survival
    pr(:,1)= exp(b(1)+b(2)*lev)./(1+exp(b(1)+b(2)*lev));                           %banks evaluate bankruptcy probability of each firm
    pr_k(:,1)=exp(b_k(1)+b_k(2)*lev_k)./(1+exp(b_k(1)+b_k(2)*lev_k));                                                                                              
                                                           
    Xi=(1-(1-theta).^(1+1./pr))./(1-(1-theta));                         
    Xi_k=(1-(1-theta).^(1+1./pr_k))./(1-(1-theta));
    %proposed rate depends on the estimated bankruptcy probability
    proposed_rate=mu*((1+r_f/theta)./Xi - theta)';
    proposed_rate_k=mu*((1+r_f/theta)./Xi_k - theta)';
    

    %for each firm the bank computes the maximum loan and gives loans up to the maximum amount 
    for i=loan_applications
        credit=B(i);
        %the bank gives a maximum credit depending
        %on  maximum expected loss
        
        maxL = (phi*totE(t)-pr(i)*deb(i))/pr(i);
        maxL=max(0,maxL); %maxL never negative
        credit=min(credit,maxL); %%credit given to firm i           

        
        deb_0=deb(i);
        deb(i)=deb(i)+credit;                                           %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot(i)=credit;                                                 %record flow of new credit for firm i
        liquidity(i)=liquidity(i)+Ftot(i);
        %compute new average interest rate
        if deb(i)>0
            interest_r(i)=(deb_0*interest_r(i)+proposed_rate(i)*credit)/deb(i);
        end
    end       
    %weighted average interest rate
     average_interest_rate(t)=mean(interest_r);

     
    %mutatis mutandis for capital firms
    for i=loan_applications_k
        credit=B_k(i);
        maxL = (phi*totE(t)-pr_k(i)*deb_k(i))/pr_k(i);
        maxL=max(0,maxL);
        credit=min(credit,maxL);
         
        deb_0=deb_k(i);
        deb_k(i)=deb_k(i)+credit;                                       %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot_k(i)=credit;                                               %record flow of new credit for firm i
        liquidity_k(i)=liquidity_k(i)+Ftot_k(i);
        if deb_k(i)>0
            interest_r_k(i)=(deb_0*interest_r_k(i)+proposed_rate_k(i)*credit)/deb_k(i);
        end
    end 
    average_interest_rate_k(t)=mean(interest_r_k);

    %%CREDIT MARKET CLOSES       
    
    
    liquidity_transfer(:)=0;
    liquidity_transfer_k(:)=0;
    
    if liquidity_transfers==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    B=wages+K_dem*p_k-liquidity;                                        %financial gap          
    B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
    gap=find(B>0);
    gap_k=find(B_k>0);
    
    for i=gap
        liquidity_transfer(i)=liquidity_transfer(i)+B(i);
        liquidity(i)=liquidity(i)+liquidity_transfer(i);
        A(i)=liquidity(i)+capital_value(i)-deb(i)-deb_g(i);
    end
    for i=gap_k
        liquidity_transfer_k(i)=liquidity_transfer_k(i)+B_k(i);
        liquidity_k(i)=liquidity_k(i)+liquidity_transfer_k(i);
        A_k(i)=liquidity_k(i)+Y_k(i)*P_k(i)-deb_k(i)-deb_g_k(i);
    end
    end
%     if liquidity_transfers==1 && t>=lock_start && lock_start>0 && t<=lock_start+1
%      liquidity(:)=liquidity(:)+10;
%      liquidity_k(:)=liquidity_k(:)+10;
%      liquidity_transfer_k(:)=liquidity_transfer_k(:)+10;
%      liquidity_transfer(:)=liquidity_transfer(:)+10;
%      A_k=liquidity_k+Y_k.*P_k-deb_k-deb_g_k;
%      A=liquidity+capital_value-deb-deb_g;
%     end

    
    if credit_guarantees==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    B=wages+K_dem*p_k-liquidity;                                        %financial gap          
    B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
    gap=find(B>0);
    gap_k=find(B_k>0);
    
    for i=gap
        credit_g=B(i);       
        deb_g(i)=deb_g(i)+credit_g;                                           %update firm's debt stock
        loans_g=loans_g+credit_g;                                             %update bank's credit stock
        liquidity(i)=liquidity(i)+credit_g;
    end
    for i=gap_k
        credit_g_k=B_k(i);       
        deb_g_k(i)=deb_g_k(i)+credit_g_k;                                           %update firm's debt stock
        loans_g=loans_g+credit_g_k;                                             %update bank's credit stock
        liquidity_k(i)=liquidity_k(i)+credit_g_k;
    end
    end
    
    
    %% JOB MARKET OPENS
     
    %determine desired labour and vacancies given available liquidity
    Ld_k=min(Ld_k, (liquidity_k)/wb);
    Ld_k=floor(Ld_k);
    Ld_k(Ld_k<1)=1;
    vacancies_k=Ld_k-Leff_k;
    surplus_k=find(vacancies_k<0);                                          %firms with too many workers
    
    %%CONSUPTION GOOD

    %%re-define labour demand given available liquidity
    Ld=min(Ld,(liquidity)/wb);                                     %demand (stock)     
    Ld=floor(Ld);
    Ld(Ld<1)=1; %%since there is the capital the firms can have positive equity and negative liquidity, in this latter case Ld would be negative, which is impossible
    vacancies=Ld-Leff;                                                 %definitive labour demand (flow)    
    
    %%JOB MARKET OPENS
    surplus=find(vacancies<0);                                          %firms with too many workers
    
    if CIGS==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    else
    for i=surplus_k
        workforce_k=find(Oc==F+i);
        rng(seeds_surplusk(t,i))
        f_k=randperm(length(workforce_k));
        f_k=f_k(1:-vacancies_k(i));                                     %take randomly "-vacancies(i)" workers and fire them
        fired_k=workforce_k(f_k);  
        Oc(fired_k)=0;
        w(fired_k)=0;
        Leff_k(i)=Ld_k(i);                                              %update no. of workers
    end
    end


 %firms with excess workforce fire
    if CIGS==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
    else
    for i=surplus
        workforce=find(Oc==i);
        rng(seeds_surplus(t,i))
        f=randperm(length(workforce));
        f=f(1:-vacancies(i));                                           %take randomly "-vacancies(i)" workers and fire them
        fired=workforce(f);  
        Oc(fired)=0;
        w(fired)=0;
        Leff(i)=Ld(i);                                                  %update no. of workers
    end
    end
    
    
    
%% UNEMPLOYED WORKERS LOOK FOR A JOB
    
    %only active unemployed participate in labour market
    unemployed=find(Oc==0 & active==1);
    
    rng(seeds_unemployed(t))
    vec = randperm(length(unemployed));
    for un=vec      
        j=unemployed(un);                                               %randomly pick an unemployed worker
        %inactives do not participate in labour market
        Z_e = permutations_un(t,j,:);
        flag=1;        
        while (Oc(j)==0 && flag<=z_e)                              %continue searching until you are unemployed and didn't search at all available firms
            f=Z_e(flag);                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            if f>F %selected firm is a capital firm
                if vacancies_k(f-F)>0  && active_k(f-F)==1                                  %if the selected firm has an open vacancy, take the job
                Oc(j)=f;                                                 %update employed status
                w(j)=wb;                                                 %salary   
                Leff_k(f-F)=Leff_k(f-F)+1;                               %firm's workforce   
                vacancies_k(f-F)=vacancies_k(f-F)-1;                     %firm's vacancies   
                end
            else %selected firm is a consuption firm
                if vacancies(f)>0 && active_f(f)==1                                       %if the selected firm has an open vacancy, take the job
                Oc(j)=f;
                w(j)=wb;
                Leff(f)=Leff(f)+1;
                vacancies(f)=vacancies(f)-1;
                end                
            end
            flag=flag+1;                                                %increase counter
        end
        
    end  
    
    %%JOB MARKET CLOSES
    %encounters in the workplace
    if epidemic==1 && t>=648 && sum(infected)>0 
    for j=1:(F)
        if active_f(j)==1 && homeoffice(j)==0
        colleagues=find(Oc==j);
        if length(colleagues)>1
        work_cons=nchoosek(colleagues,2);
        for l=1:size(work_cons,1)
            con=work_cons(l,:);
            if con_mat(con(1),con(2))==0    
               con_mat(con(1),con(2))=2;    %workplace connections are marked with a "2"
               con_mat(con(2),con(1))=2;
            end
        end
        end
        end
    end
    for j=1:(N)
        if active_k(j)==1
        colleagues=find(Oc==F+j);
        if length(colleagues)>1
        work_cons=nchoosek(colleagues,2);
        for l=1:size(work_cons,1)
            con=work_cons(l,:);
            if con_mat(con(1),con(2))==0
               con_mat(con(1),con(2))=2;
               con_mat(con(2),con(1))=2;
            end
        end
        end
        end
    end
    end
    
     
    %% production
    %produce capital
    Y_k=min(De_k,constraint_k.*Leff_k*alpha);                              
    Y_prev_k=Y_k;
    Y_k = Y_k+inventory_k;                                   %Y_k is increased by the inventories (capital good is durable)   
    %produce consuption
    Yp = min(constraint_f.*Leff*alpha, K*k);                               %production frontier given available resources (Leontieff)
    Y=min(De,Yp);                                            %actual production
    Y_prev=Y;
    
    %here we compute the average production cost to impose the minimum
    %price
    %compute interests to be paid at the end of the period
    %here we compute the average production cost to impose the minimum
    %price
    %compute interests to be paid at the end of the period
    interests=interest_r.*deb+r_f.*deb_g; 
    interests_k=interest_r_k.*deb_k+r_f.*deb_g_k;
    %total wages paid by firms
    wages=wb*Leff;
    wages_k=wb*Leff_k;
    liquidity=liquidity-wages;
    liquidity_k=liquidity_k-wages_k;

    
    %minimum prices
    
    Pl = 1.01*(wages-labour_transfer+interests+Y/k*eta*p_k+theta*deb+theta*deb_g)./Y_prev;
    Pl(isinf(Pl)) = 0;
    P(P<Pl) = Pl(P<Pl);

    
    %minimum prices capital
    Pl_k = 1.01*(wages_k-labour_transfer_k+interests_k+theta*deb_k+theta*deb_g_k)./(Y_prev_k);
    Pl_k(isinf(Pl_k)) = 0;
    P_k(P_k<Pl_k) = Pl_k(P_k<Pl_k);
    

  
    %%CAPITAL GOODS MARKET OPENS
    %%CONSUMPTION FIRMS BUY THE CAPITAL THEY NEED
    % The capital is bought now and will be used for production in the next
    % period. The capital market is a search and matching as the consumption
    % market

    capital_budget=liquidity;                              %amount of liquidity available to buy capital goods                                  
    capital_demanders=1:F;
    investment(:)=0;
    value_investments(:)=0;
    Q_k(:)=0;
    Y_kd=zeros(1,N);  %%record demand
    
    %Government buys goods from K-sector for healthcare; share based on
    %previous revenue (see below)
    health_demand_k=zeros(1,N);
    hdk=quota_health(F+1:F+N,t)*healthcare(t);
    while sum(hdk)>0 && sum(Y_k)>0
    for j=1:N
        Y_kd(j) = Y_kd(j) + hdk(j);
        if Y_k(j)>= hdk(j)
           Y_k(j)=Y_k(j)-hdk(j);   
           Q_k(j)= Q_k(j)+hdk(j);
           health_demand_k(j)=health_demand_k(j)+hdk(j);
           hdk(j)=0;
        else 
           Q_k(j)=Q_k(j)+Y_k(j);
           health_demand_k(j)=health_demand_k(j)+Y_k(j);              
           hdk(j)=hdk(j)-Y_k(j);
           Y_k(j)= 0;
        end
    end
    hresid=sum(hdk);
    quota_new=Y_k/sum(Y_k);
    hdk=quota_new*hresid;
    end
    Healthcare_k= sum(health_demand_k.*P_k);
    healthcare_k=sum(health_demand_k);
    
    
    %the capital market works exactly as the consumption market (below).
    %Search and matching!
    %take largest firm from previous period
    l_k=largest_k(capital_demanders);
    rng(seeds_capital(t))
    for firm=randperm(length(capital_demanders))
        j=capital_demanders(firm);
        %randomly chosen firms
        Z = permutations_capital(t,j,:);
        %if random firms include already the largest one from last period
        if sum(Z==l_k(firm))>0
           duplicate=Z==l_k(firm);
           PZ_k=P_k(Z);
           PZ_k(duplicate==1)=P_k(l_k(firm));
        else
        PZ_k=P_k(Z);
        Z(1,1,z_k+1)=l_k(firm);
        PZ_k(z_k+1)=P_k(l_k(firm));
        end
        [~,order]=sort(PZ_k);
        flag=1;
        visited_k=zeros(1,z_k+1);
        while (capital_budget(j)>0 && K_dem(j)>0 && flag<=z_k+1) 
            best=Z(order(flag));
            Y_kd(best)=Y_kd(best)+min(capital_budget(j)/P_k(best), K_dem(j));
                if Y_k(best)>0                                                %buy if 'best' firm has still positive stocks 
                    pk=P_k(best);
                    budget=min(capital_budget(j), K_dem(j)*pk);
                    if Y_k(best) > budget/pk           
                        Y_k(best)=Y_k(best)-budget/pk;                   %reduce stocks
                        Q_k(best)=Q_k(best)+budget/pk;                   %update sales
                        K_dem(j)=K_dem(j)-budget/pk;
                        capital_budget(j) = capital_budget(j)-budget;
                        liquidity(j)=liquidity(j)-budget;
                        investment(j)=investment(j)+budget/pk;
                        value_investments(j)= value_investments(j)+budget;
                        visited_k(order(flag))=budget/pk;                                   %j spends all its budget  
                    elseif  Y_k(best) <= budget/pk  
                        K_dem(j)=K_dem(j)-Y_k(best);
                        capital_budget(j)=capital_budget(j) - Y_k(best)*pk;
                        liquidity(j)=liquidity(j)- Y_k(best)*pk;
                        Q_k(best)=Q_k(best)+Y_k(best);    
                        investment(j)=investment(j)+Y_k(best);
                        value_investments(j)= value_investments(j)+Y_k(best)*pk;
                        visited_k(order(flag))=Y_k(best);
                        Y_k(best)=0;    
                    end    
                end
                flag=flag+1;                                                %increase counter
        end
        if sum(visited_k)>0
           largest_kfirm=find(visited_k==max(visited_k));
           largest_kfirm=largest_kfirm(1);
           largest_kfirm=Z(largest_kfirm);
           largest_k(j)=largest_kfirm;
       end
    end


    %%CAPITAL GOOD MARKET CLOSES
 
    %%CONSUMPTION GOOD MARKET OPENS
    w(w>0) = wb;
    w(Oc==0)=0;
    w(dead==1)=0;
    PA(dead==1)=0;
    %taxes on wages
    wn=w*(1-tax_rate);  
    TA(t) = TA(t)+sum(w*tax_rate);
    %workers receive wages
    interest_workers=r_d*PA(1:W);
    interest_deposits=sum(r_d*PA(1:W));
    PA(1:W)=PA(1:W)+r_d*PA(1:W);
    PA(1:W)=PA(1:W)+wn;
    

    unemployment_subsidy = unemployment_subsidy_init;    
    
    %pay unemployement subsidy to unemployed, sickpay to sick and inactive income to
    %inactive
    oPA=sum(PA(Oc==0 & active==1 & dead==0));
    nPA=sum(PA(Oc==0 & active==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate));
    dPA=nPA-oPA;
    G(t)=dPA;
    PA(Oc==0 & active==1 & dead==0)=PA(Oc==0 & active==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate);
    oPA=sum(PA(active==0 & sickpay==0 & dead==0));
    nPA=sum(PA(active==0 & sickpay==0 & dead==0)+inactive_sub*unemployment_subsidy*wb*(1-tax_rate));
    dPA=nPA-oPA;
    G(t)=G(t)+dPA;
    PA(active==0 & sickpay==0 & dead==0)=PA(active==0 & sickpay==0 & dead==0)+inactive_sub*unemployment_subsidy*wb*(1-tax_rate);
    oPA=sum(PA(active==0 & sickpay==1 & dead==0));
    nPA=sum(PA(active==0 & sickpay==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate));      %assume sickpay=unemployment benefit
    PA(active==0 & sickpay==1 & dead==0)=PA(active==0 & sickpay==1 & dead==0)+unemployment_subsidy*wb*(1-tax_rate); 
    dPA=nPA-oPA;
    G(t)=G(t)+dPA;
    if income_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+5
       PA(Oc>=0)=PA(Oc>=0)+0.5*unemployment_subsidy*wb*(1-tax_rate);
       G(t)=G(t)+unemployment_subsidy*wb*(1-tax_rate)*length(PA(Oc>=0));
    end
    
    workers_income = wn;
    workers_income(Oc==0 & active==1) = unemployment_subsidy*wb*(1-tax_rate);
    workers_income(active==0 & sickpay==0) = inactive_sub*unemployment_subsidy*wb*(1-tax_rate);
    workers_income(active==0 & sickpay==1) = unemployment_subsidy*wb*(1-tax_rate);
    if income_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+5
    workers_income = workers_income+0.5*unemployment_subsidy*wb*(1-tax_rate);
    end
    workers_income(dead==1)=0;
    workers_income=workers_income+interest_workers;
   
    
    %consumers compute their consumption budgets
    income =  ([workers_income,dividends_income,dividends_income_k]')./price(t); 

    %%what about taking the rolling mean of incomes?
    %standard deviation of consumptio = standard deviation of wealth
    permanent_income = permanent_income*xi + (1-xi)*income;

    %http://books.google.it/books?id=DWnhhYdmAucC&pg=PA93&lpg=PA93&dq=permanent+income+rules+in+agent+based+modeling&source=bl&ots=psxGRexvx2&sig=Qo2QXcsxtd2K0sQlT2rOC7iSWk8&hl=en&sa=X&ei=DLLfUpDNBcXWtAbmk4DACA&ved=0CFAQ6AEwBDgK#v=onepage&q&f=false
    target = 1*permanent_income' + chi*PA./price(t) ; %0.05
    target_b=F_b/F*target.*CPI_l(t-1)/CPI_b(t-1);
    target_l=target-target_b;
    %compute a budget for basic goods as fraction of total cons. budget.
    %Here this is also a function of relative prices CPI_l and CPI_b
    %lexicographic ordering (only spend on luxury if sufficient money
    %available)
    cons_budget_b=min(PA,target_b.*CPI_b(t-1));
    cons_budget_l=max(0,min(PA-cons_budget_b,target_l.*CPI_l(t-1)));
    if t>=648 && epidemic==1
    for i=1:W
    if distancing(i)==1
    cons_budget_b(i)=min(PA(i),c_shock_b(i)*target_b(i)*CPI_b(t-1));
    cons_budget_l(i)=max(0,min(PA(i)-cons_budget_b(i),target_l(i)*CPI_l(t-1)));   
    cons_budget_l(i)=c_shock_l(i)*cons_budget_l(i);
    c_shock_b(i)=max(1,(1-post_lock_adjustment)*c_shock_b(i)+post_lock_adjustment);
    c_shock_l(i)=min(1,(1-post_lock_adjustment)*c_shock_l(i)+post_lock_adjustment);
    end
    end
    end
    PA=PA-cons_budget_b-cons_budget_l;        
    consumers=1:(W+F+N);
    
    DC(t) = sum(cons_budget_b)/CPI_b(t-1)+sum(cons_budget_l)/CPI_l(t-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INTRO SPESA PUBBLICA
    Q(:)=0;
    Yd=zeros(1,F);
    
    %government buys from all C-firms (including "luxury" --> can view this
    %as consumption by doctors??)
    health_demand_f=zeros(1,F);
    hd=quota_health(1:F,t)*healthcare(t);
    while sum(hd)>0 && sum(Y)>0
    for i=1:F
        Yd(i) = Yd(i) + hd(i);
        if Y(i)>= hd(i)
           Y(i)=Y(i)-hd(i);   
           Q(i)= Q(i)+hd(i);
           health_demand_f(i)=health_demand_f(i)+hd(i);
           hd(i)=0;
        else 
           Q(i)=Q(i)+Y(i);
           health_demand_f(i)=health_demand_f(i)+Y(i);
           hd(i)=hd(i)-Y(i);
           Y(i)= 0;      
        end
    end
    hresid=sum(hd);
    quota_new=Y/sum(Y);
    hd=quota_new*hresid;
    end
    Healthcare_f= sum(health_demand_f.*P);
    healthcare_f=sum(health_demand_f);
    
    Healthcare(t)=Healthcare_k+Healthcare_f;
    healthcare(t)=healthcare_k+healthcare_f;
   %search and matching starts
    
    C = zeros(1,W+F+N);
    
    %same as with capital goods market
    l_b=largest_b(consumers);
    l_l=largest_l(consumers);
    rng(seeds_consumption(t))
    vec = randperm(W+F+N);
    for wor=vec     
        j=consumers(wor);                                               %randomly pick a consumer
        if j<=W && dead(j)==1
        else
        Z = permutations_consumption_b(t,j,:);
        if sum(Z==l_b(wor))>0
           duplicate=Z==l_b(wor);
           PZ=P(Z);
           PZ(duplicate==1)=P(l_b(wor));
        else
        PZ=P(Z);
        Z(1,1,z_c+1)=l_b(wor);
        PZ(z_c+1)=P(l_b(wor));
        end
        [~,order]=sort(PZ);                                            %sort prices in ascending order
        flag=1;
        visited=zeros(1,z_c+1);
        %first shop with basic firms, luxury below
        while (cons_budget_b(j)>0 && flag<=z_c+1)                              %continue buying till budget is positive and there are firms available
            best=Z(order(flag));                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            Yd(best)=Yd(best)+cons_budget_b(j)/P(best);
            if Y(best)>0                                                %buy if 'best' firm has still positive stocks 
                p=P(best);
                if Y(best) > cons_budget_b(j)/p           
                    Y(best)=Y(best)-cons_budget_b(j)/p;                   %reduce stocks
                    Q(best)=Q(best)+cons_budget_b(j)/p; 
                    C(j) = C(j)+cons_budget_b(j)/p; %update sales
                    consumption(t)=consumption(t)+cons_budget_b(j)/p;     
                    visited(order(flag))=cons_budget_b(j)/p;
                    cons_budget_b(j)=0;                                   %j spends all its budget  
                elseif  Y(best) <= cons_budget_b(j)/p  
                    cons_budget_b(j)=cons_budget_b(j)- Y(best)*p;   
                    Q(best)=Q(best)+Y(best);
                    C(j) = C(j)+Y(best); 
                    consumption(t)=consumption(t)+Y(best);
                    visited(order(flag))=Y(best);
                    Y(best)=0;
                end    
            end
            flag=flag+1;                                                %increase counter
        end 
       if sum(visited)>0
           largest_bf=find(visited==max(visited));
           largest_b(j)=Z(largest_bf(1));
       end
        Z = permutations_consumption_l(t,j,:);
        if sum(Z==l_l(wor))>0
           duplicate=Z==l_l(wor);
           PZ=P(Z);
           PZ(duplicate==1)=P(l_l(wor));
        else
        PZ=P(Z);
        Z(1,1,z_c+1)=l_l(wor);
        PZ(z_c+1)=P(l_l(wor));
        end
        [~,order]=sort(PZ);                                            %sort prices in ascending order
        flag=1;
        visited=zeros(1,z_c+1);
        while (cons_budget_l(j)>0 && flag<=z_c+1)                              %continue buying till budget is positive and there are firms available
            best=Z(order(flag));                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            if active_f(best)==1
            Yd(best)=Yd(best)+cons_budget_l(j)/P(best);
            if Y(best)>0                                                %buy if 'best' firm has still positive stocks 
                p=P(best);
                if Y(best) > cons_budget_l(j)/p           
                    Y(best)=Y(best)-cons_budget_l(j)/p;                   %reduce stocks
                    Q(best)=Q(best)+cons_budget_l(j)/p; 
                    C(j) = C(j)+cons_budget_l(j)/p; %update sales
                    consumption(t)=consumption(t)+cons_budget_l(j)/p;     
                    visited(order(flag))=cons_budget_l(j)/p;
                    cons_budget_l(j)=0;                                   %j spends all its budget  
                elseif  Y(best) <= cons_budget_l(j)/p  
                    cons_budget_l(j)=cons_budget_l(j)- Y(best)*p;   
                    Q(best)=Q(best)+Y(best);
                    C(j) = C(j)+Y(best); 
                    consumption(t)=consumption(t)+Y(best);
                    visited(order(flag))=Y(best);
                    Y(best)=0;
                end    
            end
            end
            flag=flag+1;                                                %increase counter
        end
           if sum(visited)>0
               largest_lf=find(visited==max(visited));
               largest_l(j)=Z(largest_lf(1));
           end
       end
    end    
    
    
    %search and matching ends
    %unvoluntary savings are added to the vuluntary ones
    unsatisfiedDemand(t) = sum(cons_budget_b)+sum(cons_budget_l);            
    PA=PA+cons_budget_b+cons_budget_l;
    
    % dopo vende al settore pubblico, se ha prodotto abbastanza
    %this is discretionary gov. expenditure
    pub_exp_c(t)= 1*EXP(t);
    public_dem_c(:,t)=quota_exp_c(:,t).*(pub_exp_c(t)./P');  % dove pub_exp è quella EROGABILE
    for i=1:F
        Yd(i) = Yd(i) + public_dem_c(i,t);
        if Y(i)>= public_dem_c(i,t)
           Y(i)=Y(i)-public_dem_c(i,t);     % riduco lo stock di offerta
           Q(i)= Q(i)+public_dem_c(i,t);    % ho venduto l'ammontare "public dem"
        else 
           Q(i)=Q(i)+Y(i);
           public_dem_c(i,t)=Y(i);              % così riflette il fatto che l'impresa potrebbe non avere
           Y(i)= 0;      
        end                                     % abbastanza offerta per soddisfare la domanda pubblica
    end
    exp_c(:,t)= public_dem_c(:,t).*P';     % per aggiornamento EXP effettivamente erogate complessivam
    pub_exp_c(t)=sum(exp_c(:,t));
    EXP(t)=pub_exp_c(t); % qui sempre per la storia che se faccio come lo volevo fare mi si impalla
    
%%CONSUMPTION GOOD MARKET CLOSES
    %Capital price index
    if sum(Q_k)>0 
        share_k=(Q_k)/sum(Q_k);
    else
        share_k=ones(1,N)/N;
    end
    
    price_k(t+1)=P_k*share_k';

%%ACCOUNTING

    %capital capacity update
    barYK = delta*barYK + (1-delta)*Y_prev/k;
    %capital depreciation --> this needs to be changed!!!
    dep = eta*Y_prev/k; 
    
    %capital reduced by depreciation
    K = K - dep;
    %capital increased by bought capital 
    K = K + investment;
    I(t) = sum(investment);

    %update capital value in the book
    depreciation_value =  dep.*capital_value./(K+dep-investment);
    capital_value = capital_value - depreciation_value + value_investments;
    
    %firm revenues
    RIC=P.*Q;   
    %firm uodate their liquidity, pay wages, interests and installments
    %consumption firm profits
    pi = RIC-wages-interests - depreciation_value;
    liquidity=liquidity+RIC-interests-theta*deb-theta*deb_g; %%investments already paid
    %loans are updated (the installments have been paid
    loans=loans-sum(deb)*theta;
    loans_g=loans_g-sum(deb_g)*theta;
    deb=(1-theta)*deb;
    deb_g=(1-theta)*deb_g;
    
    
    %equity law of motion!
    A = A+pi;
    
    %Capital producer accounting
    RIC_k = P_k.*Q_k;
    %update liquidity
    invent_end=Y_k.*P_k;
    pi_k=RIC_k+(invent_end-invent_start)-wages_k-interests_k-inv_dep;
    liquidity_k = liquidity_k + RIC_k-interests_k-theta*deb_k-theta*deb_g_k;
    %update loand
    loans=loans-sum(deb_k)*theta;
    loans_g=loans_g-sum(deb_g_k)*theta;
    deb_k=(1-theta)*deb_k;
    deb_g_k=(1-theta)*deb_g_k;
    
    %dividends
    div_g(:)=0;
    div_g_k(:)=0;
    pospi=find(pi>0 & liquidity>0);                                                   %pick firms with positive profits
    dividends_income(:)=0;
    for i=pospi
        TA(t)=TA(t)+pi(i)*tax_rate_f;
        tax_pr=pi(i)*tax_rate_f;
        pi(i)=(1-tax_rate_f)*pi(i);
        liquidity(i)=liquidity(i)-tax_pr;
        A(i)=A(i)-tax_pr;
        di=max(0,min(div*pi(i),liquidity(i)));  %%dividends                                              %compute dividends paid by firm i
        div_g(i)=gov_ownership(i)/(gov_ownership(i)+A_init(i))*di;
        divi=di-div_g(i); %dividends after taxes
        PA(W+i)=PA(W+i)+divi;                                          %dividends paid to firm i owner
        dividends_income(i) = divi;
        liquidity(i)=liquidity(i)-di;
        A(i)=A(i)-di;
        dividends(t)=dividends(t)+di; %lordi
        pi(i)=pi(i)-di;
       
    end
    dividends_income=dividends_income+r_d*PA(W+1:W+F);
    interest_deposits=interest_deposits+sum(r_d*PA(W+1:W+F));
    PA(W+1:W+F)=PA(W+1:W+F)+r_d*PA(W+1:W+F);
    
    A_control=A_control+pi_k;
    pospi_k=find(pi_k>0 & liquidity_k>0); 
    dividends_income_k(:)=0;%pick firms with positive profits
    for i=pospi_k
        TA(t)=TA(t)+pi_k(i)*tax_rate_k;
        tax_pr_k=pi_k(i)*tax_rate_k;
        pi_k(i)=(1-tax_rate_k)*pi_k(i);
        liquidity_k(i)=liquidity_k(i)-tax_pr_k;
        A_control(i)=A_control(i)-tax_pr_k;
        di=max(0,min((div)*pi_k(i),liquidity_k(i)));
        div_g_k(i)=gov_ownership_k(i)/(gov_ownership_k(i)+A_init_k(i))*di;
        divi=di-div_g_k(i); %dividends after taxes
        PA(W+F+i)=PA(W+F+i)+divi;                                          %dividends paid to firm i owner
        dividends_income_k(i)=divi;
        A_control(i)=A_control(i)-di;
        liquidity_k(i)=liquidity_k(i)-di;
        dividends(t)=dividends(t)+di;
        pi_k(i)=pi_k(i)-di;
       
    end
    dividends_income_k=dividends_income_k+r_d*PA(W+F+1:W+F+N);
    interest_deposits=interest_deposits+sum(r_d*PA(W+F+1:W+F+N));
    PA(W+F+1:W+F+N)=PA(W+F+1:W+F+N)+r_d*PA(W+F+1:W+F+N);
   
    A_k=liquidity_k+Y_k.*P_k-deb_k-deb_g_k;

    
    %replacement of bankrupted consumption firms 
    piB=0;                                                              %reset bank's profits
       
    
    %% time series (before bankruptcies)
    %inflation rate
    if sum(Q)>0 
        share=(Q)/sum(Q);
    else
        share=ones(1,F)/F;
    end
    RPI=P*share';                                                       %retail price index
    price(t+1)=RPI;
    P_b=P(1:F_b);
    Q_b=Q(1:F_b);
    if sum(Q_b)>0
    CPI_b(t)=sum(P_b.*(Q_b./(sum(Q_b))));
    else
    CPI_b(t)=sum(P_b.*(ones(1,F_b)./F_b));
    end
    P_l=P(F_b+1:F);
    Q_l=Q(F_b+1:F);
    if sum(Q_l)>0
    CPI_l(t)=sum(P_l.*(Q_l./(sum(Q_l))));
    else
    CPI_l(t)=sum(P_l.*(ones(1,F_l)./F_l));   
    end
    Un(t)=(sum(active)-length(Oc(Oc>0)))/sum(active);

    
    Y_nominal_k(t) = sum(Y_prev_k)*price_k(t);
    Y_nominal_c(t) = sum(Y_prev)*price(t);
    Y_nominal_tot(t)= Y_nominal_k(t)+Y_nominal_c(t);
    Y_real(t) = sum(Y_prev)*price(1) + sum(Y_prev_k)*price_k(1);
    
    gdp_deflator(t) = Y_nominal_tot(t)/ Y_real(t);
     
    Investment(t)=I(t)*price_k(1)+sum(Y_k-inventory_k)*price_k(1)+price(1)*sum(Y); %total investment is investment plus inventory variation
    
    Lev_af(t)=sum((deb+deb_g)./(A+deb+deb_g).*(A/sum(A)));
    Lev_ab(t)=sum((deb(1:F_b)+deb_g(1:F_b))./(A(1:F_b)+deb(1:F_b)+deb_g(1:F_b)).*(A(1:F_b)/sum(A(1:F_b))));
    Lev_al(t)=sum((deb((F_b+1):F_l)+deb_g((F_b+1):F_l))./(A((F_b+1):F_l)+deb((F_b+1):F_l)+deb_g((F_b+1):F_l)).*(A((F_b+1):F_l)/sum(A((F_b+1):F_l))));
    Lev_ak(t)=sum((deb_k+deb_g_k)./(A_k+deb_k+deb_g_k).*(A_k/sum(A_k)));
    debtot=[deb+deb_g,deb_k+deb_g_k];
    atot=[A,A_k];
    Lev_a(t)=sum(debtot./(atot+debtot).*atot/(sum(atot)));
 
    
    defaults_gov=0;
    negcash_k=find(A_k<=0|liquidity_k<0);
    negcash=find(A<=0|liquidity<0);
    
    
    NetEq = liquidity-deb-deb_g;
    Y_prevp=Y_prev(NetEq>0);  
    bankruptcy_rate(t) = (length(negcash_k) + length(negcash))/(F+N);
    
    equity_injection(:)=0;
    equity_injection_k(:)=0;
    meanA=mean(A(A>0));
    meanA_k=mean(A_k(A_k>0));

    
    %update bankrupted firms!
    for i=negcash                                                       %pick sequentially failed firms
        if liquidity(i)<0 && A(i)>0
           liquidity(i)=liquidity(i)+PA(W+i);
           PA(W+i)=0;
            if liquidity(i)<0
               piB=piB+liquidity(i);
               liquidity(i)=0;
            end
           A(i)=liquidity(i)+capital_value(i)-deb(i)-deb_g(i);
        else
            if equity_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
               A(i)=A(i)+PA(W+i);
               liquidity(i)=liquidity(i)+PA(W+i);
               PA(W+i)=0;
               if A(i)<0
               equity_injection(i)=-A(i)+meanA;
               else
               equity_injection(i)=meanA-A(i);
               end
               liquidity(i)=liquidity(i)+equity_injection(i);
               A(i)=A(i)+equity_injection(i);
               gov_ownership(i)=gov_ownership(i)+equity_injection(i);
               bankruptcy_rate(t) = bankruptcy_rate(t)- 1/(F+N);
            else
                defaults(t)=defaults(t)+deb(i);
                zzz=deb(i);                                                     %take residual debts     
                xxx=deb_g(i);
                piB=piB+(liquidity(i)-deb(i)-deb_g(i));                                               %account for bad debts
                loans=loans-zzz;
                loans_g=loans_g-xxx;
                defaults_gov=defaults_gov+xxx;
                A(i)=PA(W+i)+K(i)*price_k(t+1);                                                   %initialize new firm 
                capital_value(i)=K(i)*price_k(t+1);
                PA(W+i)=0;
                liquidity(i)=A(i)-K(i)*price_k(t+1);
                deb(i)=0;
                deb_g(i)=0;
                P(i)=mean(P);
                targetLev=0.2;
                mxY=((A(i)+targetLev*A(i)/(1-targetLev))*1/wb)*alpha;
                Y_prev(i)=min(trimmean(Y_prevp,10),mxY);
                Yd(i)=Y_prev(i);
                x(i)=Y_prev(i)/k/K(i);
                barK(i)=K(i);
                barYK(i)=Y_prev(i)/k;
                Y(i)=0;
                stock(i)=0;
                interest_r(i)=r_f;
                %%fire workers
                workforce=find(Oc==i);                                          %pick all firm i's workers
                fired=workforce;  
                Oc(fired)=0;
                Leff(i)=0;
            end
        end  
    end
    
    NetEq_k = A_k;
    if isempty(find(NetEq_k>0, 1))
          warning('all capital firms are bankrupted')        
    else
        Y_prevp_k=Y_prev_k(NetEq_k>0);  
    end
                  
    initialA_k(t)=sum(PA(W+F+negcash_k))/length(negcash_k);
    for i=negcash_k                                                       %pick sequentially failed firms
        if liquidity_k(i)<0 && A_k(i)>0
           liquidity_k(i)=liquidity_k(i)+PA(W+F+i);
           PA(W+F+i)=0;
            if liquidity_k(i)<0
               piB=piB+liquidity_k(i);
               liquidity_k(i)=0;
            end
           A_k(i)=liquidity_k(i)+Y_k(i)*P_k(i)-deb_k(i)-deb_g_k(i);
           A_control(i)=A_k(i);
        else
            if equity_support==1 && t>=lock_start && lock_start>0 && t<=lock_start+11
                A_k(i)=A_k(i)+PA(W+F+i);
                liquidity_k(i)=liquidity_k(i)+PA(W+F+i);
                PA(W+F+i)=0;
                if A_k(i)<0
                equity_injection_k(i)=-A_k(i)+meanA_k;
                else
                equity_injection_k(i)=meanA_k-A_k(i);    
                end
                liquidity_k(i)=liquidity_k(i)+equity_injection_k(i);
                gov_ownership_k(i)=gov_ownership_k(i)+equity_injection_k(i);
                A_k(i)=liquidity_k(i)+Y_k(i)*price_k(t)-deb_k(i)-deb_g_k(i);
            else
                defaults_k(t)=defaults_k(t)+deb_k(i);
                zzz=deb_k(i);                                                     %take residual debts   
                xxx=deb_g_k(i);
                piB=piB+(liquidity_k(i)-deb_k(i)-deb_g_k(i));                                               %account for bad debts
                loans=loans-zzz;
                loans_g=loans_g-xxx;
                defaults_gov=defaults_gov+xxx;
                A_k(i)=PA(W+F+i);                                                   %initialize new firm 
                A_control(i)=A_k(i);
                PA(W+F+i)=0;
                liquidity_k(i)=A_k(i);
                deb_k(i)=0;
                deb_g_k(i)=0;
                P_k(i)=mean(P_k);
                %maximum initial productin is given by the leverage
                targetLev=0.2;
                mxY=((A_k(i)+targetLev*A_k(i)/(1-targetLev))*1/wb)*alpha;
                Y_prev_k(i)=min(trimmean(Y_prevp_k,10),mxY);
                Y_kd(i)=Y_prev_k(i);
                Y_k(i)=0;
                stock_k(i)=0;
                interest_r_k(i)=r_f;
                workforce_k=find(Oc==F+i);                 %pick all firm i's workers
                fired_k=workforce_k;  
                Oc(fired_k)=0;
                Leff_k(i)=0; 
            end
        end
    end
    
    
    %% bank accounting
    
    piB=piB+defaults_gov+sum(interests)+sum(interests_k) + bond_interest_rate*bonds(t-1)-interest_deposits+dead_assets;                                             %bank profits  

    if piB>0 && totE(t)>0
        TA(t)=TA(t)+piB*tax_rate_b;
        piB=piB*(1-tax_rate_b);
        dividendsB(t)=div_B*piB;
        piB=(1-div_B)*piB;
        PA(W+1:W+F+N)=PA(W+1:W+F+N)+dividendsB(t)/(F+N);
    else 
        dividendsB(t) = 0;
    end    
    E=E+piB;                                                            %update bank capital

    %%add bank's dividends to income of capitalists
    dividends_income_k = dividends_income_k + dividendsB(t)/(F+N);
    dividends_income = dividends_income + dividendsB(t)/(F+N);

    
    profitsB(t)=piB;

    totE(t+1)=E;
  
  
 GB(t)= TA(t)+sum(div_g)+sum(div_g_k) - G(t)-defaults_gov- EXP(t)- bond_interest_rate*bonds(t-1)-Healthcare(t)-sum(labour_transfer)-sum(labour_transfer_k)-sum(liquidity_transfer)-sum(liquidity_transfer_k)-sum(equity_injection)-sum(equity_injection_k); %JAKOB bonds(t-1) mi impalla tutto ma non potevo

 stock_bonds(t) =sum(-GB(1:t));
 bonds(t) = max(0,stock_bonds(t)); 
 
 debGDP(t)=stock_bonds(t)/Y_nominal_tot(t);
 defGDP(t)=GB(t)/Y_nominal_tot(t);
 
money(t) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)-sum(deb_g)-sum(deb_g_k)+E;

if u_target - Un(t)>0
    wb = wb * (1+ wage_update_up * (u_target - Un(t))); 
else
    wb = wb * (1+ wage_update_down * (u_target - Un(t))); 
end

wages_t(t) = wb;

dead_assets=0;

for tt=((t-1)*4+1):((t-1)*4+4)

Hdemand(tt)=sum(health_demanded)+sum(health_demanded2);
    
health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));


if social_distancing==1
if lockdown==1
   dis_cost=dis_cost_lockdown;
end

if lock_ended==1
   dis_cost=min(dis_cost_init,dis_cost+post_lock_adjustment);
end

dis_vec=[sum(infected==1 & detected==1)-dis_threshold; sum(distancing==1)/W-sum(distancing==0)/W; -dis_cost];

if epidemic==1 && t>=648 && sum(infected)>0
for i=1:W
dis_index(i)=distancing_persistence*dis_index(i)+(1-distancing_persistence)*intensity*dis_vec; 
dis_prob(i)=1/(1+exp(-dis_index(i)));
if shocks_distancing(tt,i)<dis_prob(i) || detected(i)==1
   distancing(i)=1;
else
   distancing(i)=0;
end
if distancing(i)==1
   distanced(i)=distanced(i)+1;
end
if distanced(i)==1
   c_shock_b(i)=shock_b;
   c_shock_l(i)=shock_l;
end
end
elseif epidemic==1 && t>=648 && sum(infected)==0
for i=1:W
dis_index(i)=dis_index(i)-post_lock_adjustment;   
dis_prob(i)=1/(1+exp(-dis_index(i)));
if shocks_distancing(tt,i)<dis_prob(i)
   distancing(i)=1;
else
   distancing(i)=0;
end
end
end
end

%%here the epidemic portion of the model starts
%this part only runs if there are any infected in the system (saves time)
if epidemic==1 && t>=648 && sum(infected)>0

%update duration of disease
duration(infected==1)=duration(infected==1)+1;

inf_cum(tt)=inf_cum(tt)+inf_cum(tt-1);
serious_share(tt)=serious_share(tt)+serious_share(tt-1);

if lockdown==0 && lock_ended==0
   reduce_cons=1;
elseif lockdown==1
   reduce_cons=reduce_cons_lockdown;
   detectprob=min(0.1,detectprob+detect_adjustment*(tt-lock_start2));
elseif lock_ended==1
   reduce_cons=min(1,reduce_cons+post_lock_adjustment);
   detectprob=min(0.1,detectprob+detect_adjustment*(tt-lock_start2));
end

dprob(tt)=detectprob;

m_rand=round(reduce_cons*m_randcons);
sd_rand=round(reduce_cons*sd_randcons);

%add random connections to the network
rng(seeds_connections(tt))
new_cons=max(0,round(normrnd(m_rand,sd_rand)));
shocks_connections=rand(new_cons,1);
%new_cons=0;
adds=0;
while adds<new_cons
    nodes=randsample(W,2,true);
    conprob=1;
    if distancing(nodes(1))==1 || distancing(nodes(2))==1
       conprob=0.5;
    end
    %make sure connection does not exist yet and both agents are alive
    if con_mat(nodes(1),nodes(2))==0 && dead(nodes(1))==0 && dead(nodes(2))==0 || con_mat(nodes(2),nodes(1))==0 && dead(nodes(1))==0 && dead(nodes(2))==0
        if shocks_connections(adds+1)<=conprob
        con_mat(nodes(1),nodes(2))=1;  %temporary connections are marked with a "1"
        con_mat(nodes(2),nodes(1))=1;
        end
        adds=adds+1;
    end
end



%take upper triangular part of connection matrix to avoid double-counting connections
con_mat=triu(con_mat,1);

%get all temporary connections (workplace & random)
[row, col]=find(con_mat==1);
r_cons=[row col];
[row, col]=find(con_mat==2);
w_cons=[row col];


%make a list of all connections, permanent and temporary
%during lockdown, only a fraction of permanent connections encounter each
%other
if reduce_cons==1
    fixed_cons=fix_cons;
else
    rng(seeds_consampling(tt))
    fixed_cons=datasample(fix_cons,round(reduce_cons*size(fix_cons,1)),'Replace',false);
end
connections=[r_cons;w_cons; fixed_cons];
con_counter=0;
cons=[];
weights=[];
    for l=1:size(connections,1)
        con=connections(l,:);
        if susceptible(con(1))==1 && infected(con(2))==1 && duration(con(2))<=3 && detected(con(1))==0 && detected(con(2))==0 || susceptible(con(2))==1 && infected(con(1))==1 && duration(con(1))<=3 && detected(con(1))==0 && detected(con(2))==0
            con_counter=con_counter+1;
            cons=[cons;con];
            weights=[weights;con_mat(con(1),con(2))];
        end
    end
    new_infections=round(transmission_rate*con_counter);
    inf_counter=0;
   %iterate over all connections; disease spreads randomly
   if con_counter<new_infections
      new_infections=con_counter;
   end
   tosample=[1:size(cons,1)];
   rng(seeds_consampling2(tt))
   con_sample=datasample(tosample,new_infections,'Replace',false,'Weights',weights);
   shocks_infections=rand(new_infections,1);
   for i = 1:new_infections
       con=cons(con_sample(i),:);
       %here already assumed that detected cases isolate and hence do not
       %have contacts
       if susceptible(con(1))==1 && infected(con(2))==1 && detected(con(1))==0 && detected(con(2))==0 || susceptible(con(2))==1 && infected(con(1))==1 && detected(con(1))==0 && detected(con(2))==0 
              if shocks_infections(i)<(1-distancing_effect*distancing(con(1))-distancing_effect*distancing(con(2)))
              if susceptible(con(1))==1 && infected(con(2))==1 && serious(con(1))==1
                  serious_share(tt)=serious_share(tt)+1;
              end
              if susceptible(con(2))==1 && infected(con(1))==1 && serious(con(2))==1
                  serious_share(tt)=serious_share(tt)+1;
              end
              infected(con(1))=1;
              infected(con(2))=1;
              inf_cum(tt)=inf_cum(tt)+1;
              inf_new(tt)=inf_new(tt)+1;
              susceptible(con(1))=0;
              susceptible(con(2))=0;
              inf_counter=inf_counter+1;
              end
       end
   end
   rng(seeds_hdem(tt))
   for i=randperm(W)
       %if currently infected
       if infected(i)==1 && duration(i)<=max_dur(i)
          %random detection of disease
          if shocks_detection(tt,i)<detectprob && detected(i)==0 
             detected(i)=1;
             detect_new(tt)=detect_new(tt)+1;
             %detected cases become inactive; firm loses worker
             if active(i)==1
             active(i)=0;
             sickpay(i)=1;  %only previously active agents receive sick pay; others just continue to receive pension
             end
             Firm=Oc(i);
               if Firm>0
                  if Firm<=F
                     Leff(Firm)=Leff(Firm)-1;
                  end
                  if Firm>F
                     Leff_k(Firm-F)=Leff_k(Firm-F)-1;
                  end
               end
             Oc(i)=0;
          end
          %all serious cases are detected & leave work
          if serious(i)==1 && detected(i)==0
              detected(i)=1;
              detect_new(tt)=detect_new(tt)+1;
              if active(i)==1
              active(i)=0;
              sickpay(i)=1;  %only previously active agents receive sick pay; others just continue to receive pension
              end
              Firm=Oc(i);
               if Firm>0
                  if Firm<=F
                     Leff(Firm)=Leff(Firm)-1;
                  end
                  if Firm>F
                     Leff_k(Firm-F)=Leff_k(Firm-F)-1;
                  end
               end
              Oc(i)=0;
          end
          if serious(i)==1 && dead(i)==0
               if health_demanded(i)==0
                  health_demanded(i)=h_1*age(i)+shocks_healthcare(tt,i)*h_2;
               end
               if health_demanded(i)>0 && health_supplied(i)<health_demanded(i)
                  hdem=health_demanded(i)-health_supplied(i);
                  health_supplied(i)=health_supplied(i)+min(health_supply,hdem);
                  health_supply=max(0,health_supply-hdem);
               end
               if shocks_death(tt,i)<dieprob*age(i)^2+h_3*(health_demanded(i)-health_supplied(i))
                   infected(i)=0;
                   dead(i)=1;
                   if detected(i)==0
                   detected(i)=1;
                   detect_new(tt)=detect_new(tt)+1;
                   end
                   Firm=Oc(i);
                   if Firm>0
                      if Firm<=F
                         Leff(Firm)=Leff(Firm)-1;
                      end
                      if Firm>F
                         Leff_k(Firm-F)=Leff_k(Firm-F)-1;
                      end
                   end
                   Oc(i)=-1;
                   health_demanded(i)=0;
                   health_supplied(i)=0;
                   health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));
                   dead_assets=dead_assets+PA(i);
                   PA(i)=0;
                   active(i)=0;
                   sickpay(i)=0;
               end             
          end
       end
      if infected(i)==1 && duration(i)>max_dur(i) && dead(i)==0
      infected(i)=0;
      recovered(i)=1;
      health_demanded(i)=0;
      health_supplied(i)=0;
      health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));
      if age(i)<3
         active(i)=1;
      end
      sickpay(i)=0;
      end
   end

   %reset adjacency matrix for next period, deleting temporary ones
   con_mat=triu(con_mat,1);
   con_mat=triu(con_mat)+triu(con_mat,1)';
   con_mat(con_mat==1)=0;
   if tt==((t-1)*4+4)
   con_mat(con_mat==2)=0;    
   end
else
    health_demanded(:)=0;
    health_supplied(:)=0;
end

%calculate aggregate stats
sus_share(tt)=sum(susceptible)/W;
inf_share(tt)=sum(infected)/W;
inf_detected(tt)=sum(detected==1);
ddetect(tt)=inf_detected(tt)-inf_detected(tt-1);
rec_share(tt)=sum(recovered)/W;
dead_share(tt)=sum(dead)/W;

health_supply=max(0,healthcare(t)-sum(health_supplied)-sum(health_supplied2));

%%"normal" disease
%increment duration for infected
duration2(infected2==1)=duration2(infected2==1)+1;
%agents who have exceeded the max. duration of the disease recover
infected2(duration2>=max_dur2)=0;
recovered2(duration2>=max_dur2)=1;
duration2(recovered2==1)=0;
health_demanded2(recovered2==1)=0;
health_supplied2(recovered2==1)=0;
%recovered agents become active again
active(recovered2==1 & age<3)=1;
sickpay(recovered2==1)=0;
%find susceptible agents
sus=find(susceptible2==1 & dead==0);
%randomly infect some susceptible agents
newinf2=round(infprob2*sum(dead==0));
rng(seeds_infection2(tt))
infect=randsample(sus,newinf2,false);
infected2(infect)=1;
susceptible2(infect)=0;
sick=find(infected2==1);
rng(seeds_hdem2(tt))
svec=randperm(length(sick));
for si=svec
       i=sick(si);
       if health_demanded2(i)==0
          health_demanded2(i)=h_1*age(i)+shocks_healthcare2(tt,i)*h_2;
       end
       if health_demanded2(i)>0 && health_supplied2(i)<health_demanded2(i)
          hdem=health_demanded2(i)-health_supplied2(i);
          health_supplied2(i)=health_supplied2(i)+min(health_supply,hdem);
          health_supply=max(0,health_supply-hdem);
       end
       if active(i)==1
       active(i)=0;
       sickpay(i)=1;  %only previously active agents receive sick pay; others just continue to receive pension
       Firm=Oc(i);
        if Firm>0
          if Firm<=F
             Leff(Firm)=Leff(Firm)-1;
          end
          if Firm>F
             Leff_k(Firm-F)=Leff_k(Firm-F)-1;
          end
        end
       end
       Oc(i)=0;       
end

%a share of recovered become susceptible again
rec=find(recovered2==1 & dead==0);
rng(seeds_susceptible(tt))
new_sus=rand(1,length(rec));
susceptible2(rec(new_sus<susprob))=1;
recovered2(susceptible2==1)=0;

inf_share2(tt)=sum(infected2)/W;
sus_share2(tt)=sum(susceptible2)/W;
rec_share2(tt)=sum(recovered2)/W;
distance(tt)=sum(distancing)/W;

if sum(detected)>=lockdown_threshold && lockdown_exp==1 && lockdown==0 
    lockdown=1;
    active_f(deactivate_f)=0;   %if instead want to deactivate only a subset of L-firms (work from home for rest)
    homeoffice(F_b+1:F)=1;
    homeoffice(deactivate_f)=0;
    constraint_k(:)=lock_constraint;          %if instead we want K-firms to reduce production!
    constraint_f(homeoffice==1)=lock_constraint;
    lock_start=t;
    lock_start2=tt;
    lock_ended=0;
end

end

av_inf_new(t)=mean(detect_new(tt-3:tt));
labforce1(t)=sum(active==1)/sum(dead==0);
labforce2(t)=sum(active==1)/sum(age<3 & dead==0);

 end
et=toc;
