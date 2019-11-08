%% Description 
%This script is to for a three-state jump system with states {-1,0,1} and tests its 
% behaviour for three scenarios : detailed balance observed, detailed balance broken 
% but can be restored in the presence of a rotating frame and detailed
% balance broken but due to history dependence of the particle's transition
% rates cannot be restored. 

%% Setting random number generator

rng(0);

%% numerical
epsi = 0; %energy scale for correlation functions
eta = 4; %energy scale to evoke response functions
time = 19999; h = 0.01; % # of time steps and scale
num_part = 50; % particle ensemble
alpha = zeros(time-1,num_part); % particle's transition rate values
alpha = [ones(1,num_part);alpha]; 
rate = 0.1; lambda = 0.1; %strength and period of history dependence respectively

%probabilities of being in one of the three states {1: -1, 2: 0, 3: 1}
energies_1_cum = zeros(time+1,3); energies_2_cum = zeros(time+1,3); energies_3_cum = zeros(time+1,3);

%probabilities for the vanilla case of when all states are equivalent
energies_0 = [exp(epsi),exp(0),exp(-epsi)]/(exp(epsi) + exp(0) + exp(-epsi));

state = [ones(1,num_part);zeros(time,num_part)]; %states of 50 particles at all time points
traj_prob = [ones(1,num_part);zeros(time,num_part)]; %probability distribution of all trajectories

for i = 1:num_part
     probs = rand(time,1);
     for j = 2:time+1
        if (mod(j+2,1000) == 0 || mod(j+1,1000) == 0) %force ping for response function
            eta = 0;
        else
            eta = 0;
        end
        energies_1_cum(j,:) = [1-h*(1 + alpha(j-1,i))*exp(eta),1-h*(alpha(j-1,i) + alpha_delta)*exp(eta),1];
        energies_2_cum(j,:) = [exp(epsi-eta)*(alpha(j-1,i) + alpha_delta)*h,1-exp(epsi)*(1 - alpha_delta)*h,1];
        energies_3_cum(j,:) = [exp(epsi-eta)*(1 - alpha_delta)*h/exp(-epsi), exp(epsi-eta)*(1 - alpha_delta)*h/exp(-epsi) + exp(epsi)*h*(alpha(j-1,i) + alpha_delta)/exp(-epsi),1];
        
        energies_1 = [1-h*(1 + alpha(j-1,i))*exp(eta), h*(1 - alpha_delta)*exp(eta), h*(alpha(j-1,i) + alpha_delta)*exp(eta)];
        energies_2 = [exp(epsi-eta)*(alpha(j-1,i) + alpha_delta)*h, 1-(exp(epsi)*(1 - alpha_delta)*h)-(exp(epsi-eta)*(alpha(j-1,i) + alpha_delta)*h), exp(epsi)*(1 - alpha_delta)*h];
        energies_3 = [exp(epsi-eta)*(1 - alpha_delta)*h/exp(-epsi),exp(epsi)*h*alpha(j-1,i)/exp(-epsi), 1 - exp(epsi-eta)*(1 - alpha_delta)*h/exp(-epsi) - exp(epsi)*h*alpha(j-1,i)/exp(-epsi)];
        
        %particle transition depends on its previous state in the broken
        %detailed balance cases
        if state(j-1,i) == -1
            state(j,i) = find(probs(j-1)<energies_1_cum(j,:),1)-2;
            traj_prob(j,i) = traj_prob(j-1,i) * energies_1(find(probs(j-1)<energies_1_cum(j,:),1));
        elseif state(j-1,i) == 0
            state(j,i) = find(probs(j-1)<energies_2_cum(j,:),1)-2;
            traj_prob(j,i) = traj_prob(j-1,i) * energies_2(find(probs(j-1)<energies_2_cum(j,:),1));            
        else
            state(j,i) = find(probs(j-1)<energies_3_cum(j,:),1)-2;
            traj_prob(j,i) = traj_prob(j-1,i) * energies_3(find(probs(j-1)<energies_3_cum(j,:),1));            
        end
        part_state = state(1:j,i);
        %part_state(part_state == -1) = pi; part_state(part_state == 0) = pi/3; part_state(part_state == 1) = -pi/3;
        alpha(j,i) = alpha(1,i) + sum(rate*exp((-j+1:0)*lambda)'.*(part_state));
     end    
end

%Ensemble probabilities of the particles in the three states and over all
%time steps
prob_particle = zeros(num_part,3);
for pp = 1:num_part
    prob_particle(pp,1) = sum(state(:,pp) == -1)/(time+1);
    prob_particle(pp,2) = sum(state(:,pp) == 0)/(time+1);
    prob_particle(pp,3) = sum(state(:,pp) == 1)/(time+1);
end

prob_time = zeros(time+1,3);
for pt = 1:time+1
    prob_time(pt,1) = sum(state(pt,:) == -1)/num_part;
    prob_time(pt,2) = sum(state(pt,:) == 0)/num_part;
    prob_time(pt,3) = sum(state(pt,:) == 1)/num_part;
end

%% probabilities numerical
figure; hold on;
plot(prob_time(1:time+1,:));
xlabel('Time');
ylabel('Probabilities');

%% Calculation of state response function by averaging over all the force pings

resp_avg = zeros(1000,3);
for t = 0:18
    resp_avg = resp_avg + prob_time(t*1000+700:(t+1)*1000+699,:);
end    
resp_avg = resp_avg/(t+1);
figure; plot(resp_avg);
resp_avg_angle = sum(resp_avg .*[-1,0,1],2);
figure; plot(resp_avg_angle); %- 1/3*(pi+pi/3-pi/3));

%% Scalar autocorrelation function of the average state
%state(state == -1) = pi; state(state == 0) = pi/3; state(state == 1) = -pi/3;
avg_m = mean(mean(state));
avg_corr = zeros(time+1,1);

for ct = 0:19999
    num_corr_ind = 1:floor((time - ct)/1000) + 1;
    corr_ind = (num_corr_ind-1)*1000+1;
    state_0 = state(corr_ind,:);
    state_ct = state(corr_ind+ct,:);
    state_corr = state_0.*state_ct;
    avg_corr(ct+1) = mean(mean(state_corr));
end

figure;
plot(avg_corr(1:1000) - avg_m^2);
corr_fn = zeros(2*length(avg_corr)-1,1);
corr_fn(time+1:end,1) = avg_corr - avg_m^2;
corr_fn(1:time,1) = flipud(avg_corr(2:end) - avg_m^2);
%}
%% Tensor crosscorrelation functions of the three states
C_00 = zeros(time*2+1,1); C_01 = zeros(time*2+1,1); C_0neg1 = zeros(time*2+1,1);
C_10 = zeros(time*2+1,1); C_11 = zeros(time*2+1,1); C_1neg1 = zeros(time*2+1,1);
C_neg10 = zeros(time*2+1,1); C_neg11 = zeros(time*2+1,1); C_neg1neg1 = zeros(time*2+1,1);

for n = 1:num_part
    sigma_0 = state(:,n) == 0;
    sigma_1 = state(:,n) == 1;
    sigma_neg1 = state(:,n) == -1;
    C_00 = C_00 + xcorr(sigma_0, sigma_0);
    C_01 = C_01 + xcorr(sigma_0, sigma_1);
    C_0neg1 = C_0neg1 + xcorr(sigma_0, sigma_neg1);
    C_10 = C_10 + xcorr(sigma_1, sigma_0);
    C_11 = C_11 + xcorr(sigma_1, sigma_1);
    C_1neg1 = C_1neg1 + xcorr(sigma_1, sigma_neg1);
    C_neg10 = C_neg10 + xcorr(sigma_neg1, sigma_0);
    C_neg11 = C_neg11 + xcorr(sigma_neg1, sigma_1);
    C_neg1neg1 = C_neg1neg1 + xcorr(sigma_neg1, sigma_neg1);
end


%% State correlations to check if leaving 0 is equivalent to entering it
% This block and the next are sanity checks to understand the particle's trajectory
% better.

cnt_01 = zeros(num_part,1);
cnt_0neg1 = zeros(num_part,1);
cnt_00_leave = zeros(num_part,1);

for i = 1:num_part
    for j = 1:time
        if state(j,i) == 0 
            if state(j+1,i) == 1
                cnt_01(i) = cnt_01(i) + 1;
            elseif state(j+1,i) == -1
                cnt_0neg1(i) = cnt_0neg1(i) + 1;
                else
                cnt_00_leave(i) = cnt_00_leave(i) + 1;
            end
        end
    end
end

%% State correlations to check if entering 0 is equivalent to leaving it

cnt_10 = zeros(num_part,1);
cnt_1neg1 = zeros(num_part,1);
cnt_neg10 = zeros(num_part,1);
cnt_00 = zeros(num_part,1);

for i = 1:num_part
    for j = 2:time+1
        if state(j,i) == 0 
            if state(j-1,i) == 1
                cnt_10(i) = cnt_10(i) + 1;
            elseif state(j-1,i) == -1
                cnt_neg10(i) = cnt_neg10(i) + 1;
            else
                cnt_00(i) = cnt_00(i) + 1;
            end
        end
    end
end

for i = 1:num_part
    for j = 1:time
        if state(j,i) == 1
            if state(j+1,i) == -1
                cnt_1neg1(i) = cnt_1neg1(i) + 1;
            end
        end
    end
end

%% Analytical calculation of the above
%{
Pa = zeros(500,1); Pb = zeros(500,1); Pc = zeros(500,1);
Pa(1) = 1;
for i = 1:time
    if (i<1000 && i>997) || (i<2000 && i>1997) || (i<3000 && i>2997) || (i<4000 && i>3997) ...
        || (i<5000 && i>4997) || (i<6000 && i>5997) || (i<7000 && i>6997) || (i<8000 && i>7997) ...
        || (i<9000 && i>8997) || (i<10000 && i>9997) || (i<11000 && i>10997) || (i<12000 && i>11997) ...
        || (i<13000 && i>12997) || (i<14000 && i>13997) || (i<15000 && i>14997) || (i<16000 && i>15997)...
        || (i<17000 && i>16997) || (i<18000 && i>17997) || (i<19000 && i>18997) || (i<20000 && i>19997)
        eta = 0;
    else
        eta = 0;
    end
    rba = 1; rcb = exp(epsi+eta); rac = exp(epsi+eta)/exp(-epsi); % a = -1, b = 0, c = 1
    rab = alpha(i)*exp(epsi+eta); rbc = alpha(i)*exp(epsi+eta)/exp(-epsi); rca = alpha(i);
    Pa(i+1) = Pa(i) + h*(-rba*Pa(i) -rca*Pa(i) + rab*Pb(i) + rac*Pc(i));
    Pb(i+1) = Pb(i) + h*(-rab*Pb(i) -rcb*Pb(i) + rba*Pa(i) + rbc*Pc(i));
    Pc(i+1) = Pc(i) + h*(-rac*Pc(i) -rbc*Pc(i) + rca*Pa(i) + rcb*Pb(i));
end
time_analy = h*(1:5000);
   %}

%% Probabilities analytical
%figure;
%hold on;
%plot(Pa,'black');
%plot(Pb,'black');
%plot(Pc,'black');
%xlabel('Time(s)');
%ylabel('Probabilities');

