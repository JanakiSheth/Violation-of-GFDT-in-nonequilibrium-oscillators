%% Description 
% This script is to for a three-state jump system with states {-1,0,1} and tests its 
% behaviour for three scenarios : detailed balance obeyed, detailed balance broken 
% but can be restored in the presence of a rotating frame and detailed
% balance broken but due to history dependence of the particle's transition
% rates. The second three-state system breaks fdt but obeys gfdt, the third
% system violates both theorems.
% To reproduce results in https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.023150 paper
% compare the state correlation and response functions. The tuple
% correlation and response functions are for current work to understand
% restoration of gfdt in history-dependent systems. 


for random_num = 8
%% Setting random number generator

rng(random_num);

%% numerical
time = 199999; h = 0.01; % # of time steps and scale
num_part = 100; % particle ensemble
p = 2;%2; % particle's transition rate values
vel_array = zeros(time,num_part);
ping = 1;

%probabilities of being in one of the six tuples {1: {-1,1}, 2: {-1,0}, 3: {0,1}, 4: {0,0}, 5: {1,1}, 6: {1,0}}
energies_1_cum = zeros(time+1,3); 
energies_2_cum = zeros(time+1,3); 
energies_3_cum = zeros(time+1,3);

state = [ones(1,num_part)*-1; zeros(time,num_part)]; %states of all particles at all time points
drive = zeros(time,num_part);
interimState = zeros(time,num_part);
finalState = zeros(time,num_part);

for i = 1:num_part
     probs = rand(time,1);
     vel = 1;
     
     for j = 2:time+1       
        vel_array(j-1,i) = vel;
        if (mod(j,3000) == 0) 
            %Force ping for response function, eta \neq 0. For correlation
            %functions, eta = 0
            eta = 0;
        else
            eta = 0;
        end
        energies_1_cum(j,:) = [1 - h*(p+p)*ping*exp(eta), 1 - h*p*ping*exp(eta), 1];
        energies_2_cum(j,:) = [exp(-eta)*ping*h*p, 1 - h*p, 1];
        energies_3_cum(j,:) = [exp(-eta)*ping*h*p, exp(-eta)*ping*h*p + h*p, 1];        

        if vel == 0
            interimState(j-1,i) = state(j-1,i);
        else
            if vel == 1
                if state(j-1,i) == -1
                    interimState(j-1,i) = 0;
                elseif state(j-1,i) == 0
                    interimState(j-1,i) = 1;
                else
                    interimState(j-1,i) = -1;
                end
            elseif vel == 2
                if state(j-1,i) == -1
                    interimState(j-1,i) = 1;
                elseif state(j-1,i) == 0
                    interimState(j-1,i) = -1;
                else
                    interimState(j-1,i) = 0;
                end
            end
        end
        
        if interimState(j-1,i) == -1 
            finalState(j-1,i) = find(probs(j-1)<energies_1_cum(j,:),1)-2;
            state(j,i) = finalState(j-1,i);
        elseif interimState(j-1,i) == 0 
            finalState(j-1,i) = find(probs(j-1)<energies_2_cum(j,:),1)-2;
            state(j,i) = finalState(j-1,i);   
        elseif interimState(j-1,i) == 1 
            finalState(j-1,i) = find(probs(j-1)<energies_3_cum(j,:),1)-2;
            state(j,i) = finalState(j-1,i);  
        end      
        
        if finalState(j-1,i) == -1
            drive(j-1,i) = 1;
        elseif finalState(j-1,i) == 0
            drive(j-1,i) = 2;
        else
            drive(j-1,i) = 1;
        end
                        
        part_state = state(j,i);
        vel = drive(j-1,i);    
        
    end    
end

%% Ensemble probabilities of the particles in the three states and over all time steps

prob_time = zeros(time+1,3);
prob_time(:,1) = sum(state == -1,2)/num_part;
prob_time(:,2) = sum(state == 0,2)/num_part;
prob_time(:,3) = sum(state == 1,2)/num_part;

finalTuple_time = zeros(time,9);
% has template {final state, drive}; remember that both of these variables are from time (j-1).

for pt = 1:time
    finalTuple_time(pt,1) = sum(finalState(pt,:) == -1 & drive(pt,:) == 1)/num_part;
    finalTuple_time(pt,2) = sum(finalState(pt,:) == -1 & drive(pt,:) == 0)/num_part;
    finalTuple_time(pt,3) = sum(finalState(pt,:) == -1 & drive(pt,:) == 2)/num_part;
    finalTuple_time(pt,4) = sum(finalState(pt,:) == 0 & drive(pt,:) == 1)/num_part;
    finalTuple_time(pt,5) = sum(finalState(pt,:) == 0 & drive(pt,:) == 0)/num_part;
    finalTuple_time(pt,6) = sum(finalState(pt,:) == 0 & drive(pt,:) == 2)/num_part;
    finalTuple_time(pt,7) = sum(finalState(pt,:) == 1 & drive(pt,:) == 1)/num_part;
    finalTuple_time(pt,8) = sum(finalState(pt,:) == 1 & drive(pt,:) == 0)/num_part;
    finalTuple_time(pt,9) = sum(finalState(pt,:) == 1 & drive(pt,:) == 2)/num_part;
end

%% plotting state and tuple probabilities 
figure; hold on;
plot(prob_time(1:time+1,:));
colormap winter
xlabel('Time');
ylabel('State Probabilities');

figure; hold on;
plot(finalTuple_time(1:time,:));
xlabel('Time');
ylabel('Final tuple Probabilities');
%% Calculation of state response function by averaging over all the force pings

resp_avg = zeros(3000,65);
for t = 0:64
    resp_avg(:,t+1) = prob_time(t*3000+701:(t+1)*3000+700,3); % so response at tau=0 is effectively prob_time(2299).
end    
figure; plot(mean(resp_avg,2));

%% Probabilities matched to each state of trajectory

prob_state = zeros(size(state));
prob_tuple_state = zeros(size(finalState));
for i = 1:num_part
    for j = 1:time+1
        prob_state(j,i) = prob_time(j,state(j,i)+2);   
    end
end

for i = 1:num_part
    for j = 1:time
        if finalState(j,i)==-1 && drive(j,i)==1
            prob_tuple_state(j,i) = finalTuple_time(j,1);
        elseif finalState(j,i)==-1 && drive(j,i)==0
            prob_tuple_state(j,i) = finalTuple_time(j,2);
        elseif finalState(j,i)==-1 && drive(j,i)==2
            prob_tuple_state(j,i) = finalTuple_time(j,3);
        elseif finalState(j,i)==0 && drive(j,i)==1
            prob_tuple_state(j,i) = finalTuple_time(j,4);
        elseif finalState(j,i)==0 && drive(j,i)==0
            prob_tuple_state(j,i) = finalTuple_time(j,5);
        elseif finalState(j,i)==0 && drive(j,i)==2
            prob_tuple_state(j,i) = finalTuple_time(j,6);    
        elseif finalState(j,i)==1 && drive(j,i)==1
            prob_tuple_state(j,i) = finalTuple_time(j,7);
        elseif finalState(j,i)==1 && drive(j,i)==0
            prob_tuple_state(j,i) = finalTuple_time(j,8);
        elseif finalState(j,i)==1 && drive(j,i)==2
            prob_tuple_state(j,i) = finalTuple_time(j,9);    
        end
    end
end

%% Calculation of tuple response function by averaging over all the force pings

resp_finalTuple_1 = zeros(3000,65); resp_finalTuple_2 = zeros(3000,65); resp_finalTuple_3 = zeros(3000,65); 
resp_finalTuple_4 = zeros(3000,65); resp_finalTuple_5 = zeros(3000,65); resp_finalTuple_6 = zeros(3000,65);
resp_finalTuple_7 = zeros(3000,65); resp_finalTuple_8 = zeros(3000,65); resp_finalTuple_9 = zeros(3000,65);

for t = 0:64
    resp_finalTuple_1(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,1);
    resp_finalTuple_2(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,2);
    resp_finalTuple_3(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,3);
    resp_finalTuple_4(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,4);
    resp_finalTuple_5(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,5);
    resp_finalTuple_6(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,6);
    resp_finalTuple_7(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,7);
    resp_finalTuple_8(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,8);
    resp_finalTuple_9(:,t+1) = finalTuple_time(t*3000+701:(t+1)*3000+700,9);
end    

figure; hold on;
plot(mean(resp_finalTuple_1,2));
plot(mean(resp_finalTuple_2,2));
plot(mean(resp_finalTuple_3,2));
plot(mean(resp_finalTuple_4,2));
plot(mean(resp_finalTuple_5,2));
plot(mean(resp_finalTuple_6,2));
plot(mean(resp_finalTuple_7,2));
plot(mean(resp_finalTuple_8,2));
plot(mean(resp_finalTuple_9,2));

%% Tensor crosscorrelation functions of the three states

C_00 = zeros(time*2+1,num_part); C_01 = zeros(time*2+1,num_part); C_0neg1 = zeros(time*2+1,num_part);
C_10 = zeros(time*2+1,num_part); C_11 = zeros(time*2+1,num_part); C_1neg1 = zeros(time*2+1,num_part);
C_neg10 = zeros(time*2+1,num_part); C_neg11 = zeros(time*2+1,num_part); C_neg1neg1 = zeros(time*2+1,num_part);

for n = 1:num_part
    sigma_0 = state(:,n) == 0;
    sigma_1 = state(:,n) == 1;
    sigma_neg1 = state(:,n) == -1;
    C_00(:,n) = xcorr(sigma_0 - mean(sigma_0), sigma_0 - mean(sigma_0), 'coeff');
    C_01(:,n) = xcorr(sigma_0 - mean(sigma_0), sigma_1 - mean(sigma_1), 'coeff');
    C_0neg1(:,n) = xcorr(sigma_0 - mean(sigma_0), sigma_neg1 - mean(sigma_neg1), 'coeff');
    C_10(:,n) = xcorr(sigma_1 - mean(sigma_1), sigma_0 - mean(sigma_0), 'coeff');
    C_11(:,n) = xcorr(sigma_1 - mean(sigma_1), sigma_1 - mean(sigma_1), 'coeff');
    C_1neg1(:,n) = xcorr(sigma_1 - mean(sigma_1), sigma_neg1 - mean(sigma_neg1), 'coeff');
    C_neg10(:,n) = xcorr(sigma_neg1 - mean(sigma_neg1), sigma_0 - mean(sigma_0), 'coeff');
    C_neg11(:,n) = xcorr(sigma_neg1 - mean(sigma_neg1), sigma_1 - mean(sigma_1), 'coeff');
    C_neg1neg1(:,n) = xcorr(sigma_neg1 - mean(sigma_neg1), sigma_neg1 - mean(sigma_neg1), 'coeff');
end

%% Correlations of B agarwal and indicator functions

% To remember : the entropy correlation functions must be taken with
% respect to the h = 0 sigma indicator functions. This was a mistake I made
% earlier. It is not <A \dot{s}>|h = 0.1 - <A \dot{s}>|h=0 , rather it is
% <A (\dot{s}|h=0.1 - \dot{s}|h=0) >.

%C_ba0 = zeros(time*2-1,num_part); C_ba1 = zeros(time*2-1,num_part); C_baneg1 = zeros(time*2-1,num_part);
%C_0ba = zeros(time*2-1,num_part); C_1ba = zeros(time*2-1,num_part); C_neg1ba = zeros(time*2-1,num_part);

%C_be0 = zeros(time*2-1,num_part); C_be1 = zeros(time*2-1,num_part); C_beneg1 = zeros(time*2-1,num_part);
%C_0be = zeros(time*2-1,num_part); C_1be = zeros(time*2-1,num_part); C_neg1be = zeros(time*2-1,num_part);

C_neg10beTuple = zeros(time*2-3,num_part); C_neg1neg1beTuple = zeros(time*2-3,num_part); C_00beTuple = zeros(time*2-3,num_part);
C_0neg1beTuple = zeros(time*2-3,num_part); C_10beTuple = zeros(time*2-3,num_part); C_1neg1beTuple = zeros(time*2-3,num_part);

C_neg10bmem = zeros(time*2-1,num_part); C_neg11bmem = zeros(time*2-1,num_part); C_neg12bmem = zeros(time*2-1,num_part); 
C_00bmem = zeros(time*2-1,num_part); C_01bmem = zeros(time*2-1,num_part); C_02bmem = zeros(time*2-1,num_part); 
C_10bmem = zeros(time*2-1,num_part); C_11bmem = zeros(time*2-1,num_part);C_12bmem = zeros(time*2-1,num_part); 

prob_ss = mean(prob_time);
prob_tuple_ss = mean(finalTuple_time);

vel_0 = [-2,1,1,-1,1,0,-1,0,1]*p*h;
vel_1 = [-1,1,0,-1,0,1,-2,1,1]*p*h; % (markovian p (clockwise) and p1 (counterclockwise) with drift = 1)
% vel_1 = [-p*h,p*h,0,-p1*h,0,p1*h,-(p+p1)*h,p*h,p1*h]; % introducing non-markovian nature for p and p1 with drift 1
vel_2 = [-1,0,1,-2,1,1,-1,1,0]*p*h;

% figure out a way to make the following arrays generic
vel_0_t = [-2,0,0,1,0,0,1,0,0, -1,0,0,1,0,0,0,0,0, -1,0,0,0,0,0,1,0,0]*p*h; %velocity from tuple {state,velocity}
vel_1_t = [-1,0,1,0,0,0, 0,-1,0,0,0,1, -2,0,1,0,1,0]*p*h;
vel_2_t = [-1,0,0,0,1,0,0,-2,0,1,0,1,-1,0,1,0,0,0]*p*h;

for n = 1:num_part
    sigma_0 = state(2:end,n) == 0;
    sigma_1 = state(2:end,n) == 1;
    sigma_neg1 = state(2:end,n) == -1;
    
    sigma_01 = (finalState(:,n) == 0 & drive(:,n) == 1);
    sigma_00 = (finalState(:,n) == 0 & drive(:,n) == 0);
    sigma_02 = (finalState(:,n) == 0 & drive(:,n) == 2);
    sigma_neg11 = (finalState(:,n) == -1 & drive(:,n) == 1);
    sigma_neg10 = (finalState(:,n) == -1 & drive(:,n) == 0); 
    sigma_neg12 = (finalState(:,n) == -1 & drive(:,n) == 2);
    sigma_11 = (finalState(:,n) == 1 & drive(:,n) == 1); 
    sigma_10 = (finalState(:,n) == 1 & drive(:,n) == 0);
    sigma_12 = (finalState(:,n) == 1 & drive(:,n) == 2); 
    
    bA = zeros(time,1); %bMem = zeros(time,1);
    for m = 1:time
        if vel_array(m,n) == 0
            omegas = vel_0;
            wt = vel_0_t;
        elseif vel_array(m,n) == 1
            omegas = vel_1;
            wt = vel_1_t;
        else
            omegas = vel_2;
            wt = vel_2_t;
        end  
    end
    
    
    beTuple = (prob_tuple_state(2:end,n)-prob_tuple_state(1:end-1,n))./(h*prob_tuple_state(2:end,n));
    
    % The relevant bmem to be computed depends on the drive in the systems.
    % Comment out the others. 
    % for drive 0 while in any of the states
    
    bMem = sigma_neg10*(prob_tuple_ss(1)*0/prob_tuple_ss(2) + ...
            prob_tuple_ss(3)*0/prob_tuple_ss(2) + prob_tuple_ss(4)*0/prob_tuple_ss(2) + ...
            prob_tuple_ss(5)*(-p*h)/prob_tuple_ss(2) + prob_tuple_ss(6)*0/prob_tuple_ss(2) + ...
            prob_tuple_ss(7)*0/prob_tuple_ss(2) + prob_tuple_ss(8)*(-p*h)/prob_tuple_ss(2) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(2) - (0 + 0 + 0 + p*h + 0 + 0 + p*h + 0)) + ...
            sigma_00*(prob_tuple_ss(1)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(2)*(p*h)/prob_tuple_ss(5) + prob_tuple_ss(3)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(5) + prob_tuple_ss(6)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(7)*0/prob_tuple_ss(5) + prob_tuple_ss(8)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(5) - (0 -p*h + 0 + 0 + 0 + 0 + 0 + 0)) + ...
            sigma_10*(prob_tuple_ss(1)*0/prob_tuple_ss(8) + ...
            prob_tuple_ss(2)*(p*h)/prob_tuple_ss(8) + prob_tuple_ss(3)*0/prob_tuple_ss(8) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(8) + prob_tuple_ss(5)*0/prob_tuple_ss(8) + ...
            prob_tuple_ss(6)*0/prob_tuple_ss(8) + prob_tuple_ss(7)*0/prob_tuple_ss(8) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(8) - (0 -p*h + 0 + 0 + 0 + 0 + 0 + 0));    
    
    % for drive 1 while in any of the states
    bMem = sigma_neg11*(prob_tuple_ss(2)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(3)*0/prob_tuple_ss(1) + prob_tuple_ss(4)*(-p*h)/prob_tuple_ss(1) + ...
            prob_tuple_ss(5)*0/prob_tuple_ss(1) + prob_tuple_ss(6)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(7)*(-2*p*h)/prob_tuple_ss(1) + prob_tuple_ss(8)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(1) - (0 + 0 + p*h + 0 + 0 + 0 + 0 + 0)) + ...
            sigma_01*(prob_tuple_ss(1)*(p*h)/prob_tuple_ss(4) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(4) + prob_tuple_ss(3)*0/prob_tuple_ss(4) + ...
            prob_tuple_ss(5)*0/prob_tuple_ss(4) + prob_tuple_ss(6)*0/prob_tuple_ss(4) + ...
            prob_tuple_ss(7)*(p*h)/prob_tuple_ss(4) + prob_tuple_ss(8)*0/prob_tuple_ss(4) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(4) - (-p*h + 0 + 0 + 0 + 0 + p*h + 0 + 0)) + ...
            sigma_11*(prob_tuple_ss(1)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(7) + prob_tuple_ss(3)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(4)*(p*h)/prob_tuple_ss(7) + prob_tuple_ss(5)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(6)*0/prob_tuple_ss(7) + prob_tuple_ss(8)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(7) - (-2*p*h + 0 + 0 + p*h + 0 + 0 + 0 + 0)); 
    
    % for drive 2 while in any of the states
    bMem = sigma_neg12*(prob_tuple_ss(1)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(3) + prob_tuple_ss(4)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(5)*0/prob_tuple_ss(3) + prob_tuple_ss(6)*(-2*p*h)/prob_tuple_ss(3) + ...
            prob_tuple_ss(7)*0/prob_tuple_ss(3) + prob_tuple_ss(8)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(9)*(-p*h)/prob_tuple_ss(3) - (p*h)) + ...
            sigma_02*(prob_tuple_ss(1)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(6) + prob_tuple_ss(3)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(6) + prob_tuple_ss(5)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(7)*0/prob_tuple_ss(6) + prob_tuple_ss(8)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(9)*(p*h)/prob_tuple_ss(6) - (-2*p*h + p*h)) + ...
            sigma_12*(prob_tuple_ss(1)*0/prob_tuple_ss(9) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(9) + prob_tuple_ss(3)*(p*h)/prob_tuple_ss(9) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(9) + prob_tuple_ss(5)*0/prob_tuple_ss(9) + ...
            prob_tuple_ss(6)*(p*h)/prob_tuple_ss(9) + prob_tuple_ss(7)*0/prob_tuple_ss(9) + ...
            prob_tuple_ss(8)*0/prob_tuple_ss(9) - (-p*h + p*h)); 
        
    % for drive 1 while in any of the states but p and p1 as two diffusion
    % parameters
    bMem = sigma_neg11*(prob_tuple_ss(2)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(3)*0/prob_tuple_ss(1) + prob_tuple_ss(4)*(-p1*h)/prob_tuple_ss(1) + ...
            prob_tuple_ss(5)*0/prob_tuple_ss(1) + prob_tuple_ss(6)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(7)*(-p*h-p1*h)/prob_tuple_ss(1) + prob_tuple_ss(8)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(1) - (p*h)) + ...
            sigma_01*(prob_tuple_ss(1)*(p*h)/prob_tuple_ss(4) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(4) + prob_tuple_ss(3)*0/prob_tuple_ss(4) + ...
            prob_tuple_ss(5)*0/prob_tuple_ss(4) + prob_tuple_ss(6)*0/prob_tuple_ss(4) + ...
            prob_tuple_ss(7)*(p*h)/prob_tuple_ss(4) + prob_tuple_ss(8)*0/prob_tuple_ss(4) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(4) - (-p1*h + p*h)) + ...
            sigma_11*(prob_tuple_ss(1)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(7) + prob_tuple_ss(3)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(4)*(p1*h)/prob_tuple_ss(7) + prob_tuple_ss(5)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(6)*0/prob_tuple_ss(7) + prob_tuple_ss(8)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(7) - (-p*h-p1*h + p1*h));     
    
    % for drive 1 while in state 1 and drive 2 in state -1
    bMem = sigma_neg12*(prob_tuple_ss(1)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(3) + prob_tuple_ss(4)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(5)*(-p*h)/prob_tuple_ss(3) + prob_tuple_ss(6)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(7)*(-2*p*h)/prob_tuple_ss(3) + prob_tuple_ss(8)*0/prob_tuple_ss(3) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(3) - (0 + 0 + 0 + 0 + 0 + p*h + 0 + 0)) + ...
            sigma_00*(prob_tuple_ss(1)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(5) + prob_tuple_ss(3)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(5) + prob_tuple_ss(6)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(7)*(p*h)/prob_tuple_ss(5) + prob_tuple_ss(8)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(5) - (-p*h + 0 + 0 + 0 + 0 + 0 + 0 + 0)) + ...
            sigma_11*(prob_tuple_ss(1)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(7) + prob_tuple_ss(3)*(p*h)/prob_tuple_ss(7) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(7) + prob_tuple_ss(5)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(6)*0/prob_tuple_ss(7) + prob_tuple_ss(8)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(7) - (-2*p*h + 0 + 0 + 0 + p*h + 0 + 0 + 0));    
    
    %for drive 0 in state 0
    bMem = sigma_neg11*(prob_tuple_ss(2)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(3)*0/prob_tuple_ss(1) + prob_tuple_ss(4)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(5)*(-p*h)/prob_tuple_ss(1) + prob_tuple_ss(6)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(7)*(-2*p*h)/prob_tuple_ss(1) + prob_tuple_ss(8)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(1) - (p*h)) + ...
            sigma_00*(prob_tuple_ss(1)*(p*h)/prob_tuple_ss(5) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(5) + prob_tuple_ss(3)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(5) + prob_tuple_ss(6)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(7)*(p*h)/prob_tuple_ss(5) + prob_tuple_ss(8)*0/prob_tuple_ss(5) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(5) - (-p*h)) + ...
            sigma_11*(prob_tuple_ss(1)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(7) + prob_tuple_ss(3)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(7) + prob_tuple_ss(5)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(6)*0/prob_tuple_ss(7) + prob_tuple_ss(8)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(7) - (-2*p*h+p*h));  

    %for drive 2 in state 0
    bMem = sigma_neg11*(prob_tuple_ss(2)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(3)*0/prob_tuple_ss(1) + prob_tuple_ss(4)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(5)*0/prob_tuple_ss(1) + prob_tuple_ss(6)*(-2*p*h)/prob_tuple_ss(1) + ...
            prob_tuple_ss(7)*(-2*p*h)/prob_tuple_ss(1) + prob_tuple_ss(8)*0/prob_tuple_ss(1) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(1) - (p*h)) + ...
            sigma_02*(prob_tuple_ss(1)*(p*h)/prob_tuple_ss(6) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(6) + prob_tuple_ss(3)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(6) + prob_tuple_ss(5)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(7)*(p*h)/prob_tuple_ss(6) + prob_tuple_ss(8)*0/prob_tuple_ss(6) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(6) - (-2*p*h + p*h)) + ...
            sigma_11*(prob_tuple_ss(1)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(2)*0/prob_tuple_ss(7) + prob_tuple_ss(3)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(4)*0/prob_tuple_ss(7) + prob_tuple_ss(5)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(6)*(p*h)/prob_tuple_ss(7) + prob_tuple_ss(8)*0/prob_tuple_ss(7) + ...
            prob_tuple_ss(9)*0/prob_tuple_ss(7) - (-2*p*h + p*h)); 

    C_01bmem(:,n) = xcorr(sigma_01 - mean(sigma_01), bMem - mean(bMem),'coeff');
    C_00bmem(:,n) = xcorr(sigma_00 - mean(sigma_00), bMem - mean(bMem), 'coeff');
    C_02bmem(:,n) = xcorr(sigma_02 - mean(sigma_02), bMem - mean(bMem),'coeff');
    C_neg11bmem(:,n) = xcorr(sigma_neg11 - mean(sigma_neg11), bMem - mean(bMem), 'coeff');
    C_neg10bmem(:,n) = xcorr(sigma_neg10 - mean(sigma_neg10), bMem - mean(bMem), 'coeff');
    C_neg12bmem(:,n) = xcorr(sigma_neg12 - mean(sigma_neg12), bMem - mean(bMem), 'coeff');
    C_11bmem(:,n) = xcorr(sigma_11 - mean(sigma_11), bMem - mean(bMem), 'coeff');
    C_10bmem(:,n) = xcorr(sigma_10 - mean(sigma_10), bMem - mean(bMem), 'coeff');
    C_12bmem(:,n) = xcorr(sigma_12 - mean(sigma_12), bMem - mean(bMem), 'coeff');

end

%% Saving each random number generator iteration

save(strcat('test_state0_drive0_correlations_eta0pt2_',num2str(random_num),'.mat'),'-v7.3');
end
