
clear variables;

time_start = datestr(now, 'yymmdd_HHMMSS');         %stores time of execution of program in a string 
save_directory = ['greensfunction_with_limitcycle_sim', time_start]; % name of directory 
mkdir(save_directory); %makes directory

mu = 40;
w0 = 10;
Kb = 1.38*10^-23;
T = 300;
b = 2;
b1 = 0;

h = 0.00005;
t_final = 300000;
num = 50000;

X = zeros(t_final +1,1);
Y = zeros(t_final +1,1);
time = zeros(t_final +1,1);
X(1) = 1;
Y(1) = 1;

noise_X = randn([length(time)-1,1]);
noise_Y = randn([length(time)-1, 1]);
noise_X_amp = 0;
noise_Y_amp = 0;

Fn_1 = @(X,Y) mu*X - w0*Y - b*(X^2 + Y^2)*X - b1*(X^2 + Y^2)*Y;                   
Fn_2 = @(X,Y) mu*Y + w0*X - b*(X^2 + Y^2)*Y + b1*(X^2 + Y^2)*X;

for i = 1:num
    
    k_1 = Fn_1(X(i),Y(i));
    m_1 = Fn_2(X(i), Y(i));
    k_2 = Fn_1(X(i) + 0.5*h*k_1 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_1 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    m_2 = Fn_2(X(i) + 0.5*h*k_1 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_1 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    k_3 = Fn_1(X(i) + 0.5*h*k_2 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_2 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    m_3 = Fn_2(X(i) + 0.5*h*k_2 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_2 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    k_4 = Fn_1(X(i) + h*k_3 + sqrt(noise_X_amp*h)*noise_X(i), Y(i) + h*m_3 + sqrt(noise_Y_amp*h)*noise_Y(i));
    m_4 = Fn_2(X(i) + h*k_3 + sqrt(noise_X_amp*h)*noise_X(i), Y(i) + h*m_3 + sqrt(noise_Y_amp*h)*noise_Y(i));

    X(i+1) = X(i) + (1/6)*(k_1+ 2*k_2 + 2*k_3 + k_4)*h + sqrt(noise_X_amp*h)*noise_X(i);  
    Y(i+1) = Y(i) + (1/6)*(m_1 +2*m_2 + 2*m_3 + m_4)*h + sqrt(noise_Y_amp*h)*noise_Y(i);
    
end

[psd_X, freq_X] = pwelch(X(1000:num), [],[],length(X(1000:num))*4,1/h);

%plotting psd and calculating natural frequencies
[~, index] = max(psd_X);
natural_freq = freq_X(index);

limit_cycle = [X(num - floor(1/(h*natural_freq)+6) :num), Y(num - floor(1/(h*natural_freq)+6):num)];
phase = atan2(limit_cycle(:,2), limit_cycle(:,1));
phase_limit_cycle = unwrap(phase);
P = polyfit(h*(1:length(limit_cycle))', phase_limit_cycle,1);
T = zeros(length(limit_cycle)-1,2);
N = zeros(length(limit_cycle)-1,2);

for t_id = 1:length(T)
 T(t_id,:) =  limit_cycle(t_id+1,:) - limit_cycle(t_id,:);
 T(t_id,:) = T(t_id,:)./norm(T(t_id,:));
end

for n_id = 2:length(N)
 N(n_id,:) = T(n_id,:) - T(n_id-1,:);
 N(n_id,:) = N(n_id,:)./norm(N(n_id,:));
end

N(1,:) = T(1,:) - T(end,:);
N(1,:) = N(1,:)./norm(N(1,:));

cnt = 0;
ping_num = 0;

devi_amp_array = [];
phase_dev_array = [];
ping = [];

for i= num:(length(time)-1)
    ph = atan2(Y(i), X(i));

    if(mod(floor((i - num)*natural_freq*h), 6) == 0 && (0 < ph) &&(ph <natural_freq*2*pi*h - 0) && ping_num==0)
        
        [~, force_id] = min((X(i) - limit_cycle(:,1)).^2+ (Y(i)-limit_cycle(:,2)).^2);
        force = 1000;
        force_X = force*N(force_id,2);
        force_Y = -force*N(force_id,1);
        
        k_1 = Fn_1(X(i),Y(i));
        m_1 = Fn_2(X(i), Y(i));
        k_2 = Fn_1(X(i) + 0.5*h*k_1 +force_X*h, Y(i) + 0.5*h*m_1 + force_Y*h);
        m_2 = Fn_2(X(i) + 0.5*h*k_1 + force_X*h, Y(i) + 0.5*h*m_1 + force_Y*h);
        k_3 = Fn_1(X(i) + 0.5*h*k_2+ force_X*h, Y(i) + 0.5*h*m_2 + force_Y*h);
        m_3 = Fn_2(X(i) + 0.5*h*k_2 + force_X*h, Y(i) + 0.5*h*m_2+ force_Y*h);
        k_4 = Fn_1(X(i) + h*k_3 + force_X*h, Y(i) + h*m_3 + force_Y*h);
        m_4 = Fn_2(X(i) + h*k_3+ force_X*h , Y(i) + h*m_3 + force_Y*h);

        X(i+1) = X(i) + (1/6)*(k_1+ 2*k_2 + 2*k_3 + k_4)*h + force_X*h;  
        Y(i+1) = Y(i) + (1/6)*(m_1 +2*m_2 + 2*m_3 + m_4)*h + force_Y*h;
        
        ping_num = i;
        ping = [ping,ping_num];
        i = i+1;
        cnt = cnt+1;
        
    end  
    
    noise_X_amp = 0;
    noise_Y_amp = 0;
        
    k_1 = Fn_1(X(i),Y(i));
    m_1 = Fn_2(X(i), Y(i));
    k_2 = Fn_1(X(i) + 0.5*h*k_1 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_1 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    m_2 = Fn_2(X(i) + 0.5*h*k_1 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_1 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    k_3 = Fn_1(X(i) + 0.5*h*k_2 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_2 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    m_3 = Fn_2(X(i) + 0.5*h*k_2 + 0.5*sqrt(noise_X_amp*h)*noise_X(i), Y(i) + 0.5*h*m_2 + 0.5*sqrt(noise_Y_amp*h)*noise_Y(i));
    k_4 = Fn_1(X(i) + h*k_3 + sqrt(noise_X_amp*h)*noise_X(i), Y(i) + h*m_3 + sqrt(noise_Y_amp*h)*noise_Y(i));
    m_4 = Fn_2(X(i) + h*k_3 + sqrt(noise_X_amp*h)*noise_X(i), Y(i) + h*m_3 + sqrt(noise_Y_amp*h)*noise_Y(i));

    X(i+1) = X(i) + (1/6)*(k_1+ 2*k_2 + 2*k_3 + k_4)*h + sqrt(noise_X_amp*h)*noise_X(i);  
    Y(i+1) = Y(i) + (1/6)*(m_1 +2*m_2 + 2*m_3 + m_4)*h + sqrt(noise_Y_amp*h)*noise_Y(i);
    

    if (i == ping_num + 3000 && ping_num~=0)
        nearest_id = zeros(3000,1);
        devi_amp = zeros(3000,1);
        for j = 1:3000
            [~, nearest_id(j)] = min((X(ping_num+j) - limit_cycle(:,1)).^2+ (Y(ping_num+j)-limit_cycle(:,2)).^2);
            if nearest_id(j) ~= length(limit_cycle)
                devi_amp(j) = dot([X(ping_num+j),Y(ping_num+j)] - limit_cycle(nearest_id(j),:), N(nearest_id(j),:));        
            else
                devi_amp(j) = dot([X(ping_num+j),Y(ping_num+j)] - limit_cycle(nearest_id(j),:), N(1,:));                
            end
        end
        
        devi_amp_array = [devi_amp_array, devi_amp];
        [f_r, psd_r] = spectral_density(devi_amp, devi_amp, 1/h);
        figure; 
        plot(log10(f_r), log10(psd_r));
        rawdata = [X(ping_num:ping_num+j), Y(ping_num:ping_num+j)];
        phase = atan2(rawdata(:,2), rawdata(:,1));
        phase_rawdata = unwrap(phase);
        phase_diff = (phase_rawdata(1) +  h*P(1)*(0:length(rawdata)-1)) - phase_rawdata';
        phase_dev = (phase_diff(3:end) - phase_diff(2:end-1))/h;
        [f_phdev, psd_phdev] = spectral_density(phase_dev, phase_dev, 1/h);
        phase_dev_array = [phase_dev_array, phase_dev];
        figure; 
        plot(log10(f_phdev),log10(psd_phdev.*(f_phdev.^2)));
        ping_num = 0;
    end
end

save([save_directory, filesep, 'data.mat']);  %saves workspace
zip([save_directory, filesep, 'code_snapshot.zip'] , {'*.m'});  %saves all the m files in working directory into a zip file


