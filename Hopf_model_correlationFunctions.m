close all;
clear variables;

% If getting data from a saved figure
%{
open('/Users/janaki/Dropbox/MATLAB/smoothening data and fdt analysis/raw_data/probes_control/probe7x(t).fig');
h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
ydata = get(dataObjs, 'YData');
xdata = get(dataObjs, 'XData');
%}

time_start = datestr(now, 'yymmdd_HHMMSS');         %stores time of execution of program in a string 
save_directory = ['noise_analysis_', time_start]; % name of directory 
mkdir(save_directory); %makes directory

% If getting data from a save mat file
rawdata_matrix = load('/Users/janaki/Dropbox/MATLAB/smoothening data and fdt analysis/raw_data/justin_spikes/quiescent.mat');
rawdata = -quiescent;
sampling_time = 0.001;
time = (sampling_time:sampling_time:length(rawdata)*sampling_time)';
data = [time,rawdata(1:end)];

figure;
plot(time,data(:,2));
baseline = smooth(rawdata,0.5, 'loess');
nodrift_data = data(:,2) - baseline;
[final_data] = nodrift_data;

f1 = figure();
ax = axes(f1);
plot (ax, time, final_data);
title('Raw data');
xlabel('Time(s)');
ylabel('Position(nm)');
savefig(f1, [save_directory, filesep, 'Position of hair bundle.fig']);
print(f1, '-dpng', '-r300',[save_directory, filesep, 'Position of hair bundle.png']);

[psd_filtdt,freqdt] = pwelch(final_data, [],[],length(final_data)*4,1/sampling_time);
freqdt = 2*pi*freqdt;
[max_psd, index] = max(psd_filtdt);
natural_freq = freqdt(index);

f2 = figure();
ax = axes(f2);
plot(ax,freqdt, psd_filtdt);
title('PSD in linear axes');
xlabel('\omega(Hz)');
ylabel('PSD(nm^{2}/Hz)');
savefig(f2, [save_directory, filesep, 'PSD of hair bundle.fig']);
print(f2, '-dpng', '-r300',[save_directory, filesep, 'PSD of hair bundle.png']);

f3 = figure();
ax = axes(f3);
plot(ax, log10(freqdt), log10(psd_filtdt));
title('PSD in linear axes');
xlabel('log_{10}(\omega)(Hz)');
ylabel('log_{10}PSD(nm^{2}/Hz)');
savefig(f3, [save_directory, filesep, 'Log log PSD of hair bundle.fig']);
print(f3, '-dpng', '-r300',[save_directory, filesep, 'Log log PSD of hair bundle.png']);

y_rect = final_data;

%calculation of phase angle wrapped
hilt = hilbert(y_rect, 2.^ceil(log2(length(rawdata))));
phase_array = [y_rect, imag(hilt(1:length(rawdata))) ];
phase_array = phase_array(floor(0.2*length(rawdata)):floor(0.8*length(rawdata)),:);
phase = atan2( phase_array(:,2),  phase_array(:,1));
time = time(floor(0.2*length(rawdata)):floor(0.8*length(rawdata)))- time(floor(0.2*length(rawdata))); 

f = figure();
ax = axes(f);
plot (ax, time, phase);
title('Phase of the limit cycle using atan');
xlabel('Time(s)');
ylabel('Phase(rad)');
savefig(f, [save_directory, filesep, 'Phase of the limit cycle.fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Phase of the limit cycle.png']);

phase_un = unwrap(phase);

%fitting a straight line to derive noise in phase angle
P = polyfit(time, phase_un - phase_un(1), 1);
yfit = P(1)*time + P(2);

f = figure();
ax = axes(f);
plot (ax,time, phase_un - phase_un(1));
hold on;
plot(ax, time, yfit);
title('Unwrapped phase using acos and adding angles');
xlabel('Time(s)');
ylabel('Phase(rad)');
savefig(f, [save_directory, filesep, 'Unwrapped phase using acos and adding angles.fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Unwrapped phase using acos and adding angles.png']);

%deriving diffusion coefficient
diff_over_time = zeros(length(time),1);
corr = zeros(length(time),1);

for i = 1:1000
    count = 0;
    for j = 1:length(time)
        if(j-i >0)
            diff_over_time(i) = diff_over_time(i) + (phase_un(j) - phase_un(j-i)).^2;
            count = count + 1;
        end
       corr(i) = diff_over_time(i)/count;
    end
end

f = figure();
ax = axes(f);
l = length(find(corr > 1e-10));
plot (ax,time(1:l), corr(1:l));
title('Diffusion of total phase');
xlabel('Time(s)');
ylabel('Mean phase displacement (rad^2)');
savefig(f, [save_directory, filesep, 'Diffusion of total phase.fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Diffusion of total phase.png']);

phase_noise = phase_un - phase_un(1) - yfit;

diffsq_over_time = zeros(length(time),1);
diff_over_time = zeros(length(time),1);
corr = zeros(length(time),1);

for i = 1:1000
    count = 0;
    for j = 1:length(time)
        if(j-i >0)
            diffsq_over_time(i) = diffsq_over_time(i) + (phase_noise(j) - phase_noise(j-i)).^2;
            diff_over_time(i) = diff_over_time(i) + (phase_noise(j) - phase_noise(j-i));
            count = count + 1;
        end
       corr(i) = diffsq_over_time(i)/count - (diff_over_time(i)/count).^2;
    end
end

f = figure();
ax = axes(f);
l = length(find(corr > 1e-10));
plot (ax,time(1:l), corr(1:l));
title('Diffusion of differential phase');
xlabel('Time(s)');
ylabel('Mean phase displacement (rad^2)');
savefig(f, [save_directory, filesep, 'Diffusion of differential phase.fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Diffusion of differential phase.png']);


baseline_ph_noise = smooth(phase_noise, 0.5, 'loess');
final_phase_noise = phase_noise - baseline_ph_noise;
[psd_filtph, freqph] = pwelch(final_phase_noise, [],[],length(final_phase_noise)*4, 1/sampling_time);
freqph = 2*pi*freqph;

f = figure();
ax = axes(f);
loglog (ax,freqph, psd_filtph);
title('PSD of differential phase');
xlabel('Log(Frequency)');
ylabel('Log(PSD)');
savefig(f, [save_directory, filesep, 'Loglog plot of psd of differential phase .fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Loglog plot of psd of differential phase.png']);

%Calculating mean amplitude of the limit cycle
phase_sh = phase(3:200);

div = (-pi:2*pi/200:pi)';
sum_amp = zeros(length(div)-1,2);
norm_amp = zeros(length(div)-1,1);
devi_amp = zeros(length(phase)-1,1);
devi_phase = zeros(length(phase)-1,1);
restricted_devi_amp = zeros(length(phase)-1,1);
restricted_phase_array = zeros(length(phase) - 1,2);
norm_array = zeros(length(phase) - 1, 1);

cnt = 0;
for i = 1:length(div)-1
    ind = find(phase >= div(i) & phase < div(i+1));
    for j = 1:length(ind)
        sum_amp(i,:) = sum_amp(i,:) + phase_array(ind(j),:)/length(ind);
        norm_amp(i) = norm_amp(i) + norm(phase_array(ind(j),:))/length(ind);
    end  
   for ph = 1:length(phase)
       if(any(abs(ph-ind)<1e-10))
        devi_amp(ph) = norm(phase_array(ph,:)) - norm_amp(i);
        restricted_phase_array(ph, :) = sum_amp(i,:);
        norm_array(ph) = norm_amp(i);
       end
   end
end

%[freqdv,psddv, psd_filtdv] = spectral_density(devi_amp, devi_amp, 1/sampling_time, 2);
baseline_amp_devi = smooth(devi_amp, 0.02, 'loess');
amp_devi = devi_amp - baseline_amp_devi;
[psd_filtdv, freqdv] = pwelch(devi_amp, [],[],length(final_data)*4, 1/sampling_time);
freqdv = 2*pi*freqdv;

f = figure();
ax = axes(f);
plot (ax,phase_array(:,1), phase_array(:,2));
hold on;
plot(ax,sum_amp(:,1), sum_amp(:,2));
title('Noisy and approximated limit cycle ');
xlabel('x(\phi)');
ylabel('x(\phi + \pi/2)');
savefig(f, [save_directory, filesep, 'Limit cycles.fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Limit cycles.png']);

f = figure();
ax = axes(f);
loglog(ax,freqdv, psd_filtdv);
title('PSD of differential amplitude');
xlabel('Log(frequency)');
ylabel('Log(PSD)');
savefig(f, [save_directory, filesep, 'Loglog plot of PSD of differential amplitudes.fig']);
print(f, '-dpng', '-r300',[save_directory, filesep, 'Loglog plot of PSD of differential amplitudes.png']);

save([save_directory, filesep, 'data.mat']);  %saves workspace
zip([save_directory, filesep, 'code_snapshot.zip'] , {'*.m'});  %saves all the m files in working directory into a zip file
