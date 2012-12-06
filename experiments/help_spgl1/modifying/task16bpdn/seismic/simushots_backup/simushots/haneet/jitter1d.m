function idx = jitter1d(width, spacing)

% sapcing = 5;
% width   = 256;
% regular_samples  = 3:spacing:253;
% jittered_samples = zeros(1,255);

% regular_samples  = (spacing + 1) : spacing : (width - 1);
regular_samples  = (spacing + 1) : spacing : (width - 1);
jittered_samples = zeros(1, width);

for i = 1:size(regular_samples,2)
    
    perturb = rand()*spacing - spacing/2;
    
    if (perturb < 0)
        perturb = floor(perturb - 0.5);
    else
        perturb = floor(perturb + 0.5);
    end
    
    jittered_samples(floor(regular_samples(i) + perturb)) = 1;
    
end

% Collect indices of jittered locations
idx = find(jittered_samples);

end

% figure
% jittered_spectrum = abs(fft(jittered_samples));
% % freq_axis = [-127:127]*2*pi/length(jittered_samples);
% % plot(freq_axis, fftshift(jittered_spectrum));
% plot(jittered_spectrum);
% axis tight; title('Jittered spectrum');
% xlabel('Frequency'); ylabel('Amplitude')