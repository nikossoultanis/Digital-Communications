function [quantized_signal, centers, distortion,Ps,Pn,sqnr] = MaxLloydQuant(signal, NoB, min_value, max_value)
tic
ps = zeros(30,1);
pn = zeros(30,1);
error = 10^(-6); % least error for our calculations
signal_len = length(signal); % input signal length
stages = 2^NoB; % NoB = number of bits.
range = max_value - min_value; % range from top to bottom of the signal
max_error = range / stages; % maximum error we can achieve
distortion = [0 max_error]; % array that keeps the total distortion
delta = range/stages; %delta step
centers = (min_value:delta:max_value); % centers similar to uniform quantizer
k = 2; % first iteration

T = zeros(length(centers)); % quantization zones
T(1) = min_value; % first quant. zone
T(length(centers)) = max_value; % last quant. zone
quantized_signal = zeros(signal_len, 1); % initialize
counter = 1;
while abs(distortion(k) - distortion(k-1)) >= error % until we reach the error
    
    entries   = zeros(length(centers),1);
    summary_values = zeros(length(centers),1);
    calc_error = 0; % calculation error from input and quantized values
    
    for i=2:(length(centers)-1)
        T(i) = (centers(i) + centers(i+1))/2; % filling the rest quant. zones
    end
    
    for i=1:signal_len % for each value of the signal
        for j=1:(length(T)-1) % for each quantizing zone
            if signal(i) > T(j) && signal(i) <= T(j+1) % if value is between T(j) and T(j+1)
                quantized_signal(i) = centers(j); % quantized value = center value
                calc_error = calc_error + (abs(centers(j+1) - signal(i))); % calc. error
                summary_values(j) = summary_values(j) + signal(i); % values of the signal that went to T(j)
                entries(j)   = entries(j) + 1; % entries on the same quantizing zone
            end
        end
    end
    distortion = [distortion calc_error/signal_len]; % distortion appended on each iteration
    k = k + 1; % next iteration
    
    for j=1:(length(centers)) % calculating the new centers
        if entries(j) ~= 0 % if there are entries for this center
            centers(j+1) = summary_values(j)/entries(j); % create the new one according the the mean value
        end
    end
    Ps(counter) = mean(signal.^2);
    Pn(counter) = mean((signal-quantized_signal).^2);
    sqnr(counter) = Ps/Pn;
    counter = counter + 1;
end
toc
fprintf("Successfully exited the quantizer after %d iterations \n", length(distortion)-2);
fprintf("And the distortion is %d. \n", distortion(length(distortion)));


end