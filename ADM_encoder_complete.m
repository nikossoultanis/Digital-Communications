function [bn,xq] = ADM_encoder_complete(signal)
signal_len = length(signal);

delta = zeros(signal_len, 1);
delta(1) = 0.5;

K = 1.5;

xq = zeros(signal_len, 1);

% signal = interp(signal, 2);

en = zeros(signal_len, 1);
bn = zeros(signal_len, 1);
errorq = zeros(signal_len, 1);


if(signal(1)>0)
    bn(1) = 1;
else
    bn(1) = -1;    
end
for n=2:signal_len
    
    en(n) = signal(n) - xq(n-1);
    
    if(en(n) > 0)
        bn(n) = 1;
    else
        bn(n) = -1;
    end
    
    if(bn(n) == bn(n-1))
        delta(n) = delta(n-1)*K;
    else
        delta(n) = delta(n-1)/K;
    end
    errorq(n) = bn(n)*delta(n);
    xq(n) = xq(n-1) + errorq(n);
end
end