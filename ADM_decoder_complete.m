function [xq] = ADM_decoder_complete(bn)
bn_len = length(bn);
delta = zeros(bn_len, 1);

delta_step = 0.5;
delta(1) = delta_step;

errorq = zeros(bn_len, 1);
xq = zeros(bn_len, 1);

K = 1.5;

for n=2:bn_len
    
    if(bn(n) == bn(n-1))
        delta(n) = delta(n-1)*K;
    elseif (bn(n) ~= bn(n-1))
        delta(n) = delta(n-1)/K;
    end
    errorq(n) = bn(n)*delta(n);
    xq(n) = errorq(n) + xq(n-1);
    
end
xq = xq(1:2:end);
end