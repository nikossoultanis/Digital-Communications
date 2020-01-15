Lb = randi([10000 100000],1,1); % random Lb sequence
bergr=zeros(11,1);
sergr=zeros(11,1);
M = 4;
if (M == 4)
for SNR=0:2:20
    bersgr=0;
    sers=0;
    Lb = 10002;
    symbol_bits = log2(M); %kalo
    Lbnew= Lb + symbol_bits - (mod(Lb, symbol_bits));
    binary_series = randsrc (Lb, 1 ,[0 1]);
    bin2symbol = zeros(Lbnew / symbol_bits, 1);
    
    for i=1:Lb/log2(M)
        for j=1:symbol_bits
            bin2symbol(i) = bin2symbol(i) + binary_series(i*log2(M)-j+1)*2^(j-1); %mapping to 0 1 2 3
        end
    end
    s = zeros(length(bin2symbol),1);
    for i=1:length(bin2symbol)
        s(i) = 2*bin2symbol(i)-3;
    end
    
    constellation = zeros((log2(M)), 1);
    for i=1:M
        constellation(i)=2*i-1-M;%dianysma me ta simeia twn asterismwn
    end
    
    Tsymbol = 4;
    fc = 2.5;
    Tsample = 0.1;
    gt = sqrt(1/Tsymbol);
    % Diamorfwsi
    diamorfosi = zeros(4*length(s), 1);
    counter = 0;
    for i=1:4:length(diamorfosi)
        counter = counter+1;
        for t=0:3
            diamorfosi(i+t) = s(counter)*gt*cos(2*pi*fc*t);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       mple = repelem(mple , 4);
    
    Eb=((M^2-1)/3)/log2(M);
    stet=Eb/(2*(10^(SNR/10)));
    noise=sqrt(stet)*randn(size(diamorfosi));
    sg = diamorfosi + noise;
    distance=zeros(M, 1);
    fd=zeros(Lb/log2(M), 1);
    avga = zeros(length(sg),1);
    %Apodiamorfwtis
    for i=1:4:length(sg)
        for t=0:3
            avga(i) = avga(i)+ sg(i+t)*gt*cos(2*pi*fc*t);
        end
    end
    avga = avga(1:4:end);
    
    d=zeros(M,1);
    fd=zeros(Lb/log2(M),1);
    for i=1:size(avga)
        for j=1:M
            d(j)=sqrt((avga(i)-constellation(j))^2);
        end
        [fx, ifx]=min(d);
        fd(i,1)=constellation(ifx);
        outg(i,1)=ifx-1;
    end
    outs=dec2bin(outg);
    otp=zeros(Lb,1);
    for i=1:Lb/log2(M)
        for j=1:log2(M)
            otp((i-1)*log2(M)+j)=str2double(outs(i,j));
        end
    end
    %euresi ser
    for i=1:1:Lb/log2(M)
        if s(i)~=fd(i)
            sers=sers+1;
        end
    end
    sergr(SNR/2+1)=sers/(Lb/log2(M));
    %euresi ber
    for i=1:1:Lb
        if binary_series(i)~=otp(i)
            bersgr=bersgr+1;
        end
    end
    bergr(SNR/2+1)=bersgr/Lb;
    %aparaitites synthikes gia xrisi tis semilogy
    if bergr(SNR/2+1,1)<=0
        bergr(SNR/2+1,1)=eps;
    end
    if sergr(SNR/2+1,1)<=0
        sergr(SNR/2+1,1)=eps;
    end
    
end
end
%     aparaitites synthikes gia xrisi tis semilogy
semilogy(0:2:20,bergr(1:SNR/2+1));
hold on;
semilogy(0:2:20,sergr(1:SNR/2+1));
hold on;

