% Lb = randi([10000 100000],1,1); % random Lb sequence
ber=zeros(11,1);
ser=zeros(11,1);
M=8;
for SNR=0:2:20
    bers=0;
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
    grey_series = bin2gray(bin2symbol,'pam', M);
    for i=1:size(bin2symbol)
        s(i,1)=2*grey_series(i)+1-M;
    end
    
    constellation = zeros((log2(M)), 1);
    for i=1:M
        val(i)=bin2gray(i-1,'pam',M);
        %     va=dec2base(val,2,log2(M));%text gia asterismous
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
    out=gray2bin(outg,'pam',M);
    outs=dec2bin(out);
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
    ser(SNR/2+1)=sers/(Lb/log2(M));
    %euresi ber
    for i=1:1:Lb
        if binary_series(i)~=otp(i)
            bers=bers+1;
        end
    end
    ber(SNR/2+1)=bers/Lb;
    %aparaitites synthikes gia xrisi tis semilogy
    if ber(SNR/2+1,1)<=0
        ber(SNR/2+1,1)=eps;
    end
    if ser(SNR/2+1,1)<=0
        ser(SNR/2+1,1)=eps;
    end
    
end
semilogy(0:2:20,ber(1:SNR/2+1));
hold on;
semilogy(0:2:20,ser(1:SNR/2+1));
hold on;
constellation=zeros(log2(M),1);
for i=1:M
    constellation(i)=2*i-1-M;
end
yy=figure();
ff=zeros(M,1);
plot(constellation,ff,'-ro');
