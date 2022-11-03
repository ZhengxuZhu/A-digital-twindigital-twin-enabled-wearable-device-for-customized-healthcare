fs = 6400;      %sampling rate
N = 76800;      %data length
t = (0:N-1)/fs;
delta_f = 1*fs/N;
f = (-N/2:N/2-1)*delta_f;
wp=[20*2/fs];       %set passband digital angular frequency
ws=[35*2/fs];       %set stopband digital angular frequency
Rp = 3;     %maximum available attenuation of passband
Rs = 20;    %minimum available attenuation of stopband
[ N3,wn ] = buttord( wp , ws , Rp , Rs);        %obtain order and cut-off frequency
[ b,a ] = butter(N3,wn,'low');      %obtain transfer function coefficients
%filtering
x = VarName2;
filter_bp_s = filter(b,a,x);
X_bp_s = fftshift(abs(fft(filter_bp_s)))/N;
t1 = 3.1; t2 = t1+3;  %time domain segment selection, here take the data of the tenth group of rehabilitation training as an example11
N1 = fix(t1*fs); N2 = fix(t2*fs);
Max = max(filter_bp_s(N1:N2)); Min = min(filter_bp_s(N1:N2)); Gap = abs(abs(Min)-abs(Max))/2;
Y = zeros(1, N2-N1+1);  
for i=1:N2-N1+1
    Y(i) = filter_bp_s(N1 - 1 + i) + Gap;       %set baseline to be 0
end
Z = zeros(1, length(Y));
O = zeros(1, length(Y));
L = asin(Y(1)/max(Y));      %initial phase
cal=0;
for w1=70:1:85
    for w2=70:1:85
        for w3 = 50:1:70
            for w4 = 20:1:40
                cal = cal+1;
                for i=1:N2-N1+1
                    %expression of mathematical model
                    Z(i) = w3*sin(w1*(i-1)/fs + L)- w4*sin(w2*(i-1)/fs + L);
                    YZ(i) = Y(i)*Z(i);
                    YY(i) = Y(i)*Y(i);
                    ZZ(i) = Z(i)*Z(i);
                end
                A = trapz(YZ);B = trapz(YY);C = trapz(ZZ);
                pyz(cal) = A/sqrt(B*C);
                if abs(pyz(cal))> max(abs(pyz(1:cal-1)))
                    %output optimal parameters
                    P = pyz(cal);
                    WA = w1;
                    WB = w2;
                    WC = w3;
                    WD = w4;
                end
            end
        end
    end
end
for i=1:N2-N1+1
    %expression of optimal model
    O(i) = WC*sin(WA*(i-1)/fs + L)-WD*sin(WB*(i-1)/fs + L);
end
T = (N1/fs:1/fs:N2/fs);
plot(T,Y);hold;plot(T,O)
axis([t1 t2 -150 150])
xlabel('Time (s)');
ylabel('Velocity (mm/s)');



