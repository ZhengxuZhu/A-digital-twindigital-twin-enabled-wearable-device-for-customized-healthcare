fs = 6400;      %sampling rate
N = 76800;      %data length
t = (0:N-1)/fs;  
delta_f = 1*fs/N;
f = (-N/2:N/2-1)*delta_f;
wp=[3*2/fs];        %set passband digital angular frequency
ws=[5*2/fs];        %set stopband digital angular frequency
Rp = 5;     %maximum available attenuation of passband
Rs = 20;    %minimum available attenuation of stopband
[ N3,wn ] = buttord( wp , ws , Rp , Rs);        %obtain order and cut-off frequency
[ b,a ] = butter(N3,wn,'low');      %obtain transfer function coefficients
%filtering
x = VarName2;
filter_bp_s = filter(b,a,x);
X_bp_s = fftshift(abs(fft(filter_bp_s)))/N;
t1 = 3.278; t2 = 3.921;     %time domain segment selection, here take the data of the second group of rehabilitation training as an example
N1 = fix(t1*fs); N2 = fix(t2*fs);
Max = max(filter_bp_s(N1:N2)); Min = min(filter_bp_s(N1:N2)); Gap = abs(abs(Min)-abs(Max))/2;  
for i=1:N2-N1+1
    Y(i) = filter_bp_s(N1 - 1 + i) + Gap;  %set baseline to be 0
    T(i) = t1 + (i-1)/fs;
end
L = asin(Y(1)/max(Y));      %initial phase
cal = 0;
for A1 = 6500:10:8500
    for W1 = 5:0.1:15
        cal = cal+1;
        for i=1:N2-N1+1
            %expression of mathematical model
            Z(i) = A1*sin(W1*(i-1)/fs + L);
            YZ(i) = Y(i)*Z(i);
            YY(i) = Y(i)*Y(i);
            ZZ(i) = Z(i)*Z(i);
        end
        A = trapz(YZ);B = trapz(YY);C = trapz(ZZ);
        pyz(cal) = A/sqrt(B*C);
        if abs(pyz(cal))> max(abs(pyz(1:cal-1)))
            %output optimal parameters
            P = pyz(cal);
            WA = A1;
            WB = W1;
        end
    end
end
for i=1:N2-N1+1
    %expression of optimal model
    O(i) = WA*sin(WB*(i-1)/fs + L);
end
plot(T,Y);hold;plot(T,O);
xlabel('Time (s)');
ylabel('Acceleration (mm/s^2)');


