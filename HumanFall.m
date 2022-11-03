fs = 3200;      %sampling rate
N = 25600;      %data length
G = 10000;
t = (0:N-1)/fs;
x1 = VarName2;
x2 = VarName3;
x3 = VarName4;
wp = 2/(fs/2);      %set passband digital angular frequency
ws = 3/(fs/2);      %set stopband digital angular frequency
alpha_p = 3;        %maximum available attenuation of passband
alpha_s = 20;       %minimum available attenuation of stopband
%obtain order and cut-off frequency
[ N2,wc2 ] = buttord( wp , ws , alpha_p , alpha_s);
%obtain transfer function coefficients
[ b,a ] = butter(3,wc2,'high');
%filtering
X1 = filter(b,a,x1);
X2 = filter(b,a,x2);
X3 = filter(b,a,x3);
for i=1:N
    %expression of parameters
    SV(i) = sqrt(x1(i)^2 + x2(i)^2 + x3(i)^2);
    SVD(i) = sqrt(X1(i)^2 + X2(i)^2 + X3(i)^2);
    BVA(i) = (SV(i)^2 - SVD(i)^2 - G^2)/(2*G);
end
v1 = max(SV)/G;
v2 = max(SVD)/G;
v3 = max(BVA)/G;
figure(1);
subplot(3,1,1);
plot(t,SV/G);
xlabel('Time (s)');
ylabel('Acceleration (g)');
subplot(3,1,2);
plot(t,SVD/G);
xlabel('Time (s)');
ylabel('Acceleration (g)');
subplot(3,1,3);
plot(t,BVA/G);
xlabel('Time (s)');
ylabel('Acceleration (g)');




