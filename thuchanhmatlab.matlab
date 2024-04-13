% Cau 14 phan 1
t = 0 : 0.1 : 6*pi;
a = sin(t);
b = cos(t);
plot(t, a)
hold on 
plot(t, b)
xlabel('Thoi gian')
ylabel('gia tri')
grid on 
legend('sinx', 'cosx')

%---------------------------------------------
% Cau 15 phan 1
T = 4;
t = -T:0.01:T;
x = zeros(length(t));
ind = find((t >= -T/2)&(t <= T/2));
x(ind) = 1;
plot(t, x);
grid on
axis([-4 4 0 2]);

%---------------------------------------------
% Cau 16 phan 1
T = 2;
t = -T : 0.01 : T;
x = zeros(length(t));
N = 2 * T / 0.01;
for i = 1 : N
    x(i) = 1 - abs(t(i)) / T;
end
plot(t,x)
grid on

%---------------------------------------------
% Cau 17 
function[day] = thuchanh4(mon,year)
  switch mon
      case 1
          day = 31;
      case 3
          day = 31;
      case 5
          day = 31;
      case 7
          day = 31;
      case 8
          day = 31;
      case 10
          day = 31;
      case 12
          day = 31;
      case 2
          if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
              day = 29;
          else
              day = 28;
          end
      otherwise
          day = 30;
  end
end

%---------------------------------------------
% Cau 18
function[kq] = thuchanh5(N)
    kq = 1;
    if N == 0
        kq = 1;
    else
        for k = 1:N
            kq = k * kq; 
        end
    end
end

%---------------------------------------------
% Cau 19
% Cong thuc sum = n * (n+1) * (2n+1) / 6
function[max] = thuchanh6(N)
    M = N * 100;
    k = 1;
    while k^2 <= M
        M = M - k^2;
        max = k;
        k = k + 1;
    end
end

%---------------------------------------------
% cau 20
function[sofibo] = thuchanh7(N)
   x = zeros(1, N);
   x(1) = 0;
   x(2) = 1;
   for i = 3 : N
       x(i) = x(i - 1) + x(i - 2);
   end
   sofibo = x(N);
end

%---------------------------------------------
% Cau 21
%function[sum1] = thuchanh8(a, N)
   % sum1 = 0;
   % for k = 1 : N
   %    sum1 = sum1 + ((-2)^a) / (exp(-k));
   % end
   % sum1 =  sprintf('%.7f',sum1);
%end

% Cau 22 
function[multi] = thuchanh8(a , N)
    multi = 1;
    for k = 1: N
        multi = multi * (k^2 + sqrt((a * k + 1)/2));
    end
    multi = sprintf('%.7f',multi);
end

%---------------------------------------------
% Cau 23
% function = exp(x) * atan(x^2) / cos(x);
% dientich = thuchanh9(4,7,28,@(x) exp(x) * atan(x^2) / cos(x))
function[S] = thuchanh9(a,b,N,f)
    N = N * 100;
    h = (b - a)/N;
    S = f(a + h/2);
    for k = 2 : N
        xk = a + (k - 1/2) * h;
        S = S + f(xk);
    end
    S = S * h;
    S = sprintf('%.5f', S);
end

%---------------------------------------------
% Cau 24
% f = 4*x^3 - 13*x^2 + 13*x - 10
function[x] = thuchanh10(a, b, f)
    fa = f(a);
    fb = f(b);
    while b - a > eps*b
        if sign(fa) ~= sign(fb)
            x = (b+a)/2;
            fx = f(x);
            if sign(fx) == sign(fb)
                b = x;
                fb = fx;
            else
                a = x;
                fa = fx;
            end
        else
            disp('Phuong trinh vo nghiem trong khoang tren')
            break;
        end
    end 
end

%---------------------------------------------
% Cau 26
% f = 4*x^2 - 7*y
% co the khong can lap mang x
function[y] = thuchanh11(xmax,y0,h,f)
    x = 0 : h :xmax;
    y = zeros(1, length(x));
    y(1) = y0;
    for i = 2 : length(x)
        y(i) = y(i -1) + h* f(x(i-1), y(i-1));
    end
end

%---------------------------------------------
% Cau 27
N = 9; % so chuoi xung chu nhat
T = 4;
t = 0 : 0.01 : 2*T;
y = zeros(1,length(t));
y(t<=T) = 1;
y1 = y;
t1 = t;
for i = 1 : N-1
    t1 = [t1 t+i*2*T];
    y1 = [y1 y];
end
plot(t1, y1);
grid on
axis([0 80 0 2])

%---------------------------------------------
% Cau 28 
N = 3;
T = 30;
t = 0: 0.01 : T;
y = zeros(1,length(t));
n = T/0.01;
for k = 1 : n
    if t(k) <= T/2
        y(k) = 2*t(k)/T;
    else
        y(k) = 2 - 2*t(k)/T;
    end
end
y1 = y;
t1 = t;
for k = 1 :N-1
    t1 =[t1 t+k*T];
    y1 =[y1 y];
end
plot(t1,y1)
grid on 

%---------------------------------------------
% Cau 29
N = 9;
T = 1;
t = 0: 0.01 : T;
y = zeros(1, length(t));
n= T/0.01;
for k = 1 : n+1
    y(k) = exp(t(k)^2);
end
t1 = t;
y1 = y;
for k = 1: N-1
    t1 = [t1 t+k*T];
    y1 = [y1 y];
end
plot(t1,y1)
grid on
hold on
z = (0.5*ones(1,length(t1)))';
plot(t1, z)

%---------------------------------------------
% Cau 31
if length(x)>length(h)
   temp_var = x;
   x = h;
   h = temp_var;
end
 
y = zeros(1,(length(x)+length(h)-1));
min_val = [length(x),length(h)];
r = min(min_val)-1;  %minimum of x and h -1
hi = [h,zeros(1,r)];
 
 
for i = 1:length(x)
    temp = x(i).*hi;
   % disp(temp)
    if i>=2
        tempi = [zeros(1,i-1),temp(1:end-i+1)];
    else
        tempi = temp;
    end
    
   y = y + tempi;
end
disp('The convolution of x and h is:' );
disp(y)

%---------------------------------------------
% Cau 1 phan 2
% T = 4, A = 1.5*4 = 6
T = 4;
A = 6;
fn = 2/T;
fm = 32*fn;
t = 0: 1/fm : 8;
y = A * sin(2 * pi * fn * t);
yu = interp(y,8);

figure;
subplot(2,1,1);
stem(t,y);
title('Tín hiệu gốc');

subplot(2,1,2);
stem(yu);  % ham stem cua interp chi can 1 bien
title('tin hieu sau khi lay mau');

%---------------------------------------------
% Cau 2 phan 2
function[code,yq,sqnr]=cau2(a,M)
x = -a :0.1: a;
y = exp(x);
Nb = log2(M);
Amax = max(abs(y));
delta = (2 * Amax)/(M - 1);
Mq = -Amax:delta:Amax;
MI = 0:M-1;

yq = zeros(size(y));
ycode = yq;
for k = 1:M
    ind = find(y > Mq(k)-delta/2 & y <= Mq(k)+delta/2);
    yq(ind)=Mq(k);
    ycode(ind) = MI(k);
end
sqnr = 20*log10(norm(y)/norm(y-yq));
code = dec2bin(ycode,Nb); % theo chuan left-msb
plot(x,y);
hold on;
plot(x,yq);
end

%---------------------------------------------
% cau 5 a phan 2
% chuoi bit a = [1 1 0 0 1 0 1 1 0 0 1 0 1 0 1 1 1 1 1 1 1 0 0 1 0 1 1 0]
% Ma hoa chuoi bit theo dang NRZ
% [t,y,code] = cau5a(a,1e6,256)
function[t,y,code]= cau5a(a,R,Ns,type)
    Tb = 1/R;
    Nb = length(a);
    time = Tb * Nb;
    ts = time/(Ns-1);
    t = 0: ts:time;
    y = zeros(size(t));
    code=[];
    if nargin <= 3
        type = 'unipol';
    end
    for k = 1 : Ns
        n = fix(t(k)/Tb) + 1;    %Ham fix lay phan nguyen, bo phan sau day phay
        if n >= Nb
            n = Nb;
        end
        switch(type)
            case 'unipol'
                y(k) = a(n);
                code(n) = a(n);
            case 'pol'
                y(k) = 2*a(n)-1;
                code(n) = 2*a(n)-1;
        end
    end
    plot(t,y)
end

%---------------------------------------------
% cau 5 b phan 2
% chuoi bit a = [1 1 0 0 1 0 1 1 0 0 1 0 1 0 1 1 1 1 1 1 1 0 0 1 0 1 1 0]
% Ma hoa chuoi bit theo dang RZ
% [t,y,code] = cau5b(a,1e6,256)
function[t,y,code]= cau5b(a,R,Ns)
    Tb = 1/R;
    Nb = length(a);
    time = Tb * Nb;
    ts = time/(Ns-1);
    t = 0: ts:time;
    y = zeros(size(t));
    code=[];
    for k = 1 : Ns
        n = fix(t(k)/Tb) + 1;    %Ham fix lay phan nguyen, bo phan sau day phay
        if n >= Nb
            n = Nb;
        end
        if a(n) == 0
            code(n) = 0;
            if t(k) < Tb*n-Tb/2
                y(k) = -1;
            else
                y(k) = 0;
            end
        else
            code(n) = 1;
            if t(k) < Tb*n-Tb/2
                y(k) = 1;
            else
                y(k) = 0;
            end 
        end    
    end
    plot(t,y)
    grid on
end

%---------------------------------------------
% cau 5 c phan 2
% chuoi bit a = [1 1 0 0 1 0 1 1 0 0 1 0 1 0 1 1 1 1 1 1 1 0 0 1 0 1 1 0]
% Ma hoa chuoi bit theo dang AMI
% [t,y,code] = cau5c(a,1e6,256)
function[t,y,code]= cau5c(a,R,Ns)
    Tb = 1/R;
    Nb = length(a);
    time = Tb * Nb;
    ts = time/(Ns-1);
    t = 0: ts:time;
    y = zeros(size(t));
    code=[];
    s = 1;
    index = 0;
    for k = 1 : Ns
        n = fix(t(k)/Tb) + 1;    %Ham fix lay phan nguyen, bo phan sau day phay
        if n >= Nb
            n = Nb;
        end
        if a(n) == 0
            y(k) = 0;
            code(n) = 0;
        else
            if index ~= n
                index = n;
                s = s + 1;
                if mod(s,2) == 0
                    code(n) = 1;
                    y(k) = 1;
                else
                    code(n) = -1;
                    y(k) = -1;
                end
            else
                y(k) = code(n);
            end
        end    
    end
    plot(t,y)
    grid on
end