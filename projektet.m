

omega = 19;
n = 200;
xs = 0.2;
ys = 0.2;
aa = 1;
noiselevel = 0.2;                % I denna sektion följer jag stegen i uppgift 3.1-3.5. 
peta = trapets2d(n, omega);       % Att faktiskt hitta källor gör jag i nästa sektion. 
S = @(x,y) aa*S0(x-xs,y-ys);     % Koden som faktiskt är relevant är alltså längre ner.

eta = trapets2d(n, omega); 

[Bound,Sol]=hhsolver(omega,S,n);
N = 10;
i = 1:N;
alfa = 2*i*pi/N;

g = Bound.un;

noise = max(abs(g))*randn(size(g))*noiselevel;
for i = 1:N
    intvector(i) = numeriskcos(Bound, omega, alfa(i), noise);
end






plot(alfa, intvector)
hold on
plot(alfa, eta*vc(xs, ys, alfa, omega))
%

delta = 1;

numsinpihalv = numerisksin(Bound, omega, pi/2, noise);
numsinnoll = numerisksin(Bound, omega, 0, noise);
numcosnoll = numeriskcos(Bound, omega, 0, noise);

h = 0.0001;
derivIcpihalv = (numeriskcos(Bound, omega, pi/2+h, noise)-numeriskcos(Bound, omega, pi/2-h, noise))/(2*h);
derivIcnoll = (numeriskcos(Bound, omega, h, noise)-numeriskcos(Bound, omega, -h, noise))/(2*h);

x0 = derivIcpihalv/(omega*numsinpihalv);
y0 = -derivIcnoll/(omega*numsinnoll);

aa = 1/eta * sqrt(numcosnoll^2+numsinnoll^2);


while max(abs(delta)) > 10^(-10)

    delta = Jacobi(aa, x0, y0, alfa', omega, eta) \ -Ic(aa, x0, y0, alfa', omega, eta, intvector');
    aa = aa + delta(1);
    x0 = x0 + delta(2);
    y0 = y0 + delta(3);
    
end

disp("aa: " + aa)
disp("x0: " + x0)    
disp("y0: " + y0)

%%

omega = 19;
n = 200;
xs = 0.2;
ys = 0.2;
aa = 1;
                
eta = trapets2d(n, omega);       
S = @(x,y) aa*S0(x-xs,y-ys);     


[Bound,Sol]=hhsolver(omega,S,n);
N = 10;
i = 1:N;
alfa = 2*i*pi/N;

h = 0.0001;

g = Bound.un;
basenoise = max(abs(g))*randn(size(g));
for noiselevel = 0.01:0.01:0.5

    
    noise = basenoise*noiselevel;

    numsinpihalv = numerisksin(Bound, omega, pi/2, noise);
    numsinnoll = numerisksin(Bound, omega, 0, noise);
    numcosnoll = numeriskcos(Bound, omega, 0, noise);
    
    derivIcpihalv = (numeriskcos(Bound, omega, pi/2+h, noise)-numeriskcos(Bound, omega, pi/2-h, noise))/(2*h);
    derivIcnoll = (numeriskcos(Bound, omega, h, noise)-numeriskcos(Bound, omega, -h, noise))/(2*h);
    

    for i = 1:N
        intvector(i) = numeriskcos(Bound, omega, alfa(i), noise);
    end

    x0 = derivIcpihalv/(omega*numsinpihalv);
    y0 = -derivIcnoll/(omega*numsinnoll);
    aa = 1/eta * sqrt(numcosnoll^2+numsinnoll^2);
    
    delta = 1;
    while max(abs(delta)) > 10^(-10)
    
        delta = Jacobi(aa, x0, y0, alfa', omega, eta) \ -Ic(aa, x0, y0, alfa', omega, eta, intvector');
        aa = aa + delta(1);
        x0 = x0 + delta(2);
        y0 = y0 + delta(3);
        
    end
    
    noroundoff = round(noiselevel*100);

    x0brastartg(noroundoff) = x0;
    y0brastartg(noroundoff) = y0;
    aabrastartg(noroundoff) = aa;

    
    x0 = 0.3;
    y0 = 0.1;
    aa = 1;
    delta = 1;
    while max(abs(delta)) > 10^(-10)
    
        delta = Jacobi(aa, x0, y0, alfa', omega, eta) \ -Ic(aa, x0, y0, alfa', omega, eta, intvector');
        aa = aa + delta(1);
        x0 = x0 + delta(2);
        y0 = y0 + delta(3);
        
    end

    x0vanligstartg(noroundoff) = x0;
    y0vanligstartg(noroundoff) = y0;
    aavanligstartg(noroundoff) = aa;
    
end

figure(1)
plot(0.01:0.01:0.5, x0brastartg)


figure(2)
plot(0.01:0.01:0.5, y0brastartg)


figure(3)
plot(0.01:0.01:0.5, aabrastartg)


figure(4)
plot(0.01:0.01:0.5, x0vanligstartg)


figure(5)
plot(0.01:0.01:0.5, y0vanligstartg)


figure(6)
plot(0.01:0.01:0.5, aavanligstartg)






%%

n = 200;
noiselevel = 0;
format long
h = 0.000001;
N = 100;
i = 1:N;
alfa = 2*i*pi/N;    % Vi har en massa vinklar som vi vill testa mot

for k = 1:5
    load("source"+k+".mat")   % Vi försöker hitta alla fem källor.
    
  
    eta = trapets2d(n, omega);
    
   
    numsinpihalv = numerisksin(B, omega, pi/2, 0);      % Enskilda funktioner för att integrera över v_c*g och v_s*g
    numsinnoll = numerisksin(B, omega, 0, 0);
    
    derivIcpihalv = (numeriskcos(B, omega, pi/2+h, 0)-numeriskcos(B, omega, pi/2-h, 0))/(2*h);
    derivIcnoll = (numeriskcos(B, omega, h, 0)-numeriskcos(B, omega, -h, 0))/(2*h);
    
    x0 = derivIcpihalv/(omega*numsinpihalv);
    y0 = -derivIcnoll/(omega*numsinnoll);
      %Vi tar startgissning för vår startgissning i princip
    
    
    
    intCvektor = numeriskcos(B, omega, alfa, noiselevel);
    intSvektor = numerisksin(B, omega, alfa, noiselevel);
    Cderiv = (numeriskcos(B, omega, alfa+h, 0)-numeriskcos(B, omega, alfa-h, 0))/(2*h);   %C och S i namnen understryker skillnad mellan vc och vs.
    
    delta = 1;
    
    while max(delta) > 10^(-10)
    
        delta = battreJ(alfa', intSvektor', omega) \ -battrestartF(x0, y0, alfa', intSvektor', Cderiv', omega);
        x0 = x0 + delta(1);
        y0 = y0 + delta(2);
                                        % Här optimerar vi startgissningen
                                        % för gauss-newton av
                                        % cosinus-planvågen.
    end
    
    aa = 1/eta * sum(sqrt(intSvektor.^2+intSvektor.^2))/N;    %Medelvärde
    
   
    
    for i = 1:N
        intvector(i) = numeriskcos(B, omega, alfa(i), noiselevel);   %Numerisk integration med simpsons för att kunna anpassa x, y och a-värden
    end
    
    delta = 1;
    while max(abs(delta)) > 10^(-10)
        
        delta = Jacobi(aa, x0, y0, alfa', omega, eta) \ -Ic(aa, x0, y0, alfa', omega, eta, intvector');
        aa = aa + delta(1);
        x0 = x0 + delta(2);
        y0 = y0 + delta(3);
                                        % Gauss-Newton av cosinusplanvågen.
    end
    disp("Källa " + k + ":")
    disp("x = " + x0)
    disp("y = " + y0)
    disp("a = " + aa)

    figure(k)
    plot(alfa, intvector)
    hold on
    plot(alfa, aa*eta*vc(x0, y0, alfa, omega))  
    figure(k+5)
    plot(alfa, intvector-aa*eta*vc(x0, y0, alfa, omega))
    
    
end

%%



omega = 30;      %I denna sektion hittar vi den optimala placeringen av TVn
n = 1000;
N = 100;
i = 1:N;
xs = 0.6+0.4*i/N;
ys = 0.6;


for j = 1:N
    S = @(x,y) S0(x-xs(j),y-ys);
    [Bound,Sol]=hhsolver(omega,S,200);
    
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    A(j) = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));    %Här bestämmervi ungefär var mimimum finns. (Typ mellan 0.65 ovh 0.7)
end

plot(xs, A)

[varde, index] = min(A);

bestx = 0.6 + 0.4*index/N;
b = ceil(bestx*20)/20;
a = floor(bestx*20)/20;

alfa = (sqrt(5)-1)/2;
x1 = b - (b-a)*alfa;
x2 = a + (b-a)*alfa;

S = @(x,y) S0(x-x1,y-ys);
[Bound,Sol]=hhsolver(omega,S,n);
    
w = find(Sol.x<=0.25 & Sol.y>=0.5);
f1 = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));

S = @(x,y) S0(x-x2,y-ys);
[Bound,Sol]=hhsolver(omega,S,n);
    
w = find(Sol.x<=0.25 & Sol.y>=0.5);
f2 = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));

while abs(b-a) > 10^(-5)   %För att vara extra säker. Det är ända bara ett par extra iterationer.
    if f1 < f2
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = b - (b-a)*alfa;

        S = @(x,y) S0(x-x1,y-ys);
        [Bound,Sol]=hhsolver(omega,S,n);
    
        w = find(Sol.x<=0.25 & Sol.y>=0.5);
        f1 = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
    else
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = a + (b-a)*alfa;

        S = @(x,y) S0(x-x2,y-ys);
        [Bound,Sol]=hhsolver(omega,S,n);
    
        w = find(Sol.x<=0.25 & Sol.y>=0.5);
        f2 = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
    end
end

minimum = (a + b)/2;
disp(minimum)


S = @(x,y) S0(x-minimum,y-ys);      %Vi kör solver igen bara för att visa resultatet.
[Bound,Sol]=hhsolver(omega,S,1000);
w = find(Sol.x<=0.25 & Sol.y>=0.5);
minimalA = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));

disp(minimalA)

figure(2)
mesh(Sol.x,Sol.y,Sol.u)



S = @(x,y) S0(x-0.8,y-ys);  %Plot för en godtycklig källa vid väggen.
[B,Sol]=hhsolver(omega,S,n);


figure(3)
mesh(Sol.x,Sol.y,Sol.u)


options = optimset('Display', 'iter');

for xkoordinater = 0.6:0.05:1
    for ykoordinater = 0.4:0.05:0.6
        [vektor, varde] = fminsearch(@FUN, [xkoordinater, ykoordinater], options);
        if vektor(2)>0.6 || vektor(1) > 1
            varde = 1;      %Ibland erhålls en minpunkt utanför rummet. Dessa räknas bort med denna rad.
        end
        matrisfunktionsvarden(round(xkoordinater*20)-11, round(ykoordinater*20)-7) = varde;
    end
end

%Det finns många lokala minimum i två variabler, därför testar vi många
%olika startgissningar för att se vilket som är bäst. 


[minstavarde, linjarIndex] = min(matrisfunktionsvarden(:));

[radIndex, colIndex] = ind2sub(size(matrisfunktionsvarden), linjarIndex);
radvarde = (radIndex+11)/20;
colvarde = (colIndex+7)/20;
[vektor, varde] = fminsearch(@FUN2, [radvarde, colvarde]);
%Vi plockar fram bästa minimum och sen gör optimerar vi med den punkten som
%startgissning med fun2 som använder n = 1000 i solver.

disp("x-värde: " + vektor(1))
disp("y-värde: " + vektor(2))
disp("Ljudnivå: " + minstavarde)




S = @(x,y) S0(x-vektor(1),y-vektor(2));   %Vi anropar igen med lägre n för att plottningen inte ska ta för lång tid. 
[B,Sol]=hhsolver(omega,S,n);

figure(4)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off

figure(5)
mesh(Sol.x,Sol.y,Sol.u)






%%

function  fun = FUN(X0)
    omega = 30;
    n = 200;
    S = @(x,y) S0(x-X0(1),y-X0(2));
    [Bound,Sol]=hhsolver(omega,S,n);
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    minimalA = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));      %Funktion som beräknar kvoten A givet en källas position

    fun = minimalA;
end

function  fun = FUN2(X0)
    omega = 30;
    n = 1000;
    S = @(x,y) S0(x-X0(1),y-X0(2));                  %Samma med större n.
    [Bound,Sol]=hhsolver(omega,S,n);
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    minimalA = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));

    fun = minimalA;
end




function F = battrestartF(x0, y0, alfa, intS, Cderiv, omega)     %Implementation av Gauss-Newton för att få 
    F = omega*(x0.*sin(alfa)-y0.*cos(alfa)).*intS-Cderiv;        %bättre startgissningar till Gauss-Newton sen.

end

function J = battreJ(alfa, intS, omega)
    J = [omega.*sin(alfa).*intS, -omega.*cos(alfa).*intS];    %Vi har deriverat map x0 och y0.
end


function F = Ic(a, x0, y0, alfa, omega, eta, intvector)
    F = a * eta * cos(omega*(x0*cos(alfa)+y0*sin(alfa))) - intvector;  %Vi utnyttjar Ic = a\eta *vc. Denna funktion är noll vid rätt x, y och a
end

function J = Jacobi(a, x0, y0, alfa, omega, eta)
    J = [eta .* cos(omega*(x0*cos(alfa)+y0*sin(alfa))), -omega .* cos(alfa) .* a .* eta .* sin(omega*(x0*cos(alfa)+y0*sin(alfa))), -omega .* sin(alfa) .* a .* eta .* sin(omega*(x0*cos(alfa)+y0*sin(alfa)))];

    %Derivering med avseende på a, x och y.
end


function cosint = numeriskcos(Bound, omega, alfa, noise)
    g = Bound.un + noise;
   
    
    sum = vc(Bound.x(1), Bound.y(1), alfa, omega)*g(1) + vc(Bound.x(end), Bound.y(end), alfa, omega)*g(end);
    for j = 2:length(Bound.s)-1
        if mod(j, 2) == 0  
            sum = sum + 2*vc(Bound.x(j), Bound.y(j), alfa, omega)*g(j);
        else
            sum = sum + 4*vc(Bound.x(j), Bound.y(j), alfa, omega)*g(j);
        end
    end
    cosint = sum * (Bound.s(2)-Bound.s(1))/3;
   
end                                                                    % I båda dessa integrerar vi med hjälp av simpsons formel.

function sinint = numerisksin(Bound, omega, alfa, noise)
    g = Bound.un + noise;
   
    
    sum = vs(Bound.x(1), Bound.y(1), alfa, omega)*g(1) + vs(Bound.x(end), Bound.y(end), alfa, omega)*g(end);
    for j = 2:length(Bound.s)-1
        if mod(j, 2) == 0  
            sum = sum + 2*vs(Bound.x(j), Bound.y(j), alfa, omega)*g(j);
        else
            sum = sum + 4*vs(Bound.x(j), Bound.y(j), alfa, omega)*g(j);
        end
    end
    sinint = sum * (Bound.s(2)-Bound.s(1))/3;
   
end



function kallstyrka = S0(x, y)
    kallstyrka = cos(20.*sqrt(x.^2+y.^2)).*exp(-1000*(x.^2+y.^2));
end

function vc = vc(x, y, alfa, omega)
    vc = cos(omega*(x*cos(alfa)+y*sin(alfa)));
end

function vs = vs(x, y, alfa, omega)
    vs = sin(omega*(x*cos(alfa)+y*sin(alfa)));
end



function eta = trapets2d(n, omega)
    edge = 0.2;     %Nästan noll vid 0.2

    sum = S0(edge, edge) * cos(omega*edge);    %Hörnen, *4/4

    for i = 1:n-1
        sum = sum + S0(edge, -edge + 2*edge*i/n)*cos(omega*edge);   
        sum = sum + S0(-edge + 2*edge*i/n, edge)*cos(omega*(-edge + 2*edge*i/n));    %Kanter
    end


    for i = 1:n-1
        for j = 1:n-1
            sum = sum + S0(-edge + 2*edge*i/n, -edge + 2*edge*j/n)*cos(omega*(-edge + 2*edge*i/n));   %Resten
        end
    end

    eta = sum*(2*edge/n)^2;    %2*edge/n blir motsvarande h.
end


