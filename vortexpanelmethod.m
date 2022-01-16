%NACA 4 airfoil ;
%linearly varying vortex strength on each panel no sources.



%% defining the airfoil
a = input('Enter the angle of attack in degrees ');
a = a*pi/180;
uf = input('Enter the magnitude of free stream velocity ');
N = input('Enter the no of panels ');
c = input('Enter the chord length ');
M = input('Enter the first digit of four digit naca no ');
P = input('Enter second digit of four digit naca  ');
TT = input('Enter the last two digits of four digit naca ');
syms x t
tmax = (TT*c)/100;
t  = tmax*(2.969 * ((x./c).^(0.5)) - 1.260*(x./c)- 3.516*((x./c).^(2)) + 2.843*((x./c).^(3)) - 1.036*((x./c).^(4)));
%o =  angle of mean chamber line;
m  = M/100;
p = P/10;
z = piecewise(0<=x<=p*c,(m/(p^2))*x.*( (- x./c)+ 2*p ),p*c <= x <= c, (m*(c)/(1-p)^2)*(1-2*p + 2*p*(x/c)- (x/c).^2));
%fplot(x,z);
%axis ([-1 c 0 c]);
diffz = (diff(z,x));
o = atan(diffz);
xu = x - t.*sin(o)./2;
zu =(z) + t.*cos(o)./2;
xl = x + t.*sin(o)./2;
zl =(z) - t.*cos(o)./2;
od = linspace(0,2*pi,N+1);
for i = 1:N+1
    x(i) = (c/2)*(1+cos(od(i)));
end
zun = double(subs(zu));
zln = double(subs(zl));
xln = double(subs(xl));
xun = double(subs(xu));
xln(1) = c; xln(N/2 +1) = 0; xln(N + 1) = c;
xun(1) = c; xun(N/2 +1) = 0; xun(N + 1) = c;
zln(1) = 0; zln(N/2 +1) = 0; zln(N + 1) =0;
zun(1) = 0; zun(N/2 +1) = 0; zun(N + 1) =0;
for i = 1:(N/2 + 1)
    zf(i) = zln(i);
    zf(i + N/2 ) = zun(i + N/2);
end
for i = 1:(N/2 + 1)
    xf(i) = xln(i);
    xf(i + N/2 ) = xun(i + N/2);
end
figure(1)
plot(xf,zf);
daspect([1 1 1]);
axis([-1 2 -1.5 1.5]);
hold on;



%% airfoil coordinates
% xf=naca2412data(:,1);
% fake=naca2412data(:,2);
% for i = 1:201
%    zf(i) = fake(202-i);
% end

figure(2)
fill(xf,zf,'y');
axis([-0.5 c+0.5 -c c]);
daspect([1 1 1]);
hold on;





%% generation of control points at the middle of each panel
for i = 1:N
    xc(i)= ( xf(i) + xf(i+1) )/2;
end
for i = 1:N
    zc(i) = (zf(i) + zf(i+1))/2;
end




%% defining length of each panels
for i = 1:N
    oi(i,1)    = atan2( (zf(i+1)-zf(i)) , (xf(i+1)-xf(i)) );
    RHS(i,1)   = sin(oi(i) - a);
    oj(i,1)    = atan2( (zf(i+1)-zf(i)) , (xf(i+1)-xf(i)) );
    S(i)       = sqrt((zf(i+1)-zf(i))^2 + (xf(i+1)-xf(i))^2);
end




%% defining the matrix to solve for source and vortex strength
for i = 1:N
    for j = 1:N
        if (i == j)
            CN1(i,j) = -1;
            CN2(i,j) = 1 ;
            CT1(i,j) = 0.5*pi;
            CT2(i,j) = 0.5*pi;
        else
            A = -(xc(i) - xf(j))*cos(oj(j)) - (zc(i)- zf(j))*sin(oj(j));
            B = (xc(i) - xf(j))^2 + (zc(i) - zf(j))^2;
            C = sin(oi(i)-oj(j));
            D = cos(oi(i) - oj(j));
            F = log(1 + (S(j)^2 +2*A*S(j))/B);
            E = (xc(i) - xf(j))*sin(oj(j)) - (zc(i)- zf(j))*cos(oj(j));
            G = atan2((E*S(j)),(A*S(j) + B));
            P = ((xc(i) - xf(j)) * sin(oi(i) - 2*oj(j))) + ((zc(i) - zf(j)) * cos(oi(i) - 2*oj(j)));
            Q = ((xc(i) - xf(j)) * cos(oi(i) - 2*oj(j))) - ((zc(i) - zf(j)) * sin(oi(i) - 2*oj(j)));
            
            
            CN2(i,j) = D + ((0.5*Q*F)/S(j)) - ((A*C + D*E)*(G/S(j)));
            CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
            CT2(i,j) = C + ((0.5*P*F)/S(j)) + ((A*D - C*E)*(G/S(j)));
            CT1(i,j) = 0.5*C*F - D*G - CT2(i,j);
        end
    end
end




%% Computation of Influence Coefficients
for i = 1:N
    AN(i,1) = CN1(i,1); % AN = the normal velocity influence coefficients representing the influences of Gama(i) on the normal velocity at the ith control point. 
    AN(i,N+1) = CN2(i,N);
    AT(i,1) = CT1(i,1);
    AT(i,N+1) = CT2(i,N);
    for j = 2:N
        AN(i,j) = CN1(i,j) + CN2(i,j-1);
        AT(i,j) = CT1(i,j) + CT2(i,j-1);
    end
end
AN(N+1,1) = 1;
AN(N+1,N+1) = 1;
for j = 2:N
    AN(N+1,j) = 0;
end
RHS(N+1) = 0;





%% Solve for Gamma and velocity and coefficient of pressure
Gamma = AN\RHS;                                 % Solving for a system of linear equations

for i = 1:N
    Vst(i) = cos(oi(i)-a);
    for j = 1:N+1
        Vst(i) = Vst(i) + AT(i,j)*Gamma(j); %this is the dimensional velocity tangential to the surface of the airfoil  ie Vt/uf 
        Cp(i) = 1 - (Vst(i))^2;
    end
end

figure(5)
grid on;
plot(xc./c,Cp,'or');
set(gca, 'YDir','reverse');
axis([0 c -4 2]);
xlabel('x/c');
ylabel('Coefficient of Pressure,Cp');
legend('Vortexpanel','darmofalfig9.6' );
hold on;


%% Calculation of centre of pressure
delp = zeros(N/2 , 1);
pinf = 101325; %atmospheric pressure
rho  = 1.0;  %density of air at sea level and temp 15C
for i = 1 : N/2 
    p(i)       = (pinf + Cp(i)*0.5*rho*(uf^2))*abs(xf(i+1) - xf(i));
    p(i + N/2) = (pinf + Cp(i  + N/2)*0.5*rho*(uf^2))*abs(xf(i + N/2)-xf(i + 1  + N/2));
end
for i = 1:N/2
    delp(i,1) = p(i) - p(N + 1 - i);
end 
copx = 0; copy = 0;
for i = 1:N/2
    copx = copx + xc(i)*delp(i);
    copy = copy + zc(i)*delp(i);
end

copx = copx /sum(delp);
copy = copy /sum(delp);

xs     = 0.5; %stiffness centre
moment = (xs-copx)*sum(delp);
Cl     = sum(delp)/(0.5*rho*(uf^2)*c);


%% plotting the streamline
[X,Y] = meshgrid(-1:0.01:2,-1.5:0.01:1.5);
U  = uf*cos(a)*ones(301,301);
V  = uf*sin(a)*ones(301,301);
for k = 1:N  %no of panels
    for i = 1:301
        for j = 1:301
            vorx(i,j) = uf*((Gamma(k)+Gamma(k + 1))/2)*S(k)*(((Y(i,j) - zf(k))/((X(i,j) - xf(k))^2 +(Y(i,j) - zf(k))^2)));
            vory(i,j) =-uf*((Gamma(k)+Gamma(k + 1))/2)*S(k)*(((X(i,j) - xf(k))/((X(i,j) - xf(k))^2 +(Y(i,j) - zf(k))^2)));
        end
    end
    U = U + vorx ;
    V = V + vory ;
end

figure(2)
starty = -1.5:0.05:1;
startx = -1*ones(size(starty));
streamline(X,Y,U,V,startx,starty); 

XY = stream2(X,Y,U,V,startx,starty);   
streamparticles(XY,20,'Animate',30,'ParticleAlignment','on','MarkerSize',5);




%% finding the centre of pressure


