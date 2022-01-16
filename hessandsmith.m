%NACA 4 airfoil for even panels;
%hess and smith vortex source method;
%incompressible flow


%% defining the airfoil
a  = input('Enter the angle of attack in degrees ');
a  = a*pi/180;
uf = input('Enter the magnitude of free stream velocity ');
N  = input('Enter the no of panels ');
c  = input('Enter the chord length ');
M  = input('Enter the first digit of four digit naca no ');
P  = input('Enter second digit of four digit naca  ');
TT = input('Enter the last two digits of four digit naca ');
syms x t
tmax = (TT*c)/100;
t    = tmax*(2.969 * ((x./c).^(0.5)) - 1.260*(x./c)- 3.516*((x./c).^(2)) + 2.843*((x./c).^(3)) - 1.036*((x./c).^(4)));
%o   = angle of mean chamber line;
m    = M/100;
p    = P/10;
z    = piecewise(0<=x<=p*c,(m/(p^2))*x.*( (- x./c)+ 2*p ),p*c <= x <= c, (m*(c)/(1-p)^2)*(1-2*p + 2*p*(x/c)- (x/c).^2));
%fplot(x,z);
%axis ([-1 c 0 c]);
diffz = (diff(z,x));
o     = atan(diffz);
xu    = x - t.*sin(o)./2;
zu    = z + t.*cos(o)./2;
xl    = x + t.*sin(o)./2;
zl    = z - t.*cos(o)./2;
x     = zeros(N+1);
od    = linspace(0,2*pi,N+1);
for i = 1:N+1
    x(i) = (c/2)*(1+cos(od(i)));
end
zun = double(subs(zu));
zln = double(subs(zl));
xln = double(subs(xl));
xun = double(subs(xu));
xln(1) = c; xln(N/2 +1) = 0; xln(N + 1) = c;
xun(1) = c; xun(N/2 +1) = 0; xun(N + 1) = c;
zln(1) = 0; zln(N/2 +1) = 0; zln(N + 1) = 0;
zun(1) = 0; zun(N/2 +1) = 0; zun(N + 1) = 0;



%% Combining the upper and lower matrices
xf = zeros(1,N+1); zf = zeros(1,N+1); 
for i = 1:(N/2 + 1)
    zf(i)        = zln(i);
    zf(i + N/2 ) = zun(i + N/2);
end
for i = 1:(N/2 + 1)
    xf(i)        = xln(i);
    xf(i + N/2 ) = xun(i + N/2);
end
figure(1)
grid on;
plot(xf,zf);
daspect([1 1 1]);
axis([-1 2 -1.5 1.5]);
xlabel('x-axis');
ylabel('y-axis');
title('Airfoil Geometry');
hold on;





%% plotting airfoil coordinates


figure(2)
fill(xf,zf,'m');
axis([-0.5 c+0.5 -c c]);
daspect([1 1 1]);
hold on;







%% generation of control points at the middle of each panel
xc = zeros(1,N); zc = zeros(1,N);
for i = 1:N
    xc(i) = ( xf(i) + xf(i+1) )/2;
end
for i = 1:N
    zc(i) = (zf(i) + zf(i+1))/2;
end






%% defining length and angle of orientation of each panels
oi = zeros(1,N);
S  = zeros(1,N);
for i = 1:N
    oi(i,1)    = atan2( (zf(i+1)-zf(i)) , (xf(i+1)-xf(i)) );
    S(i)       = sqrt((zf(i+1)-zf(i))^2 + (xf(i+1)-xf(i))^2);
end







%% normal and tangential vectors
% figure(2)
% quiver(xc,zc,-sin(oi),cos(oi));
% hold on
% quiver(xc,zc,cos(oi),sin(oi));







%% Set variables for the use of making the coefficient matrix.
   A  = zeros(N+1);
   b  = zeros(N+1,1);
   Bn = zeros(N);
   An = zeros(N);
   dx = xf(2:N + 1)-xf(1:N );
   dy = zf(2:N + 1)-zf(1:N );
%    sintheta= dy./sqrt(dx.*dx+dy.*dy);
%    costheta= dx./sqrt(dx.*dx+dy.*dy);
  
   
   dxj  = xc-xf(1:N);
   dxjp = xc-xf(2:N + 1);
   dyj  = zc-zf(1:N);
   dyjp = zc-zf(2:N + 1);
   flog = 0.5*log((dxjp.*dxjp+dyjp.*dyjp)./(dxj.*dxj+dyj.*dyj));
   ftan = atan2(dyjp.*dxj-dxjp.*dyj,dxjp.*dxj+dyjp.*dyj);
   
  
   
   
   
 %% Create influence matrix
 
   for i=1:N
      for j=1:N
      	flog = 0;
      	ftan = pi;
      	if j~=i
		      dxj  = xc(i)-xf(j);
		      dxjp = xc(i)-xf(j+1);
		      dyj  = zc(i)-zf(j);
		      dyjp = zc(i)-zf(j+1);
		      flog = 0.5*log((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj));
      		  ftan = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj);
      	end
      	ctimtj    = cos(oi(i))*cos(oi(j))+sin(oi(i))*sin(oi(j));
      	stimtj    = sin(oi(i))*cos(oi(j))-sin(oi(j))*cos(oi(i));
      	An(i,j)   = (0.5/pi)*(ftan*ctimtj+flog*stimtj);
      	Bn(i,j)   = (0.5/pi)*(flog*ctimtj-ftan*stimtj);
        A(i,j)    = An(i,j);
        A(i,N + 1)= A(i,N + 1)+Bn(i,j);
      	if i == 1||i == N
			%If ith panel touches trailing edge, add contribution to kutta condition
			A(N + 1,j)     = A(N + 1,j) - Bn(i,j);
   	   	    A(N + 1,N + 1) = A(N + 1,N + 1) + A(i,j);
        end
       end
		%Fill in known sides
		b(i) = sin(oi(i))*cos(a)-cos(oi(i))*sin(a);
    end
	b(N + 1)=  -(cos(oi(1))+cos(oi(N)))*cos(a)-(sin(oi(1))+sin(oi(N)))*sin(a);
    Ts      =  A\b;
    At      =  -Bn;
    Bt      =  An;
    
    
    
    
    
%% Use solution to find coefficient of pressure
	Cp  = zeros(N,1);
    vst = zeros(N,1);
for i=1:N
    vst(i)=cos(a)*cos(oi(i))+sin(a)*sin(oi(i));
        for j=1:N
            vst(i) = vst(i)-Bn(i,j)*Ts(j)+Ts(N+1)*A(i,j);
        end
    Cp(i)=1-vst(i).^2;
end
%xflr5naca4412aoa5 = table2array(xflr5naca4412cpaoa5pan204080160);
figure(3)
grid on;
plot(xc./c,Cp,'or');%xflr5naca4412aoa5(:,9),xflr5naca4412aoa5(:,10),'-k');
set(gca, 'YDir','reverse');
axis([0 c -4 2]);
xlabel('x/c');
ylabel('Coefficient of Pressure,Cp');
legend('Hess and smith','xflr5-xfoil' );
hold on;






%% Surface tangential velocities 
vn1(N) = 0;vn2(N) = 0;vt1(N) = 0; vt2(N) = 0;
vn(N) = 0; vt(N) = 0;
for i = 1:N
    for j= 1:N
        vn1(i) = vn1(i) + An(i,j)*Ts(j)*uf;
        vn2(i) = vn2(i) + Bn(i,j)*Ts(N+1)*uf;
        vt1(i) = vt1(i) + At(i,j)*Ts(j)*uf;
        vt2(i) = vt2(i) + Bt(i,j)*Ts(N+1)*uf;
    end
end
for i = 1:N
    vn(i) = vn1(i) +vn2(i) +uf*(-sin(oi(i))*cos(a)+cos(oi(i))*sin(a));
    vt(i) = vt1(i) +vt2(i) +uf*(cos(a)*cos(oi(i))+sin(a)*sin(oi(i)));
end




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
U     = uf*cos(a)*ones(301,301);
V     = uf*sin(a)*ones(301,301);
vorx  = zeros(301,301); vory = zeros(301,301); soux = zeros(301,301); souy = zeros(301,301);
for k = 1:N  %no of panels
    for i = 1:301
        for j = 1:301
            vorx(i,j) = (uf*S(k)/(2*pi))*(Ts(N+1))*(((Y(i,j) - zc(k))/((X(i,j) - xc(k))^2 +(Y(i,j) - zc(k))^2)));
            vory(i,j) =-(uf*S(k)/(2*pi))*(Ts(N+1))*(((X(i,j) - xc(k))/((X(i,j) - xc(k))^2 +(Y(i,j) - zc(k))^2)));
            soux(i,j) = (uf*S(k)/(2*pi))*(Ts(k))*(((X(i,j) - xc(k))/((X(i,j) - xc(k))^2 +(Y(i,j) - zc(k))^2)));
            souy(i,j) = (uf*S(k)/(2*pi))*(Ts(k))*(((Y(i,j) - zc(k))/((X(i,j) - xc(k))^2 +(Y(i,j) - zc(k))^2)));
        end
    end
    U = U + vorx + soux ;
    V = V + vory + souy ;
end

figure(2)
starty = -1.5:0.05:1;
startx = -1*ones(size(starty));
streamline(X,Y,U,V,startx,starty); 

XY = stream2(X,Y,U,V,startx,starty);   
streamparticles(XY,20,'Animate',30,'ParticleAlignment','on','MarkerSize',5,'r');














