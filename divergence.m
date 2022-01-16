%NACA 4 airfoil for even panels;
%hess and smith vortex source method;
clear variables;


%% defining the airfoil
ra  = input('Enter the rigid angle in degrees ');
ra  = ra*pi/180;
uf = input('Enter the magnitude of free stream velocity ');
N  = input('Enter the no of panels ');
c  = input('Enter the chord length ');
M  = input('Enter the first digit of four digit naca no ');
P  = input('Enter second digit of four digit naca  ');
TT = input('Enter the last two digits of four digit naca ');
syms x t
tmax = (TT*c)/100;
t    = tmax*(2.969 * ((x./c).^(0.5)) - 1.260*(x./c)- 3.516*((x./c).^(2)) + 2.843*((x./c).^(3)) - 1.036*((x./c).^(4)));
%o   =  angle of mean chamber line;
m    = M/100;
p    = P/10;
z    = piecewise(0<=x<=p*c,(m/(p^2))*x.*( (- x./c)+ 2*p ),p*c <= x <= c, (m*(c)/(1-p)^2)*(1-2*p + 2*p*(x/c)- (x/c).^2));
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
zln(1) = 0; zln(N/2 +1) = 0; zln(N + 1) =0;
zun(1) = 0; zun(N/2 +1) = 0; zun(N + 1) =0;





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




%% THEORETICAL DIVERGENCE SPEED CALCULATION
K     = 3000; % assumption for  the stiffness coefficient
slope = 2*pi;
copx  = 0;
xs    = c/2;
% Cl    = zeros(1,12);
% for k = -1:10
%     [Cl(k+2),moment,copx] = influmat(uf,k,c,N,xf,zf);  
% end
% a4cl  = -1:10;
% for i = 1:11
%     slope = slope + Cl(i+1) - Cl(i) ;
% end
ec    = xs - copx ;
rho   = 1.0 ; 
qdiv  = K/((ec^2)*slope);
syms vel ta
xd    = 0.5*rho*(vel^2)/qdiv;
ta    = (xd*ra)/(1-xd) -ra;
vel   = 0:2:200;
ta    = double(subs(ta));

%% TEST FOR THE CONVERGENCE
theta     = ra;
for i = 1:60
    [Cl,moment(i),copx] = influmat(uf,theta,c,N,xf,zf);
    angleiter(i)        = 180*theta/pi;
    theta               = ra + moment(i)/K;
end

figure(2)
plot(angleiter,'or');
xlabel('No of iterations');ylabel('Total Twist + rigid angle (deg)');title('Divergence');
%For different freestream velocities
theta = ra;
for i = 0:100
    for j = 1:60
        [Cl,moment(i+1),copx] = influmat(i*2,theta,c,N,xf,zf);
        angleiter           = 180*theta/pi;
        theta               = ra + moment(i+1)/K;
    end
    anglefree(i+1)   = angleiter;
end
figure(1)
plot(vel,anglefree,'o',vel,(ta + 2*ra)*(180/pi));
xlabel('Freestream velocity (ms^-^1)');
ylabel('Total Twist + rigid angle (deg)');
title('Divergence');
legend('Computational Divergence','Theoretical Divergence');
axis([0 60 0 60]);
daspect([1 1 1]);
hold on;


%% function returns the Cl for each alpha.
function [Cl,moment,copx] = influmat(uf,a,c,N,xf,zf) 


%% generation of control points at the middle of each panel

xc = zeros(1,N); zc = zeros(1,N);
for i = 1:N
    xc(i) = (xf(i) + xf(i+1))/2;
end
for i = 1:N
    zc(i) = (zf(i) + zf(i+1))/2;
end




%% defining length and angle of orientation of each panels

oi = zeros(N,1);
S  = zeros(1,N);
for i = 1:N
    oi(i,1)    = atan2( (zf(i+1)-zf(i)) , (xf(i+1)-xf(i)) );
    S(i)       = sqrt((zf(i+1)-zf(i))^2 + (xf(i+1)-xf(i))^2);
end



%% Set variables for the use of making the coefficient matrix.

   A  = zeros(N+1);
   b  = zeros(N+1,1);
   Bn = zeros(N);
   An = zeros(N);
 
  
   
   dxj  = xc-xf(1:N);
   dxjp = xc-xf(2:N + 1);
   dyj  = zc-zf(1:N);
   dyjp = zc-zf(2:N + 1);
   flog = 0.5*log((dxjp.*dxjp+dyjp.*dyjp)./(dxj.*dxj+dyj.*dyj));
   ftan = atan2(dyjp.*dxj-dxjp.*dyj,dxjp.*dxj+dyjp.*dyj);
   
    
   
 %% Creating influence matrix
 
   for i = 1:N
      for j = 1:N
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
    
        
%% Use solution to find coefficient of pressure,drag and lift


	Cp  = zeros(N,1);
    vst = zeros(N,1);
for i=1:N
    vst(i)=cos(a)*cos(oi(i))+sin(a)*sin(oi(i));
        for j=1:N
            vst(i)=vst(i)-Bn(i,j)*Ts(j)+Ts(N+1)*A(i,j);
        end
    Cp(i)=1-vst(i).^2;
end


%% Calculation of centre of pressure


delp = zeros(N/2 , 1);
pr   = zeros(1,N);
pinf = 101325; %atmospheric pressure
rho  = 1.0;  %density of air
for i = 1 : N/2 
    pr(i)       = (pinf + Cp(i)*0.5*rho*(uf^2))*abs(xf(i+1) - xf(i));
    pr(i + N/2) = (pinf + Cp(i  + N/2)*0.5*rho*(uf^2))*abs(xf(i + N/2)-xf(i + 1  + N/2));
end
for i = 1:N/2
    delp(i,1) = pr(i) - pr(N + 1 - i);
end 
copx = 0; copy = 0;
for i = 1:N/2
    copx = copx + xc(i)*delp(i);
    copy = copy + zc(i)*delp(i);
end

copx = copx / sum(delp);
copy = copy / sum(delp);

xs     = c/2; %stiffness centre
moment = (xs-copx)*sum(delp);


%% calculation of coefficient of lift
 Cl = sum(delp)/(0.5*rho*(uf^2)*c);
end

