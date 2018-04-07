clear all;
%frequency
f1 = 329.63;
%mass per unit length
m = 0.00008;
%string length
l = 1;
%corresponfing tension
t = m*((2*l*f1)^2);
%string mesh
j = 101;
%damping factor
r = 0.00000008;
%steps in string
dx = l/(j-1);
%steps in recording dispalcement
nskip = ceil((2*f1*(j-1))/8192);
nskip = nskip*6;
%time step
dt = 1/(nskip*8192) ;
%total time
time = 2;
%number of loops
clockmax = ceil (time/dt);
%bridgemass
mb = 0.0009;
%bridge velocity
vb = 0;
%spring constant for bridge
k = 100;
%bridge movement
hb = 0;

% ------Bridge Design ---------
bridgey=0.1;
bridgebr= 0.06;
for i = 1:j
    if i*dx<= (bridgebr)
        bval(i) = -((0.05*bridgey)/(bridgebr^2))*((i*dx)^2)+bridgey;
        maxnode=i;
    end
end



% ----- Initialization --------

%height at j points
h = zeros(1, j);
%velocity at j points
v = zeros(1, j);
%plucking distance
xp = bridgebr+0.04;
%plucking height
hp = bridgey+0.035;
%string initialization
for i = 1:j
    x = (i-1)*dx;
    if x<xp
        h(i)= (hp-bridgey)*(x/xp)+bridgey;
    else
        h(i)= hp*((l-x)/(l-xp));
    end
end

%--------- process -------------

%visualization set up
set(gcf, 'double', 'on')
hplbr = plot(dx*(0:maxnode-1),bval);
hold on
hplbr1 = plot(dx*(0:maxnode-1),bval);
hpl = plot(0:dx:l, h);
axis([-l/10, l + l/10,-1.5*hp,1.5*hp]);
axis manual;
hold off

%process setup
jj = 2:(j-1);
s = zeros(1, ceil(clockmax/nskip));
count = 0;
%process

for clock = 1:clockmax
    %sound generation
    
    v(jj) = v(jj) + (dt/(dx)^2)*(t/m)*(h(jj+1)-2*h(jj)+h(jj-1)) + ( (dt*r)/(m*(dx^2)) )*(v(jj+1)-2*v(jj)+v(jj-1));
    h(jj) = h(jj) + dt*v(jj);
    
    vb = vb + (dt/mb)*(-k*hb);
    hb = hb +dt*vb;
    
    for i = 1:maxnode
        if h(i)<=(bval(i)+hb)
            %tempv = v(i);
            
            %v(i) = v(i)*((m*dx-mb)/(m*dx+mb)) + vb*((2*mb)/(m*dx+mb));
            %vb = vb*((mb-m*dx)/(m*dx+mb)) + tempv*((2*m*dx)/(m*dx+mb));
            
            h(i) = bval(i);
            v(i) = 0;
            
            %for simulation, to be set to normal later
            
            
        end
    end
    

    set(hpl,'ydata',h);
    set(hplbr, 'ydata', bval+hb);
    set(hplbr1, 'ydata', bval);
    %drawnow
    
    

    %recoring
    if mod(clock, nskip) == 0
        count = count + 1;
        s(count) = h(maxnode+5);
    end
    
   
    
end


save playvalue;
soundsc(s);
plot(s);






        
    
