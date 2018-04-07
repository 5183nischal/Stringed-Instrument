clear all;
%frequency
f1 = 329.63;
%mass per unit length
m = 0.00008;
%string length
l = 0.66;
%corresponfing tension
t = m*((2*l*f1)^2);
%string mesh
j = 801;
%damping factor
r = 0.000000001;
%steps in string
dx = l/(j-1);
%steps in recording dispalcement
nskip = ceil((2*f1*(j-1))/8192);
%time step
dt = 1/(nskip*8192) ;
%total time
time = 4;
%number of loops
clockmax = ceil (time/dt);
%bridge velocity
vb = 0;
%bridge height
hb = 0;
%bridge mass
mb = 0.005;
%spring constant for bridge
k = 2;


% ------Bridge Design ---------
bnodes = 0.007*j;
bmax = 0;
bmin = 0.003;
for i = 1:bnodes
    bval(i) = - ((bmin)/((bnodes-1)^2))*((i-1)^2) ;
end



% ----- Initialization --------

%height at j points
h = zeros(1, j);
h1 = zeros(1, j);
%velocity at j points
v = zeros(1, j);
v1 = zeros(1, j);
%plucking distance
xp = 0.04;
%plucking height
hp = 0.035;
%string initialization
for i = 1:j
    x = (i-1)*dx;
    if x<xp
        h(i)= hp*(x/xp);
    else
        h(i)= hp*((l-x)/(l-xp));
    end
end

%--------- process -------------

%visualization set up
set(gcf, 'double', 'on')
hplbr = plot(dx*(0:bnodes-1),bval);
hold on
hpl(1) = plot(0:dx:l, h);
hpl(2) = plot(0:dx:l, h);
axis([-l/10, l + l/10,-1.5*hp,1.5*hp]);
axis manual;
hold off

%process setup
jj = 2:(j-1);
s = zeros(1, ceil(clockmax/nskip));
count = 0;
rec = 1;
%process

for clock = 1:clockmax
    %sound generation
    
    
    
    v(jj) = v(jj) + (dt/(dx)^2)*(t/m)*(h(jj+1)-2*h(jj)+h(jj-1)) + ( (dt*r)/(m*(dx^2)) )*(v(jj+1)-2*v(jj)+v(jj-1));
    v1(jj) = v1(jj) + (dt/(dx)^2)*(t/m)*(h1(jj+1)-2*h1(jj)+h1(jj-1)) + ( (dt*r)/(m*(dx^2)) )*(v1(jj+1)-2*v1(jj)+v1(jj-1));
    
    vb = vb + (dt/mb)*( -k*hb + t*((h(rec+1)-h(rec))/dx) );
    hb = hb +dt*vb;
    
    h(jj) = h(jj) + dt*v(jj);
    h1(jj) = h1(jj) + 1.2*dt*v1(jj);
    
    
    h1(1) = hb;
    for i = 1:rec
        h(i) = h(i) + hb;
    end
    
    for i = 1:bnodes
        bval(i) = bval(i)+hb;
        if h(i)<=(bval(i))
            
            %h(i) = bval(i) + (bval(i)-h(i));
            %v(i) = -v(i);
            
            h(i) = bval(i);
            rec = i;
            
        end
    end
    
%     for i = 1:bnodes
%         if (h(i)<bval(i))
%             error('h is below the bridge')
%         end
%     end
        
    set(hpl(1),'ydata',h);
    set(hpl(2),'ydata',h1);
    set(hplbr,'ydata',bval) 
    drawnow
    
    
    %recoring
    if mod(clock, nskip) == 0
        count = count + 1;
        s(count) = h(20) + h1(2);
    end
    
    %restoring values
    
    for i = 1:rec
        h(i) = h(i) - hb;
    end
    
    for i = 1:bnodes
        bval(i) = bval(i)-hb;
    end
    
end


save playvalue;
soundsc(s);
plot(s);






        
    
