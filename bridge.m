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
j = 201;
%damping factor
r = 0.000000001;
%steps in string
dx = l/(j-1);
%steps in recording dispalcement
nskip = ceil((2*f1*(j-1))/8192);
nskip = nskip*5;
%time step
dt = 1/(nskip*8192) ;
%total time
time = 2;
%number of loops
clockmax = ceil (time/dt);

% ------Bridge Design ---------
bnodes = 0.022*j;
bmax = 0;
bmin = 0.003;
for i = 1:bnodes
    bval(i) = - ((bmin)/((bnodes-1)^2))*((i-1)^2) ;
    maxnode = i;
end



% ----- Initialization --------

%height at j points
h = zeros(1, j);
%velocity at j points
v = zeros(1, j);
%plucking distance
xp = 0.04;
%plucking height
hp = 0.04;
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
plot(dx*(0:bnodes-1),bval);
hold on
hpl = plot(0:dx:l, h);
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
    h(jj) = h(jj) + dt*v(jj);
    
    for i = 1:bnodes
        if h(i)<=(bval(i))
            
            %h(i) = bval(i) + (bval(i)-h(i));
            %v(i) = -v(i);
            
            h(i) = bval(i);
            v(i) = 0;
            
            
            
        end
    end
    
%     for i = 1:bnodes
%         if (h(i)<bval(i))
%             error('h is below the bridge')
%         end
%     end
        
    set(hpl(1),'ydata',h);
    %drawnow
    
    
    %recoring
    if mod(clock, nskip) == 0
        count = count + 1;
        s(count) = h(maxnode+1);
    end
    
   
    
end


save playvalue;
soundsc(s);
plot(s);






        
    
