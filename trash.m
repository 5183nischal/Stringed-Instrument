%frequency
f1 = 329.63;
%mass per unit length
m = 0.00008;
%string length
l = 0.66;
%corresponfing tension
t = m*((2*l*f1)^2);
%string mesh
j = 101;
%damping factor
r = 0.0000002;
%resonance factor
k = 23000;
%steps in string
dx = l/(j-1);
%steps in recording dispalcement
nskip = ceil((2*f1*(j-1))/8192);
%time step
dt = 1/(2*nskip*8192) ;
%total time
time = 1;
%number of loops
clockmax = ceil (time/dt);
%number of strings
nstrings = 6;

% ----- Initialization --------

%height at j points
h = zeros(1, j);
%velocity at j points
v = zeros(1, j);
%plucking distance
xp = 0.01;
%plucking height
hp = 0.001;
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
hpl = plot(0:dx:l, h);
axis([-l/10, l + l/10,-0.15*hp,0.15*hp]);
axis manual;

%process setup
jj = 2:(j-1);
s = zeros(1, ceil(clockmax/nskip));
count = 0;

%process

for clock = 1:clockmax
    %sound generation
    
    
    
    v(jj) = v(jj) + (dt/(dx)^2)*(t/m)*(h(jj+1)-2*h(jj)+h(jj-1)); 
    h(jj) = h(jj) + dt*v(jj);
    
    %set(hpl,'ydata',h);
    %drawnow
    
    %recoring
    if mod(clock, nskip) == 0
        count = count + 1;
        s(count) = h(2);
    end
end


soundsc(s);
%srate = (clock/nskip)/(time*60);
filename = 'firste.wav';
s = s*500;
audiowrite(filename,s,9700);
plot(s);





        
    
