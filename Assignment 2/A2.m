clear all
r=input('Enter Link Lengths as matrix i.e. [r1 r2 r3 r4]:'); 
% link lengths to be given as an array with [r1 r2 r3 r4] as a specified order
a=r(1);
b=r(2);
c=r(3);
d=r(4);
S=sort(r);
s=S(1);
l=S(4);
p=S(2);
q=S(3);
T24=[];
T23=[];% all possible values of theta2 in degrees
k1=d/a;
k2=d/c;
k3=((a^2)-(b^2)+(c^2)+(d^2))/(2*a*c);
k4=d/b;
k5=((c^2)-(d^2)-(a^2)-(b^2))/(2*a*b);
T4_1=[];
T4_2=[];
T3_1=[];
T3_2=[];

if s+l<=p+q                    
    % Grashof Linkage Condition
    g='Grashof Linkage';
    % For differentiating between different types of Grashof Linkage Mechanisms
    if a==s                    
        gty='Double Crank Mechanism';
    end
    if c==s
        gty='Double Rocker Mechanism'; 
    end
    if b==s || d==s
        gty='Crank Rocker Mechanism';
    end
    if s+l==p+q
        gty='Double Crank/Crank Rocker Mechanism';
    end
else
    g='Non-Grashof Linkage';
    gty='Triple Rocker Mechanism';
end

for t2=0:360
    % Computing different parameters for each value of theta2 in order to
    % determine the positions of different end points of 4 bar linkage
    % mechanisms by calculating theta3 and theta4
    A=cos(degtorad(t2))-k1-(k2*cos(degtorad(t2)))+k3;
    B=-2*sin(degtorad(t2));
    C=k1-((k2+1)*cos(degtorad(t2)))+k3;
    if B^2-(4*A*C)>=0
        % As the equation to determine theta4 is a quadratic equation, there
        % will be two solutions which are t4_1 and t4_2 respectively
        t4_1=2*radtodeg(atan2((-1*B)+sqrt(B^2-(4*A*C)),(2*A)));
        T4_1=[T4_1,t4_1];
        t4_2=2*radtodeg(atan2((-1*B)-sqrt(B^2-(4*A*C)),(2*A)));
        T4_2=[T4_2,t4_2];
        T24=[T24,t2];
    end
    
    D=cos(degtorad(t2))-k1-(k4*cos(degtorad(t2)))+k5;
    E=-2*sin(degtorad(t2));
    F=k1+((k4-1)*cos(degtorad(t2)))+k5;
    
    if E^2-(4*D*F)>=0    
        % As the equation to determine theta3 is a quadratic equation, there
        % will be two solutions which are t3_1 and t3_2 respectively
        t3_1=2*radtodeg(atan2((-1*E)+sqrt(E^2-(4*D*F)),(2*D)));
        T3_1=[T3_1,t3_1];
        t3_2=2*radtodeg(atan2((-1*E)-sqrt(E^2-(4*D*F)),(2*D)));
        T3_2=[T3_2,t3_2];
        T23=[T23,t2];
    end
    
end
% plotting the graphs as requested
figure
plot(T24,T4_1)
title(strcat(g,': ',gty))
xlabel('\theta_2')
ylabel('\theta_4(1)')

figure
plot(T24,T4_2)
title(strcat(g,': ',gty))
xlabel('\theta_2')
ylabel('\theta_4(2)')

figure
plot(T23,T3_1)
title(strcat(g,': ',gty))
xlabel('\theta_2')
ylabel('\theta_3(1)')

figure
plot(T23,T3_2)
title(strcat(g,': ',gty))
xlabel('\theta_2')
ylabel('\theta_3(2)')
