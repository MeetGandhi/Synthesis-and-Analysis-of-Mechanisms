aa=input('Input Array [R rho1 rho2 rho3 alpha2 alpha3] anles in deg:');
R=aa(1);
rho1=deg2rad(aa(2));
rho2=deg2rad(aa(3));
rho3=deg2rad(aa(4));
alpha2=deg2rad(aa(5));
alpha3=deg2rad(aa(6));

A=exp(1i*alpha2)*R*exp(1i*rho3)-exp(1i*alpha3)*R*exp(1i*rho2);
B=exp(1i*alpha3)*R*exp(1i*rho1)-R*exp(1i*rho3);
C=-exp(1i*alpha2)*R*exp(1i*rho1)+R*exp(1i*rho2);

a=conj(A)*B;
b=conj(A)*A+conj(B)*B-conj(C)*C;
c=A*conj(B);

t1=(-b+sqrt(b^2-4*a*c))/2*a;
t2=(-b-sqrt(b^2-4*a*c))/2*a;

s1=(-1*A+(B*t1))/C;
s2=(-1*A+(B*t2))/C;

beta2_1=rad2deg(angle(t1));
beta2_2=rad2deg(angle(t2));

beta3_1=rad2deg(angle(s1));
beta3_2=rad2deg(angle(s2));

z_1=abs(((R*exp(1i*rho2))-(R*exp(1i*rho1)*exp(1i*beta2_1)))/(exp(1i*alpha2)-(exp(1i*beta2_1))));
phi_1=rad2deg(angle(((R*exp(1i*rho2))-(R*exp(1i*rho1)*exp(1i*beta2_1)))/(exp(1i*alpha2)-(exp(1i*beta2_1)))));

w_1=abs(R*exp(1i*rho1)-z_1*exp(1i*phi_1));
theta_1=rad2deg(angle(R*exp(1i*rho1)-z_1*exp(1i*phi_1)));

z_2=abs(((R*exp(1i*rho2))-(R*exp(1i*rho1)*exp(1i*beta2_1)))/(exp(1i*alpha2)-(exp(1i*beta2_2))));
phi_2=rad2deg(angle(((R*exp(1i*rho2))-(R*exp(1i*rho1)*exp(1i*beta2_1)))/(exp(1i*alpha2)-(exp(1i*beta2_2)))));

w_2=abs(R*exp(1i*rho1)-z_2*exp(1i*phi_2));
theta_2=rad2deg(angle(R*exp(1i*rho1)-z_2*exp(1i*phi_2)));