
m=0.11;
g=-9.8;
d=0.03;
L=1;
J=9.99*10^(-6);
R=0.015;
A=[0 1;0 0];
B=[0;-(m*g*d)/(L*((J/R^2)+m))];
C=[1 0];
D=[0];
sistema=ss(A,B,C,D)
[num,den]=ss2tf (A,B,C,D);
H=tf(num,den)
t=0.01;
Hd=c2d(H,t)
[G,HD,Cd,Dd]=c2dm(A,B,C,D,t,'zoh');
G
HD
Cd
Dd
figure(1)
step(Hd,'--')
axis([0 50 0 300])
figure(2)

step(H,'g')

figure (3)
step(Hd,'--',H,'g')
axis([0 50 0 300])

figure(4)
step(sistema)
