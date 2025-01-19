data=readmatrix('scope_1.csv', 'NumHeaderLines',2)
t=data(:,1)
u=data(:,2)
y1=data(:,3)
y2=data(:,4)
figure
plot(t, [u, y1,y2]), grid, shg
legend('Semnal de intrare u', 'Semnal de iesire y1', 'Semnal de iesire y2');
xlabel('Timp [s]');
ylabel('Amplitudinea(u, y1, y2)');
title('Date măsurate');

%% Metoda neparametrica pt y1

K=mean(y1)/mean(u) %factor propotionalitate 
K1=(mean(y1(1:995))/(mean(u(1:995))))
u_max_idx=194;
u_min_idx=214;
y1_max_idx=201;
y1_min_idx=221;

A_iesire=y1(y1_max_idx)-y1(y1_min_idx);
A_intare=u(u_max_idx)-u(u_min_idx);
Mr=A_iesire/A_intare  %Factorul de amplificare in rezonanta

zeta = sqrt((Mr-sqrt(Mr^2-1))/(2*Mr))

Trez=2*(t(y1_min_idx)-t(y1_max_idx)) %perioada de rezonanta
wr=(2*pi)/Trez    %pulsatia de rezonanta 

wn=wr/sqrt(1-2*zeta^2)  %pulsatia naturala

Num=K*wn^2;
Den=[1, 2*zeta*wn, wn^2];
H=tf(Num, Den)

A=[0 1; -wn^2 -2*zeta*wn];
B=[0; K*wn^2];
C=[1 0];
D=0;

sys=ss(A,B,C,D)
figure
ysim=lsim(sys, u, t, [y1(1), (y1(2)-y1(1))/(t(2)-t(1))]);  %valori initiale
hold on
plot(t, [y1 ysim]), grid
title('Suprapunerea dintre semnalul masurat si semnalul identificat')
xlabel('Timp[s]')
ylabel('Amplitudine')
legend('Iesire', 'Iesire identificata')

%eroarea medie patratica
J=norm(y1-ysim)/sqrt(length(y1))
Empn=norm(y1-ysim)/norm(y1-mean(y1))*100

%% Metode parametrice y1  %obtiunute 2 modele pt y1 unul validat prin autocorelatie unul validat prin intercorelatie 

dt=t(2)-t(1)  % pasul de esantionare
date_y1=iddata(y1, u, dt)

%validare orin autocorelatie 
armax_y1=armax(date_y1, [2 1 2 1])

M_armax_y1=idpoly(armax_y1)
Hz_y1armax=tf(M_armax_y1.B, M_armax_y1.A, dt, 'variable', 'z^-1')
Hs_y1armax=minreal(zpk(d2c(Hz_y1armax, 'zoh')))
[A,B,C,D]=tf2ss(M_armax_y1.B, M_armax_y1.A)
Y1c=dlsim(A,B,C,D, u, [y1(1), y1(2)-A(1,1)*y1(1)-C(1,1)*u(1)])
plot(t, [y1 Y1c])

figure; resid(armax_y1, date_y1, 5)
title('Validare autocorelatie pentru y1')
figure; compare(date_y1, armax_y1)
title('Comparația între semnalul y1 și modelul ARMAX')

%validare intercorelatie 
oe1=oe(date_y1, [1 2 1]);
Moe=idpoly(oe1)
Hz_y1oe=tf(Moe.B, Moe.F, dt, 'variable', 'z^-1')  %in z 
Hs_y1oe=minreal(zpk(d2c(Hz_y1oe, 'zoh')))

figure; resid(date_y1, oe1, 5);
title('Validare intercorelatie pentru y1')
figure; compare(date_y1, oe1);
title('Comparația între semnalul y1 și modelul oe')

%% Metoda paramaterica pt y2
date_y2=iddata(y2, u, dt)

armax_y2=armax(date_y2, [2 2 2 1])

Marmax_y2=idpoly(armax_y2)
Hz_y2=tf(Marmax_y2.B, Marmax_y2.A, dt, 'variable', 'z^-1')
Hs_y2=minreal(zpk(d2c(Hz_y2, 'zoh')))

figure; resid(armax_y2, date_y2, 5)
title('Validare autocorelatie pentru y2')
figure; compare(date_y2, armax_y2)
title('Comparația între semnalul y2 și modelul ARMAX')

 % validare intercorelatie 
iv2=iv4(date_y2, [2 2 0]);

Miv=idpoly(iv2)
Hz_y2iv=tf(Miv.B, Miv.A, dt,'variable', 'z^-1')  %in z 
Hs_y2iv=d2c(Hz_y2iv, 'zoh')

figure; resid(date_y2, iv2, 5);
title('Validare intercorelatie pentru y2')
figure; compare(date_y2, iv2);
title('Comparația între semnalul y2 și modelul iv')
%% Estimarea raspunsului in frecventa
w1=pi/(t(121)-t(91))
w2=pi/(t(148)-t(121))
w3=pi/(t(172)-t(148))

m1=(y1(96)-y1(127))/(u(91)-u(122))
m2=(y1(127)-y1(154))/(u(121)-u(148)) 
m3=(y1(154)-y1(178))/(u(148)-u(172))

ph1=(t(91)-t(96))*w1*180/pi
ph2=(t(121)-t(127))*w2*180/pi
ph3=(t(148)-t(154))*w3*180/pi

ph_r=(t(194)-t(201))*wr*180/pi

w6=pi/(t(495)-t(482))
m6=(y1(491)-y1(502))/(u(482)-u(495))
ph6=(t(482)-t(491))*w6*180/pi

w7=pi/(t(811)-t(801)) 
w8=pi/(t(821)-t(811)) 

m7=(y1(809)-y1(819))/(u(801)-u(811))
m8=(y1(819)-y1(828))/(u(811)-u(821)) 

ph7=(t(801)-t(809))*w7*180/pi
ph8=(t(811)-t(819))*w8*180/pi

w=logspace(3, 4)
[M,ph]=bode(H,w)
figure
subplot(211)
semilogx([w1 w2 w3 wr w6 w7 w8], 20*log10([m1 m2 m3 Mr m6 m7 m8]), 'o'); grid, shg
hold on 
semilogx(w, squeeze(20*log10(M)))
title('Diagrama de modul')

subplot(212)
semilogx([w1 w2 w3 wr w6 w7 w8],[ph1 ph2 ph3 ph_r ph6 ph7 ph8], 'o'); grid, shg
hold on 
semilogx(w, squeeze(ph))
title('Diagrama de faza')

panta=20*log10(m7/m6)*10*w6/w7
