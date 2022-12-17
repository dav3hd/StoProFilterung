% Stochastische Prozesse - Beleg 3 - Filterung

clear; close all; clc; format long

y = readmatrix("aufgabe9.txt"); % Input Zeitreihe

% Aufgabe 1

m1 = 32;
m2 = 36;
% m3 = 3; % Soll nur zur Anschauung dienen
pas = pasdrei(m2);
b1 = pas(m1+1,:);

b2 = pas(m2+1,:);

% pas3 = pasdrei(m3);
% b3 = pas3(m3+1,:);

clearvars pas pas3

%
%%
% Durchlasscharakteristik (Aufgabe 3)

% bei Abtastweite dt = 2[h] --> Frequenz domega = 1/dt => 1/2[1/h]
% (180/2=90[°/h] == pi/2[rad/h])
% --> Nyquist-Frq.: 1/(2*dt) => 1/4[1/h]
omega_g = 90;
domega = 360/2000;

n(:,1) = 1:round(omega_g/domega);
mu(:,1) = n.*domega;

G_m1(:,1) = (cos((pi/2)*(mu/omega_g))).^m1;
G_m2(:,1) = (cos((pi/2)*(mu/omega_g))).^m2;
% G_m3(:,1) = (cos((pi/2)*(mu/omega_g))).^m3;

figure(1)
hold on
title("Durchlasscharakteristiken Binomialfilter (Tiefpass)")
plot(mu,G_m1,'b-')
plot(mu,G_m2,'r-')
% plot(mu,G_m3,'g-')
legend("m1","m2","m3")
grid on
ylabel("Durchlasscharakteristik G")
xlabel("Frequenz [°/h]")
hold off

saveas(1,"images/1_G_Tiefpass.png")

G_hp1 = 1-G_m1;
G_hp2 = 1-G_m2;
% G_hp3 = 1-G_m3;

figure(2)
hold on
title("Durchlasscharakteristik Hochpassfilter")
plot(mu,G_hp1,'b-')
plot(mu,G_hp2,'r-')
% plot(mu,G_hp3,'g-')
legend("m1","m2","m3")
grid on
ylabel("Durchlasscharakteristik G")
xlabel("Frequenz [°/h]")
hold off

saveas(2,"images/2_G_Hochpass.png")

G_bp1 = G_m1 - G_m2;
% G_bp2 = G_m1 - G_m3;
% G_bp3 = G_m3 - G_m2;

figure(3)
hold on
title("Durchlasscharakteristik Bandpassfilter")
plot(mu,G_bp1,'b-')
% plot(mu,G_bp2,'r-')
% plot(mu,G_bp3,'g-')
legend("m1-m2","m1-m3","m3-m2")
grid on
ylabel("Durchlasscharakteristik G")
xlabel("Frequenz [°/h]")
hold off

saveas(3,"images/3_G_Bandpass.png")


%
%%
% (a) Tiefpassfilter der Ordnung m = 32


y_m1 = [];


n = 1;
for r = m1:length(y)-m1
    y_m1(n,1) = 0;
    for k = 1:m1+1

        y_m1(n,1) = y_m1(n,1) + b1(k) * y(r+k-m1/2,3);

    end
    n = n+1;
end

y_m1(:,1) = y_m1(:,1) * 1/(2^m1);


%
%%
% Hochpassfilter der Ordnung m = 36

y_m2 = [];


n = 1;

for r = 36:length(y)-36
    y_m2(n,1) = 0;
    for k = 1:m2+1

        y_m2(n,1) = y_m2(n,1) + b2(k) * y(r+k-m2/2,3);

    end
    n = n+1;
end

y_m2(:,1) = y(1:length(y_m2),3) - (y_m2(:,1) * 1/(2^m2));

%
%%
% Bandpassfilter mit Ordnungen m1 = 32 und m2 = 36

y_m3 = [];


n = 1;
for r = 32:length(y)-32
    zw1 = 0;
    zw2 = 0;
    for k = 1:m2+1

        zw1 = zw1 + b1(k) * y(r+k-m1/2,3);
        zw2 = zw2 + b2(k) * y(r+k-m2/2,3);

    end
    y_m3(n,1) = (zw1 * 1/(2^m1)) - (zw2 * 1/(2^m2));
    n = n+1;
end


%

%
%%
% Plot Aufgabe 1

f4 = figure;
f4.Position(3:4) = [800 800];
figure(f4)
% title("Vergleich verschiedener Glättungsfilter mit der Ausgangszeitreihe")
hold on
subplot(4,1,1), plot(y(:,2),y(:,3)), subtitle("Ausgangszeitreihe"), ylabel("x[t]"), xlabel("Zeit [t]")
subplot(4,1,2), plot(y(1:length(y_m1),2),y_m1(:,1)), subtitle("Binomiales Glättungsfilter: Tiefpassfilter"), ylabel("x[t]"), xlabel("Zeit [t]")
subplot(4,1,3), plot(y(1:length(y_m2),2),y_m2(:,1)), subtitle("Hochpassfilter"), ylabel("x[t]"), xlabel("Zeit [t]")
subplot(4,1,4), plot(y(1:length(y_m3),2),y_m3(:,1)), subtitle("Bandpassfilter"), ylabel("x[t]"), xlabel("Zeit [t]")
hold off

saveas(f4,"images/4_Filter.png")

%
%%
% Aufgabe 2

clc

k = length(y)-length(y_m1);
fprintf("Kürzung bei Tiefpassfilter um %.0f Stellen.\n",k)
k = length(y)-length(y_m2);
fprintf("Kürzung bei Hochpassfilter um %.0f Stellen.\n",k)
k = length(y)-length(y_m3);
fprintf("Kürzung bei Bandpassfilter um %.0f Stellen.\n",k)
fprintf("----------------\n\n")

%
%%
% Aufgabe 4

[maxGB, k] = max(G_bp1);
fprintf("Maximaler Wert Durchlasscharakteristik: %.6f\n",maxGB)
fprintf("Frequenz bei Maximum: %.4f\n",mu(k))
fprintf("----------------\n\n")

w = mu(k)*(omega_g)^-1;





%
%%
% Aufgabe 5

% *************************************************************************
% Ausgangszeitreihe
% *************************************************************************

I = length(y);

X = fft(y(:,3));
c = 2*sqrt(X.*conj(X))/I;
c = c(2:I/2);

n = 2:I/2;
mu = n.*domega;

figure(5)
hold on
title("Amplitudenspektrum der Ausgangszeitreihe")
ylabel("Amplitude")
xlabel("Frequenz [°/h]")
plot(mu,c)
hold off

saveas(5,"images/5_AmpSpek_Ausgang.png")

% *************************************************************************
% Tiefpassfilter
% *************************************************************************

I = length(y_m1);

Y1 = fft(y_m1);
c1 = 2*sqrt(Y1.*conj(Y1))/I;
c1 = c1(2:I/2);

n = 2:I/2;
mu = n.*domega;

figure(6)
hold on
title("Amplitudenspektrum der Tiefpass-gefilterten Zeitreihe")
ylabel("Amplitude")
xlabel("Frequenz [°/h]")
plot(mu,c1)
hold off

saveas(6,"images/6_AmpSpek_Tiefpass.png")

% *************************************************************************
% Hochpassfilter
% *************************************************************************

I = length(y_m2);

Y2 = fft(y_m2);
c2 = 2*sqrt(Y2.*conj(Y2))/I;
c2 = c2(2:I/2);

n = 2:I/2;
mu = n.*domega;

figure(7)
hold on
title("Amplitudenspektrum der Hochpass-gefilterten Zeitreihe")
ylabel("Amplitude")
xlabel("Frequenz [°/h]")
plot(mu,c2)
hold off

saveas(7,"images/7_AmpSpek_Hochpass.png")

% *************************************************************************
% Bandpassfilter
% *************************************************************************

I = length(y_m3);

Y3 = fft(y_m3);
c3 = 2*sqrt(Y3.*conj(Y3))/I;
c3 = c3(2:I/2);

n = 2:I/2;
mu = n.*domega;

figure(8)
hold on
title("Amplitudenspektrum der Bandpass-gefilterten Zeitreihe")
ylabel("Amplitude")
xlabel("Frequenz [°/h]")
plot(mu,c3)
hold off

saveas(8,"images/8_AmpSpek_Bandpass.png")


close all
fprintf("================\n")



