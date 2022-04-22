clear all
clc
format long g

%% NO.1  Import data
% Data StarPack StarPack.ascii
[long_gps,lat_gps,yrs_gps,jul_gps,hrs_gps,min_gps,sec_gps,hgt_gps,knt_gps] = textread('StarPack.ascii','%f %f %f %f %f %f %f %f %f', 'headerlines',1);
%Trans EM710 EM710Attitude.ascii
[jul_em,hrs_em,min_em,sec_em,fracsec_em,GPSsec_em,heading_em,roll_em,pitch_em,heave_em] = textread('EM710Attitude.ascii','%f %f %f %f %f %f %f %f %f %f', 'headerlines',1);
%Trans Knudsen KnudsenAttitude.ascii
[jul_k,hrs_k,min_k,sec_k,fracsec_k,GPSsec_k,heading_k,roll_k,pitch_k,heave_k] = textread('KnudsenAttitude.ascii','%f %f %f %f %f %f %f %f %f %f', 'headerlines',1);
% Read tide file, StJohnsPredictedTide.ascii
[jul_tide, hrs_tide, min_tide, tide] = textread('StJohnsPredictedTide.ascii','%f %f %f %f');

%% NO.2 Mengubah Waktu Ke Jam Desimal
%GPS
HrsGPS = hrs_gps+(min_gps/60)+(sec_gps/3600);
num2str(HrsGPS,'%7.6f');
GPS= table(long_gps,lat_gps,yrs_gps,HrsGPS,hgt_gps,knt_gps)

%%Knudsen Attitude
HrsKnud= hrs_k+2+((min_k+30)/60)+(sec_k/3600);
num2str(HrsKnud,'%7.6f');
KnudsenAtt= table(HrsKnud,heading_k,roll_k,pitch_k,heave_k)

%EM710
HrsEM710= (hrs_em+(min_em/60)+((sec_em+fracsec_em)/3600));
num2str(HrsEM710,'%7.6f');
EM710= table (HrsEM710,heading_em,roll_em,pitch_em,heave_em)

%Tide
HrsTide= (hrs_tide+(min_tide/60));
num2str(HrsTide,'%7.6f');
Tide= table(HrsTide, tide)

%% No.3 Interpolasi longitude dan Latitute (GPS)
%GPS-Knudsen Attitude
KnLong= interp1(HrsGPS,long_gps,HrsKnud,'linear');
KnLat= interp1(HrsGPS,lat_gps,HrsKnud,'linear');
KnudLongLat = table(HrsKnud,KnLong,KnLat);

%GPS-EM710
EMLong= interp1(HrsGPS,long_gps,HrsEM710,'linear');
EMLat= interp1(HrsGPS,lat_gps,HrsEM710,'linear');
EMLongLat= table (HrsEM710,EMLong,EMLat);

%Plot position
figure(20);
lat = mean(EMLat);
long = mean(EMLong);
geoplot(lat, long, 'r.', 'MarkerSize',15);
% geobasemap streets;
geobasemap colorterrain;
legend('Lokasi Survey');


%Plot Jalur EM710 & Knudsen
figure(19);
lat = (EMLat);
long = (EMLong);
geoplot(lat, long, 'r*');
hold on
lat = (KnLat);
long = (KnLong);
geoplot(lat, long, 'b*');
title('GeoPlot of EM710 & Knudsen');
legend('Em710', 'Knudsen');
geobasemap colorterrain;

%% NO.4 Menggunakan tide pada jule_tide 233, jam_tide 13, min_tide 0~15
File_row49_low= Tide(49,:)
File_row49_high= Tide(50,:)
Tide_file3= [File_row49_low; File_row49_high];
jul_tide3= [233;233];
hour_tide3= [12;12.25];
tide3= [0.6;0.59];
Tide_file4= table(jul_tide3, hour_tide3, tide3);

%Knudsen ke Tide 
Knud_Tides= interp1(hour_tide3,tide3,HrsKnud,'linear');
Height_Knud= table(HrsKnud,Knud_Tides);
num2str(Knud_Tides,'%7.6f');

%EM710 ke Tides
EM_Tides= interp1(hour_tide3,tide3,HrsEM710,'linear');
Height_EM710= table(HrsEM710,EM_Tides);
num2str(EM_Tides,'%7.6f');

%Hitung dlong,dlat, dan dheight dengan matrik rotasi
%hitung cell pertama dengan cos sin
roll_k_cos=cos(roll_k);
roll_k_sin=sin(roll_k);
roll_k_negsin=-sin(roll_k);
pitch_k_cos=cos(pitch_k);
pitch_k_sin=sin(pitch_k);
pitch_k_negsin=-sin(pitch_k);
heading_k_cos=cos(heading_k);
heading_k_sin=sin(heading_k);
heading_k_negsi=-sin(heading_k);

%% No.5 
% Menggunakan offser bersama dengan heading pitch dan roll untuk mendapat 
% koordinat masing-masing sounder dan antenna GPS di bingkai navigasi 
%(terapkantiga rotasi dalam urutan yang benar)

%Diketahui (X,Y,Z)
%Reference Point (RP)	0.0	0.0	0.0 (1.310 m above water line)
%33 kHz Transducer	    7.670	-3.240	3.420
%200 kHz Transducer	    5.140	-3.240	3.560
%EM710 Tx Transducer	6.904	-3.158	3.477
%EM710 Rx Transducer	5.990	-3.240	3.554
%Starpack antenna	   -5.220	-0.139	-6.629;

%Knudsen
rl=length(roll_k)
ph=length(pitch_k)
hd=length(heading_k)

%Roll
for roll=1:rl
    Rot_roll=[1 0 0; 0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
end 

%Pitch 
for pitch1=1:ph
    Rot_pitch=[cos(pitch1) 0 sin(pitch1); 0 1 0; -sin(pitch1) 0 cos(pitch1)];
end

%Heading
for heading=1:hd
    Rot_head=[cos(heading) -sin(heading) 0; sin(heading) cos(heading) 0; 0 0 1];
end
%rotasi 
R=Rot_head*Rot_pitch*Rot_roll
%koordinat 
X_knud33= R*[7.670;-3.240;3.420]
X_knud200= R*[5.140;-3.240;3.560];

%Matriks Rotasi Em710
roll = length(roll_em)
pitch1= length(pitch_em)
hd_EM= length(heading_em)
for roll=1:roll
    Rot_roll2= [1 0 0;...
        0 cos(roll) -sin(roll);...
        0 sin(roll) cos(roll)];
end
for pitch1=1:pitch1
        Rot_pitch2= [cos(pitch1) 0 sin(pitch1);...
           0 1 0;...
            -sin(pitch1) 0 cos(pitch1)];
end
for heading= 1:hd_EM
            Rot_head2= [cos(heading) -sin(heading) 0;...
                sin(heading) cos(heading) 0;...
                0 0 1];
end
%Rotasi 
R2=Rot_head2*Rot_pitch2*Rot_roll2

%Koordinat Knudsen 
X_EM710tx = R2*[6.904;-3.158;3.477]
X_EM710rx= R2*[5.990;-3.240;3.554]

%% NO.6
%Konvenrsi ke sistem koordinat lintang dan bujur 

%Parameter datum referensi (wgs84)
a=6378137;                      %semi major axis
b=6356752.314;                  %semi mainor axis
f=(a-b)/a;                      %flattening 
Eksentrisitas=(a^2-b^2)/(a^2);  %nilai eksentrisitas 1

%Hitung dLat dan dLon dari EM710
EM710_M=a*(1-Eksentrisitas)./(1-(Eksentrisitas.*((sind(EMLat)).^2))).^(3/2); %radius kelengkungan 
EM710_N = a./(1-(Eksentrisitas.*((sind(EMLat)).^2))).^(1/2); %jari-jari kelengkungan vertikal prima
%% No.7


%% No.8



%% Plot
% NO 1. Lat vs. long of the sounders and GPS antenna on the same plot.

[P,L]=size(lat_gps)
[C,R]=size(long_gps)
lat_gps(P,:) = [];
long_gps(P,:) = [];

% %Knudsen
% figure (1);
% plot(lat_gps, long_gps,'g',KnLat, KnLong);
% title ('Grafik Latitude Longitude GPS dan Knudsen');
% xlabel('Latitude');
% ylabel('longitude');
% grid on
% 
% %Em710
% figure (2);
% plot(lat_gps, long_gps,'g',EMLat, EMLong);
% title ('Grafik Latitude Longitude GPS dan EM710');
% xlabel('Latitude');
% ylabel('longitude');
% grid on

%GeoPlot
%Plot Knudsen
figure(18);
geoplot(lat_gps, long_gps, 'r*');
hold on
geoplot(KnLat, KnLong, 'b.');
title('GeoPlot of sounders and GPS antenna');
legend('GPS', 'Knudsen');
%Plot Knudsen
figure(17);
geoplot(lat_gps, long_gps, 'r*');
hold on
geoplot(EMLat, EMLong, 'y.');
title('GeoPlot of sounders and GPS antenna');
legend('GPS', 'Em710');

% No. 2 Plot of the Pitch vs  time
%Knudsen
figure(3);
plot(HrsKnud,pitch_k);
title('Grafik Pitch Knudsen');
xlabel('Time');
ylabel('Pitch Knudsen');
grid on

%Em710
figure(4);
plot(HrsEM710,pitch_em);
title('Grafik Pitch Em710');
xlabel(' Time');
ylabel('Pitch Em710');
grid on

% No. 3 Plot of the Roll vs  time
%Knudsen
figure(5);
plot(HrsKnud,roll_k);
title('Grafik Roll Knudsen');
xlabel('Time');
ylabel('Roll Knudsen');
grid on

%Em710
figure(6)
plot(HrsEM710,roll_em);
title('Grafik Roll Em710');
xlabel(' Time');
ylabel('Roll Em710');
grid on

% No. 4 Plot of the Heading vs UTC time
%Knudsen
figure(7);
plot(HrsKnud,heading_k);
title('Grafik Heading Knudsen');
xlabel('Time');
ylabel('Heading');
grid on

%Em710
figure(8);
plot(HrsEM710,heading_em);
title('Grafik Heading Em710');
xlabel('Time');
ylabel('Heading');
grid on

% No. 5 Heave vs UTC Time and induced heave vs. UTC time on same plot

% No. 6 Height variation of a sounder from GPS

% No. 7 Height Above chart datum vs UTC time each (which is RP above 1.310)
y = repmat(1.310,1,length(tide));
figure (9);
plot(HrsTide,tide,'b',HrsTide,y);

% No. 8 EM710 heave vs UTC time and Knudsen total heave vs UTC time on the same plot
%Knudsen
figure(10);
plot(HrsKnud,heave_k);
title('Grafik Heave');
xlabel(' Time');
ylabel('Heave');
grid on
hold on
%Em710
plot(HrsEM710,heave_em);
grid on
legend('Grafik Heave Knudsen','Grafik Heave Em710');