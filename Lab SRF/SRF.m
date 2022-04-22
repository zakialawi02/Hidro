clear all
clc
format long g
%% load data
EM710Att = load('EM710Attitude.ascii', '-ascii');
KnudsenAtt = load('KnudsenAttitude.ascii', '-ascii');
StJohnsPredic = load('StJohnsPredictedTide.ascii', '-ascii');
% StarP = load('StarPack.ascii', '-ascii');
file_id=fopen('StarPack.ascii','r');
cek=str2num(fgets(file_id))
[CC,RR]=size(cek);
filename = 'StarPack.ascii';
fid = fopen(filename);
%skip headers, read and throw away remainder of last header line
% textscan(fid, '%[^\n]', 1, 'CommentStyle', {'#', '~A'});
StarP = cell2mat( textscan(fid, repmat('%f', 1, RR), 'CollectOutput', true) );
fclose(fid);

[P,L]=size(EM710Att);
[C,R]=size(KnudsenAtt);
[B,S]=size(StJohnsPredic);
%% Remap data/Classify data/memisahkan data
Tide= StJohnsPredic(:,4);
tTide= [StJohnsPredic(:,1) StJohnsPredic(:,2) StJohnsPredic(:,3)];

X= [EM710Att(:,1) EM710Att(:,2) EM710Att(:,3) EM710Att(:,4) EM710Att(:,5) EM710Att(:,6)];
XX= [EM710Att(:,7) EM710Att(:,8) EM710Att(:,9) EM710Att(:,10)];

Y= [KnudsenAtt(:,1) KnudsenAtt(:,2) KnudsenAtt(:,3) KnudsenAtt(:,4) KnudsenAtt(:,5) KnudsenAtt(:,6)];
YY= [KnudsenAtt(:,7) KnudsenAtt(:,8) KnudsenAtt(:,9) KnudsenAtt(:,10)];

CoordX= StarP(:,1); [CC,RR]=size(CoordX); CoordX(CC,:) = [];
CoordY= StarP(:,2); [CC,RR]=size(CoordY); CoordY(CC,:) = [];
Coord= [CoordX CoordY];
WWW= [StarP(:,3) StarP(:,4) StarP(:,5) StarP(:,6) StarP(:,7)]; [CC,RR]=size(WWW); WWW(CC,:) = [];
QQQ= StarP(:,8); [CC,RR]=size(QQQ); QQQ(CC,:) = [];
EEE= StarP(:,9); [CC,RR]=size(EEE); EEE(CC,:) = [];

%% ?

% Rn=a/sqrt(1-e^2*sind(Lintang)*sind(Lintang));
% Rm=(a*(1-e^2))/(1-e^2*sind(Lintang)*sind(Lintang))^(3/2);
% du= dN/Rm
% dy= dE/Rn*cosd(Lintang)











