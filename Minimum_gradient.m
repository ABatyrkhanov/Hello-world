% Code for collecting the ARPES data and applying minimum gradient method
% to see the peak spectral features
tic
% Reading files
Ef=22.05; %MUST input known fermilevel here!
theta0=-10;  %MUST input the first angle in the series RELATIVE TO SURFACE NORMAL 
dtheta=2;     %MUST input a constant angle increment
bias=5;      %MUST manually input bias (don't include negative sign)
filemarker='*_AR*'; %marker to identify ARPES data
files=dir(filemarker); %collecting files in the directory with the marker
N=length(files); %number of files
format1=['%f %*f']; %format to skip the second column while reading files
format2=['%*f %f']; %format to skip the first column while reading files

%store one filename in order to get its length
%because they all have the same length
%more time efficient
filenames={};
filenames=[filenames,files(1).name]; 
length_check=[]; %not allocated list
length_check=[length_check,importdata(filenames{1})]; %fill the list from file
M=length(length_check); %this will be the length of the 

%datafiles=zeros(1,N);
KE=zeros(M,N);
BE=zeros(M,N);
k=zeros(M,N);
I=zeros(M,N);

KE(:,1)=length_check(:,1); %use the length_check array to fill the first row of KE
const=sqrt(2*9.109*1.602*10^-50)/(1.055*10^-34)*10^-10;

for i=1:N
   fileID=fopen(files(i).name,'r');
   I(:,i)=fscanf(fileID,format2);
   fclose(fileID);
   KE(:,i)=KE(:,1); %energies are the same over the angles
   k(:,i)=const*sqrt(KE(:,1))*sin((theta0+(i-1)*dtheta)*pi/180);
end
BE=KE-Ef;
KE=KE-bias;

%Minimum Gradient Method
G=zeros(M,N);
Ginv=zeros(M,N);
dBE=BE(2,1)-BE(1,1); %step size along energy axis
g_E=zeros(N-2);
g_W=zeros(N-2);

%G(1,1)=sqrt();
for j=2:(N-1)
    for i=2:(M-1)
        g_N=(I(i+1,j)-I(i,j))/dBE;
        g_S=(I(i-1,j)-I(i,j))/dBE;
        g_NE=(I(i+1,j+1)-I(i,j))/sqrt(dBE^2+(k(i,j+1)-k(i,j))^2);
        g_NW=(I(i+1,j-1)-I(i,j))/sqrt(dBE^2+(k(i,j-1)-k(i,j))^2);
        g_SE=(I(i-1,j+1)-I(i,j))/sqrt(dBE^2+(k(i,j+1)-k(i,j))^2);
        g_SW=(I(i-1,j-1)-I(i,j))/sqrt(dBE^2+(k(i,j-1)-k(i,j))^2);
        g_E=(I(i,j+1)-I(i,j))/(k(i,j+1)-k(i,j));
        g_W=(I(i,j-1)-I(i,j))/(k(i,j-1)-k(i,j));
        G(i,j)=sqrt(g_N^2+g_S^2+g_E^2+g_W^2+g_NE^2+g_NW^2+g_SE^2+g_SW^2);
        %G(i,j)=sqrt(g_N^2+g_S^2+g_E^2+g_W^2);
        Ginv(i,j)=1/G(i,j);
%        if (I(i,j)<63690) && (I(i,j)>63670)
%            i,j
%        end
    end
end

%{
I(377,12)
N=I(378,12)-I(377,12)
S=I(376,12)-I(377,12)
W=I(377,11)-I(377,12)
E=I(377,13)-I(377,12)
NE=I(378,13)-I(377,12)
NW=I(378,11)-I(377,12)
SE=I(376,13)-I(377,12)
SW=I(376,11)-I(377,12)
%}
I(1:3,1:11)
G(1:3,1:11)

I(375:380,7:15)
G(375:380,7:15)

%I=I./G;

%Just copied this part from the curvature code
%====================================================================
shadecontrol=0; %'1'=manual control, '0'=auto control
cmax=1E6; %color limit maximum value. Ignore if shadecontrol=0
cmin=0.001;    %color limit minimum value. Ignore if shadecontrol=0

brightcontrol=0; %'1'=manual control, '0'=auto control
brightness=-0.7; %Adjust between -1 and 1 to change brightness of colormap. Ignore if brightcontrol=0
%Kpar vs Angle
xaxis=1; %'1'=Kpar plot, '0'=Angle plot 
% Background. (See 'Background Removal' block in script for details)
bcoeff=1;  %This coefficient scales the background to be subtracted.

cMap = hot(200); %choose a color scheme
dataMax = 1E9;
dataMin = 0;
% the following 3 parameters are the ones we should tune for
% individual graphs
centerPoint =1E1; %max variation around this point
scalingIntensity = 5; %scaling factor
g=2.7; %coefficient for exp function
cs = 1:length(cMap); 
cs = cs - (centerPoint-dataMin)*length(cs)/(dataMax-dataMin);
cs = scalingIntensity * cs/max(abs(cs));
cs = sign(cs).* exp(abs(cs))*g;
newMap = interp1(cs, cMap, 1:400);

figure('Position',[100 0 700 600])
surf(k,BE,I)
colormap (newMap)
shading interp
if shadecontrol==1
    caxis([cmin,cmax])
end
if brightcontrol==1
    brighten(brightness)
end
view(0,89.9)
colorbar('westoutside')
axis([ min(min(k(:,:))) max(max(k(:,:))) min(BE(:,1)) max(BE(:,1)) min(min(I(:,:))) max(max(I(:,:)))])
ylabel('Binding Energy (eV)','FontSize',16)
xlabel('KPar','FontSize',16,'Units','Normalized','Position',[0.5 -0.06])
zlabel('Intensity','FontSize',16,'Units','Normalized','Position',[0 0 0])
title 'Raw data ARPES (No BG Subtraction)';
%title 'Normalized I/|G| (No BG Subtraction)';
%title 'Gradient of I (No BG Subtraction)';

figure('Position',[100 0 700 600])
surf(k,BE,G)
colormap (newMap)
shading interp
if shadecontrol==1
    caxis([cmin,cmax])
end
if brightcontrol==1
    brighten(brightness)
end
view(0,89.9)
colorbar('westoutside')
axis([ min(min(k(:,:))) max(max(k(:,:))) min(BE(:,1)) max(BE(:,1)) min(min(G(:,:))) max(max(G(:,:)))])
ylabel('Binding Energy (eV)','FontSize',16)
xlabel('KPar','FontSize',16,'Units','Normalized','Position',[0.5 -0.06])
zlabel('Intensity','FontSize',16,'Units','Normalized','Position',[0 0 0])
title 'Gradient of I (No BG Subtraction)';


%======================================================================

toc