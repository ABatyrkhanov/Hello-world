% Code for collecting the ARPES data and applying minimum gradient method
% to see the peak spectral features
tic
% Reading files
Ef=20.33; %MUST input known fermilevel here!
theta0=-10;  %MUST input the first angle in the series RELATIVE TO SURFACE NORMAL 
dtheta=2;     %MUST input a constant angle increment
bias=3;      %MUST manually input bias (don't include negative sign)
filemarker='*_AR*'; %marker to identify ARPES data
files=dir(filemarker); %collecting files in the directory with the marker
N=length(files); %number of files
format1=['%f %*f']; %format to skip the second column while reading files
format2=['%*f %f']; %format to skip the first column while reading files

%store one filename in order to get its length
%because they all have the same length
%more time efficient
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
KE=KE-bias;

for i=1:N
   fileID=fopen(files(i).name,'r');
   I(:,i)=fscanf(fileID,format2);
   fclose(fileID);
   KE(:,i)=KE(:,1); %energies are the same over the angles
   k(:,i)=const*sqrt(KE(:,1))*sin((theta0+(i-1)*dtheta)*pi/180);
end
BE=KE-Ef;


%Just copied this part from the curvature code
%====================================================================
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
axis([ min(min(k(:,:))) max(max(k(:,:))) min(BE(:,1)) max(BE(:,1)) min(min(I(:,:))) max(max(I(:,:)))+max(max(I(:,1))*0.5)])
ylabel('Binding Energy (eV)','FontSize',16)
xlabel('KPar','FontSize',16,'Units','Normalized','Position',[0.5 -0.06])
zlabel('Intensity','FontSize',16,'Units','Normalized','Position',[0 0 0])
title 'Raw data ARPES (No BG Subtraction)';
%======================================================================
toc