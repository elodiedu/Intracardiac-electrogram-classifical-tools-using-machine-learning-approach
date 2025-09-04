%Step 1&2
%Step 1 2 3 5 6 in this file/Step 4 in hallmark.mlx-Step 8 exists only for PVI
%
clc
clear;
 close all
 clear all;

load('Matfile/2-LA_SINUS.mat')
%%Step 2
load('Matfile/c2.mat')
load('Matfile/vertices.mat')
%%
dir_new='C:\01_Study\03_2022_summer\01_Lab\New_data\ExportData25_01_96 09_35_16\Patient 2022_03_18\Study 1\Study2024\P1'
filename="1-RA.mat"
filepath = fullfile(dir_new, filename);
data=load(filepath)
RA_data=data.data
filename="2-LA.mat"
filepath = fullfile(dir_new, filename);
data=load(filepath)
LA_data=data.data
%%

 DirLocation = 'C:\01_Study\03_2022_summer\01_Lab\New_data\ExportData25_01_96 09_35_16\Patient 2022_03_18\Study 1\Export_Study-1-01_25_2096-09-21-42'; 
mainfo = dir(fullfile(DirLocation,'*.mesh')); 


K=3;  %%%%LA Fib Mesh
thisfilename = mainfo(K).name;
[vertices,faces] = convertmesh(DirLocation,thisfilename,0);
%%
RA_data = RA_data(4:end, :);
[~, numColumns] = size(RA_data);
RA_DF_result = zeros(1, numColumns);
RA_MSF_result = zeros(1, numColumns);
RA_Kt_result = zeros(1, numColumns);
RA_MSE_result = zeros(1, numColumns);
for i=4:numColumns
    RA_DF_result(i)=dominant_freq(RA_data(:, i),1000);
    RA_MSF_result(i)=MSF_1D(transpose(RA_data(:, i)),1000,RA_DF_result(i));
    RA_Kt_result(i)=kurtosis(RA_data(:, i));
    RA_MSE_result(i)=SampEn(RA_data(:, i),5);
end
LA_data = LA_data(4:end, :);
[~, numColumns] = size(LA_data);
LA_DF_result = zeros(1, numColumns);
LA_MSF_result = zeros(1, numColumns);
LA_Kt_result = zeros(1, numColumns);
LA_MSE_result = zeros(1, numColumns);
for i=1:numColumns
    LA_DF_result(i)=dominant_freq(LA_data(:, i),1000);
    LA_MSF_result(i)=MSF_1D(transpose(LA_data(:, i)),1000,LA_DF_result(i));
    LA_Kt_result(i)=kurtosis(LA_data(:, i));
    LA_MSE_result(i)=SampEn(LA_data(:, i),5);
  end
%%
%%

%%Step 3 Calculate DF/MSF/MSE/Kt/SE or MI
%For comparing accuracy, we need to test MSF KT SE /MSE KT SE/MSF MSE SE
for i=1:numColumns
 

    signal_all=cell2mat(UniqueBiEGMs(i));
    signal=signal_all(:,j)
    PVI_SIN_DF{i}(j)=dominant_freq(signal,1000);
    PVI_SIN_MSF{i}(j)=MSF_1D(transpose(signal),1000,PVI_SIN_DF{i}(j));
    PVI_SIN_Kt{i}(j)=kurtosis(signal);
    PVI_SIN_MSE{i}(j)=SampEn(signal,5);
    PVI_SIN_SE{i}(j)=Shn_ent(signal);
    %PVI_MI{i}=findmicorr(signal_all);
end

%%Step 3 finish-> savefiles by hand

%% 
%%Step 6-option 1 , convert 3 cells to a 3-dimensional array, MSF, Kt, MSE
% option -2: MSE, SE, MSF
SIN_MSF_KT_MSE_0=zeros(length(UniqueBiEGMs),15,3);
for i =1:length(UniqueUiEGMs)
    for j =1:15
        %original with 4 parameter version: MSF,KT,MSE,SE
        %test MSF KT SE  filename:Sin_3
        %test MSF MSE SE  filename:Sin_4
        SIN_MSF_KT_SE_0(i,j,1)=PVI_SIN_MSF{i}(j);
        SIN_MSF_KT_SE_0(i,j,2)=PVI_SIN_Kt{i}(j);
%AF_MSF_KT_MSE_0(i,j,2)=PVI_SIN_MSE{i}(j);
        SIN_MSF_KT_SE_0(i,j,3)=PVI_SIN_SE{i}(j);
    end
end

%save file to mat
%
%% 
%%Step 8-only for PVI label to match data：transfer Idx to label data
for j=1:66
    for k =1:22
        vertice_index=Idx(j,k);
        Idx(j,k)=Vla_Label(vertice_index);

    end
end

%save Idx data as name-matfile

%% Step 5 
for i=1:length(UniqueUiEGMs)
    
for j =1:20
    signal_all=cell2mat(UniqueUiEGMs(i));
    signal=signal_all(:,j)
    PVI_DF{i}(j)=dominant_freq(signal,1000);
    PVI_MSF{i}(j)=MSF_1D(transpose(signal),1000,PVI_DF{i}(j));
    PVI_Kt{i}(j)=kurtosis(signal);
    PVI_MSE{i}(j)=SampEn(signal,5);
    PVI_SE{i}(j)=Shn_ent(signal);
    %PVI_MI{i}=findmicorr(signal_all);
end
end




%%
for i=1:1:length(ElectrodePosition)

Pentaryneartestpointonshell=ElectrodePosition{i};
BipolariEGMs=BiEGMs{i};
UipolariEGMs=UiEGMs{i};
%BiMi{i}=findmicorr(BipolariEGMs);
 for j=1:1:15             %%%%%C alculate MI here


% BiSampEn{i}(j)=SampEn(BipolariEGMs(:,j),5);
%  BiSampSE{i}(j)=Shn_ent(BipolariEGMs(:,j));
%  BiSampKt{i}(j)=kurtosis(BipolariEGMs(:,j));
  BiSampdf{i}(j)=dominant_freq(BipolariEGMs(:,j),1000);
%  BiSampmsf{i}(j)=MSF_1D(transpose(BipolariEGMs(:,j)),1000,BiSampdf{i}(j));

% BiSampEn{i}(j)=SampEn(BipolariEGMs(:,j),1);
% % UiSampEn{i}(j)=SampEn(UipolariEGMs(:,j),1);
% end
%tic% #Measures the MI runtime
%MI{i}=findmicorr(UipolariEGMs);

%toc
%for j=1:1:20   for Unipolar       %%%Calculate MI here

 
 %UiMi{i}(j)=SampEn(UipolariEGMs(:,j),5);
%UiSampEn{i}(j)=SampEn(UipolariEGMs(:,j),5);
 %UniSampSE{i}(j)=Shn_ent(UipolariEGMs(:,j));
 %UniSampKt{i}(j)=kurtosis(UipolariEGMs(:,j));
 %UniSampdf{i}(j)=dominant_freq(UipolariEGMs(:,j),1000);
 %UniSampmsf{i}(j)=MSF_1D(transpose(UipolariEGMs(:,j)),1000,UniSampdf{i}(j));

end

end 

% save MSEres BiSampEn UiSampEn
%  firstUniqueMSE=[];
% % 
% for i=1:1:20
% 
% temp=UiSampEn{i};
% firstUniqueMSE=[firstUniqueMSE;temp];
% 
% end









% 
%% 
%% 
% BiSampEn{i}(j)=SampEn(BipolariEGMs(:,j),5);
%  BiSampSE{i}(j)=Shn_ent(BipolariEGMs(:,j));
%  BiSampKt{i}(j)=kurtosis(BipolariEGMs(:,j));
for i =1:227
    Bi_PVI_sample_KT{i}=kurtosis(Final_point(:,i));
    Bi_PVI_sample_MSE{i}=SampEn(Final_point(:,i),5);
    Bi_PVI_sample_SE{i}=Shn_ent(Final_point(:,i));
end
%%Assumption: MI/MSE high EMD. MSF/Kt low EMD


%%

for i=1:1:length(UniqueLocationsLocs)
 c1 = 0.01.*ones(length(vertices),1);
 Pentaryneartestpointonshell=ElectrodePosition{UniqueLocationsLocs(i)};
 Pentaryneartestpointonshell=Pentaryneartestpointonshell(3:22,:);
   %MSE=UiSampEn{UniqueLocationsLocs(i)};
   %SE=UniSampSE{UniqueLocationsLocs(i)};
   KT=BiSampdf{UniqueLocationsLocs(i)};
%   %DF=UniSampdf{UniqueLocationsLocs(i)};
   %MI_unique=transpose(MI{UniqueLocationsLocs(i)});
 Name='KT'%change each time
 % MSF=UniSampmsf{UniqueLocationsLocs(i)};



% 
Idx= knnsearch(vertices,Pentaryneartestpointonshell,'K',50);
[m,n]=size(Idx)
CoveredPointId=reshape(Idx,[m*n,1]);
for ii=1:1:length(CoveredPointId)
c1(CoveredPointId(ii),:)=0.015;
end


Idx= knnsearch(vertices,Pentaryneartestpointonshell,'K',36);
[m,n]=size(Idx);

 for k=1:1:14
     Value=KT(k);%change each time
   
 c1(Idx(k,1:12))=Value;
 c1(Idx(k,13:24))=Value;
% 
 c1(Idx(k,25:36))=Value;
% 
 end


map = [166 166 166;
       62 38 168;
       69 55 213;
       70 60 222;
       71 65 229;
       66 105 254;
       62 111 254;
       56 117 254;
       50 123 252;
       47 129 250;
       45 140 243;
       43 145 238;
       39 150 235;
       35 160 229;
       24 173 219;
       17 177 214;
       7 181 208;
       1 184 202;
       2 186 195;
       11 189 188;
       24 191 182;
       36 193 174;
       44 195 167;
       55 200 151;
       63 202 142;
       74 203 132;
       171 199 57;
       185 196 49;
       197 194 42;
       209 191 39;
       220 189 41;
       230 187 45;
       239 186 53;
       248 186 61;
       254 189 60;
       245 227 39;
       247 245 27;
       249 251 20]/255;
c1_KT=c1;%change each time
f=figure('position',[0 0 2000 1000])

plot3(ElectrodePosition{UniqueLocationsLocs(i)}(1,1),ElectrodePosition{UniqueLocationsLocs(i)}(2,2),ElectrodePosition{UniqueLocationsLocs(i)}(1,3),'','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')
% view([60,90,350]);
% hold on;

fc = patch('Faces',faces,'Vertices',vertices,'FaceColor','interp','EdgeColor','black','edgealpha',0.5,'FaceVertexCData',c1);

 colormap(map);
 colorbar;
title('Catheter unique positions',i)
saveas(f,Name,'jpg')
saveas(f,Name,'fig')
F(i)=getframe(gcf)
 writerObj = VideoWriter('MSE_new_changedview_3.avi');
  writerObj.FrameRate = 1;

open(writerObj);
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);









 
 




end




map = [166 166 166;
       62 38 168;
       69 55 213;
       70 60 222;
       71 65 229;
       66 105 254;
       62 111 254;
       56 117 254;
       50 123 252;
       47 129 250;
       45 140 243;
       43 145 238;
       39 150 235;
       35 160 229;
       24 173 219;
       17 177 214;
       7 181 208;
       1 184 202;
       2 186 195;
       11 189 188;
       24 191 182;
       36 193 174;
       44 195 167;
       55 200 151;
       63 202 142;
       74 203 132;
       171 199 57;
       185 196 49;
       197 194 42;
       209 191 39;
       220 189 41;
       230 187 45;
       239 186 53;
       248 186 61;
       254 189 60;
       245 227 39;
       247 245 27;
       249 251 20]/255;

% figure
% fc = patch('Faces',faces,'Vertices',vertices,'FaceColor','interp','EdgeColor','none','FaceVertexCData',c1);
% colormap(map);


















% for i = 1:length(C_pos_lateral)
%     c1(C_pos_lateral(i),1) = Samp_En_2(i,10);
% end
% 
% % Nearest neighbour interpolation
% 
% for j = 1 : length(C_pos_lateral)
%     d1=[];
%     for i=1:length(Vla_lateral)
%         d1(i,2)=sqrt((Vla(C_pos_lateral(j),1)-Vla_lateral(i,1))^2 +(Vla(C_pos_lateral(j),2)-Vla_lateral(i,2))^2 +(Vla(C_pos_lateral(j),3)-Vla_lateral(i,3))^2);
%         d1(i,1)=Vla_lateral(i,4);
%     end
%     [~,indx] = sort(d1(:,2),'ascend');
%     sorted_d1 = d1(indx,:);
%     val = 0;
%     for id = 1:2
%         c1(sorted_d1(id,1))= c1(C_pos_lateral(j),1);
%     end
% 
%     for id = 3:4
%         c1(sorted_d1(id,1))= 0.80*c1(C_pos_lateral(j),1);
%     end
% 
%     for id = 5:6
%         c1(sorted_d1(id,1))= 0.6*c1(C_pos_lateral(j),1);
%     end
% 
% end
% 
% %/////////////////////////////////////////////////////
% 
% Vla_filled_lateral = [];
% pts = [];
% p=1;
% k=1;
% for i = 1: length(Vla_lateral)    
%     if c1 (Vla_lateral(i,4)) == 0
%         pts(k) = Vla_lateral (i,4);
%         k = k + 1;
%     else
%         Vla_filled_inf(p) = Vla_lateral (i,4);
%         p = p + 1;
%     end
% end
% 
% for j = 1: length(pts)
% 
%     for i = 1 : length(Vla_filled_inf)
% 
%         dist(i,2) = sqrt((Vla(pts(j),1)-Vla(Vla_filled_inf(i),1))^2 +(Vla(pts(j),2)-Vla(Vla_filled_inf(i),2))^2 +(Vla(pts(j),3)-Vla(Vla_filled_inf(i),3))^2);
% 
%         dist(i,1) = Vla_filled_inf(i);
% 
%     end
%     [~,indx] = sort(dist(:,2),'ascend');
%     sorted_dist = dist(indx,:);
% 
%     id = 1:6; arr = c1(sorted_dist(id,1),1);
% 
%     c1(pts(j))= mean(arr);
% end
% 
% % to fill out the grey spots in septal wall
% for i = 1 : length(Vla_lateral)  
%     if c1(Vla_lateral(i,4)) < 0.05
%         c1(Vla_lateral(i,4)) = 0.055;
%     end
% end
% 





% %%
% clc
% close all
% clear all
% DirLocation = 'C:\Users\kongx\Downloads\ExportData25_01_96 09_35_16\Patient 2022_03_18\Study 1\Export_Study-1-01_25_2096-09-21-42'; 
% mainfo = dir(fullfile(DirLocation,'*.mesh')); 
% K=6;  %%%%Fib Mesh
% thisfilename = mainfo(K).name;
% [vertices,faces] = convertmesh(DirLocation,thisfilename,0);
% 
% searchECGfilename=[thisfilename(1:length(thisfilename)-5) '*' 'ECG_Export.txt'];
% subfo=dir(fullfile(DirLocation,searchECGfilename));
% [~,idx] = sort([subfo.datenum]);
% subfo = subfo(idx);
%  for J=1:1:length(subfo)
%     ElectrodePosition{J}=[];
% 
% ECGfileName = subfo(J).name;
% [Pentary_Eleclectrode_Positions,Pentary_Eleclectrode_PositionsOn,Pentary_Eleclectrode_PositionsBi,Pentary_Eleclectrode_PositionsBiOn] = SearchMAGNETIC_EleclectrodeLocation(DirLocation,ECGfileName);
% ElectrodePosition{J}=Pentary_Eleclectrode_PositionsOn(:,3:5);
% 
% [UiEGMs{J} BiEGMs{J} ] = extractECG(DirLocation,ECGfileName);
% 
% 
% 
% 
% 
%  end
% 
%  save iEGMSinformation.mat UiEGMs BiEGMs ElectrodePosition UniqueLocationsLocs UniqueLocationspeaks
% 
% 
%%
%Route B, matfile: B_AF_MSFMSEKT.mat + B_Sinus_MSFMSEKT.mat
%Route A matfile: msf-kt-mse.mat + Frequencyoutput/AF/AF_MSF_KT_MSE.mat
% 
BASE_Af=load('Study2\AF_MSF_KT_MSE.mat');
PVI_AF=load('Study2\PVI_AF_MSF_KT.mat');
PVI_SIN=load('Study2\PVI_SIN_MSF_KT.mat')
%%

figure;

for j=1:92

for k =1:20
    
scatter3(PVI_SIN_MSF_KT_MSE(j,k,2),PVI_SIN_MSF_KT_MSE(j,k,3),PVI_SIN_MSF_KT_MSE(j,k,1),'MarkerFaceColor',[1 0 0])
hold on;
%scatter3(La_Sinus_Kt{j},La_Sinus_MSE{j},La_Sinus_MSF{j},'MarkerFaceColor',[0.1 .05 .35])%\
if j<63
   scatter3(PVI_AF_MSF_KT_MSE(j,k,2),PVI_AF_MSF_KT_MSE(j,k,3),PVI_AF_MSF_KT_MSE(j,k,1),'MarkerFaceColor',[0.1 1 0.05])
if j<=38
   scatter3(PVI_AF_MSF_KT_MSE_0(j,k,2),PVI_AF_MSF_KT_MSE_0(j,k,3),PVI_AF_MSF_KT_MSE_0(j,k,1),'MarkerFaceColor',[0.06 .02 .96])
   hold on;
end
end
end
end
hold on;
xlabel('Kt','FontSize', 10)
ylabel('MSE','FontSize', 10)
zlabel('MSF','FontSize', 10)
legend('Sinus','BASELINE AF',"PVI AF")
title('Sinus vs. AF','FontSize',25)
%%

figure;
scatter3(LA_Kt_result(:), LA_MSE_result(:), LA_MSF_result(:), 'b', 'filled'); % Plot matrix A in red
hold on; % Hold the plot to add the second matrix
scatter3(RA_Kt_result(:), RA_MSE_result(:), RA_MSF_result(:),  'MarkerFaceColor',  [1, 0.5, 0]); % Plot matrix B in blue
hold off;

xlabel('Kt','FontSize', 10)
ylabel('MSE','FontSize', 10)
zlabel('MSF','FontSize', 10)
legend('LA','RA')
%%
