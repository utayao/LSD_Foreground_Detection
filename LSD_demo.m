% % **************************************************************************************************************
% % This is a code for Low-rank and Structured sparse matrix Decomposition (LSD)
% % If you happen to use this source code, please cite our papers:
% %
% % [1] "Background Subtraction Based on Low-rank and Structured Sparse Decomposition"
% % Xin Liu, Guoying Zhao, Jiawen Yao, Chun Qi.
% % IEEE Transactions on Image Processing, Vol. 24, No. 8, pp. 2502 - 2514, 2015.
% %
% % [2] "Foreground detection using Low rank and Structured sparsity"
% % Jiawen Yao, Xin Liu, Chun Qi.
% % IEEE ICME 2014. pp. 1 - 6, 2014 
% % 
% % Please note this code is only for LSD, the first-pass of background subtraction as described in our paper. 
% % The motion sailency check and Group-sparse RPCA (second-pass) are not included. 
% % 
% % For problems about our source code, please email Xin: linuxsino@gmail.com  or  Jiawen: yjiaweneecs@gmail.com 
% % **************************************************************************************************************

clc;
clear all

addpath(genpath('spams-matlab')); %add path for SPAMS tools

frame_idx=[
    1,48;
    ];%% can be tuned


% frame_idx=[
%     1,200;
%     201,400;
%     401,600;
%     ];

dataName = 'WaterSurface';

load(['data\' dataName ],'ImData'); 

num_task = size(frame_idx,1);
ImData0 = ImData;

for task = 1:num_task
    
frame_st = frame_idx(task,1); frame_ed = frame_idx(task,2);
ImData = ImData0(:,:,frame_st:frame_ed); 



%% down sampling
ratio = 4; % can be tuned, ratio = 4 in papers
ImData2 = downsample(ImData,ratio);
[M,N,T]=size(ImData);

%% first-pass: LSD 
sizeImg = [size(ImData2,1),size(ImData2,2)];
numImg = size(ImData2,3);
D2 = mat2gray(ImData2); % 0~1
ImMean = mean(D2(:)); 
D2 = D2 - ImMean; % subtract mean is recommended
D3 = reshape(D2,prod(sizeImg),numImg);   
disp('first-pass LSD');

graph = getGraphSPAMS(sizeImg,[3,3]);% for  lsd
[L_d S_d iter] = inexact_alm_lsd(D3,graph);   % lsd without L_21 TIP 2015
%[L_d S_d iter] = rpca_lsd(D3,graph);  %lsd with L_21 ICME 2014

%% using a very smple method to binarize foreground 
S = ForegroundMask(S_d,D3,L_d,0); 
S = reshape(S,[sizeImg,numImg]);
S_f = zeros(size(ImData));
L_d = uint8(mat2gray(reshape(L_d,size(ImData2))+ImMean)*256); %Background

for i = 1:T
    S_f(:,:,i) = imresize(S(:,:,i),[M,N],'box');  % back to origin size
    background(:,:,i) = imresize(L_d(:,:,i),[M,N],'box');  
end

foreground = uint8(256.*abs(S_f));


for i = 1:size(ImData,3)
    figure(1); clf;
     subplot(1,3,1);
     imshow(ImData(:,:,i)), axis off, colormap gray; axis off;
     title('Frames','fontsize',12);
     outtext = ['results\Img\', num2str(2000+frame_st+i-1),'.png'];
     imwrite(ImData(:,:,i),outtext);  
     
     subplot(1,3,2);
     imshow(foreground(:,:,i)), axis off,colormap gray; axis off;   
     title('Foreground','fontsize',12);
     outtext = ['results\FG\', num2str(2000+frame_st+i-1),'.png'];
     imwrite(foreground(:,:,i),outtext);      
 
     subplot(1,3,3);
     imshow(background(:,:,i)), axis off,colormap gray; axis off;   
     title('Background','fontsize',12);
     outtext = ['results\BG\', num2str(2000+frame_st+i-1),'.png'];
     imwrite(background(:,:,i),outtext);  
     drawnow; 
end  
end




