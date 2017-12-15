%%downsample function
function ImData2 = downsample(ImData,ratio)
%% Read image from file 
[M,N,T] = size(ImData);
for i = 1:T
    ImData2(:,:,i) = imresize(ImData(:,:,i),[ceil(M/ratio),ceil(N/ratio)]);
end