clear all
clc
%%Reading image
img_org=imread('pencere.jpg'); %we read our image
img_grey=rgb2gray(img_org);%we convert it from RGB to grey for  reducing the number of channel
img_double=im2double(img_grey); 
im_resize= imresize(img_double, 2,'bicubic');
%Padding image and creating mask
img_pad=padarray(im_resize,[16 16],'symmetric','both'); %we padded our image because we use 8 by 8 patch and we are searching its surroundings
patch_big=[40 40];
patch_small=[8 8];
    mask=ones(patch_big(1,1),patch_big(1,2));%we create a mask for parts that we dont want to see in our offset
    mask(12:28,12:28)=nan;
    mask(16:24,16:24)=nan;
    col=im2col(mask,patch_small,'sliding');
    V=sum(col,'omitnan');
    offset=[];
    dif_main=[];
%%
for j=1:size(img_pad,2)-patch_big(1,2)
    for i=1:size(img_pad,1)-patch_big(1,1)
    mat_40=img_pad(i:i+39,j:j+39);%this is our patch which is go through the image
    mat_8=im2col(mat_40,patch_small,'sliding'); %this is the small patch which is go through the big patch (patch_40)
    middle=mat_8(:,ceil(size(mat_8,2)/2));%middle one is the part that we are looking for the similer of it
    mat_8(:,V~=64)=nan;% we eliminate the unwanted part
    ind_middle=ceil(size(mat_8,2)/2);%we find the index of the patch  
    repeat=repmat(middle,1,size(mat_8,2));
    subs=repeat-mat_8; %we subtract mat_8 from repmat to find the most similar patch
    summ=sum(subs.^2);
    summ(:,ceil(size(mat_8,2)/2))=max(summ);%we equalize the patches itself to max because it is the closest one 
    A=ceil(ind_middle/sqrt(size(mat_8,2))); %%Location of middle in 40*40 patch
    B=ceil(ind_middle/sqrt(size(mat_8,2))); %%Location of middle in 40*40 patch
    pixel_A=A+i-1; %pixel value of the patch's itself in the picture
    pixel_B=B+j-1; %pixel value of the patch's itself in the picture
     
    [mins,indx]=find(summ == min(summ));% finding minimums and their locations in summ
    J=[]; 
    pixel_row=[];
    I=[]; 
    pixel_column=[];
    dif=[];
    J=ceil(indx/sqrt(size(mat_8,2))); %find minimum's location in the image
    I=indx-(J-1)*sqrt(size(mat_8,2));%find minimum's location in the image
    pixel_row=I+i-1;
    pixel_column=J+j-1; 
    % In below part we repmat the original pixels index and subtruct it
    % other patches indexes to find minimum distance
    pixel_AA=repmat(pixel_A,1,size(pixel_row,2)); %repmat the ooriginal index 
    pixel_BB=repmat(pixel_B,1,size(pixel_column,2));
   
    dif=sqrt((pixel_AA-pixel_row).^2+ (pixel_BB-pixel_column).^2);%we find the distance of the minimum ones and compare it with each other
    ind=indx(find(dif==min(dif),1,'first'));% we select the pne that havethe minimum distances
    %after find the correct minimum we calculate their index in original
    %image
    J=ceil(ind/sqrt(size(mat_8,2))); 
    I=ind-(J-1)*sqrt(size(mat_8,2));
    pixel_row=I+i-1;
    pixel_column=J+j-1;
    X=[pixel_A-pixel_row pixel_B-pixel_column];%we subtract to find the vector of the most similar one
    offset=vertcat(offset,X); 
    end
end
%%
figure
ctrs{1}=-50:50; %we determine an interval for histogram
ctrs{2}=-50:50;
hist3(offset,ctrs)%we create a histogram
%%Task3 padding and scroll the image up to size of offset value
abs_off=abs(offset);
max_off=[max(abs_off(:,1)), max(abs_off(:,2))]; %gives me the max dimentions to extend the image.

fake=zeros(size(img_double,1)*2,size(img_double,2)*2);%we create a fake matrix and we will fill it the extended version of image
%nancol=nan(size(fake,1),1);%matlab indexes begin with 1 so we create a nan column to fix this problem
for i=1:size(img_double,1)
    for j=1:size(img_double,2)
    fake(2*i-1,2*j-1)=img_double(i,j); %we extend our image we add zeros near the each pixels
    end
end
%fake=horzcat(nancol, fake);%matlab indexes begin with 1 so we create a nan column to fix this problem
pad2ext=padarray(fake,[max_off(1),max_off(2)],'symmetric','both');
%to make all shifted versions equal size we padd extended image by the size of maximum offset
zeross=zeros(size(pad2ext,1),size(pad2ext,2));%makes the remaining parts 0 after shift the image
img_strt=[max_off(1)+1,max_off(2)+1]; %where fake image starts


%%
%%WEIGHTED AVERAGING
%Offset Calculation
[Au,~,ic] = unique(offset,'rows'); %Au means offset values with deleted duplicates
k = accumarray(ic,1); %k gives us how many times these elements occured in offset
fprintf('\n\tPoints\tFrequency\n') %To tabulate
fprintf('\t%2d %2d\t\t%2d\n', [Au k]') %Also to tabulate
mask2=pad2ext;
mask2(pad2ext~=0)= NaN;
mask2(pad2ext==0)= 1;
%Averaging
weight_channel=[];
sum_mask=[];


for on=1:size(Au,1) %iterates until it reaches the number of elements in off vector
    ind_r=img_strt(1,1)+Au(on,1); %where shifted image starts as row
    ind_c=img_strt(1,2)+Au(on,2); %where shifted image starts as column
    zeross(ind_r:ind_r+size(fake,1)-1,ind_c:ind_c+size(fake,2)-1)=fake; %%SHIFTED IMAGES
    shifted_img=zeross.*mask2;%%SHIFTED IMAGES sadece boþ olan yerlerdeki dolma miktarýný görüyoruz
    mask3=shifted_img;
    weight_mask=mask3*k(on);
    weight_mask(isnan(weight_mask))=0;
    weight_channel(:,:,on)=weight_mask;
    mask3(isnan(mask3))=0;
    mask3(mask3~=0)=k(on);
    sum_mask(:,:,on)=mask3;
    zeross=zeros(size(pad2ext,1),size(pad2ext,2));  
end

%Without Normalizing
tot_sum=0;
img_weigts=0;
for on=1:size(Au,1)
    tot_sum=tot_sum+sum_mask(:,:,on);
    img_weigts=img_weigts+weight_channel(:,:,on); %Sumpix colects numerator of the weighted averaging calculation(preform before dividing total number of offset)
end
% img=rgb2gray(im_resize);%we convert it from RGB to grey for  reducing the number of channel
% img_pad=padarray(img,[16 16],'symmetric','both'); %we padded our image because we use 8 by 8 patch and we are searching its surroundings

bolum=img_weigts./tot_sum;
bolum(isnan(bolum))=0;
weighted_image=pad2ext+bolum;
% weighted_image=weighted_image+weighted_image;
% croped=weighted_image(17:474,17:456);
subplot(1,3,1)
imshow(weighted_image,[])
title('weighted_image')
subplot(1,3,2)
imshow(img_grey,[])
title('original')
subplot(1,3,3)
imshow(img_pad,[])
title('expanded')

%  figure(2)
% subplot(1,3,1)
% uint8_w=im2uint8(weighted_image);
% imshow(uint8_w)
% title('uint8 of weighted result')
% subplot(1,3,2)
% imshow(weighted_image)
% title('double of weighted result')
% subplot(1,3,3)
% imshow(weighted_image,[])
% title('[] of weighted result')
% % 
MSE=sum(sum(sqrt((weighted_image-img_pad).^2)))
percent_MSE=MSE/(size(im_resize,1)*size(im_resize,2))*100
%%