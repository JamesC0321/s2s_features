function feature = SoS(img)

scalenum = 3;
im2color=img;

R = 1; P = 8;  
lbp_type = { 'ri' 'u2' 'riu2' }; 
y = 3;
mtype = lbp_type{y};
MAPPING = getmapping( P, mtype ); 
feat = [];
imdist=double(rgb2gray(im2color));
for itr_scale = 1:scalenum

    [map] = MLVMap(imdist);
    LBPMap = lbp_new(map,R,P,MAPPING,'x');
    wLBPHist = [];
    weightmap = map;
    wintesity = weightmap(2:end-1, 2:end-1);
    wintensity = abs(wintesity);
    for k = 1:max(LBPMap(:))+1 
        idx = LBPMap == k-1; 
        kval = sum(wintensity(idx));
        wLBPHist = [wLBPHist kval]; 
    end
    wLBPHist = wLBPHist/sum(wLBPHist);
    feat = [feat wLBPHist];


     imdist                   = imresize(imdist,0.5);
end    
%%  Gabor

im=rgb2gray(img);
num_scale =3 ;   
num_orien =4;   
minWaveLength =3;   
mult = 2;       
sigmaOnf = 0.65;     
dThetaOnSigma =1.5;  

gf = gaborconvolve(im,  num_scale, num_orien, minWaveLength, mult, sigmaOnf, dThetaOnSigma);

A=[];
C=[];
D=[];
for i = 1:num_scale    
    for j = 1:num_orien    
           
          nrf=10.*abs(gf{i,j});
          ss=nrf.*log2(nrf);
          ss=-ss;
          C  = [C mean(ss(:))];   

          sub_im = log(abs(gf{i,j})+.1);  
          win = [1 0 -1; 0 0 0; -1 0 1;];
          deri_im_7 = filter2(win, sub_im, 'same');
          [Gx ,Gy] = gradient(deri_im_7);
          [sigma_1, alpha_1] = gaussian_para_esti(Gx(:));
          [sigma_2, alpha_2] = gaussian_para_esti(Gy(:));
          GM2 = (abs(Gx) + abs(Gy))/2;
          phat1  = wblfit(GM2(:));
          A=[A sigma_1 alpha_1 sigma_2 alpha_2 phat1 ];
    end
end
 

feature = [feat  A  C];

