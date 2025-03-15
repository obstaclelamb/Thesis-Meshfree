% pathMediaGen = "C:\Users\admin\Downloads\程式類\MATLAB\媒體生成"; addpath(pathMediaGen);
close all; clear;

% kernel = 'Heaviside';
kernel = 'Tent';
% kernel = 'Quadratic B';
% kernel = 'Cubic B';
% kernel = 'Quartic';
% kernel = 'Quintic B';


% kernel={'Heaviside', 'Tent', 'Quadratic B', 'Cubic B', 'Quartic', 'Quintic B'};

figure; hold on;
% for i = 1:6
  ZL=0:0.01:2;
  [PHIZL,DPHIZL,DDPHIZL]=getPhiDD(ZL,kernel);  
  ZL=[-flip(ZL),ZL];
  PHIZL=[flip(PHIZL),PHIZL];
  plot(ZL,PHIZL,'LineWidth',5)
% end
axis([-2,2,0,1.1]);
title(['Kernel: ',kernel]);
% title('Various Kernel');
xlabel('z')
ylabel('\phi_a')
% legend(kernel)

outpath = "C:\Users\admin\Downloads\成大\碩一\無網格法-林冠中\RKPM2D open source\out";
% pngmaker(append(outpath,'\kernel5-',kernel));
pngmaker(append(outpath,'\kernel1-Tent'));