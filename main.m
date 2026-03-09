clc
clear all
close all

T2 = [5, 10, 15, 20];       % T2* relaxation times
TE = [1:1.375:16.5]';       % Time to eecho
s0 = [155, 255, 355, 455];  % Initial signal inensity values
sigma2 = 0;                 % noise variance
Nrow = 32;                  % Phantom image size
Ncol = 32;                  % Phantom image size

%% Create phantom
for uu = 1:length(T2)
    Phantom_WO_NoiseTemp{uu} = createPhantoms('exp', TE, T2(uu), s0(uu), Nrow, Ncol);
end
Phantom_WO_Noise = [Phantom_WO_NoiseTemp{1}, Phantom_WO_NoiseTemp{2}; Phantom_WO_NoiseTemp{3}, Phantom_WO_NoiseTemp{4}];
[Nrow_, Ncol_, bands] = size(Phantom_WO_Noise);

%% Add noise if required
Y = Phantom_WO_Noise + sqrt(sigma2) * randn(Nrow_, Ncol_, bands);
yReshaped = reshape(Y, Nrow_*Ncol_, bands)';

%% Run T2* estimation algorithm
lambdaA = 1e-5;
lambdaR = 1e-5;
[a, r, g, f] = relaxationEst(yReshaped, TE, Nrow_, Ncol_, lambdaA, lambdaR);

%% Plot stuff - MODIFY AS REQUIRED

%% -------- Figure 1 --------
figure(1)
set(gcf,'Units','normalized','Position',[0.05 0.3 0.9 0.4])

tiledlayout(1,4,'TileSpacing','compact','Padding','compact')

% Ground truth s0
nexttile
S0_Image = [s0(1)*ones(Nrow, Ncol), s0(2)*ones(Nrow, Ncol); ...
            s0(3)*ones(Nrow, Ncol), s0(4)*ones(Nrow, Ncol)];

imagesc(S0_Image)
axis image off
caxis([min(Y(:,:,1),[],'all') max([a(:); s0(:)])])
colorbar
title('a_0','FontSize',14)

% Observed signal
nexttile
imagesc(Y(:,:,1))
axis image off
caxis([min(Y(:,:,1),[],'all') max([a(:); s0(:)])])
colorbar
title('a_1','FontSize',14)

% Estimated a
nexttile
imagesc(reshape(a(:,:,1),Nrow_,Ncol_))
axis image off
caxis([min(Y(:,:,1),[],'all') max([a(:); s0(:)])])
colorbar
title('Estimated a_0','FontSize',14)

% Mean estimated a
a_reshaped = reshape(a,Nrow_,Ncol_);

a_mean = [mean(a_reshaped(1:32,1:32),'all')*ones(Nrow,Ncol), ...
          mean(a_reshaped(1:32,33:end),'all')*ones(Nrow,Ncol); ...
          mean(a_reshaped(33:end,1:32),'all')*ones(Nrow,Ncol), ...
          mean(a_reshaped(33:end,33:end),'all')*ones(Nrow,Ncol)];

nexttile
imagesc(a_mean)
axis image off
caxis([min(Y(:,:,1),[],'all') max([a_mean(:); s0(:)])])
colorbar
title('Mean of Estimated a_0','FontSize',14)


%% -------- Figure 2 --------
figure(2)
set(gcf,'Units','normalized','Position',[0.05 0.05 0.9 0.4])

tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

% Estimated T2*
nexttile
imagesc(reshape(1./r(:,:,1),Nrow_,Ncol_))
axis image off
caxis([0 max([1./r(:); T2(:)])])
colorbar
colormap hsv
title('Estimated T2*','FontSize',14)

% Mean estimated T2*
r_reshaped = reshape(1./r,Nrow_,Ncol_);

r_mean = [mean(r_reshaped(1:32,1:32),'all')*ones(Nrow,Ncol), ...
          mean(r_reshaped(1:32,33:end),'all')*ones(Nrow,Ncol); ...
          mean(r_reshaped(33:end,1:32),'all')*ones(Nrow,Ncol), ...
          mean(r_reshaped(33:end,33:end),'all')*ones(Nrow,Ncol)];

nexttile
imagesc(r_mean)
axis image off
caxis([0 max([r_mean(:); T2(:)])])
colorbar
colormap hsv
title('Mean Estimated T2*','FontSize',14)