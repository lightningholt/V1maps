%Function to find the linear selectivities in 2D to the correlation matrices

%The equation for the difference in Won and Woff growth is the following
%(note W = Won - Woff

%tau_w dW/dt = K_s (Cnn-Cnf) W + K_f(Cnn_fast + Cnf_fast) W 

%signs are correct.
% so W should become selective to the features in Cnn - Cnf and
% Cnn_fast-Cnf_fast but only one should develop a map

%IMPORTANTLY: This is assuming that Cnn = Cff and Cnf = Cfn
function [wp, CL] = linear2Dselectivity(Inp, Won, Woff)
%Extract relevant variables
Cnn = Inp.Cnn;
Cnf = -0.5*Cnn;

dC_slow = Cnn - Cnf;

Cnn_fast = Inp.Cnn_fast;
Cnf_fast = -0.5*Cnn_fast;

dC_fast = Cnn_fast - Cnf_fast; 

W = Won - Woff; 

%Find the principal components of dC's 

if length(dC_slow) ~= Inp.nV1
    dC_slow = globalDist(dC_slow, Inp.nDimV1, Inp.d);
end

if length(dC_fast) ~= Inp.nV1
    dC_fast = globalDist(dC_fast, Inp.nDimV1, Inp.d);
end

[e1, v1] = svds(dC_slow,1);
[e2, v2] = svds(dC_fast,1);

e3 = randn(length(e1),1);
e3 = e3 - (e2'*e3)*e2 - (e1'*e3)*e1;
e3 = e3/(e3'*e3);

lambda1 = Inp.lambda1;
lambda2 = Inp.lambda2;

if e2'*e1 > 0.01
    warning('e1 and e2 are not orthogonal yo')
    e2'*e1
elseif e2'*e3 > 0.01
    warning('e2 and e3 could be more orthogonal')
    e2'*e3
elseif e1'*e3 > 0.01
    warning('e1 and e3 could be more orthogonal')
    e1'*e3
end

if v1/v2 > 1.2
    warning('Eigenvalue 1 is too large')
    v1/v2
elseif v1/v2 <0.8
    warning('Eigenvalue 1 too small')
    v1/v2
end
   
%wpX = w projected on either slow (1) or fast (2) principal axis
wp1 = W*e1;
slow_mean = mean(wp1.^2);
wp1 = reshape(wp1, Inp.nDimV1, Inp.nDimV1);
wp1_roughness = nearNeighbor2D(wp1, Inp)/slow_mean;
CL1 = correlationLengthCalc(wp1, Inp);

wp2 = W*e2;
fast_mean = mean(wp2.^2);
wp2 = reshape(wp2, Inp.nDimV1, Inp.nDimV1);
wp2_roughness = nearNeighbor2D(wp2, Inp)/fast_mean;
CL2 = correlationLengthCalc(wp2, Inp);

wp3 = W*e3;
neither_mean = mean(wp3.^2);
wp3 = reshape(wp3, Inp.nDimV1, Inp.nDimV1);
wp3_roughness = nearNeighbor2D(wp3, Inp)/neither_mean;
CL3 = correlationLengthCalc(wp3, Inp);

%find the min and max of the color scheme
cmax = max([max(wp1(:)), max(wp2(:)), max(wp3(:))]);
cmin = min([min(wp1(:)), min(wp2(:)), min(wp3(:))]);

expectedScaledDifference = (lambda1*v1*max(abs(e1)))/(lambda2*v2*max(abs(e2)))
actualScaledDifference = max(abs(wp1(:)))/max(abs(wp2(:)))


Wi = Inp.Won - Inp.Woff;
wp1i = Wi*e1;
slow_mean_i = mean(wp1i.^2);
wp1i = reshape(wp1i, Inp.nDimV1, Inp.nDimV1);
wp1i_roughness = nearNeighbor2D(wp1i, Inp)/slow_mean;

wp2i = Wi*e2;
fast_mean_i = mean(wp2i.^2);
wp2i = reshape(wp2i, Inp.nDimV1, Inp.nDimV1);
wp2i_roughness = nearNeighbor2D(wp2i, Inp)/slow_mean;

wp3i = Wi*e3;
neither_mean_i = mean(wp3i.^2);
wp3i = reshape(wp3i, Inp.nDimV1, Inp.nDimV1);
wp3i_roughness = nearNeighbor2D(wp3i, Inp)/slow_mean;

figure(65)
plot(1:length(e1), e1, 1:length(e1),e2, 1:length(e1), e3)
legend('PC1', 'PC2', 'PC3')
xlabel(num2str(max(abs(e1))/max(abs(e2))))


figure(66);
subplot(1,3,1)
imagesc(wp1)
title("Selectivity to slow C's")
botstr = strcat('Mean Selectivity =', num2str(slow_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp1_roughness)];
xlabel(botstr);
caxis([cmin, cmax])
colorbar

subplot(1,3,2)
imagesc(wp2)
title("Selectivity to fast C's")
botstr = strcat('Mean Selectivity =', num2str(fast_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp2_roughness)];
xlabel(botstr);
caxis([cmin, cmax])
colorbar

subplot(1,3,3)
imagesc(wp3)
title("Selectivity to neither")
botstr = strcat('Mean Selectivity =', num2str(neither_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp3_roughness)];
xlabel(botstr);
caxis([cmin, cmax])
colorbar

figure(67);
subplot(1,3,1)
imagesc(wp1)
title("Selectivity to slow C's")
botstr = strcat('Mean Selectivity =', num2str(slow_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp1_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(1,3,2)
imagesc(wp2)
title("Selectivity to fast C's")
botstr = strcat('Mean Selectivity =', num2str(fast_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp2_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(1,3,3)
imagesc(wp3)
title("Selectivity to neither")
botstr = strcat('Mean Selectivity =', num2str(neither_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp3_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

figure(68)
n_space = 1:length(CL1);

plot(n_space, CL1, n_space, CL2, n_space, CL3, 'LineWidth', 2.5);
% hold on
% plot(n_space, dC_slow(1,n_space)/max(dC_slow(1,n_space)),'k', n_space, dC_fast(1,n_space)/max(dC_fast(1,n_space)),'k*', 'LineWidth',1)
% hold off
legend('slow feature', 'fast feature', 'neither', 'Slow Corr','Fast Corr')
title('Correlation lengths')
xlabel('Neural distance')

figure(69)

subplot(3,3,1)
imagesc(wp1)
title("Selectivity to slow C's")
botstr = strcat('Mean Selectivity =', num2str(slow_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp1_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(3,3,2)
imagesc(wp2)
title("Selectivity to fast C's")
botstr = strcat('Mean Selectivity =', num2str(fast_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp2_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(3,3,3)
imagesc(wp3)
title("Selectivity to neither")
botstr = strcat('Mean Selectivity =', num2str(neither_mean));
botstr = [botstr newline ' Roughness = ' num2str(wp3_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar


subplot(3,3,4)
imagesc(wp1i)
title("Init Selectivity to slow C's")
botstr = strcat('Mean Selectivity =', num2str(slow_mean_i));
botstr = [botstr newline ' Roughness = ' num2str(wp1_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(3,3,5)
imagesc(wp2i)
title("Init Selectivity to fast C's")
botstr = strcat('Mean Selectivity =', num2str(fast_mean_i));
botstr = [botstr newline ' Roughness = ' num2str(wp2i_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(3,3,6)
imagesc(wp3i)
title("Init Selectivity to neither")
botstr = strcat('Mean Selectivity =', num2str(neither_mean_i));
botstr = [botstr newline ' Roughness = ' num2str(wp3i_roughness)];
xlabel(botstr);
% caxis([cmin, cmax])
colorbar

subplot(3,3,7)
imagesc(wp1- wp1i)
title("Growth of Selectivity to slow C's")
xlabel('Final proj - Initial')
% caxis([cmin, cmax])
colorbar

subplot(3,3,8)
imagesc(wp2 - wp2i)
title("Growth of Selectivity to fast C's")
xlabel('Final proj - Initial')
% caxis([cmin, cmax])
colorbar

subplot(3,3,9)
imagesc(wp3 - wp3i)
title("Growth of Selectivity to neither")
xlabel('Final proj - Initial')
% caxis([cmin, cmax])
colorbar

wp(:,:, 1) = wp1;
wp(:,:,2) = wp2;
wp(:,:,3) = wp3;

CL = [CL1, CL2, CL3];
