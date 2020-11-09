% Make Arbor fig


A = fftshift(reshape(Inp.A(1,:), Inp.nDimV1,Inp.nDimLGN));

figure(87)
subplot(1,2,1)
axis square
midline = ceil(Inp.nDim/2)+1;
shift = 4
xx = (-shift-Inp.R_arbor):(shift+Inp.R_arbor);
xx_inds = midline + xx
plot(xx, A(midline, xx_inds), 'LineWidth',3, 'Color', [0.4940, 0.1840, 0.5560])
%set(gca, 'Visible', 'off');
%axes('XAxisLocation', 'bottom', 'YAxisLocation','left')
yticks([0, 1])
ylabel('Arbor strength')
xticks(xx(1:2:end));
xlabel('Neural distance (grid spacing)')
box off

subplot(1,2,2)
imagesc(A)
axis square
colorbar
colormap('gray')
axis off 
rectangle('Position', [0, midline, Inp.nDim, 1], 'EdgeColor',[0.4940, 0.1840, 0.5560], 'LineWidth', 2)
