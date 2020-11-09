

function higherOrderPC_selectivities(Big_Results, Inp, order)


[e_slow, s_slow] = eigs(globalDist(Inp.Cnn, Inp.nDim,Inp.d) ,order);
[e_fast, s_fast] = eigs(globalDist(Inp.Cnn_fast, Inp.nDim,Inp.d) ,order);

e_slow = e_slow(:, 1:order);
e_fast = e_fast(:, 1:order);

s_slow = s_slow(1:order, 1:order)
s_fast = s_fast(1:order, 1:order)

Wdyn = Big_Results.Res.Won - Big_Results.Res.Woff;
Wmill = Big_Results.Mill_Res.Won - Big_Results.Mill_Res.Woff;

dyn_lin_slow = Wdyn*e_slow;
dyn_lin_fast = Wdyn*e_fast;

mill_lin_slow = Wmill*e_slow;
mill_lin_fast = Wmill*e_fast;

CL_slow_dyn = [];
CL_fast_dyn = [];
CL_slow_mill = [];
CL_fast_mill = [];

for oo = 1:order
    wp_dyn_slow = reshape(dyn_lin_slow(:,oo), Inp.nDimLGN, Inp.nDimV1);
    wp_dyn_fast = reshape(dyn_lin_fast(:,oo), Inp.nDimLGN, Inp.nDimV1);
    
    wp_mill_slow = reshape(mill_lin_slow(:,oo), Inp.nDimLGN, Inp.nDimV1);
    wp_mill_fast = reshape(mill_lin_fast(:,oo), Inp.nDimLGN, Inp.nDimV1);
    
    CL_temp = correlationLengthCalc(wp_dyn_slow, Inp);
    CL_slow_dyn = [CL_slow_dyn, CL_temp];
    
    CL_temp = correlationLengthCalc(wp_dyn_fast, Inp);
    CL_fast_dyn = [CL_fast_dyn, CL_temp];
    CL_temp = correlationLengthCalc(wp_mill_slow, Inp);
    CL_slow_mill = [CL_slow_mill, CL_temp];
    CL_temp = correlationLengthCalc(wp_mill_fast, Inp);
    CL_fast_mill = [CL_fast_mill, CL_temp];
end
    
xx = (1:length(CL_temp))-1;


h = figure(12345);
fontsize = 14
set(gca, 'Colormap', redblue)
axis equal

nrows = order;
ncols = 5;
slow_plots_dyn = 1:ncols:ncols*order-4;
fast_plots_dyn = 2:ncols:ncols*order-3;
slow_plots_mill = 3:ncols:ncols*order-2;
fast_plots_mill = 4:ncols:ncols*order-1;

cl_plots = 5:ncols:ncols*order;

if max(abs(dyn_lin_slow(:))) > max(abs(dyn_lin_fast(:)))
    cmax1 = max(dyn_lin_slow(:));
    cmin1 = min(dyn_lin_slow(:));
else
    cmax1 = max(dyn_lin_fast(:));
    cmin1 = min(dyn_lin_fast(:));
end


if max(abs(mill_lin_slow(:))) > max(abs(mill_lin_fast(:)))
    cmax2 = max(mill_lin_slow(:));
    cmin2 = min(mill_lin_slow(:));
else
    cmax2 = max(mill_lin_fast(:));
    cmin2 = min(mill_lin_fast(:));
end


if abs(cmax1) > abs(cmin1)
    cmin1 = -cmax1;
else
    cmax1 = abs(cmin1);
end

if abs(cmax2) > abs(cmin2)
    cmin2 = -cmax2;
else
    cmax2 = abs(cmin2);
end


for ii= 1:order
    subplot(nrows,ncols,slow_plots_dyn(ii))
    imagesc(reshape(dyn_lin_slow(:,ii), Inp.nDimV1, Inp.nDimV1))
    c = colorbar('TickDirection','out');
    c.Label.String = 'Selectivity';
    caxis([cmin1, cmax1])
    colormap(redblue)
    % axis off
    xticks([])
    yticks([])
    xticklabels([])
    yticklabels([])
    box on
    xstr = ['Dyn Slow PC ',num2str(ii)];
    xlabel(xstr)
    %text(-0.15, 1, 'A', 'Units','normalized', 'FontSize', fontsize)

    subplot(nrows, ncols, fast_plots_dyn(ii))
    imagesc(reshape(dyn_lin_fast(:,ii),Inp.nDimV1, Inp.nDimV1));
    xstr = ['Dyn Fast PC ',num2str(ii)];
    xlabel(xstr)
    c = colorbar('TickDirection','out');
    c.Label.String = 'Selectivity';
    colormap(redblue)
    caxis([cmin1, cmax1])
    box on
    xticklabels([])
    yticklabels([])
    xticks([])
    yticks([])
    %text(-0.15, 1, 'B', 'Units','normalized', 'FontSize', fontsize)
    
    subplot(nrows,ncols,slow_plots_mill(ii))
    imagesc(reshape(mill_lin_slow(:,ii), Inp.nDimV1, Inp.nDimV1))
    c = colorbar('TickDirection','out');
    c.Label.String = 'Selectivity';
    colormap(redblue)
    caxis([cmin2, cmax2])
    % axis off
    xticks([])
    yticks([])
    xticklabels([])
    yticklabels([])
    box on
    xstr = ['Mill Slow PC',num2str(ii)];
    xlabel(xstr)
    %text(-0.15, 1, 'A', 'Units','normalized', 'FontSize', fontsize)

    subplot(nrows, ncols, fast_plots_mill(ii))
    imagesc(reshape(mill_lin_fast(:,ii),Inp.nDimV1, Inp.nDimV1));
    xstr = ['Mill Fast PC', num2str(ii)];
    xlabel(xstr)
    c = colorbar('TickDirection','out');
    c.Label.String = 'Selectivity';
    colormap(redblue)
    caxis([cmin2, cmax2])
%     caxis([cmin1, cmax1])
    box on
    xticklabels([])
    yticklabels([])
    xticks([])
    yticks([])
    
    subplot(nrows, ncols, cl_plots(ii))
    plot(xx, CL_slow_dyn(:,ii))
    hold on 
    plot(xx, CL_fast_dyn(:,ii))
    plot(xx, CL_slow_mill(:,ii))
    plot(xx, CL_fast_mill(:,ii))
    hold off
    
    if ii == 1
        legend({'Dyn slow', 'Dyn fast', 'Mill Slow', 'Mill Fast'},'Location','eastoutside')
    end
    
end
%{
subplot(2,2,3)
imagesc(lin_control);
xlabel('Control Feature')
c = colorbar('TickDirection','out');
c.Label.String = 'Selectivity';
colormap(redblue)
caxis([cmin1, cmax1])
% axis off
xticks([])
yticks([])
xticklabels([])
yticklabels([])
box on
text(-0.15, 1, 'C', 'Units','normalized', 'FontSize', fontsize)

subplot(2,2,4)
plot(x, CL(:, 1:2), 'LineWidth',3);
hold on
plot(x, zeros(length(x),1), '--k', 'LineWidth', 3)
hold off
legend({'Slow feature', 'Fast feature'})
legend('boxoff')
xlabel('Neural distance (grid spacing)')
ylabel('Correlation')
text(-0.15, 1, 'D', 'Units','normalized', 'FontSize', fontsize)
%}