%Creates pictures similar to the ones seen in Miller 94
function MapsMiller(Inp, WRFon, WRFoff, varargin)

nDimV1 = Inp.nDimV1;
nDimLGN = Inp.nDimLGN;
Nx_arbor = Inp.Nx_arbor; 

WRFploton = reshape(WRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFploton = permute(WRFploton, [3,1,4,2]);
WRFploton = reshape(WRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

WRFplotoff = reshape(WRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFplotoff = permute(WRFplotoff, [3,1,4,2]);
WRFplotoff = reshape(WRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

%definition of huemap. 
hmap(1:256,1) = linspace(0,1,256);
hmap(:,[2 3]) = 0.7; %brightness
huemap = hsv2rgb(hmap);
colormap(huemap)

figure(1)
imagesc(WRFploton)
colorbar
title('On Cells WRF')

figure(2)
imagesc(WRFplotoff)
colorbar
title('Off Cells WRF')

figure(3)
imagesc(WRFploton-WRFplotoff)
colorbar
title('Difference between On & Off')
colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
colormap('gray')

%{
figure(4);
subplot(2,2,1)
imagesc(fftshift(fft2(Inp.Cnn)))
colorbar
title('fft2 Cnn')

subplot(2,2,2)
imagesc(fftshift(fft2(Inp.Cnn_fast)))
colorbar
title('fft2 Cnn-fast')

subplot(2,2,3)
imagesc(fftshift(fft2(Inp.G1)))
colorbar
title('fft2 G1')

subplot(2,2,4)
imagesc(fftshift(fft2(Inp.G2)))
colorbar
title('fft2 G2')

figure(5);
subplot(1,2,1)
imagesc(fftshift(fft2(Inp.G1).*fft2(Inp.Cnn)))
colorbar
title('GC1')

subplot(1,2,2)
imagesc(fftshift(fft2(Inp.G2).*fft2(Inp.Cnn_fast)))
colorbar
title('GC2')
%}
if nargin > 3
    Results = varargin{1};
    [PrefTheta, PrefPhi, zThetaPref, zPhiPref, zThetaSel, zPhiSel, spatialCorrRF, DThetaNeighbors,DThetaRandos, DPhiNeighbors, DPhiRandos] = extractResults(Results);
    
%     PrefTheta = varargin{1};
%     PrefPhi = varargin{2};
%     zThetaPref = varargin{3};
%     zPhiPref = varargin{4};
%     zThetaSel =  varargin{5};
%     zPhiSel =  varargin{6};
%     
    
    figure(111)
    % ax1 = subplot(4, 1, 1);
    % imagesc(WRFploton - WRFplotoff)
    % colormap(ax1, 'gray')
    % title('On/Off Dominance')
    % colorbar
    
    ax2 = subplot(3,2,1);
    imagesc(PrefTheta);
    colormap(ax2, jet)
    colorbar
    title('Pref Ori Map (Max Resp)')
    
    ax3 = subplot(3,2,2);
    imagesc(PrefPhi)
    colormap(ax3, jet)
    colorbar
    title('Pref Phase Map (Max Resp)')
    
    ax4 = subplot(3,2,3);
    imagesc(zThetaPref)
    colormap(ax4, jet)
    colorbar
    title('Pref Ori Map - arg(z)')
    
    ax5 = subplot(3,2,4);
    imagesc(zPhiPref)
    colormap(ax5, jet)
    colorbar
    title('Pref Phase Map - arg(z)')
    
    ax6 = subplot(3,2,5);
    imagesc(zThetaSel)
    colormap(ax6, jet)
    colorbar
    title('Ori Sel Map - abs(z)')
    
    ax7 = subplot(3,2,6);
    imagesc(zPhiSel)
    colormap(ax7, jet)
    colorbar
    title('Phase Sel Map - abs(z)')
    if nargin > 4
        edges = 0:10:90;
        
        
%         spatialCorrRF = varargin{7};
%         DThetaNeighbors = varargin{8};
%         DThetaRandos =  varargin{9};
        
        figure(23)
%         subplot(5 , 1, 5)
        histogram(spatialCorrRF, 'Normalization', 'probability')
        title(['Correlation Coefficient between Neighbors'])
        xlabel('Correlation Coefficient')
        ylabel('Percentage')
        
        figure(45)
        subplot(4, 2, 1);
        histogram(DThetaNeighbors,edges,  'Normalization', 'probability')
        title('Histogram of Preferred Ori Difference Neighbors')
        xlabel('Diff of Preferred Ori')
        ylabel('Probability')
        
        
        
        subplot(4, 2, 3)
        histogram(DThetaRandos,edges, 'Normalization', 'probability')
        title('Histogram of Preferred Ori Difference Random Pairs')
        xlabel('Diff of Preferred Ori')
        ylabel('Probability')
        
        subplot(2 ,2,2)
        imagesc(zThetaPref);
        title('Preferred Ori Map (arg z)')
        % imagesc(PrefTheta)
        % title('Preferred Ori Map (Max Resp)')
        colormap('jet');
        colorbar;
        
        edgesP = 0:10:180;
        subplot(4, 2, 5)
        histogram(DPhiNeighbors,edgesP,  'Normalization', 'probability')
        title('Histogram of Preferred Phase Difference Neighbors')
        xlabel('Diff of Preferred Phase')
        ylabel('Probability')
        
        subplot(4, 2, 7)
        histogram(DPhiRandos, edgesP, 'Normalization', 'probability')
        title('Histogram of Preferred Phase Difference Random Pairs')
        xlabel('Diff of Preferred Phase')
        ylabel('Probability')
        
        subplot(2 ,2, 4)
        imagesc(zPhiPref);
        title('Preferred Phase Map (arg z)')
        % imagesc(PrefTheta)
        % title('Preferred Ori Map (Max Resp)')
        colormap('jet');
        colorbar;
        
        if nargin > 5
            stitle = varargin{3};
            
            figure(4);
            if ispc
                suptitle(stitle);
            else
                sgtitle(stitle);
            end
            %make Rec field half
            %         ax41 = subplot(2, 1, 1);
            %         ax41 = subplot(3, 2, 1);
            ax41 = subplot(2, 2, 1);
            imagesc(WRFploton-WRFplotoff)
            colorbar
            title('Difference between On & Off')
            colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
            colormap(ax41, 'gray')
            axis('square')
            
            %make correlation quadrant
            %         subplot(2,3,6);
            subplot(2,2,2);
            %         subplot(3,2,2);
            histogram(spatialCorrRF, 'Normalization', 'probability')
            axis('square')
            title(['Correlation Coefficient between Neighbors'])
            xlabel('Correlation Coefficient')
            ylabel('Percentage')
            
            %make theta half quadrant
            %         subplot(8, 12, 49:51)
            %         subplot(6,2,5);
            %         subplot(8,4, 17);
            subplot(4,2, 5);
            histogram(DThetaNeighbors,edges,  'Normalization', 'probability')
            hold on
            histogram(DThetaRandos, edges, 'Normalization', 'probability', 'DisplayStyle','stairs','LineWidth', 2)
            legend('Neighbors', 'Random', 'Location', 'southeast');
            title('Histogram of Preferred Ori Difference')
            xlabel('Diff of Preferred Ori')
            ylabel('Probability')
            hold off
            
            %         subplot(8, 4, 21)
            %         subplot(8, 12, 61:63)
            %         subplot(6, 2, 7)
            %         subplot(8, 2, 11)
            %         histogram(DThetaRandos,edges, 'Normalization', 'probability')
            %         title('Histogram of Preferred Ori Difference Random Pairs')
            %         xlabel('Diff of Preferred Ori')
            %         ylabel('Probability')
            
            %         ax42 = subplot(4 ,12,29:32);
            %         ax42 = subplot(3,2, 4);
            ax42 = subplot(4,2, 6);
            %         ax42 = subplot(4 ,4,10);
            imagesc(zThetaPref);
            title('Preferred Ori Map (arg z)')
            % imagesc(PrefTheta)
            % title('Preferred Ori Map (Max Resp)')
            colormap(ax42, hsv);
            colorbar;
            
            %make phi map half quadrant
            %         subplot(8, 12, 73:75)
            %         subplot(6, 2, 9)
            %         subplot(8, 2, 13)
            % %         subplot(8, 4, 25)
            %         histogram(DPhiNeighbors,edgesP,  'Normalization', 'probability')
            %         title('Histogram of Preferred Phase Difference Neighbors')
            %         xlabel('Diff of Preferred Phase')
            %         ylabel('Probability')
            %
            % %         subplot(8, 12, 85:87)
            % %         subplot(6, 2, 11)
            %         subplot(8, 2, 15)
            % %         subplot(8, 4, 29)
            %         histogram(DPhiRandos, edgesP, 'Normalization', 'probability')
            %         title('Histogram of Preferred Phase Difference Random Pairs')
            %         xlabel('Diff of Preferred Phase')
            %         ylabel('Probability')
            subplot(4,2, 7);
            histogram(DPhiNeighbors,edgesP,  'Normalization', 'probability')
            hold on
            histogram(DPhiRandos, edgesP, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2)
            legend('Neighbors', 'Random', 'Location', 'southeast');
            title('Histogram of Preferred Phase Difference')
            xlabel('Diff of Preferred Phase')
            ylabel('Probability')
            hold off
            
            %         ax43 = subplot(4 ,12, 41:44);
            %         ax43 = subplot(3 ,2, 6);
            ax43 = subplot(4 ,2, 8);
            %         ax43 = subplot(4 , 4, 14);
            imagesc(zPhiPref);
            title('Preferred Phase Map (arg z)')
            % imagesc(PrefTheta)
            % title('Preferred Ori Map (Max Resp)')
            colormap(ax43, hsv);
            colorbar;
        end
    end
end