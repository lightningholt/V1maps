


figure(199);
subplot(2,2,1);
imagesc(Big_Results.Res.zAbsPhiPref);
colormap(twilight);
colorbar
caxis([0, 360])


subplot(2,2,2);
imagesc(Big_Results.Res.zPhiPref);
colormap(twilight);
colorbar
caxis([0, 360])

subplot(2,2,3);
imagesc(Big_Results.Mill_Res.zAbsPhiPref);
colormap(twilight);
colorbar
caxis([0, 360])


subplot(2,2,4);
imagesc(Big_Results.Mill_Res.zPhiPref);
colormap(twilight);
colorbar
caxis([0, 360])