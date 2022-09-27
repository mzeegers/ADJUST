%Plot from the Disks phantom
T = readmatrix('AllOtherPhantoms_MatrixTXray100chan_42_MinEn5_MaxEn35_mat_CORRECTED');

% imagesc(T);
% axis image; xlabel('Energy channel'); ylabel('Material number');
% set(gca, 'FontSize', 16);
% pbaspect([2 1 1]);
% yticklabels = 23:64;
% yticks = linspace(1, size(T, 1), numel(yticklabels));
% yticksnew = [yticks(1), yticks(3), yticks(8), yticks(13), yticks(18), yticks(23), yticks(28), yticks(33), yticks(38), yticks(42)];
% yticklabelsnew = [yticklabels(1), yticklabels(3), yticklabels(8), yticklabels(13), yticklabels(18), yticklabels(23), yticklabels(28), yticklabels(33), yticklabels(38), yticklabels(42)]
% set(gca, 'YTick', yticksnew, 'YTickLabel', yticklabelsnew)

plot([T(1,:)' T(5,:)' T(10,:)' T(18,:)'], 'LineWidth',1.5)
set(gca,'FontSize',10);
xlabel('Energy bin')
xlim([1 100])
pbaspect([2 1 1]);
ylabel('Attenuation')
legend('Vanadium', 'Cobalt', 'Arsenic', 'Zirocnium')