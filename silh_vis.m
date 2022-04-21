function [plt,cluster_s,overall_s,s] = silh_vis(X,labels)
% Input is data matrix X (rows observations, columns features)
% and vector labels that gives that contains the true labels 
% Create silhouette plot plt, also returns pointwise (cluster_s), 
% clusterwise (cluster_s) and overall (s) silhouette score
distance = 'Euclidean';
label = unique(labels);
[s, plt] = silhouette(X,labels,distance);
% ax = plt.Children;
% ax.Position  = [.1 .1 .8 .8];
% ax.Box = 'off';
% ax.XLim = [-1 1];
% ax.XTick = [-1:.2:1];
% ax.XTickLabel = [num2str([-1:.2:1]',3);'    '];
% b = ax.Children;
% lim = ax.YLim;
% 
% mini = lim(1);
% ax.YLim = lim - mini;
% maxi = ax.YLim(2);
% ax.YLim = ax.YLim/maxi;
% ax.YTick = (ax.YTick - mini)./maxi;
% b.XData = (b.XData - mini)./maxi;
% x = labels;

str = [];
cluster_s = zeros(length(label),1);
for i = 1:max(size(label))
    
    indi = find(labels==label(i));
    meani = mean(s(indi));
    cluster_s(i) = meani;
    %str = [str,[newline,'Avg. for ',' grp ',num2str(label(i)),': ',num2str(meani,3)]];
    annotation('textbox',[.68 (1-.85*i/length(label)) .4 .05],'interpreter','latex','String',[num2str(meani,3)],'FitBoxToText','on','EdgeColor','none')
end
overall_s = mean(s);
str = [['Mean silhouette index: ',num2str(overall_s,3)],newline,str];
%annotation('textbox',[.68 .7 .4 .05],'interpreter','latex','String',str,'FitBoxToText','on','EdgeColor','red')
title(['$\overline{S} =$ ',num2str(overall_s,3)],'interpreter','latex')

end

