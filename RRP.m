function plot = RRP(D1,D2)
% Create RRP plot given distance matrices D1 and D2
    fast = true;
    D1sq = squareform(D1);
    n = size(D1sq,1);
    r1 = dist_to_rank(D1);
    r2 = dist_to_rank(D2);
    score = mean(mean(abs(r1(:) - r2(:))/(n-1)));
    
    
    plot = smoothed_2d_hist([r1(:),r2(:)-r1(:)],floor(n*0.05),floor(n*0.05),fast);
    title(['$\overline{R}=$ ',num2str(score,3)],'interpreter','latex')
end

