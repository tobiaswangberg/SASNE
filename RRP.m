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

function [plt] =  smoothed_2d_hist(data,x_width,y_width,fast)

if ~exist('x_width', 'var')
    x_width = 1;
end
if ~exist('y_width', 'var')
    y_width = 1;
end
max_rank = max(data(:));

x_edges = [-.5:x_width:max_rank];
y_edges = [-(max_rank+.5):y_width:max_rank];


edges = {x_edges y_edges};

[N,c] = hist3(data,'Edges',edges);

bandwidth = x_width;

centers_x = c{1}(:);
centers_y = c{2}(:);


nbins = length(centers_x);

N = N.^2;
plt = pcolor(centers_x,centers_y,N');




dens_x = zeros(nbins,1);
dens_y = zeros(2*nbins,1);
if ~fast
    dens = zeros(nbins,2*nbins);
    for i = 1:nbins  
        dens_x = normpdf(data(:,1),centers_x(i),bandwidth);
        for j = 1:2*nbins
            dens_y = normpdf(data(:,2),centers_y(j),bandwidth);
            dens(i,j) = mean(dens_x.*dens_y);
        end
    end
    size(data,1).*dens'
    
    m = max(N(:));
    plt = pcolor(centers_x,centers_y,m.*dens'/max(dens(:)));
end
hold on

%mycolormap = customcolormap([0 .5 (1-1/max_rank) 1], [1 0 0; 0 1 0; 1/4 1/4 3.5/4; 1 1 1]);
myColorMap = jet(256);
%myColorMap(1,:) = 1;
colormap(myColorMap);



shading interp

set(plt,'edgecolor','none');

plot([0 max_rank],[(max_rank -1) 0],'r-','LineWidth',.8);
hold on
plot([0 (max_rank)],[0 -(max_rank - 1)],'r-','LineWidth',.8);
% ncx = length(centers_x);
% ncy = length(centers_y);
% for i = 1:ncx
%     for j = 1:ncy
%         curr_cx = centers_x(ncx);
%         curr_cy = centers_y(ncx);
%         if curr_cx < abs(curr_cy)
%             continue
%         end
%         disp('hej')
%         rectangle('Position',[(curr_cx - x_width/2),(curr_cy - y_width/2),x_width,y_width],'EdgeColor','none','FaceColor',[1 1 1]);
%         hold on 
%     end
% end
box off
xlim([-10,max_rank+5])
ylim([-max_rank-30,max_rank+5])
xticks([0,max_rank])
yticks([-max_rank,0,max_rank])
xticklabels([0 1])
yticklabels([-1 0 1])
xlabel('Original rank')
ylabel('Rank error')



end

function pointwise_ranks = dist_to_rank(D)
    D = squareform(D);
    n = size(D,1);
    [~,ind] = sort(D);
    pointwise_ranks = repmat([1:n-1]',1,n);
    for i = 1:n
        r = 1:n;
        r(ind(:,i)) = r;
        r = r(setdiff(1:n,i)) - 1;
        pointwise_ranks(:,i) = r;
    end
end