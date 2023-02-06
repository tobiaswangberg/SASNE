function [comps,count] = find_comps(A)
% Takes as input an adjacency matrix A
% returns the disconneted components of the graph comps and the number
% of disconnected components.
    n = size(A,1);
    unmarked = 1:n;
    count = 0;
    comps = ones(1,n);
    while ~isempty(unmarked)
        marked = [];
        v_curr = unmarked(1);
        marked = union(v_curr,marked);
        unmarked = setdiff(unmarked,v_curr);
        [comp,marked,unmarked] = dfs(v_curr,marked,unmarked,A);
        count = count + 1;
        comps(marked) = count*ones(length(marked),1);
    end

end