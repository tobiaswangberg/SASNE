function [comp,marked,unmarked] = dfs(v,marked,unmarked,A)
    NNs = A(v,:);
    NNs = find(NNs > 0);
    %disp(['knns ',int2str(length(knns)),newline])
    NNs = intersect(NNs,unmarked);
    comp = NNs;
    marked = union(marked,comp);
    unmarked = setdiff(unmarked,comp);
    if ~isempty(comp)
        n = length(comp);
        for i = 1:n
            [new_comp,marked,unmarked] = dfs(NNs(i),marked,unmarked,A);
            comp = union(comp,new_comp);
        end
    end
end