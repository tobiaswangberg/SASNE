function connected = is_connected(A)
    n = size(A,1);
    % initiliase the recursive dfs
    marked = 1; 
    unmarked = 2:n;
    v_curr = 1;

    [comp,marked,unmarked] = dfs(v_curr,marked,unmarked,A);
    connected = isempty(unmarked);
end