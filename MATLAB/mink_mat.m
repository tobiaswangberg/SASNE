function [row,col,vals] = mink_mat(M,k)
    [nrow,ncol] = size(M);
    [~,am] = mink(M(:),k);
    rem = mod(am,nrow);
    col = (am - rem)/nrow + (rem ~= 0);
    row = rem + nrow*(rem==0);
    vals = diag(M(row,col));
end