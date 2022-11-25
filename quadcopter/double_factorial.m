function nee = double_factorial( n )
% nee = double_factorial( n )
%
% Returns nee = n!!

    if numel(n) > 1
        nee = zeros(size(n));
        for i=1:size(n,1)
            for j=1:size(n,2)
                nee(i,j) = double_factorial(n(i,j));
            end
        end
    elseif numel(n) == 1
        if n < 0
            nee = 1;
        elseif rem(n,2) == 0
            nee = 2^(n/2) * factorial(n/2);
        else
            nee = factorial(n) / double_factorial(n-1);
        end
    else
        nee = [];
    end
end

