function y = ModelFunction(x,A12, A21)
    y = zeros(size(x));

    for i=1:length(x)
        gE_Mod = x(i).*(1-x(i)).*(A12.*x(i)+A21.*(1-x(i)));
        y(i) = gE_Mod;
    end
end