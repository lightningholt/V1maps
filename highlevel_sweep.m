%high level wrapper

lambda2 = 2:2:20
ii = 9;

dirname = 'data_' + num2str(ii);
mkdir dirname
for ii = 1:length(lambda2)
    temp_l2 = lambda2(ii);
    
    [Inp, Big_Res] = C10VisMaps(ii, temp_l2, dirname);
end

