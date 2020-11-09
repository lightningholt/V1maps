%2D nearest neighbor difference

function roughness = nearNeighbor2D(proj, Inp)


if length(size(proj)) == 2
    %This wants proj to be a 2D grid. NOT a vector
    R_arbor = sqrt(2);
    d = Inp.d;
    
    % [A2,~,~,~,~] = Arbor2(Inp.Xdist, R_arbor,d, N1,N2, round(N1/2), 1, 1);
    %
    % for i = 1:N1
    %     for j = 1:N2
    %         arbor_i = reshape(A2(i,:), N1, N2);
    %         temp_rough = (proj(i,j) - proj).^2;
    %         temp_rough(~arbor_i) = [];
    %         roughness(i,j) = mean(temp_rough);
    %     end
    % end
    
    r1 = diff(proj,[], 1).^2;
    r2 = diff(proj,[], 2).^2;
    
    roughness = mean(r1(:)+r2(:));
else
    time = size(proj,2);
    types = size(proj,3);
    nDim = Inp.nDim;
    
    roughness = zeros(time, types);
    
    for jj = 1:types
        for tt = 1:time
            proj2D = reshape(squeeze(proj(:,tt,jj)), nDim, nDim);
            
            r1 = diff(proj2D, [], 1).^2;
            r2 = diff(proj2D, [], 2).^2;
            
            roughness(tt, jj) = mean(r1(:)+r2(:));
        end
    end
end

            
        