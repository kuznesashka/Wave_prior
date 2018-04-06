function loc = RAP_MUSIC_scan(G2, U, thresh)
% -------------------------------------------------------
% RAP-MUSIC scan
% -------------------------------------------------------
% FORMAT:
%   loc = RAP_MUSIC_scan(G2, U, thresh)
% INPUTS:
%   G2 -- Nch x Nsources*2 forward model matrix
%   U -- Nch x k signal subspace matrix
%   thresh -- minimal considered correlation value
%
% OUTPUT:
%   loc -- indexes of the sources with spikes
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

corr = MUSIC_scan(G2, U);
Nsens = size(G2,1);
Nsites = size(G2,2)/2;
[val_max, ind_max] = max(corr);
loc = [ind_max];
if val_max>thresh
    A = G2(:,(ind_max*2-1:ind_max*2));
    P = eye(Nsens, Nsens)-A*inv(A'*A)*A';
    Upr = P*U;
    Upr = orth(Upr);
%     [u s v] = svd(Upr);
%     Upr = u(:,1:size(Upr,2));
%     Upr = bsxfun(@rdivide, Upr, sqrt(diag(Upr'*Upr))');
    G2pr = P*G2;
%     range2 = 1:2;
%     for i = 1:Nsites
%         g = G2pr(:,range2);
%         [u s v] = svd(g, 'econ');
%         G2pr(:,range2) = u(:,1:2);
%         range2 = range2+2;
%     end
%     G2pr_norm = sqrt(sum(G2pr.^2, 1));
%     G2pr = bsxfun(@rdivide, G2pr, G2pr_norm);
    G2pr = orth(G2pr);
    corr = MUSIC_scan(G2pr, Upr);
    [val_max, ind_max] = max(corr);
    if val_max>thresh
        loc = [loc, ind_max];
    end
end

end

