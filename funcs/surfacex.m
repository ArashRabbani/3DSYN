% function []=surfacex(A,Title)            
%             if nargin==1; Title='Voxel Values'; end
%             faces=6;
%             A=squeeze(double(A));
%             A=permute(A,[2,1,3]);
%             S=size(A);
%             [X,Y,Z] = meshgrid(1:S(1),1:S(2),1:S(3));
%             xslice = [1,S(1)]; yslice = [1,S(2)]; zslice = [1,S(3)];
%             h=slice(X,Y,Z,permute(A,[2 1 3]),xslice,yslice,zslice);
%             axis equal tight;
%             for I=1:faces; h(I).EdgeColor='none'; end
%             c=colorbar; c.Label.String = Title;
%             plotcube([min(X(:)),max(X(:)),min(Y(:)),max(Y(:)),min(Z(:)),max(Z(:))]);
%             view(30,40);
% end

function []=surfacex(A, Title)
    if nargin==1 
        Title='Voxel Values'; 
    end

    % Get current axes
    current_ax = gca;
    
    % Process input volume
    A = squeeze(double(A));
    A = permute(A,[2,1,3]);
    S = size(A);
    
    % Create meshgrid
    [X,Y,Z] = meshgrid(1:S(1),1:S(2),1:S(3));
    xslice = [1,S(1)]; 
    yslice = [1,S(2)]; 
    zslice = [1,S(3)];
    
    % Create slice plot and immediately set properties
    set(slice(current_ax, X,Y,Z,permute(A,[2 1 3]),xslice,yslice,zslice), 'EdgeColor', 'none');
    
    % Set axis properties
    axis(current_ax, 'equal', 'tight');
    
    % Add colorbar
    c = colorbar(current_ax);
    c.Label.String = Title;
    
    % Plot cube outline
    plotcube([min(X(:)),max(X(:)),min(Y(:)),max(Y(:)),min(Z(:)),max(Z(:))]);
    
    % Set view angle
    view(current_ax, 30, 40);
end