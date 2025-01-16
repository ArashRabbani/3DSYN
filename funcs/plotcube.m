        function []=plotcube(c,thick,color)
            if nargin==1; thick=2; color=[0,0,0]; end
            if nargin==2; color=[0,0,0]; end
            line([c(1),c(2)],[c(3),c(3)],[c(5),c(5)],'LineWidth',thick,'color',color); line([c(2),c(2)],[c(3),c(4)],[c(5),c(5)],'LineWidth',thick,'color',color); line([c(2),c(1)],[c(4),c(4)],[c(5),c(5)],'LineWidth',thick,'color',color); line([c(1),c(1)],[c(4),c(3)],[c(5),c(5)],'LineWidth',thick,'color',color);
            line([c(1),c(1)],[c(3),c(3)],[c(5),c(6)],'LineWidth',thick,'color',color); line([c(2),c(2)],[c(3),c(3)],[c(5),c(6)],'LineWidth',thick,'color',color); line([c(2),c(2)],[c(4),c(4)],[c(5),c(6)],'LineWidth',thick,'color',color); line([c(1),c(1)],[c(4),c(4)],[c(5),c(6)],'LineWidth',thick,'color',color);
            line([c(1),c(2)],[c(3),c(3)],[c(6),c(6)],'LineWidth',thick,'color',color); line([c(2),c(2)],[c(3),c(4)],[c(6),c(6)],'LineWidth',thick,'color',color); line([c(2),c(1)],[c(4),c(4)],[c(6),c(6)],'LineWidth',thick,'color',color); line([c(1),c(1)],[c(4),c(3)],[c(6),c(6)],'LineWidth',thick,'color',color);
        end