% Procedure otdm_uo_plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = otdm_uo_plot(f_fun, xk, gk, xylim)
hold off;
grid on;
box on;
set(gca,'Fontsize',14);
hold on;
% Plot area
xmin=xylim(1);
xmax=xylim(2);
ymin=xylim(3);
ymax=xylim(4);
a=.5;
max_xk1   = max(xk(1,:));
min_xk1   = min(xk(1,:));
w_xk1 = max_xk1-min_xk1;
max_xk1   = max_xk1 + a*w_xk1;
min_xk1   = min_xk1 - a*w_xk1;
w_xk1 = max_xk1-min_xk1;

max_xk2   = max(xk(2,:));
min_xk2   = min(xk(2,:));
w_xk2 = max_xk2-min_xk2;
max_xk2   = max_xk2 + a*w_xk2;
min_xk2   = min_xk2 - a*w_xk2;
w_xk2 = max_xk2-min_xk2;

difwidth = w_xk1-w_xk2;
if difwidth > 0
    max_xk2 = max_xk2 + difwidth/2;
    min_xk2 = min_xk2 - difwidth/2;
elseif difwidth < 0
    max_xk1 = max_xk1 - difwidth/2;
    min_xk1 = min_xk1 + difwidth/2;
end

if xmax~=xmin
    if xmax < max_xk1 xmax = max_xk1; end
    if xmin > min_xk1 xmin = min_xk1; end
else
    xmax = max_xk1;
    xmin = min_xk1;
end
if ymax~=ymin
    if ymax < max_xk2 ymax = max_xk2; end
    if ymin > min_xk2 ymin = min_xk2; end
else
    ymax = max_xk2;
    ymin = min_xk2;
end
%
set(gcf,'Color','w');  
set(findobj('Type','line'),'LineWidth',4.0)
x=linspace(xmin,xmax);
y=linspace(ymin,ymax);
[X,Y]=meshgrid(x,y);
[m,sx]=size(x);
[m,sy]=size(y);
Z=zeros(sy,sx);
for i=1:sx
    for j=1:sy
        Z(j,i)=f_fun([x(i);y(j)]);
    end
end
zmax=max(max(Z));
zmin=min(min(Z));
zmin=zmin-0.5*(zmax-zmin);
axis([xmin xmax ymin ymax zmin inf]);
view(45, 20);
surfc(X,Y,Z)
xlabel('x_1');ylabel('x_2');zlabel('f(x)');title('');
[n,iter] = size(xk);
zl=zlim;
z(1:iter)=zl(1);
for i=1:iter
    zf(i)=f_fun(xk(:,i))+0.0*abs(f_fun(xk(:,i)));
end
plot3( xk(1,:)',  xk(2,:)', z(:), ':ok','LineWidth',3)
plot3( xk(1,:)',  xk(2,:)', zf(:), ':oy','LineWidth',6)
quiver3( xk(1,1:iter)', xk(2,1:iter)', z(:), gk(1,1:iter)', gk(2,1:iter)', zeros(iter,1), 'LineWidth',2,'Color','r')    
%if iter > 1, line( xk(1,:)',  xk(2,:)','LineStyle',':','Color','k','LineWidth',3), end
% End Procedure otdm_uo_plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
