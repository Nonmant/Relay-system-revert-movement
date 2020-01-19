%»сследование зависимости длительности перехода в либрацию от начального
%положени€

xLims=[-pi pi];%рад
yLims=[-1 1];%рад/с

%—ортируем по возрастанию
xLims=sort(xLims);
yLims=sort(yLims);

maxSteps=50;%ќграничение по числу обработанных точек

a=1;%эффективность управл€ющего воздействи€, рад/с2
m=0.1;%эффективность гравитационного момента, рад/с2
k=0.5;%коэф. при скорости (x+ky)
alpha=deg2rad(1);%точность угловой стабилизаци, рад
h=deg2rad(0.1);%ширина петли гистерезиса, рад
tol=deg2rad(0.01);%погрешность решени€

[passivePoints, ff]=getPointsPassiveMove(tol, m, k, alpha, h, 'figure', 'no');

ppL2=[];
ppL4=[];

for i=1:length(passivePoints)
    if passivePoints(i).lNum==2
        ppL2=[ppL2 passivePoints(i)];
        continue
    end
    
    if passivePoints(i).lNum==4
        ppL4=[ppL4 passivePoints(i)];
        continue
    end 
end

if (length(ppL2)==2)
    ppL2yLims=sort([ppL2(1).y ppL2(2).y]);
    func2endL2 = @(x,y)...
        (abs( x+k*y - (alpha-h))>tol)*10*tol+...%проверка нахождени€ на линии
        (y<ppL2yLims(1))*10*tol+...%y > нижней границы
        (y>ppL2yLims(end))*10*tol;%y < верхней границы
else
    func2endL2 = @(x,y)10*tol;
end

if (length(ppL4)==2)
    ppL4yLims=sort([ppL4(1).y ppL4(2).y]);
    func2endL4 = @(x,y)...
        (abs( x+k*y + (alpha-h))>tol)*10*tol+...%проверка нахождени€ на линии
        (y<ppL4yLims(1))*10*tol+...%y > нижней границы
        (y>ppL4yLims(end))*10*tol;%y < верхней границы
else
    func2endL4 = @(x,y)10*tol;
end

func2end=@(x,y)min(func2endL4(x,y),func2endL2(x,y));

p2s=[];

xStart=0.0141:((0.0183-0.0141)/1000):0.0183;
yStart=ones(1, length(xStart))*0.0102;
stepsForXStart=zeros(length(xStart),1);

points2Test=repmat(DataPoint(0,0,1,1,uint32(0),uint32(0),"Ќачальна€ точка"), 1, length(xStart));

for i=1:length(xStart)
    
    points2Test(i).x=xStart(i);
    points2Test(i).y=yStart(i);
    
    p2s=[points2Test(i)];
    
    [~,~,stepsForXStart(i)]= forwardMove(p2s, tol, a, m, k, alpha, h, maxSteps, yLims, func2end, 'quiteMode');
end

f=figure;
f.set('Color', 'white')
axes1=findall(f,'type','axes');
if isempty(axes1)
    axes1 = axes('Parent',f);
end
hold(axes1,'on');
box(axes1,'on');
xlabel('x, рад')
ylabel('„исло шагов')
set(axes1,'FontName','Times New Roman');
hold on
grid minor

plot(xStart, stepsForXStart, 'Color', [0.49,0.18,0.56], 'LineWidth', 1.5)


axes1=findall(fig,'type','axes');
hold(axes1,'on');
plot([xStart(1) xStart(1)], [yStart(1) 15e-3], 'Color', 'black')
plot(xStart, yStart(1)+(stepsForXStart*(15e-3 - yStart(1)))/50, 'Color', [0.49,0.18,0.56])