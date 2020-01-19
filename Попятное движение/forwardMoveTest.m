xLims=[-pi pi];%рад
yLims=[-1 1];%рад/с

%—ортируем по возрастанию
xLims=sort(xLims);
yLims=sort(yLims);

maxSteps=10;%ќграничение по числу обработанных точек

a=1;%эффективность управл€ющего воздействи€, рад/с2
m=0.1;%эффективность гравитационного момента, рад/с2
k=0.5;%коэф. при скорости (x+ky)
alpha=deg2rad(1);%точность угловой стабилизаци, рад
h=deg2rad(0.1);%ширина петли гистерезиса, рад
tol=deg2rad(0.01);%погрешность решени€

[passivePoints, ff]=getPointsPassiveMove(tol, m, k, alpha, h);

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

p2s=[DataPoint(0.017,0.01,1,1,uint32(0),uint32(0),"Ќачальна€ точка")];

[pp2e, pp2s, i] = forwardMove(p2s, tol, a, m, k, alpha, h, maxSteps, yLims, func2end);

plotDataPoints([pp2e pp2s], 10*tol, a, m, k, alpha, h, ff)