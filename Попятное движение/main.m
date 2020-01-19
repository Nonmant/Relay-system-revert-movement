%Программа для определения траекторий попятного движения системы (подробное
%описание системы в документации).
%На входе задаётся массив начальных значений points2Start, состоящий из
%точек класса DataPoint. На выходе - массив points2End той же структуры, в
%который записаны все точки переключения, которые удалось найти программе.
%
%Каждая точка имеет уникальный идентификатор (DataPoint.id), а также
%идентификатор предыдущей точки (DataPoint.prevId), с помощью которых может
%быть восстановлена последовательность точек. Для начальных точек
%идентификаторы начинаются с 0 и идут по очереди.
%
%Для графического представления результатов используется функция
%plotDataPoints
%
%Линии переключения:
%\ \    \ \
% \ \    \ \
%  \ \    \ \
%  L3 L4  L2 L1

xLims=[-pi pi];%рад
yLims=[-1 1];%рад/с

%Сортируем по возрастанию
xLims=sort(xLims);
yLims=sort(yLims);

points2End=[];%Выходной массив
maxSteps=250;%Ограничение по числу обработанных точек
count=0;

a=1;%эффективность управляющего воздействия, рад/с2
m=0.1;%эффективность гравитационного момента, рад/с2
k=0.5;%коэф. при скорости (x+ky)
alpha=deg2rad(1);%точность угловой стабилизаци, рад
h=deg2rad(0.1);%ширина петли гистерезиса, рад
tol=deg2rad(0.005);%погрешность решения

%начальные точки:

%{
%Первая начальная точка
p1=DataPoint(alpha-h-k*(0.1),...%x
    0.1,...%y
    1,...%f
    2,...%l num
    uint32(0),...%id
    uint32(0),...%prev id
    'Начальная точка 1');fig=figure;
%Вторая начальная точка
%p2=DataPoint(alpha-h-k*(-0.1),-0.1,1,2,uint32(1),uint32(1), 'Начальная точка 2');
points2Start=[p1];%points2Start=[p1 p2];

%points2Start=[p1];%только для тестирования 
%}

[points2Start,fig]=getPointsPassiveMove(tol, m, k, alpha, h);

id=uint32(length(points2Start));

%Отображать ли текущие точки на графике
simplePlot=false(1);%true(1);
%Речь идёт о предварительном графике, наиболее наглядный строится функцией plotDataPoints

%системные переменные
syms x y

if(simplePlot) plot([],[]); end
if(simplePlot) hold('on'); end

while(count<maxSteps)&&(~isempty(points2Start))
          
    count=count+1;
    point = points2Start(1,1);%берём первую точку
    points2Start=removeElementByIndex(points2Start,1);%удаляем её из очереди на обработку
        
    %если точка не лежит в заданных пределах y
    if (point.y>yLims(end))||(point.y<yLims(1))
        %выкидываем точку
        warning(strcat("Точка [",num2str(point.x,4), " ", num2str(point.y,4),...
            "] лежит вне заданных диапазонов y:[",...
            num2str(yLims(1),4), ":", num2str(yLims(end),4),...
            "]\nНе обрабатываем эту точку"),...
            "")
        
        if isempty(point.comm)
            point.comm='Точка вне пределов y';
        else
            point.comm=char(string(point.comm) + newline +"Точка вне пределов y");
        end
        
        points2End=[points2End point];
        continue
    end
    
    %предотвращение повторения точек
    minLenSq=tol^2;
    removePoint=false(1);
    for i=1:length(points2End)
        %сравниваем точки с одинаковым значением f
        if point.f~=points2End(i).f
            continue
        end
        
        len=(k*(point.y-points2End(i).y))^2 + (point.x - points2End(i).x)^2;
        if len<=minLenSq
            %выкидываем точку
            removePoint=true(1);
                break
        end
    end
    
    if removePoint
        warning(strcat("Точка слишком близка к уже обработанной\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
        continue
    end
    
    if (point.lNum == 1)%начало на L1
        
        if (point.f == -1)
            %выкидываем точку
            warning(strcat("Точка попала на L1 по траектории F=-1\nЭто противоречит логике\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
        end
        
        
        %при активном управляющем воздействии
        if (point.f == 1)
            
            % ищем пересечение с верхней границей
            eqns = [y==yLims(end), y^2 + a*x==point.y^2+a*point.x];
            [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
            
            %если не нашли - выкидываем точку
            if(isempty(sol_x))
                warning(strcat("Не получилось попасть с верхней границы на L1\nЭто противоречит логике\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
            end
            
            
            while(~isempty(sol_x))%если есть решения
                
                xd=double(sol_x(end));
                yd=double(sol_y(end));
                
                sol_x=sol_x(1:end-1);
                sol_y=sol_y(1:end-1);
                
                %если точка ниже оси y
                if( yd<0 )
                    %пропускаем точку
                    continue
                end
                
                %если точка левее линии L1
                if xd+k*yd<alpha
                    %пропускаем точку
                    continue
                end
                
                %добавляем конечную точку
                %с текущими координатами, F=1, 0-я линия (не на линии), id
                %- текущее значение, prevId - id текущей точки
                p2add=DataPoint(xd,yd,1,0,id,point.id, "Пересечение с y_{max}");
                id=id+1;%увеличиваем счётчик id
                points2End=[points2End p2add];%добавляем точку сразу к обработанным
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %считаем точку обработанной
            points2End=[points2End point];
            
            continue
        end
        
        %при пассивном управляющем воздействии
        if(point.f==0)
            
            %ищем пересечение с линией выключения L4
            eqns = [x+k*y==-(alpha-h), y^2 + m*x^2==point.y^2 + m*point.x^2];
            [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
            
            %цикл по всем решениям
            while(~isempty(sol_x))%если есть решения
                
                xd=double(sol_x(end));
                yd=double(sol_y(end));
                
                sol_x=sol_x(1:end-1);
                sol_y=sol_y(1:end-1);
                
                %если решение комплексное
                if(~isreal(xd))||(~isreal(yd))
                    %пропускаем точку
                    continue
                end
                
                %если слишком близко к начальной точке
                if( (xd-point.x)^2 + (yd-point.y)^2<tol^2)
                    %пропускаем точку
                    continue
                end
                
                %если выше верхней границы
                if(yd>yLims(end))
                    %пропускаем точку
                    continue
                end
                
                %если левее левой границы - пока не используем
%                 if(xd<xLims(1))
%                     %пропускаем точку
%                     continue
%                 end
                
                %если ниже оси y
                if(yd<0)
                    %пропускаем точку
                    continue
                    %спорно, возможно, есть варианты
                end
                
                %если знак y меняется
                if sign(yd)~=sign(point.y)
                    %пропускаем точку
                    continue
                    %см. описание
                end
                
                %добавляем точку
                %с текущими координатами, F=-1, 4-я линия
                p2add=DataPoint(xd,yd,-1,4, id, point.id);
                id=id+1;%увеличиваем счётчик id
                points2Start=[points2Start p2add];%добавляем точку в очередь обработки
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %считаем точку обработанной
            points2End=[points2End point];
            
            continue
        end
    end
    
    if (point.lNum ==2)%начало на L2
        %при f!=1 управляющем воздействии
        if (point.f ~= 1)
            %выкидываем точку
            warning(strcat("Точка попала на L2 не по траектории F=1\nЭто противоречит логике\nВыкидываем точку [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        % ищем пересечение с линией включения L1
        eqns = [x+k*y==alpha, y^2 + a*x==point.y^2+a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %если не нашли - выкидываем точку
        if(isempty(sol_x))
            warning(strcat("Не получилось попасть с на L1 на L2\nЭто противоречит логике\nВыкидываем точку [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %убираем неподходящие точки
        
        iCount=length(sol_x);
        while iCount>0
            
            xd=double(sol_x(iCount));
            yd=double(sol_y(iCount));
            
            %если точка ниже начальной или выше верхнего предела
            if (yd<point.y)||(yd>yLims(end))
                %удаляем точку
                sol_x=removeElementByIndex(sol_x,iCount);
                sol_y=removeElementByIndex(sol_y,iCount);
            end
            iCount=iCount-1;
        end
        
        %если все точки оказались неподходящими
        if(isempty(sol_x))
            warning(strcat("Все переходы с на L1 на L2 неподходящие\nЭто противоречит логике\nВыкидываем точку [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %цикл по всем подходящим решениям
        while(~isempty(sol_x))%если есть решения
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %добавляем точку
            
            %с текущими координатами, F=1, 1-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,1,1,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            %с текущими координатами, F=0, 1-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,0,1,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        
        continue
    end
    
    if (point.lNum ==4)%начало на L4
        %при f!=-1 управляющем воздействии
        if (point.f ~= -1)
            %выкидываем точку
            warning(strcat("Точка попала на L4 не по траектории F=-1\nЭто противоречит логике\nВыкидываем точку [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        % ищем пересечение с линией включения L3
        eqns = [x+k*y==-alpha, y^2 - a*x==point.y^2 - a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %если не нашли - выкидываем точку
        if(isempty(sol_x))
            warning(strcat("Не получилось попасть с L3 на L4\nЭто противоречит логике\nВыкидываем точку [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %убираем неподходящие точки
        
        iCount=length(sol_x);
        while iCount>0
            
            xd=double(sol_x(iCount));
            yd=double(sol_y(iCount));
            
            %если точка выше начальной или ниже нижнего предела
            if (yd>point.y)||(yd<yLims(1))
                %удаляем точку
                sol_x=removeElementByIndex(sol_x,iCount);
                sol_y=removeElementByIndex(sol_y,iCount);
            end
            iCount=iCount-1;
        end
        
        %если все точки оказались неподходящими
        if(isempty(sol_x))
            warning(strcat("Все переходы с на L4 на L3 неподходящие\nЭто противоречит логике\nВыкидываем точку [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %цикл по всем подходящим решениям
        while(~isempty(sol_x))%если есть решения
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %добавляем точку
            
            %с текущими координатами, F=-1, 3-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,-1,3,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            %с текущими координатами, F=0, 3-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,0,3,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        
        continue
    end
    
    if (point.lNum ==3)%начало на L3
        if (point.f == 1)
            %выкидываем точку
            warning(strcat("Точка попала на L3 по траектории F=1\nЭто противоречит логике\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
        end
        
        
        %при активном управляющем воздействии
        if (point.f == -1)
            
            % ищем пересечение с нижней границей
            eqns = [y==yLims(1), y^2 - a*x==point.y^2 - a*point.x];
            [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
            
            %если не нашли - выкидываем точку
            if(isempty(sol_x))
                warning(strcat("Не получилось попасть с нижней границы на L3\nЭто противоречит логике\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
            end
            
            
            while(~isempty(sol_x))%если есть решения
                
                xd=double(sol_x(end));
                yd=double(sol_y(end));
                
                sol_x=sol_x(1:end-1);
                sol_y=sol_y(1:end-1);
                
                %если точка выше оси y
                if( yd>0 )
                    %пропускаем точку
                    continue
                end
                
                %если точка правее линии L3
                if xd+k*yd>-alpha
                    %пропускаем точку
                    continue
                end
                
                %добавляем конечную точку
                %с текущими координатами, F=-1, 0-я линия (не на линии), id
                %- текущее значение, prevId - id текущей точки
                p2add=DataPoint(xd,yd,-1,0,id,point.id,"Пересечение с y_{min}");
                id=id+1;%увеличиваем счётчик id
                points2End=[points2End p2add];%добавляем точку сразу к обработанным
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %считаем точку обработанной
            points2End=[points2End point];
            
            continue
        end
        
        %при пассивном управляющем воздействии
        if(point.f==0)
            
            %ищем пересечение с линией выключения L2
            eqns = [x+k*y==alpha-h, y^2 + m*x^2==point.y^2 + m*point.x^2];
            [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
            
            %цикл по всем решениям
            while(~isempty(sol_x))%если есть решения
                
                xd=double(sol_x(end));
                yd=double(sol_y(end));
                
                sol_x=sol_x(1:end-1);
                sol_y=sol_y(1:end-1);
                
                %если решение комплексное
                if(~isreal(xd))||(~isreal(yd))
                    %пропускаем точку
                    continue
                end
                
                %если слишком близко к начальной точке
                if( (xd-point.x)^2 + (yd-point.y)^2<tol^2)
                    %пропускаем точку
                    continue
                end
                
                %если ниже нижней границы
                if(yd<yLims(1))
                    %пропускаем точку
                    continue
                end
                
                %если правее правой границы - пока не используем
%                 if(xd>xLims(end))
%                     %пропускаем точку
%                     continue
%                 end
                
                %если выше оси y
                if(yd>0)
                    %пропускаем точку
                    continue
                    %спорно, возможно, есть варианты
                end
                
                %если знак y меняется
                if sign(yd)~=sign(point.y)
                    %пропускаем точку
                    continue
                    %см. описание
                end
                
                %добавляем точку
                %с текущими координатами, F=1, 2-я линия
                p2add=DataPoint(xd,yd,1,2, id, point.id);
                id=id+1;%увеличиваем счётчик id
                points2Start=[points2Start p2add];%добавляем точку в очередь обработки
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %считаем точку обработанной
            points2End=[points2End point];
            
            continue
        end
    end
end
if count>=maxSteps
    warning(strcat("Программа завершена по достижении ",num2str(count),...
        " шагов\nФазовый портрет может быть неполным\n"),...
                "")
end
plotDataPoints([points2End points2Start], 10*tol, a, m, k, alpha, h, fig)
