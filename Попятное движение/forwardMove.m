function [points2End, points2Start, count] = forwardMove(dataPointsIn, tol, a, m, k, alpha, h, maxSteps, yLims, varargin)
%Функция для прямого движения системы
%   dataPointsIn -  массив точек начала
%   tol - минимальное расстояние между последовательными решениями
%   a - эффективность управляющего воздействия, рад/с2
%   m - эффективность гравитационного момента, рад/с2
%   k - коэф. при скорости (x+ky)
%   alpha - точность угловой стабилизаци, рад
%   h - ширина петли гистерезиса, рад
%   maxSteps - максимальное число шагов
%
%   func2end - функция (x, y), при возврате значения <tol пропускается
%   точка
%
%   points2End - обработанные точки
%   points2Start - необработанные точки
%   count - совершённое число шагов

func2end=@(x,y)10*tol;

quiteMode=false(1);

%Ищем функцию остановки в дополнительных входных параметрах
for i=1:length(varargin)
    if isa(varargin{i},'function_handle')
        if nargin(varargin{i})==2 %если 2 входных аргумента
            func2end=varargin{i}; %заменяем стандартную функцию
        end
    end
    
    if isa(varargin{i},'string')||isa(varargin{i},'char')
        switch varargin{i}
            case 'quiteMode'
                quiteMode=true(1);
        end
    end
    
end

points2Start=dataPointsIn;
points2End=[];
count=0;

id=uint32(length(points2Start));

yLims=sort(yLims);

%системные переменные
syms x y
solFound=false(1);

%Уравнения линий переключения
LEqns=[...
    x+k*y==alpha,...%L1
    x+k*y==alpha-h,...%L2
    x+k*y==-alpha,...%L3
    x+k*y==-(alpha-h)...%L4
    ];
while(count<maxSteps)&&(~isempty(points2Start))
    count=count+1;
    point = points2Start(1,1);%берём первую точку
    points2Start=removeElementByIndex(points2Start,1);%удаляем её из очереди на обработку
    
    if func2end(point.x, point.y) < tol %условие выхода
        if ~quiteMode
            disp("Выход по функции " + func2str(func2end) + newline)
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        continue
    end
    
    %если точка не лежит в заданных пределах y
    if (point.y>yLims(end))||(point.y<yLims(1))
        %выкидываем точку
        if ~quiteMode
            warning(strcat("Точка [",num2str(point.x,4), " ", num2str(point.y,4),...
                "] лежит вне заданных диапазонов y:[",...
                num2str(yLims(1),4), ":", num2str(yLims(end),4),...
                "]\nНе обрабатываем эту точку"),...
                "")
        end
        if isempty(point.comm)
            point.comm='Точка вне пределов y';
        else
            point.comm=char(string(point.comm) + newline +"Точка вне пределов y");
        end
        
        points2End=[points2End point];
        continue
    end
    
    if (point.lNum == 1)%начало на L1
        %при активном управляющем воздействии
        
        solFound=false(1);%флаг нахождения пересечения с L2
        
        % ищем пересечение с L2
        eqns = [LEqns(2), y^2 + a*x==point.y^2+a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
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
            
            %если точка выше исходной
            if yd>point.y
                %пропускаем точку
                continue
            end
            
            %если точка ниже нижней границы
            if yd<yLims(1)
                %пропускаем точку
                continue
            end
            
            %добавляем точку
            %с текущими координатами, F=1, 2-я линия
            p2add=DataPoint(xd,yd,0,2, id, point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            solFound=true(1);
        end
        
        % если было найдено пересечение с L2 - не нужно искать пересечений
        % с нижней границей
        if(solFound)
            %считаем точку обработанной
            points2End=[points2End point];
            continue
        end
        
        % ищем пересечение с нижней границей
        eqns = [y==yLims(1), y^2 + a*x==point.y^2+a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %если не нашли - выкидываем точку
        if(isempty(sol_x))
            if ~quiteMode
                warning(strcat("Не получилось попасть с L1 на нижнюю границу\nЭто противоречит логике\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
            end
            continue
        end
        
        
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
            
            %если точка выше исходной
            if yd>point.y
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
            p2add=DataPoint(xd,yd,1,0,id,point.id, "Пересечение с y_{min}");
            id=id+1;%увеличиваем счётчик id
            points2End=[points2End p2add];%добавляем точку сразу к обработанным
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        
        continue
        
    end
    
    if (point.lNum ==2)%начало на L2
        %при пассивном управляющем воздействии
        
        solFound=false(1);%флаг нахождения пересечения с L1
        
        % ищем пересечение с линией включения L1
        eqns = [LEqns(1), y^2 + m*x^2==point.y^2 + m*point.x^2];
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
            
            %если точка выше начальной
            if yd>point.y
                %пропускаем точку
                continue
            end
            
            %добавляем точку
            
            %с текущими координатами, F=0, 1-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,1,1,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            solFound=true(1);
        end
        
        % если было найдено пересечение с L1 - не нужно искать пересечений
        % с L3
        if(solFound)
            %считаем точку обработанной
            points2End=[points2End point];
            continue
        end
        
        % ищем пересечение с линией включения L3
        eqns = [LEqns(3), y^2 + m*x^2==point.y^2 + m*point.x^2];
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
            
            %если точка выше начальной
            if yd>point.y
                %пропускаем точку
                continue
            end
            
            %добавляем точку
            
            %с текущими координатами, F=-1, 3-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,-1,3,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        
        continue
    end
    
    if (point.lNum == 3)%начало на L3
        %при активном управляющем воздействии
        
        solFound=false(1);%флаг нахождения пересечения с L4
        
        % ищем пересечение с L4
        eqns = [LEqns(4), y^2 - a*x==point.y^2 - a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
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
            
            %если точка ниже исходной
            if yd<point.y
                %пропускаем точку
                continue
            end
            
            %если точка выше верхней границы
            if yd>yLims(end)
                %пропускаем точку
                continue
            end
            
            %добавляем точку
            %с текущими координатами, F=-1, 4-я линия
            p2add=DataPoint(xd,yd,0,4, id, point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            solFound=true(1);
        end
        
        % если было найдено пересечение с L4 - не нужно искать пересечений
        % с верхней границей
        if(solFound)
            %считаем точку обработанной
            points2End=[points2End point];
            continue
        end
        
        % ищем пересечение с верхней границей
        eqns = [y==yLims(end), y^2 - a*x==point.y^2 - a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %если не нашли - выкидываем точку
        if(isempty(sol_x))
            if ~quiteMode
                warning(strcat("Не получилось попасть с L3 на верхнюю границу\nЭто противоречит логике\nВыкидываем точку [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
            end
            continue
        end
        
        
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
            
            %если точка ниже исходной
            if yd<point.y
                %пропускаем точку
                continue
            end
            
            %если точка правее линии L4
            if xd+k*yd>-(alpha-h)
                %пропускаем точку
                continue
            end
            
            %добавляем конечную точку
            %с текущими координатами, F=-1, 0-я линия (не на линии), id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,-1,0,id,point.id, "Пересечение с y_{max}");
            id=id+1;%увеличиваем счётчик id
            points2End=[points2End p2add];%добавляем точку сразу к обработанным
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        
        continue
        
    end
    
    if (point.lNum ==4)%начало на L4
        %при пассивном управляющем воздействии
        
        solFound=false(1);%флаг нахождения пересечения с L3
        
        % ищем пересечение с линией включения L3
        eqns = [LEqns(3), y^2 + m*x^2==point.y^2 + m*point.x^2];
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
            
            %если точка выше начальной
            if yd>point.y
                %пропускаем точку
                continue
            end
            
            %добавляем точку
            
            %с текущими координатами, F=0, 3-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,-1,3,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
            
            solFound=true(1);
        end
        
        % если было найдено пересечение с L3 - не нужно искать пересечений
        % с L1
        if(solFound)
            %считаем точку обработанной
            points2End=[points2End point];
            continue
        end
        
        % ищем пересечение с линией включения L1
        eqns = [LEqns(1), y^2 + m*x^2==point.y^2 + m*point.x^2];
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
            
            %если точка выше начальной
            if yd>point.y
                %пропускаем точку
                continue
            end
            
            %добавляем точку
            
            %с текущими координатами, F=0, 1-я линия, id
            %- текущее значение, prevId - id текущей точки
            p2add=DataPoint(xd,yd,1,1,id,point.id);
            id=id+1;%увеличиваем счётчик id
            points2Start=[points2Start p2add];%добавляем точку в очередь обработки
        end
        
        %считаем точку обработанной
        points2End=[points2End point];
        
        continue
    end
    
end



