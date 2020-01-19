function [pointsOut,fig] = getPointsPassiveMove(tol, m, k, alpha, h, varargin)
%�������� ����� ������� �� ����� ����������, ������������� ��������� �
%������������ ��������

createFigure=true(1);

for i=1:(length(varargin)-1)
    if (isa(varargin{i},'string')||isa(varargin{i},'char'))&&...
            (isa(varargin{i+1},'string')||isa(varargin{i+1},'char'))
        switch varargin{i}
            case 'figure'
                switch varargin{i+1}
                    case 'no'
                        createFigure=false(1);
                    case 'yes'
                        createFigure=true(1);
                end
        end
    end
end

pointsOut=[];

g=0;%���������� ������� �����������, �*�, ���� �� ������������

id=uint32(0);

syms x y

%c0 - ���������� � ��������� �������
if g~=0
    c0=m*((sign(g)*alpha - g/m)^2)/(1+m*k^2);
else
    c0=(m*alpha^2)/(1+m*k^2);
end

eqEllipse=m*(x-g/m)^2+y^2==c0;

%����������� � L2
eqns = [x+k*y==alpha-h, eqEllipse];
[sol_x, sol_y]=vpasolve(eqns, [x y]);

while(~isempty(sol_x))
    xd=double(sol_x(end));
    yd=double(sol_y(end));
    
    sol_x=sol_x(1:end-1);
    sol_y=sol_y(1:end-1);
    
    %���� ������� �����������
    if(~isreal(xd))||(~isreal(yd))
        %���������� �����
        continue
    end
    
    %��������� �����
    %� �������� ������������, F=1, 2-� �����
    p2add=DataPoint(xd,yd,1,2, id, id, '����� ����. ��������');
    id=id+1;%����������� ������� id
    pointsOut=[pointsOut p2add];%��������� ����� � ������� ���������
end

%����������� � L4
eqns = [x+k*y==-(alpha-h), eqEllipse];
[sol_x, sol_y]=vpasolve(eqns, [x y]);

while(~isempty(sol_x))
    xd=double(sol_x(end));
    yd=double(sol_y(end));
    
    sol_x=sol_x(1:end-1);
    sol_y=sol_y(1:end-1);
    
    %���� ������� �����������
    if(~isreal(xd))||(~isreal(yd))
        %���������� �����
        continue
    end
    
    %��������� �����
    %� �������� ������������, F=-1, 4-� �����
    p2add=DataPoint(xd,yd,-1,4, id, id, '����� ����. ��������');
    id=id+1;%����������� ������� id
    pointsOut=[pointsOut p2add];%��������� ����� � ������� ���������
end

if ~createFigure
    fig=[];
    return;
end

%����������
fig=figure;

xPlot=[-sqrt(c0/m):tol:sqrt(c0/m) sqrt(c0/m)]+g/m;

yPoz=sqrt(c0-m*(xPlot-g/m).^2);

%������� ����������� ��������
for i=length(xPlot):-1:1
    if ~isreal(yPoz)
        xPlot=removeElementByIndex(xPlot, i);
        yPoz=removeElementByIndex(yPoz, i);
    end
end

yNeg=-yPoz;

plot(xPlot, yPoz,'Color', [0.3010 0.7450 0.9330]); hold on;
plot(xPlot, yNeg,'Color', [0.3010 0.7450 0.9330]);
end

