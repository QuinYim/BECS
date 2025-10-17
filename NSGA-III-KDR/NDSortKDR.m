% function [FrontNo,MaxFNo] = NDSortKDR(PopObj,nSort)
% % 使用高斯核诱导距离进行非支配排序
% 
%     N      = size(PopObj,1);
%     %M      = size(PopObj,2);
%     Z_ide  = min(PopObj, [], 1);
%     sigma  = 10; % 根据需要调整 sigma 的值
% 
%     % 计算高斯核函数值
%     sq_dist = sum((PopObj - Z_ide).^2, 2);
%     K_xz = exp(-sq_dist / (2 * sigma^2));
% 
%     % 计算高斯核诱导距离
%     dk   = sqrt(2 * (1 - K_xz));
% 
%     % 计算解之间的角度
%     cosine = 1 - pdist2(PopObj,PopObj,'cosine');
%     cosine(logical(eye(length(cosine)))) = 0;
%     Angle  = acos(cosine);
% 
%     % 计算最小角度
%     temp = sort(min(Angle,[],2));
%     if length(temp) >= floor(N/2)
%         minA = temp(floor(N/2));
%     else
%         minA = temp(end);
%     end
% 
%     Theta = max(1, Angle./minA);
% 
%     % 非支配关系判断
%     dominate = false(N);
%     for i = 1 : N-1
%         for j = i+1 : N
%             if dk(i)*Theta(i,j) < dk(j)
%                 dominate(i,j) = true;
%             elseif dk(j)*Theta(j,i) < dk(i)
%                 dominate(j,i) = true;
%             end
%         end
%     end
% 
%     % 分配前沿编号
%     FrontNo = inf(1,N);
%     MaxFNo  = 0;
%     while sum(FrontNo~=inf) < min(nSort,N)
%         MaxFNo  = MaxFNo + 1;
%         current = ~any(dominate,1) & FrontNo==inf;
%         FrontNo(current)    = MaxFNo;
%         dominate(current,:) = false;
%     end
% end


function [FrontNo, MaxFNo] = NDSortKDR(PopObj, nSort)
% 使用高斯核诱导距离进行非支配排序，并对数据进行归一化

    % 数据归一化
    [PopObj_norm, Z_ide] = normalizeData(PopObj);
    
    N      = size(PopObj_norm, 1);
    % M    = size(PopObj_norm, 2); % 未使用，可根据需要启用
    sigma  = 1; % 根据需要调整 sigma 的值

    % 计算高斯核函数值
    sq_dist = sum((PopObj_norm - Z_ide).^2, 2);
    K_xz = exp(-sq_dist / (2 * sigma^2));

    % 计算高斯核诱导距离
    dk   = sqrt(2 * (1 - K_xz));

    % 计算解之间的角度
    cosine = 1 - pdist2(PopObj_norm, PopObj_norm, 'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    Angle  = acos(cosine);

    % 计算最小角度
    temp = sort(min(Angle, [], 2));
    if length(temp) >= floor(N / 2)
        minA = temp(floor(N / 2));
    else
        minA = temp(end);
    end

    % 计算 Theta 矩阵
    Theta = max(1, Angle ./ minA);

    % 非支配关系判断
    dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            if dk(i) * Theta(i, j) < dk(j)
                dominate(i, j) = true;
            elseif dk(j) * Theta(j, i) < dk(i)
                dominate(j, i) = true;
            end
        end
    end

    % 分配前沿编号
    FrontNo = inf(1, N);
    MaxFNo  = 0;
    while sum(FrontNo ~= inf) < min(nSort, N)
        MaxFNo  = MaxFNo + 1;
        current = ~any(dominate, 1) & FrontNo == inf;
        FrontNo(current)    = MaxFNo;
        dominate(current, :) = false;
    end
end

function [PopObj_norm, Z_ide, objMin, objMax] = normalizeData(PopObj)
% 对 PopObj 进行最小-最大归一化
    objMin = min(PopObj, [], 1);
    objMax = max(PopObj, [], 1);
    range  = objMax - objMin;
    
    % 防止除以零
    range(range == 0) = 1;
    
    PopObj_norm = (PopObj - objMin) ./ range;
    
    % 计算归一化后的理想点 Z_ide
    Z_ide = min(PopObj_norm, [], 1);
end