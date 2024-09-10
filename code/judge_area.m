function flag = judge_area(p,i)
% 通过把手前一时刻所处区域和当前位置判断把手当前时刻所处区域
    global p_inter2;

    flag = i;
    
    if i==1 & norm(p)<4.5
        flag = 2;
    end

    if i==2 & p(1)>p_inter2(1)
        flag = 3;
    end

    if i==3 & norm(p)>=4.5
        flag = 4;
    end

end