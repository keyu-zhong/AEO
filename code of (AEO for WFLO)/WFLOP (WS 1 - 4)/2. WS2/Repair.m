function Position = Repair(pos, Nt, Cons_L)
    Position = pos;
    dim = size(Position);
    actual_N = sum(Position);
    
    if actual_N > Nt %% 多于规定的机组数量
        while 1
            %%%  随机选择一个维度，判断是0还是1，若为1则赋值为0；若为0则继续选择其他维度
            rd_ind = randi(dim);
            if Position(rd_ind) == 1
                Position(rd_ind) = 0;
            end
            if sum(Position) == Nt
                break;
            else
                continue;
            end
        end
    elseif actual_N < Nt  %% 少于规定的机组数量
        while 1
            %%%  随机选择一个维度，判断是0还是1，若为0则赋值为1；若为1则继续选择其他维度
            rd_ind = randi(dim);
            if Position(rd_ind) == 0 && ismember(rd_ind,Cons_L) == 0
                Position(rd_ind) = 1;        
            end
            if sum(Position) == Nt
                break;
            else
                continue;
            end
        end
    end
end
