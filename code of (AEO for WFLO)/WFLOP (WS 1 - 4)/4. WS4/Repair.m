function Position = Repair(pos, Nt, Cons_L)
    Position = pos;
    dim = size(Position);
    actual_N = sum(Position);
    
    if actual_N > Nt %% ���ڹ涨�Ļ�������
        while 1
            %%%  ���ѡ��һ��ά�ȣ��ж���0����1����Ϊ1��ֵΪ0����Ϊ0�����ѡ������ά��
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
    elseif actual_N < Nt  %% ���ڹ涨�Ļ�������
        while 1
            %%%  ���ѡ��һ��ά�ȣ��ж���0����1����Ϊ0��ֵΪ1����Ϊ1�����ѡ������ά��
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
