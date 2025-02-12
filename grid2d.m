classdef grid2d
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        offset_stdpts (1,2) double = [nan nan]
        a1_stdpts (1,2) double = [nan nan]
        a2_stdpts (1,2) double = [nan nan]
        sz_stdpts (1,2) double = [nan nan]
        clr_a1 (1,3) double = [1,0,0]
        clr_a2 (1,3) double = [0,1,0]
    end
    properties (Dependent)
        augMat (3,3) table
        baryBss (2,3) table
    end

    methods
        function v = get.augMat(obj)
            v = [obj.a1_stdpts' obj.a2_stdpts' obj.offset_stdpts';0 0 1];
            v = array2table(v,"RowNames",{'x','y','aug'} ,"VariableNames",{'a1','a2','o'});
        end
        function v = get.baryBss(obj)%barycentric basis
            v = [[0;0] obj.a1_stdpts' obj.a2_stdpts']-obj.offset_stdpts';
            v = array2table(v,"RowNames",{'x','y'} ,"VariableNames",{'p0(=o)','p1(=o+v1)','p2(=o+v2)'});
        end
        function lattice_coords_in = fnd_Lattice_in_window(obj)
            I_oB__pC = table2array(obj.augMat) ;
            sz = obj.sz_stdpts;
            rect = [0 sz(1) sz(1) 0;
                0 0 sz(2) sz(2)]
            lattice_coords_in =  fnd_Lattice_in_window(I_oB__pC,rect)
        end



        function visBss(obj)
            hold on
            for ibasis = 1:2
                arr = obj.baryBss(:,[1 1+ibasis]);
                plot(arr(1,:),arr(2,:))
            end
            hold off
        end
        % function tbl = t_augMat(obj)
        %     tbl = array2table(obj.augMat,"RowNames",{'x','y','aug'} ,"VariableNames",{'a1','a2','o'});
        % end
        % function tbl = t_baryBss(obj)
        %     tbl = array2table(obj.baryBss,"RowNames",{'x','y'} ,"VariableNames",{'p0(=o)','p1(=o+v1)','p2(=o+v2)'});
        % end

    end
end