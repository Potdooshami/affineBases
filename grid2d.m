classdef grid2d
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        offset_stdpts (1,2) double = [nan nan] %[o-p]_pC
        a1_stdpts (1,2) double = [nan nan] % [a_1]_pC
        a2_stdpts (1,2) double = [nan nan] % [a_2]_pC
        sz_stdpts (1,2) double = [nan nan] % [a_3]
        clr_a1 = [1 0 0]
        clr_a2 = [0 1 0]
        clr_faceWindow = [0 1 1];
        clr_edgeWindow = [0 0 1];
    end
    properties (Dependent)
        augMat (3,3) table % [I]_pC^oB
        baryBss (2,3) table
        inLatt table
        inGrid table
    end
    properties (Dependent,Hidden)
        rect (2,4) double        
    end
    
    methods % about dependent properties
        function v = get.augMat(obj)
            v = [obj.a1_stdpts' obj.a2_stdpts' obj.offset_stdpts';0 0 1];
            v = array2table(v,"RowNames",{'x','y','aug'} ,"VariableNames",{'a1','a2','o'});
        end
        function v = get.baryBss(obj)%barycentric basis            
            v = [[0;0] obj.a1_stdpts' obj.a2_stdpts']-obj.offset_stdpts';
            v = array2table(v,"RowNames",{'x','y'} ,"VariableNames",{'p0(=o)','p1(=o+v1)','p2(=o+v2)'});
        end
        function v = get.rect(obj)
            sz = obj.sz_stdpts;
            v = [0 sz(1) sz(1) 0;
                0 0 sz(2) sz(2)];
        end
        function v = get.inLatt(obj)
            isnan =  any(isnan(table2array(obj.augMat)),"all");
            if isnan
            else
                I_oB__pC = inv(table2array(obj.augMat)) ;
            v =  fnd_Lattice_in_window(I_oB__pC,obj.rect);
            end            
        end
        function v = get.inGrid(obj)
            isnan =  any(isnan(table2array(obj.augMat)),"all");
            if isnan
            else
                I_oB__pC = inv(table2array(obj.augMat)) ;
                tbl_line = fnd_Grid_in_window(I_oB__pC,obj.rect);
            end            
        end
    end   
    methods % calculate graphic object in Window
    end
    methods % plot graphic object
        function pWindow(obj)
            fill(obj.rect(1, :), obj.rect(2, :), obj.clr_faceWindow, 'FaceAlpha', 0.3, 'EdgeColor', obj.clr_edgeWindow, 'LineWidth', 1.5);
        end
        function pBases(obj)
            O = obj.offset_stdpts;
            b1 = obj.a1_stdpts;
            b2 = obj.a2_stdpts;
            quiver(O(1), O(2), b1(1), b1(2), 0,'Color' ,obj.clr_a1, 'LineWidth', 2, 'MaxHeadSize', 0.5); % a1 벡터
            hold on
            quiver(O(1), O(2), b2(1), b2(2), 0,'Color' ,obj.clr_a2, 'LineWidth', 2, 'MaxHeadSize', 0.5); % a2 벡터
            hold off
        end
        function pLatticeWin(obj)
        end
        function pGridlineWin(obj)
        end

    end
end