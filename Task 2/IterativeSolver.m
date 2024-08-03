classdef IterativeSolver < handle
    properties (Access = public)
        mat
        vec
    end


    methods (Access = public)

        function obj = IterativeSolver(mat,vec)
            obj.init(mat,vec)
        end

        function sol = solve(obj)
            sol = pcg(obj.mat,obj.vec);
        end

    end


    methods (Access = private)

        function init(obj,mat,vec)
            obj.mat = mat;
            obj.vec = vec;
        end

    end

end