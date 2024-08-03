classdef Solver < handle
    properties (Access = public)
        mat
        vec
        n % 0: direct solver. 1: iterative solver
    end

    properties (Access = private)
        directSolver
        iterativeSolver
    end


    methods (Access = public)

        function obj = Solver(mat,vec)
            obj.init(mat,vec)
            obj.createDirectSolver();
            obj.createIterativeSolver();
        end

        function sol = solve(obj,n)
            switch n
                case 0
                    sol = obj.directSolver.solve();
                case 1
                    sol = obj.iterativeSolver.solve();
            end
        end

    end


    methods (Access = private)

        function init(obj,mat,vec)
            obj.mat = mat;
            obj.vec = vec;
        end

        function createDirectSolver(obj) % private method that creates an instance of the class DirectSolver as an attribute of this class Solver
            obj.directSolver = DirectSolver(obj.mat,obj.vec);
        end
        
        function createIterativeSolver(obj) % private method that creates an instance of the class IterativeSolver as an attribute of this class Solver
            obj.iterativeSolver = IterativeSolver(obj.mat,obj.vec);
        end
       
    end

end