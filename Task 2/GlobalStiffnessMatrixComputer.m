classdef GlobalStiffnessMatrixComputer < handle
	properties (Access = public)
		K
	end

	properties (Access = private)
		Kel
		Td
		nDof
	end

	methods (Access = public)
		
		function obj = GlobalStiffnessMatrixComputer(Kel,Td,nDof)
			obj.init(Kel,Td,nDof);
		end

		function compute(obj)
            obj.K = zeros(obj.nDof);
			nel = size(obj.Kel,3);
            for e=1:nel
      		  	obj.K(obj.Td(e,:),obj.Td(e,:)) = obj.K(obj.Td(e,:),obj.Td(e,:)) + obj.Kel(:,:,e);
            end
        end
	end
	
	methods (Access = private)
			
		function init(obj,Kel,Td,nDof)
			obj.Kel  = Kel;
			obj.Td   = Td;
			obj.nDof = nDof; 
		end
    end
end
