function formStiffness3D(local_nodecoord, D){
	K_element = math.zeros(24, 24);
	GaussPoint = [-0.577350269189626, 0.577350269189626];

	nrows = math.range(0, local_nodecoord.size()[0]);
	elenodecoordx = local_nodecoord.subset(math.index(nrows, 0));
	elenodecoordy = local_nodecoord.subset(math.index(nrows, 1));
	elenodecoordz = local_nodecoord.subset(math.index(nrows, 2));


for (p = 0; p < GaussPoint.length; p++){
    for (q = 0; q < GaussPoint.length; q++){
        for (r = 0; r < GaussPoint.length; r++){
        	//console.log("r " + r);
        	//console.log("GPl " +  GaussPoint.length);

            xi = GaussPoint[p];
            eta = GaussPoint[q];
            sigma = GaussPoint[r];
            //console.log("r Second Look " + r);
            
            // Shape Function derivatives
            dShape = math.dotMultiply((1/8), 
            		math.matrix([[-(1-eta)*(1-sigma), -(1-eta)*(1+sigma), -(1+eta)*(1+sigma), 
            					-(1+eta)*(1-sigma),	(1-eta)*(1-sigma), (1-eta)*(1+sigma), 
            					(1+eta)*(1+sigma), (1+eta)*(1-sigma)],
                            	[-(1-xi)*(1-sigma), -(1-xi)*(1+sigma), (1-xi)*(1+sigma),   
                            	(1-xi)*(1-sigma), -(1+xi)*(1-sigma), -(1+xi)*(1+sigma),  
                            	(1+xi)*(1+sigma), (1+xi)*(1-sigma)],
                                [-(1-xi)*(1-eta), (1-xi)*(1-eta), (1-xi)*(1+eta), 
                                -(1-xi)*(1+eta), -(1+xi)*(1-eta), (1+xi)*(1-eta),    
                                (1+xi)*(1+eta), -(1+xi)*(1+eta)]])); 
          // Jacobians
            nColRange = math.range(0, dShape.size()[1]);

            Jacob = [[Number((math.multiply(dShape.subset(math.index(0, nColRange)), elenodecoordx))._data), 
            		 Number(math.multiply(dShape.subset(math.index(0,nColRange)), elenodecoordy)._data), 
            		 Number(math.multiply(dShape.subset(math.index(0,nColRange)), elenodecoordz)._data)],
                     [Number(math.multiply(dShape.subset(math.index(1,nColRange)), elenodecoordx)._data), 
                     Number(math.multiply(dShape.subset(math.index(1,nColRange)), elenodecoordy)._data), 
                     Number(math.multiply(dShape.subset(math.index(1,nColRange)), elenodecoordz)._data)],
                     [Number(math.multiply(dShape.subset(math.index(2,nColRange)), elenodecoordx)._data), 
                     Number(math.multiply(dShape.subset(math.index(2,nColRange)), elenodecoordy)._data), 
                     Number(math.multiply(dShape.subset(math.index(2,nColRange)), elenodecoordz)._data)]];
            invJacobian = math.inv(math.matrix(Jacob));
            //disp(invJacobian);
            //disp(dShape);
            XYZderivatives =  math.multiply(invJacobian, dShape);
            
            // B Matrix
            nColRange = math.range(0, XYZderivatives.size()[1]);

            Brow1 = XYZderivatives.subset(math.index(0,nColRange))._data.concat(math.zeros(1,8)._data, math.zeros(1,8)._data);
            Brow2 = math.zeros(1,8)._data.concat(XYZderivatives.subset(math.index(1,nColRange))._data, math.zeros(1,8)._data);
            Brow3 = math.zeros(1,8)._data.concat(math.zeros(1,8)._data, XYZderivatives.subset(math.index(2,nColRange))._data);
            Brow4 = XYZderivatives.subset(math.index(1,nColRange))._data.concat(XYZderivatives.subset(math.index(0,nColRange))._data, math.zeros(1,8)._data);
            Brow5 = math.zeros(1,8)._data.concat(XYZderivatives.subset(math.index(2,nColRange))._data, XYZderivatives.subset(math.index(1,nColRange))._data);
            Brow6 = XYZderivatives.subset(math.index(2,nColRange))._data.concat(math.zeros(1,8)._data, XYZderivatives.subset(math.index(0,nColRange))._data); 
            B = [[].concat.apply([], Brow1),
                 [].concat.apply([], Brow2),
                 [].concat.apply([], Brow3),
                 [].concat.apply([], Brow4),
                 [].concat.apply([], Brow5),
                 [].concat.apply([], Brow6)];
            //B = math.matrix(B);
           // B = math.transpose(B);
            // Calculate the Stiffness Matrix
            K_element = math.add(K_element, math.multiply(math.transpose(B), D, B, math.det(Jacob)));                                                                           
            }
    }
}
return K_element;
}


