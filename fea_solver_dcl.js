let results = []

//Callback on result reception
function onBatchPropertiesResult (ev) {
  //results[ev.sliceNumber] = ev.result;
  //document.querySelector('pre#results').innerText = JSON.stringify(results, null, 2);

  for (i = 0; i < 3; i++){

    inp = results[i][0].input;
    fea_data = ({E : 120e9 * inp, poisson : 0.4* inp, numeleX : 2, numeleY : 2,  numeleZ : 2,
    Pz : 100e3, Py : 0, Px : 0, Lx : 1* inp, Ly : 0.1* inp, Lz : 0.1* inp});

    document.getElementById('Scen_' + inp + "_Settings").innerHTML = JSON.stringify(fea_data);
    for (j = 0; j < 8; j ++) {
    document.getElementById('Scen_' + results[i][0].input + "_Node_" + results[i][j + 1].coordinates + "_coords").innerHTML = JSON.stringify(fea_data);
    document.getElementById('Scen_' + results[i][0].input + "_Node_" + results[i][j + 1].displacement + "_Displacement").innerHTML = JSON.stringify(fea_data);
    }
  }
}

function fea_solver_dcl_Properties (inp) {
  fea_data = ({E : 120e9 * inp, poisson : 0.4* inp, numeleX : 2, numeleY : 2,  numeleZ : 2,
    Pz : 100e3, Py : 0, Px : 0, Lx : 1* inp, Ly : 0.1* inp, Lz : 0.1* inp});

  E = fea_data.E;
  poisson = fea_data.poisson;
  numeleX = fea_data.numeleX;
  numeleY = fea_data.numeleY;
  numeleZ = fea_data.numeleZ;
  Pz = fea_data.Pz;
  Py = fea_data.Py;
  Px = fea_data.Px;
  Lx = fea_data.Lx;
  Ly = fea_data.Ly;
  Lz = fea_data.Lz;

  Dcoeff = E/((1 + poisson) * (1 - 2 * poisson));
  D = math.dotMultiply(Dcoeff, 
    [[(1 - poisson), poisson, poisson, 0, 0, 0],
    [poisson, (1 - poisson), poisson, 0, 0, 0],
    [poisson, poisson, (1 - poisson), 0, 0, 0],
    [0, 0, 0, ((1 - 2 * poisson)/2), 0, 0],
    [0, 0, 0, 0, ((1 - 2 * poisson)/2), 0],
    [0, 0, 0, 0, 0, ((1 - 2 * poisson)/2)]]);

  numele = numeleX * numeleY * numeleZ;

  fixedPosnsX = [];
  fixedPosnsY = [];
  fixedPosnsZ = [];
  rightBord = [];
  // Set up the coordinates of the nodes in the Global coordinates
  nodecoord = math.zeros(numeleX * numeleY  * numeleZ, 3);
  count = 0;
  for (i = 0; i < numeleX; i ++){
    for (j= 0; j < numeleY; j++){
      for (k = 0; k < numeleZ; k++){
        nodecoord = math.subset(nodecoord, math.index(count, 0), Lx/(numeleX - 1) * i);
        nodecoord = math.subset(nodecoord, math.index(count, 1), Ly/(numeleY - 1) * j);
        nodecoord = math.subset(nodecoord, math.index(count, 2), Lz/(numeleZ - 1) * k);

        if((Lx/(numeleX - 1) * i) == 0){
          fixedPosnsX = fixedPosnsX.concat(count);
        }
        if((Ly/(numeleY - 1) * i) == 0){
          fixedPosnsY = fixedPosnsY.concat(count);
        }
        if((Lz/(numeleZ - 1) * i) == 0){
          fixedPosnsZ = fixedPosnsZ.concat(count);
        }
        if(Lx/(numeleX - 1) * (i) == Lx){
          rightBord = rightBord.concat(count);
        }
        count = count + 1;
      }
    }
  }

  //Set Degree of Freedom
  xx = math.subset(nodecoord, math.index(math.range(0, count),0));
  yy = math.subset(nodecoord, math.index(math.range(0, count),1));
  zz = math.subset(nodecoord, math.index(math.range(0, count),2));
  //numnodes = count - 1;
  GDof = 3 * numele;//numnodes; 

  // Set up the nodes that comprise each element
  elenodes = math.zeros(1, 4);//numele,4);
  count = 0;
  for (x = 1; x < numeleX; x++){
    seed = x + (3 * (x - 1)) - 1;
    elenodes = math.subset(elenodes, math.index(count, 0), seed); 
    elenodes = math.subset(elenodes, math.index(count, 1), seed + 1);
    elenodes = math.subset(elenodes, math.index(count, 2), seed + 3);
    elenodes = math.subset(elenodes, math.index(count, 3), seed + 2);
    elenodes = math.subset(elenodes, math.index(count, 4), seed + 4);
    elenodes = math.subset(elenodes, math.index(count, 5), seed + 5);
    elenodes = math.subset(elenodes, math.index(count, 6), seed + 7);
    elenodes = math.subset(elenodes, math.index(count, 7), seed + 6);
    count = count + 1;
  }


  // Calculate the Stiffness Matrix
  stiffness = math.zeros(GDof, GDof);
  for (e = 0; e < Math.cbrt(numele) - 1; e++){
    local_nodes = elenodes;//math.subset(elenodes, math.index(e,math.range(0, ));
    local_dofx = math.subtract(math.dotMultiply(3, local_nodes), 0);
    local_dofy = math.subtract(math.dotMultiply(3, local_nodes), -1);
    local_dofz = math.subtract(math.dotMultiply(3, local_nodes), -2);
    
    local_nodecoord = nodecoord;
    size = local_nodes.size();
    size = size[1];
    for (i = 0; i < size; i++){
      idx = local_nodes.subset(math.index(0, i));
      local_nodecoord = math.subset(local_nodecoord, math.index(i, math.range(0, 3)), 
        math.subset(nodecoord, math.index(idx, math.range(0, 3))));
    }

    K_element = formStiffness3D(local_nodecoord, D);
      K = math.zeros(GDof, GDof);

      K = K.subset(math.index((local_dofx._data[0]), (local_dofx._data[0])), K_element.subset(math.index(math.range(0, 8),math.range(0, 8))));
      K = K.subset(math.index((local_dofy._data[0]), (local_dofy._data[0])), K_element.subset(math.index(math.range(8, 16),math.range(8, 16))));
      K = K.subset(math.index((local_dofz._data[0]), (local_dofz._data[0])), K_element.subset(math.index(math.range(16,24),math.range(16,24))));
      K = K.subset(math.index((local_dofx._data[0]), (local_dofy._data[0])), K_element.subset(math.index(math.range(0,8),math.range(8, 16))));
      K = K.subset(math.index((local_dofx._data[0]), (local_dofz._data[0])), K_element.subset(math.index(math.range(0, 8),math.range(16,24))));
      K = K.subset(math.index((local_dofy._data[0]), (local_dofx._data[0])), K_element.subset(math.index(math.range(8, 16), math.range(0, 8))));
      K = K.subset(math.index((local_dofy._data[0]), (local_dofz._data[0])), K_element.subset(math.index(math.range(8, 16),math.range(16,24))));
      K = K.subset(math.index((local_dofz._data[0]), (local_dofx._data[0])), K_element.subset(math.index(math.range(16,24), math.range(0, 8))));
      K = K.subset(math.index((local_dofz._data[0]), (local_dofy._data[0])), K_element.subset(math.index(math.range(16,24),math.range(8, 16))));
      stiffness = math.add(stiffness, K);  
    }

    // Setting up BCs
    BC = math.zeros(GDof,1);

    fixedPosnsX0 = fixedPosnsX;
    fixedPosnsX = math.add(fixedPosnsX0, 1);
    fixedPosnsY = math.add(fixedPosnsX0, 1);
    fixedPosnsZ = math.add(fixedPosnsX0, 1);

    fixedPosnsX = math.subtract(math.multiply(math.matrix(fixedPosnsX), 3), 2);
    fixedPosnsY = math.subtract(math.multiply(math.matrix(fixedPosnsY), 3), 1);
    fixedPosnsZ = math.multiply(math.matrix(fixedPosnsZ), 3);

    fixedPosnsX = math.subtract(fixedPosnsX, 1);
    fixedPosnsY = math.subtract(fixedPosnsY, 1);
    fixedPosnsZ = math.subtract(fixedPosnsZ, 1);

    for(i = 0; i < fixedPosnsX.size(); i++){
      BC = BC.subset(math.index(fixedPosnsX.get([i]), 0), 1); 
    }
    for(i = 0; i < fixedPosnsY.size(); i++){
      BC = BC.subset(math.index(fixedPosnsY.get([i]), 0), 1); 
    }
    for(i = 0; i < fixedPosnsZ.size(); i++){
      BC = BC.subset(math.index(fixedPosnsZ.get([i]), 0), 1); 
    }
  //Force

  force = math.zeros(GDof,1);
  //rightBord = find(nodecoord(:,1) == Lx);

  rightBord1 = math.add(rightBord, 1);
  rightBord1 = math.subtract(math.multiply(rightBord1, 3), 2);
  rightBord1 = math.subtract(rightBord1, 1);
  for(i = 0; i < rightBord1.length; i++){
    force = force.subset(math.index(rightBord1[i], 0), Px/4); 
  }
  rightBord1 = math.add(rightBord, 1);
  rightBord1 = math.subtract(math.multiply(rightBord1, 3), 1);
  rightBord1 = math.subtract(rightBord1, 1);
  for(i = 0; i < rightBord1.length; i++){
    force = force.subset(math.index(rightBord1[i], 0), Py/4); 
  }
  rightBord1 = math.add(rightBord, 1);
  rightBord1 = math.multiply(rightBord1, 3);
  rightBord1 = math.subtract(rightBord1, 1);
  for(i = 0; i < rightBord1.length; i++){
    force = force.subset(math.index(rightBord1[i], 0), Pz/4); 
  }


  K_BC = [];

  for (i = 0; i < (BC.size())[0]; i++){
    if (BC.subset(math.index(i, 0))){
      K_BC = K_BC.concat(i);
    }
  }

  // Reduce and Solve for Displacement
  //K_BC = find(BC);


  stiffness = stiffness.subset(math.index(
    K_BC, math.range(0, (stiffness.size())[1])), 
  math.zeros(K_BC.length, (stiffness.size())[1])._data);  

  stiffness = stiffness.subset(math.index(
    math.range(0, (stiffness.size())[0]), K_BC), 
  math.zeros(stiffness.size()[0], K_BC.length));

  for (i = 0; i < K_BC.length; i++){
    stiffness = stiffness.subset(math.index(K_BC[i], K_BC[i]), 1); 
  }

  displacements = math.multiply(math.inv(stiffness), force); 
  nodesLocations = nodecoord._data;
1
  return [{input: inp}, {coordinates: nodesLocations[0], displacement: displacements.slice(0, 3)}, 
  {coordinates: nodesLocations[1], displacement: displacements.slice(3, 6)},
  {coordinates: nodesLocations[2], displacement: displacements.slice(6, 9)},
  {coordinates: nodesLocations[3], displacement: displacements.slice(9, 12)},
  {coordinates: nodesLocations[4], displacement: displacements.slice(12, 15)},
  {coordinates: nodesLocations[5], displacement: displacements.slice(15, 18)},
  {coordinates: nodesLocations[6], displacement: displacements.slice(18, 21)}
  ];
}



//Controller for scale
async function batchPropertyChanges() {

  let generator = compute.for(1, 3 , fea_solver_dcl_Properties);
  generator.requires('./libraries/math.js', "./formStiffness3D.js");
  generator.on('batchPropertyResult', onBatchPropertiesResult);
  let batchPropertyResults = await generator.localExec();
  document.querySelector('pre#results').innerText = JSON.stringify(batchPropertyResults, null, 2);
}
