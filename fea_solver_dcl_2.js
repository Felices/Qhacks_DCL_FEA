let results = []

//Callback on result reception
function onBatchPropertiesResult (ev) {
  //results[ev.sliceNumber] = ev.result;
  //document.querySelector('pre#results').innerText = JSON.stringify(results, null, 2);

  for (i = 0; i < 3; i++){

    inp = results[i][0].input;
    fea_data = ({E : 120e9 * inp, poisson : 0.4* inp, numeleX : 2, numeleY : 2,  numeleZ : 2,
    Pz : 100e3, Py : 0, Px : 0, Lx : 1* inp, Ly : 0.1* inp, Lz : 0.1* inp});

    /*document.getElementById('Scen_' + inp + "_Settings").innerHTML = JSON.stringify(fea_data);
    for (j = 0; j < 8; j ++) {
    document.getElementById('Scen_' + results[i][0].input + "_Node_" + results[i][j + 1].coordinates + "_coords").innerHTML = JSON.stringify(fea_data);
    document.getElementById('Scen_' + results[i][0].input + "_Node_" + results[i][j + 1].displacement + "_Displacement").innerHTML = JSON.stringify(fea_data);
    */}
  }
}

function fea_solver_dcl_Properties (inp) {
  fea_data = ({E : 120e9 * inp, poisson : 0.4* inp, numeleX : 2, numeleY : 2,  numeleZ : 2,
    Pz : 100e3, Py : 0, Px : 0, Lx : 1* inp, Ly : 0.1* inp, Lz : 0.1* inp});

  let E = fea_data.E;
  let poisson = fea_data.poisson;
  let numeleX = fea_data.numeleX;
  let numeleY = fea_data.numeleY;
  let numeleZ = fea_data.numeleZ;
  let Pz = fea_data.Pz;
  let Py = fea_data.Py;
  let Px = fea_data.Px;
  let Lx = fea_data.Lx;
  let Ly = fea_data.Ly;
  let Lz = fea_data.Lz;
  let nodecoordCounter = 0;
  let elenodesCounter = 0;

  numele = numeleX* numeleY*numeleZ;
  GDof = 3 * numele;
  
  outsourceDVals();
  outsourceNodeCoordInitialization();
  outsourceElenodes();

  while (elenodesCounter < numele);
  outsourceDofInitialization();
  while (nodecoordCounter < numele);
  

  //seed function
  function getDVals(E, poisson){
    Dcoeff = E/((1 + poisson) * (1 - 2 * poisson));
    D = math.dotMultiply(Dcoeff, 
    [[(1 - poisson), poisson, poisson, 0, 0, 0],
    [poisson, (1 - poisson), poisson, 0, 0, 0],
    [poisson, poisson, (1 - poisson), 0, 0, 0],
    [0, 0, 0, ((1 - 2 * poisson)/2), 0, 0],
    [0, 0, 0, 0, ((1 - 2 * poisson)/2), 0],
    [0, 0, 0, 0, 0, ((1 - 2 * poisson)/2)]]);

    return {D, Dcoeff};
  }
  function buildDVals(ev){
    D = ev.result.D;
    Dcoeff = ev.result.Dcoeff;
  }
  //outsource function
  async function outsourceDVals(){
    let gen_outsourceDVals = compute.for(E, poisson, getDVals(i, j));
    gen_outsourceDVals.on("outsourceDVals_Complete", buildDVals);
    let D = await gen_outsourceDVals.localExec();
    let Dcoeff = await gen_outsourceDVals.localExec();
  }

  function shardNodeCoordInitialization(xIdx, yIdx, zIdx){
    nodecoord = math.getZeros(1, 3);
    flag_fixedPosX = flag_fixedPosY = flag_fixedPosZ = flag_rightBord = false;

    nodecoord = math.subset(nodecoord, math.index(count, 0), Lx/(numeleX - 1) * xIdx);
    nodecoord = math.subset(nodecoord, math.index(count, 1), Ly/(numeleY - 1) * yIdx);
    nodecoord = math.subset(nodecoord, math.index(count, 2), Lz/(numeleZ - 1) * zIdx); 

    if((Lx/(numeleX - 1) * i) == 0){
      flag_fixedPosX = true;
    }
    if((Ly/(numeleY - 1) * i) == 0){
      flag_fixedPosY = true;
    }
    if((Lz/(numeleZ - 1) * i) == 0){
      flag_fixedPosZ = true;
    }
    if(Lx/(numeleX - 1) * (i) == Lx){
      flag_rightBord = true;
    }

    return {nodecoord:nodecoord, xIdx:xIdx, yIdx:yIdx, zIdx:zIdx, flag_fixedPosX:flag_fixedPosX, flag_fixedPosY:flag_fixedPosY, flag_fixedPosZ:flag_fixedPosZ};
  }

  async function outsourceNodeCoordInitialization(){
    let gen_outsourceNodeCoordInitialization = compute.for(numeleX, numeleY, numeleZ, shardNodeCoordInitialization(i, j, k));
    generator.on('shardNodeCoordInitialization_Complete', nodeCoordInitialization);
    let nodecoords = NaN;
    let fixedPosnsX = [];
    let fixedPosnsY = [];
    let fixedPosnsZ = [];
    let rightBord= [];
    nodecoords = await gen_outsourceNodeCoordInitialization.localExec();
  }

  function nodeCoordInitialization(){
    if isNaN(nodecoords){
      nodecoords = [];   
    }
    nodecoords.push(e.result.nodecoord);
    if(ev.result.flag_fixedPosX){
      fixedPosnsX.push(nodecoords.length);
    }
    if(ev.result.flag_fixedPosY){
      fixedPosnsY.push(nodecoords.length);
    }
    if(ev.result.flag_fixedPosZ){
      fixedPosnsZ.push(nodecoords.length);
    }
    if(ev.result.flag_rightBord){
      rightBord.push(nodecoords.length);
    }
    nodecoordCounter++;
  }
  
  
  async function buildNodeCoords(numeleX, numeleY, numeleZ){
    let xx = math.subset(nodecoord, math.index(math.range(0, count),0));
    let yy = math.subset(nodecoord, math.index(math.range(0, count),1));
    let zz = math.subset(nodecoord, math.index(math.range(0, count),2));  
    return {_numele: numele, _fixedPosnsX: fixedPosnsX, _fixedPosnsY: fixedPosnsY, _fixedPosnsZ: fixedPosnsZ, _rightBord:rightBord, _nodecoord:nodecoord};
  }
  function buildElenodes(ev){
    if isNaN(builtElenodes){
      elenodes = math.zeros(numeleX, 8);
    }
    elenodes = math.subset(builtElenodes, math.index(ev.result.x, math.range(0, 8)), ev.result.elenode.subset(math.index(ev.result.x, math.range(0, 8))));
    elenodesCounter++;
  }

  function shardElenodes(x){
    seed = x + (3 * (x - 1)) - 1;
    elenodes = math.subset(elenodes, math.index(x, 0), seed); 
    elenodes = math.subset(elenodes, math.index(x, 1), seed + 1);
    elenodes = math.subset(elenodes, math.index(x, 2), seed + 3);
    elenodes = math.subset(elenodes, math.index(x, 3), seed + 2);
    elenodes = math.subset(elenodes, math.index(x, 4), seed + 4);
    elenodes = math.subset(elenodes, math.index(x, 5), seed + 5);
    elenodes = math.subset(elenodes, math.index(x, 6), seed + 7);
    elenodes = math.subset(elenodes, math.index(x, 7), seed + 6);
    return {x: x, elenode: elenodes};
  }
  async function outsourceElenodes(){
    let gen_outsourceElenodes = compute.for(0, numeleX, shardElenodes(i));
    gen_outsourceElenodes.on('shardElenodes_Complete', buildElenodes);
    let elenodes = NaN;
    elenodes = await gen_outsourceElenodes.localExec();
  } 

  async function outsourceDofInitialization(){
    let gen_outsourceDofInitialization = compute.for(0, local_nodes.size()[0], 0, local_nodes.size()[1], shardDofInitialization(i, j));
    //let local_nodes = elenodes;
    gen_outsourceDofInitialization.on("shardDofInitialization_Complete", DofInitialization);
    local_dofx = NaN;
    local_dofx = await gen_outsourceDofInitialization.localExec();
    local_dofy = await gen_outsourceDofInitialization.localExec();
    local_dofz = await gen_outsourceDofInitialization.localExec();
  }
  function shardDofInitialization(i, j){
    element = local_nodes.subset(math.index(i, j))
    local_dofx_elem = 3*element;
    local_dofy_elem = 3*element+1;
    local_dofz_elem = 3*element+2;  
    return {local_dofx_elem: local_dofx_elem, local_dofy_elem: local_dofy_elem, local_dofz_elem: local_dofz_elem};
  }
  function DofInitialization(ev){
    if isNaN(local_dofx){
      local_dofx = [];
      local_dofy = [];
      local_dofz = [];   
    }   
    local_dofx.push(ev.result.local_dofx_elem);
    local_dofy.push(ev.result.local_dofy_elem);
    local_dofz.push(ev.result.local_dofz_elem);
  }
  
  function outsourceLocalNodecoordsInitialization(){
    let gen_LocalNodecoordsInitialization = compute.for(0, elenodes.size()[1], shardLocalNodecoordsInitialization(i));
    let local_nodecoord = nodecoords;
    gen_outsourceDofInitialization.on("shardLocalNodecoordsInitialization_Complete", DofInitialization);
    local_dofx = NaN;
    local_dofx = await gen_outsourceDofInitialization.localExec();
    ;
    
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

}

async function K_BCBuilder(stiffness){
  K_BC = [];

  for (i = 0; i < (BC.size())[0]; i++){
    if (BC.subset(math.index(i, 0))){
      K_BC.push(i);
    }
  }

  stiffness = stiffness.subset(math.index(
        K_BC, math.range(0, (stiffness.size())[1])), 
        math.zeros(K_BC.length, (stiffness.size())[1])._data); 
  stiffness = stiffness.subset(math.index(
          math.range(0, (stiffness.size())[0]), K_BC), 
        math.zeros(stiffness.size()[0], K_BC.length));
  for (i = 0; i < K_BC.length; i++){
    stiffness = stiffness.subset(math.index(K_BC[i], K_BC[i]), 1); 
  }
}


async function multiplyDisplacements(invStiffRow, force){
  return math.multiply(invStiffRow, forceEntry);
}



//Controller for scale
async function batchPropertyChanges() {

  let generator = compute.for(1, 3 , fea_solver_dcl_Properties);
  generator.requires('./libraries/math.js', "./formStiffness3D.js");
  generator.on('batchPropertyResult', onBatchPropertiesResult);
  let batchPropertyResults = await generator.localExec();
  document.querySelector('pre#results').innerText = JSON.stringify(batchPropertyResults, null, 2);
}
