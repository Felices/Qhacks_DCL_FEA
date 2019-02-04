console.log("Hello world");
const fea_data = require("./fea_job.json");
const fea_solver = require("./fea_solver.js");
const fs = require("fs");

results = fea_solver.fea_solver(fea_data);

nodesLocations = results[0];
displacements = results[1];

resultData = [{coordinates: nodesLocations[0], displacement: displacements.slice(0, 3)}, 
	{coordinates: nodesLocations[1], displacement: displacements.slice(3, 6)},
	{coordinates: nodesLocations[2], displacement: displacements.slice(6, 9)},
	{coordinates: nodesLocations[3], displacement: displacements.slice(9, 12)},
	{coordinates: nodesLocations[4], displacement: displacements.slice(12, 15)},
	{coordinates: nodesLocations[5], displacement: displacements.slice(15, 18)},
	{coordinates: nodesLocations[6], displacement: displacements.slice(18, 21)}
	];

resultDataJSON = JSON.stringify(resultData);
fs.writeFile("./fea_results.json", resultDataJSON);


const fea_solver_dcl_1 = require("./fea_solver_dcl.js");
