
#################### Input Data
#Network Node Data
nodes = pd.DataFrame({
					  'name':['N1', 'N2', 'N3','N4','So','Si'],
                      'pressure':[0,30,0,0,0,1],
					  'flow':[0,0,0,0,+10,-10],
                      'pboundary':[False,True,False,False,False,False],
                      'fboundary':[False,False,False,False,True,True]
})

#Network Edge Data 
edges = pd.DataFrame({
	'name':['PSo', 'P1', 'P2','P3','P4','P5','PSi'], 
    'from':['So','N1','N1','N2','N2','N3','N4'],
    'to':['N1','N2','N3','N3','N4','N4','Si'],
    'flow':[ 10, 5, 5, 0, 5, 5,10],
    'length':[5,10, 10, 10, 10,10,5],
    'diam':[ 0.3, 0.3, 0.1, 0.3, 0.1, 0.3, 0.3 ],
    'pumpHeadGain':[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ],
    'frictionFactor':[ 0.000001, 0.00294093238243111, 0.0147046619121555,0.00294093238243111, 0.0147046619121555, 0.00294093238243111,0.000001],
    'pipeRoughness':[ 0.001, 0.0026, 0.036,0.0026, 0.0026, 0.0026,0.0001]
})
 
#Settings
settings = {
    "deltaConv": 1e-3,
    "residualConv": -1,
    "Printout": False,
    "maxruns": 100,
    "dynamicFF": True,
    "monotonic": False,
    "outstring": "wiki_example_1",
    "visualise": True
}
