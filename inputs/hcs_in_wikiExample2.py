
#Network Node Data
nodes = pd.DataFrame({
					  'name':['N1', 'N2', 'N3','N4','So','Si'],
                      'pressure':[0,0,0,0,0,100],
					  'flow':[0,0,0,0,0,0],
                      'pboundary':[False,False,False,False,False,True],
                      'fboundary':[False,False,False,False,False,False]
})


#Load Network Edge Data 
edges = pd.DataFrame({
	'name':['PSo', 'P1', 'P2','P3','P4','P5','PipeSi','Pump'], 
    'from':['So','N1','N1','N2','N2','N3','N4','Si'],
    'to':['N1','N2','N3','N3','N4','N4','Si','So'],
    'flow':[ 10, 5, 5, 0, 5, 5,10,10],
    'length':[5,10, 10, 10, 10,10,5,10],
    'diam':[ 0.3, 0.3, 0.1, 0.3, 0.1, 0.3, 0.3,0.3 ],
    'pumpHeadGain':[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,35.0 ],
    'frictionFactor':[ 0.000001, 0.00294093238243111, 0.0147046619121555,0.00294093238243111, 0.0147046619121555, 0.00294093238243111,0.000001,0.0001],
    'pipeRoughness':[ 0.0001, 0.0026, 0.040,0.0026, 0.0026, 0.0026,0.0001,0.0001]
})
 
#Settings
settings = {
    "deltaConv": 1e-4,
    "residualConv": 1e-4,
    "Printout": False,
    "maxruns": 100000,
    "dynamicFF": False,
    "monotonic": True,
    "outstring": "wiki_example_2",
    "visualise": True,
    "accelerator" : 1.0
}
