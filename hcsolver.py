import fluidmech
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd



defaultSettings = {
    "deltaConv":1e-3,
    "residualConv":-1,
    "Printout":False,
    "maxruns": 100,
    "dynamicFF":False,
    "monotonic": False,
    "accelerator": 1.0
}



def solve(pipeNetwork,settings=defaultSettings):
    Printout = settings["Printout"]
#    dynamicFF = settings["dynamicFF"]
    deltaConv = settings["deltaConv"]
    residualConv = settings["residualConv"]
    maxruns = settings["maxruns"]
    monotonic = settings["monotonic"]

    if Printout:
        print("starting solve with settings:")
        print(settings)

    lastDelta = 99999999999
    lastResidual = 999999999999

    mce = pipeNetwork.getMCE()[0]
    if math.isclose(mce,0,abs_tol=deltaConv):
        run = True
        iterationCounter = 0
        while(run):

            if Printout:
                print("Iteration: ", iterationCounter)

            (currentDelta,currentResidual) = iterate(pipeNetwork,settings)

            if (currentDelta < deltaConv or deltaConv < 0) and (currentResidual < residualConv or residualConv < 0):
                print("convergence criterion reached")
                run = False
            if iterationCounter > maxruns:
                # Break if iterations maxes out.
                print("Error! Max Runs reached.")
                run = False
            if monotonic and (currentDelta > lastDelta):
                # break if monotonically decreasing convergence is expected
                print("Error!, convergence is no longer monotonically decreasing.")
                run = False

            lastDelta = currentDelta
            lastResidual = currentResidual

            print("Iteration, convergence criterion, Residual:", iterationCounter, lastDelta, lastResidual)
            # printout Data
            if Printout:
                for node in pipeNetwork.nodes:
                    for key, value in pipeNetwork.nodes[node].items():
                        print(node, ' - ', key, " : ", value)
                for edge in pipeNetwork.edges:
                    for key, value in pipeNetwork.edges[edge].items():
                        print(edge, ' - ', key, " : ", value)

            iterationCounter += 1
    else:
        print("Warning! mass flows are not balanced: ",mce)

def iterate(pipeNetwork,settings=defaultSettings):

    deltaFlowLoops=computedF(pipeNetwork,settings)
    updateFlows(pipeNetwork,deltaFlowLoops,settings)
    residual = computeP(pipeNetwork,settings)
    totDelta = sum(abs(delta) for delta in deltaFlowLoops)
    pipeNetwork.appendHistory(totDelta,residual)
    return (totDelta,residual)


def computeP(pipeNetwork,settings=defaultSettings):
    Printout = settings["Printout"]

    flows = np.array([[e[-1] ** 2 for e in pipeNetwork.edges(data='flow')]]).transpose()

#    print("flows1",flows)
    hydraulicResistances = np.array([[k[-1]  for k in pipeNetwork.edges(data='c')]]).transpose()
#    print("hydraulicRes",hydraulicResistances)
    pumpheadGains = np.array([[p[-1]  for p in pipeNetwork.edges(data='pumpHeadGain')]]).transpose()
#    print("pumpheadGains",pumpheadGains)

#    print("product", hydraulicResistances * pumpheadGains)
    flows = flows + hydraulicResistances * pumpheadGains

#    print("flows2",flows)

    conductanceMat = nx.incidence_matrix(pipeNetwork.digraph, oriented=True, weight='c').toarray().transpose()

    # Precondition matrix with known boundary conditions
    pressureBoundaries = np.array([[0.0] * len(pipeNetwork.edges)]).transpose()
    i = 0
    for node in pipeNetwork.nodes:
        if pipeNetwork.nodes[node]['pboundary']:
            pressureBoundaries += np.array([conductanceMat[:, i]]).transpose() * pipeNetwork.nodes[node][
                'pressure']
            conductanceMat = np.delete(conductanceMat, i, 1)
        i += 1
    if Printout:
        print("solving", conductanceMat, "\n x pressures =", flows, "-\n", pressureBoundaries)
    solution = np.linalg.lstsq(conductanceMat, (flows - pressureBoundaries), rcond=None)
    if Printout:
        print("pressure Residuals:", solution[1])
    pressures = solution[0]

    i = 0
    s = 0
    for node in pipeNetwork.nodes:
        if not pipeNetwork.nodes[node]['pboundary']:
            pipeNetwork.nodes[node]['pressure'] = pressures[i - s][0]
        else:
            s += 1
        i += 1
    return solution[1][0]

def computedF(pipeNetwork,settings=defaultSettings):
    Printout = settings["Printout"]
    dynamicFF = settings["dynamicFF"]
    deltaConv = settings["deltaConv"]
    residualConv = settings["residualConv"]
    maxruns = settings["maxruns"]
    monotonic = settings["monotonic"]



    ##STEP 3.1 for every edge, update hydraulic resistance and conductance (k and c where k=1/c)
    for edge in pipeNetwork.edges:
        flow = pipeNetwork.edges[edge]["flow"]
        length = pipeNetwork.edges[edge]["length"]
        diam = pipeNetwork.edges[edge]["diam"]
        if dynamicFF:
            pipeRoughness = pipeNetwork.edges[edge]["pipeRoughness"]
            pipeNetwork.edges[edge]["frictionFactor"] = fluidmech.calcff(pipeRoughness, diam,
                                                                         fluidmech.calcRe(diam, flow))
        ff = pipeNetwork.edges[edge]["frictionFactor"]
        k = fluidmech.calck(length, diam, ff)
        pipeNetwork.edges[edge]['k'] = k
        if ff != 0 and length != 0:
            pipeNetwork.edges[edge]['c'] = -1 / k
        else:
            pipeNetwork.edges[edge]['c'] = 1.7976931348623158e+308

    for edge in pipeNetwork.uedges:
        flow = pipeNetwork.uedges[edge]["flow"]
        length = pipeNetwork.uedges[edge]["length"]
        diam = pipeNetwork.uedges[edge]["diam"]
        if dynamicFF:
            pipeRoughness = pipeNetwork.uedges[edge]["pipeRoughness"]
            pipeNetwork.uedges[edge]["frictionFactor"] = fluidmech.calcff(pipeRoughness,
                                                                          diam, fluidmech.calcRe(diam, flow))
        ff = pipeNetwork.uedges[edge]["frictionFactor"]
        k = fluidmech.calck(length, diam, ff)
        pipeNetwork.uedges[edge]['k'] = k
        if ff != 0 and length != 0:
            pipeNetwork.uedges[edge]['c'] = -1 / k
        else:
            pipeNetwork.uedges[edge]['c'] = 1.7976931348623158e+308

    ##STEP 3.2 for each loop, determine head loss (kQ^n) in clockwise and counterclockwise directions.
    ##CW and CCW is determined as following and anti-following the direction of the graph.

    hlLoops = []
    for loop in pipeNetwork.loops:
        hlLoop = []
        for node in loop:
            thisnode = node
            nextnode = node
            if loop.index(node) + 1 < len(loop):
                nextnode = loop[loop.index(node) + 1]
            else:
                nextnode = loop[0]
            CWtpl = (thisnode, nextnode)
            CCWtpl = (nextnode, thisnode)

            if CWtpl in list(pipeNetwork.edges):
                cwflow = pipeNetwork.digraph[thisnode][nextnode]["flow"]
                pump = pipeNetwork.digraph[thisnode][nextnode]['pumpHeadGain']
                k = pipeNetwork.digraph[thisnode][nextnode]['k']
                hl = fluidmech.calchl(k, cwflow, pump)
                hlLoop.append(hl)
                pipeNetwork.digraph[thisnode][nextnode]["headloss"] = hl
            if CCWtpl in list(pipeNetwork.edges):
                ccwflow = pipeNetwork.digraph[nextnode][thisnode]["flow"]
                pump = pipeNetwork.digraph[nextnode][thisnode]['pumpHeadGain']
                k = pipeNetwork.digraph[nextnode][thisnode]['k']
                hl = fluidmech.calchl(k, ccwflow, pump)
                hlLoop.append(-hl)
                pipeNetwork.digraph[nextnode][thisnode]["headloss"] = -hl
        hlLoops.append(hlLoop)

    ##STEP 4 For each loop, find sum head losses (sum kQ^n)
    sumhlLoops = []
    for loop in hlLoops:
        sumhl = sum(loop)
        sumhlLoops.append(sumhl)

    ##STEP 5 For each loop, find sum of derivative of head loss w.r.t. flow (sum nkQ^(n-1))
    drhlLoops = []
    for loop in pipeNetwork.loops:
        drhlLoop = []
        for node in loop:
            thisnode = node
            nextnode = node
            if loop.index(node) + 1 < len(loop):
                nextnode = loop[loop.index(node) + 1]
            else:
                nextnode = loop[0]

            flow = pipeNetwork.ugraph[thisnode][nextnode]["flow"]

            if dynamicFF:
                diam = pipeNetwork.ugraph[thisnode][nextnode]["diam"]
                length = pipeNetwork.ugraph[thisnode][nextnode]["length"]
                eps = pipeNetwork.ugraph[thisnode][nextnode]["pipeRoughness"]
                pump = pipeNetwork.ugraph[thisnode][nextnode]['pumpHeadGain']
                dhl = fluidmech.numerical_dhldQ(flow,length,diam,eps,pump,997.77,1.0016E-3,1e-3)
            else:
                k = pipeNetwork.ugraph[thisnode][nextnode]["k"]
                dhl = fluidmech.calcdhl(k, flow)
            #            print(thisnode, nextnode, dhl)
            drhlLoop.append(dhl)
        drhlLoops.append(drhlLoop)

    sumdrhlLoops = []
    for loop in drhlLoops:
        sumdrhlLoops.append(sum(loop))

    ##STEP 6 Determine change in flow given as output of Step 4 over output of Step 5
    deltaFlowloop = []
    for i in range(0, len(pipeNetwork.loops)):
        delta = sumhlLoops[i] / sumdrhlLoops[i]
        deltaFlowloop.append(delta)
    if Printout:
        print(deltaFlowloop)
    return deltaFlowloop

def updateFlows(pipeNetwork,deltaFlowLoops,settings=defaultSettings):
    accelerator = settings["accelerator"]

    i = 0
    for loop in pipeNetwork.loops:
        delta = deltaFlowLoops[i]
        for node in loop:
            thisnode = node
            nextnode = node
            if loop.index(node) + 1 < len(loop):
                nextnode = loop[loop.index(node) + 1]
            else:
                nextnode = loop[0]
            CWtpl = (thisnode, nextnode)
            CCWtpl = (nextnode, thisnode)

            if CWtpl in list(pipeNetwork.edges):
                pipeNetwork.digraph[thisnode][nextnode]["flow"] -= delta * accelerator
            if CCWtpl in list(pipeNetwork.edges):
                pipeNetwork.digraph[nextnode][thisnode]["flow"] += delta * accelerator
        i += 1


class PipeNetwork:
    def __init__(self, nodes: pd.DataFrame, edges: pd.DataFrame):

        #TODO: Write input checker
        # mandatoryAttributes = ["to","from","name","flow","length","daim"]
        # mustHaveOneAttrubte = ["pipeRoughness","frictionFactor"]
        # optionalAttrubutes = ["pumpHeadGain"]

        self.digraph  = nx.from_pandas_edgelist(edges, 'from', 'to',
                                       edges.columns.values.tolist(),
                                       create_using=nx.DiGraph())

        nx.set_node_attributes(self.digraph, nodes.set_index('name').to_dict('index'))

        self.ugraph = nx.from_pandas_edgelist(edges, 'from', 'to',
                                        edges.columns.values.tolist(),
                                        create_using=nx.Graph())
        self.loops = list(nx.cycle_basis(self.ugraph))
        self.edges=self.digraph.edges
        self.nodes=self.digraph.nodes
        self.uedges=self.ugraph.edges

        self.history = []
        self.appendHistory(999999,9999999)


    def residual(self):
        return computeP(self)[1][0]

    def deltaFlow(self):
        deltaFlowLoop = computedF(self)
        totDelta = sum(abs(delta) for delta in deltaFlowLoop)
        return totDelta

    def visualise(self,filepath):
        pos = nx.kamada_kawai_layout(self.digraph, weight='length')
        edgeWidths = [self.digraph[u][v]['flow'] / (0.25 * nx.to_pandas_edgelist(self.digraph)['flow'].max()) for u, v in self.edges()]
        pressures = list(nx.get_node_attributes(self.digraph, 'pressure').values())
        nodeColours = [n / pd.DataFrame.from_dict(self.nodes, orient='index')['pressure'].max() for n in pressures]

        # Plot it, providing a continuous color scale with cmap:
        nx.draw(self.digraph, pos, with_labels=True, node_color=nodeColours, width=edgeWidths, cmap=plt.cm.Blues, arrows=True)
        plt.savefig(filepath)

    def updateHistory(self):
        self.history.append([self.deltaFlow(), self.residual(), pd.DataFrame.from_dict(self.nodes, orient='index'), nx.to_pandas_edgelist(self.digraph)])

    def appendHistory(self,deltaFlow,residual):
        self.history.append([deltaFlow, residual, pd.DataFrame.from_dict(self.nodes, orient='index'),
                             nx.to_pandas_edgelist(self.digraph)])


    def nodeCSV(self):
        nodeDF = pd.DataFrame.from_dict(self.nodes, orient='index')
        nodeDF.reset_index(inplace=True)
        nodeDF.rename(columns={'index': 'name'}, inplace=True)

        return nodeDF.to_csv()

    def edgeCSV(self):
        edgeDF = nx.to_pandas_edgelist(self.digraph)
        edgeDF.drop(columns=['source', 'target'],inplace=True)
        return edgeDF.to_csv()

    def historyCSV(self):
            dfcolumns = ["Iteration","convcrit","residual"]
            for node in self.nodes:
               dfcolumns.append(str(node) + '_P')
            for edge in self.edges:
                dfcolumns.append(str(self.edges[edge]["name"]) + '_F')

            unpackedData=[]
            i=0
            for iteration in self.history:
                unpackedIteration=[i,iteration[0],iteration[1]]
                unpackedIteration += iteration[2]["pressure"].tolist() + iteration[3]["flow"].tolist()
                #print(unpackedIteration)
                unpackedData.append(unpackedIteration)
                i+=1
            formatteddata = unpackedData
            historyDataframe = pd.DataFrame(formatteddata,columns=dfcolumns)
            return historyDataframe.to_csv()


    def getMCE(self):
        totalFlow = 0
        debugstrlst = []
        for node in self.nodes:
            if 'fboundary' not in self.nodes[node].keys():
                flowaccumulator = 0
            else:
                flowaccumulator = self.nodes[node]['flow']
            for frm, to, edat in self.digraph.in_edges(node, data=True):
                flowaccumulator += edat['flow']
            for frm, to, edat in self.digraph.out_edges(node, data=True):
                flowaccumulator -= edat['flow']
            debugstr = str(node) + ":" + str(flowaccumulator)
            debugstrlst.append(debugstr)
            totalFlow += flowaccumulator
        return (totalFlow,debugstrlst)