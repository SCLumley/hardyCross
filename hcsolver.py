import fluidmech
import math
import numpy as np
import networkx as nx
import pipeNetwork

defaultSettings = {
    "deltaConv": 1e-3,
    "residualConv": -1,
    "Printout": False,
    "maxruns": 100,
    "dynamicFF": False,
    "monotonic": False,
    "accelerator": 1.0,
    "visualise": False,
    "visDisplayLvl": 1,
    "rho": fluidmech.stprho,
    "mu": fluidmech.stpmu
}

class Solver:
    def __init__(self,config=defaultSettings):
        self.configuration = defaultSettings
        self.configuration.update(config)


    def updateSettings(self,settings):
        self.configuration.update(settings)


    def solve(self,pipeNetwork:pipeNetwork.PipeNetwork, settings= None):
        if settings == None:
            settings = self.configuration

        Printout = settings["Printout"]
        deltaConv = settings["deltaConv"]
        residualConv = settings["residualConv"]
        maxruns = settings["maxruns"]
        monotonic = settings["monotonic"]
        rho = settings["rho"]
        mu = settings["mu"]

        if Printout:
            print("starting solve with settings:")
            print(settings)

        # setting up some initial values:
        for edge in pipeNetwork.edges:
            if 'frictionFactor' not in pipeNetwork.edges[edge]:
                if Printout:
                    print("Initial friction factor not specified. Calculating from pipe characteristics and initial flow.:")
                flow = pipeNetwork.edges[edge]["flow"]
                diam = pipeNetwork.edges[edge]["diam"]
                pipeRoughness = pipeNetwork.edges[edge]["pipeRoughness"]
                ff = fluidmech.calcff(pipeRoughness, diam, fluidmech.calcRe(diam, fluidmech.calcmdot(flow, rho), mu))
                pipeNetwork.edges[edge]["frictionFactor"] = ff
                if Printout:
                    print(pipeNetwork.edges[edge])
            pipeNetwork.edges[edge]["headloss"] = 0

        for edge in pipeNetwork.uedges:
            if 'frictionFactor' not in pipeNetwork.uedges[edge]:
                flow = pipeNetwork.uedges[edge]["flow"]
                diam = pipeNetwork.uedges[edge]["diam"]
                pipeRoughness = pipeNetwork.uedges[edge]["pipeRoughness"]
                ff = fluidmech.calcff(pipeRoughness, diam, fluidmech.calcRe(diam, fluidmech.calcmdot(flow, rho), mu))
                pipeNetwork.uedges[edge]["frictionFactor"] = ff

            pipeNetwork.uedges[edge]["headloss"] = 0

        lastDelta = 99999999999
        lastResidual = 999999999999

        #Preflight Checks to make sure there are no fundamentally broken aspects of the input network.
        mce = pipeNetwork.getMCE(rho)
        if not math.isclose(mce[0], 0, abs_tol=deltaConv):
            raise ValueError("Error! mass flows are not balanced to within the convergence criteria. Check your inputs and try again.: ", mce[1])

        if not self.checkBoundaryConditions(pipeNetwork,settings):
            raise ValueError("Error! Boundary conditions are incoherent. Check your inputs and try again. ")

        run = True
        iterationCounter = 0
        while (run):

            if Printout:
                print("Iteration: ", iterationCounter)

            #Main calculation executed in each loop
            (currentDelta, currentResidual) = self.iterate(pipeNetwork, settings)

            print("Iteration, convergence criterion, Residual:", iterationCounter, currentDelta, currentResidual)
            if (currentDelta < deltaConv or deltaConv < 0) and (currentResidual < residualConv or residualConv < 0):
                print("Success! Convergence criterion reached")
                run = False
            if iterationCounter > maxruns:
                # Break if iterations maxes out.
                print("Abort! Max Runs reached.")
                run = False
            if monotonic and (currentDelta > lastDelta):
                # break if monotonically decreasing convergence is expected
                print("Abort! Convergence is no longer monotonically decreasing.")
                run = False

            lastDelta = currentDelta
            lastResidual = currentResidual

            # printout Data
            if Printout:
                for node in pipeNetwork.nodes:
                    for key, value in pipeNetwork.nodes[node].items():
                        print(node, ' - ', key, " : ", value)
                for edge in pipeNetwork.edges:
                    for key, value in pipeNetwork.edges[edge].items():
                        print(edge, ' - ', key, " : ", value)

            iterationCounter += 1



    def iterate(self,pipeNetwork:pipeNetwork.PipeNetwork, settings = None):
        if settings == None:
            settings = self.configuration
        deltaFlowLoops = self.computedF(pipeNetwork, settings)
        self.updateFlows(pipeNetwork, deltaFlowLoops, settings)
        residual = self.computeP(pipeNetwork, settings)
        totDelta = sum(abs(delta) for delta in deltaFlowLoops)
        pipeNetwork.updateHistory(totDelta, residual)

        if math.isinf(totDelta) or math.isnan(totDelta) or math.isinf(residual) or math.isnan(residual):
            raise ValueError('Error! System has become numerically unstable.')

        return (totDelta, residual)


    def computeP(self,pipeNetwork:pipeNetwork.PipeNetwork, settings=None):
        if settings == None:
            settings = self.configuration
        Printout = settings["Printout"]

        if Printout:
            for node in pipeNetwork.nodes:
                print(pipeNetwork.nodes[node])

        flows = np.array([[e[-1] ** 2 for e in pipeNetwork.edges(data='flow')]]).transpose()
        hydraulicResistances = np.array([[k[-1] for k in pipeNetwork.edges(data='c')]]).transpose()
        pumpheadGains = np.array([[p[-1] for p in pipeNetwork.edges(data='pumpHeadGain')]]).transpose()
        flows = flows + hydraulicResistances * pumpheadGains
        conductanceMat = nx.incidence_matrix(pipeNetwork.digraph, oriented=True, weight='c').toarray().transpose()

        # Precondition matrix with known boundary conditions
        pressureBoundaries = np.array([[0.0] * len(pipeNetwork.edges)]).transpose()
        i = 0
        for node in pipeNetwork.nodes:
            if 'pboundary' in pipeNetwork.nodes[node] and pipeNetwork.nodes[node]['pboundary']:
                pressureBoundaries += np.array([conductanceMat[:, i]]).transpose() * pipeNetwork.nodes[node][
                    'pressure']
                conductanceMat = np.delete(conductanceMat, i, 1)
            i += 1
        if Printout:
            print("solving", conductanceMat, "\n x pressures =", flows, "-\n", pressureBoundaries)
        solution = np.linalg.lstsq(conductanceMat, (flows - pressureBoundaries), rcond=None)
        if Printout:
            print("Solution and Residuals:", solution)
        pressures = solution[0]

        i = 0
        s = 0
        for node in pipeNetwork.nodes:
            if 'pboundary' in pipeNetwork.nodes[node] and not pipeNetwork.nodes[node]['pboundary']:
                pipeNetwork.nodes[node]['pressure'] = pressures[i - s][0]
            else:
                s += 1
            i += 1
        return solution[1][0]


    def computedF(self,pipeNetwork:pipeNetwork.PipeNetwork, settings=None):
        if settings == None:
            settings = self.configuration
        Printout = settings["Printout"]
        dynamicFF = settings["dynamicFF"]
        rho = settings["rho"]
        mu = settings["mu"]


        ##STEP 3.1 for every edge, update hydraulic resistance and conductance (k and c where k=1/c),
        # and calculate headlosses (kQ^2) and derivatives
        for edge in pipeNetwork.edges:
            flow = pipeNetwork.edges[edge]["flow"]
            length = pipeNetwork.edges[edge]["length"]
            diam = pipeNetwork.edges[edge]["diam"]
            pump = pipeNetwork.edges[edge]["pumpHeadGain"]
            pipeRoughness = pipeNetwork.edges[edge]["pipeRoughness"]
            if dynamicFF:
                pipeNetwork.edges[edge]["frictionFactor"] = fluidmech.calcff(
                        pipeRoughness,
                        diam,
                        fluidmech.calcRe(diam, fluidmech.calcmdot(flow, rho),mu)
                    )
            ff = pipeNetwork.edges[edge]["frictionFactor"]
            k = fluidmech.calck(length, diam, ff)
            pipeNetwork.edges[edge]['k'] = k
            hl = fluidmech.calchl(k, flow, pump)
            pipeNetwork.edges[edge]["headloss"] = hl
            if dynamicFF:
                dhl = fluidmech.numerical_dhldQ(flow, length, diam, pipeRoughness, pump, rho, mu, 1e-3)
            else:
                dhl = fluidmech.calcdhl(k, flow)
            pipeNetwork.edges[edge]["dhl"] = dhl

            if ff != 0 and length != 0:
                pipeNetwork.edges[edge]['c'] = -1 / k
            else:
                pipeNetwork.edges[edge]['c'] = float('inf') # zero resistance implies infinite conductance.


        ##STEP 3.2 for each loop, determine head loss in clockwise and counterclockwise directions.
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
                    hl = pipeNetwork.digraph[thisnode][nextnode]["headloss"]
                    hlLoop.append(hl)
                if CCWtpl in list(pipeNetwork.edges):
                    hl = pipeNetwork.digraph[nextnode][thisnode]["headloss"]
                    hlLoop.append(-hl)
            hlLoops.append(hlLoop)


        ##STEP 4 For each loop, sum head losses (sum kQ^n)
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
                dhl = pipeNetwork.ugraph[thisnode][nextnode]["dhl"]
                drhlLoop.append(dhl)
            drhlLoops.append(drhlLoop)

        sumdhlLoops = []

        for loop in drhlLoops:
            sumdhlLoops.append(sum(loop))

        ##STEP 6 Determine change in flow given as output of Step 4 over output of Step 5
        deltaFlowloops = []
        for i in range(0, len(pipeNetwork.loops)):
            delta = sumhlLoops[i] / sumdhlLoops[i]
            deltaFlowloops.append(delta)
        if Printout:
            print("deltaFlowLoop:",deltaFlowloops)
        return deltaFlowloops


    def updateFlows(self,pipeNetwork:pipeNetwork.PipeNetwork, deltaFlowLoops, settings=None):
        if settings == None:
            settings = self.configuration
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


    def getMCE(self,pipeNetwork:pipeNetwork.PipeNetwork,settings=None):
        if settings == None:
            settings = self.configuration
        rho = settings["rho"]
        totalFlow = 0
        debugstrlst = []
        for node in self.nodes:
            if 'fboundary' not in self.nodes[node].keys():
                flowaccumulator = 0
            else:
                flowaccumulator = fluidmech.calcmdot(self.nodes[node]['flow'],rho)
            for frm, to, edat in self.digraph.in_edges(node, data=True):
                flowaccumulator += fluidmech.calcmdot(edat['flow'],rho)
            for frm, to, edat in self.digraph.out_edges(node, data=True):
                flowaccumulator -= fluidmech.calcmdot(edat['flow'],rho)
            debugstr = str(node) + ":" + str(flowaccumulator)
            debugstrlst.append(debugstr)
            totalFlow += flowaccumulator
        return (totalFlow, debugstrlst)

    def getFlowDelta(self,pipeNetwork:pipeNetwork.PipeNetwork,settings=None):
        if settings == None:
            settings = self.configuration
        deltaFlowLoops = self.computedF(pipeNetwork, settings)
        return sum(abs(delta) for delta in deltaFlowLoops)


    def getResidual(self,pipeNetwork:pipeNetwork.PipeNetwork,settings=None):
        if settings == None:
            settings = self.configuration
        return self.computeP(pipeNetwork, settings)

    #TODO write something to check boundary conditions are coherent
    def checkBoundaryConditions(self,pipeNetwork:pipeNetwork.PipeNetwork,settings=None):
        print("Warning! There is no current checking of boundary condition coherency. This is just a placeholder")
        return True

