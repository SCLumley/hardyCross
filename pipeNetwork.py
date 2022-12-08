import fluidmech
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd


class PipeNetwork:
    """
    Description: The PipeNetwork class is predominantly a data-structure that uses a directed graph to model a network made of fluid pipes and junctions
    The class is for storing the inital data that defines a network, storing history of changes to that network, and retreving datab about the prior data.
    It is not designed to perform any serious calculation itself (other than limited input validation). That is handled by a solver class.
    """


    def __init__(self, nodes: pd.DataFrame, edges: pd.DataFrame):
        """
        @param nodes: A pandas dataframe that describes list of nodes defining the junctions of the pipe network. For example:
            nodes = pd.DataFrame({
                          'name':['N1', 'N2', 'N3','N4','So','Si'],
                          'pressure':[0,30,0,0,0,1],
                          'flow':[0,0,0,0,+10,-10],
                          'pboundary':[False,True,False,False,False,False],
                          'fboundary':[False,False,False,False,True,True]
            })
        Each node must have the following parameters defined:
            name: A short string that is an identifiable name for the junction. Names must be unique.
            pressure: The starting pressure at that node in Pa.
            flow: The starting loss or gain of water that is put in or taken out of the network at that node in m^3s^(-1).
            pboundary: A bool defining if the node is a pressure boundary condition. If True, the node pressure will be held fixed during the simulation.
            fboundary: (NOT YET FULLY IMPLIMENTED) A bool defining if the node is a flow boundary condition. If True, the node flow will be held fixed during the simulation


        @param edges: A pandas dataframe that describes list of edges defining the pipes of the pipe network. For example:
            edges = pd.DataFrame({
                'name':['PSo', 'P1', 'P2','P3','P4','P5','PSi'],
                'from':['So','N1','N1','N2','N2','N3','N4'],
                'to':['N1','N2','N3','N3','N4','N4','Si'],
                'flow':[ 10, 5, 5, 0, 5, 5,10],
                'length':[5,10, 10, 10, 10,10,5],
                'diam':[ 0.3, 0.3, 0.3, 0.3, 0.1, 0.3, 0.3 ],
                'pumpHeadGain':[ 0.0, 10000.0, 0.0, 0.0, 0.0, 0.0, 0.0 ],
                'frictionFactor':[ 0.000001, 0.00294093238243111, 0.0147046619121555,0.00294093238243111, 0.0147046619121555, 0.00294093238243111,0.000001],
                'pipeRoughness':[ 0.001, 0.0026, 0.036,0.0026, 0.0026, 0.0026,0.0001]
            })
        Each edge must have the following parameters defined:
                name: A short string that is an identifiable name for the pipe. Names must be unique.
                from:The name of the node the pipe starts at. This name MUST be identical to a name in the nodes input.
                to:The name of the node the pipe ends at. This name MUST be identical to a name in the nodes input.
                flow: The starting Flow along the pipe in m^3s^(-1), where positive numbers align with flow originating at the "from" node and going to the "to" node.
                length: The length of the pipe in meters.
                diam: The diameter of the pipe in meters.
                pumpHeadGain: The pressure supplied by any pumps that may be on this length of pipe, in Pa

        The edge must have at least one of the following parameters defined:
                frictionFactor: The moody friction factor of the pipe.
                pipeRoughness:the surface roughness of the pipe in meters.
        """



        #input checking
        mandatoryEdgeAttributes = ["to","from","name","flow","length","diam"]
        for MA in mandatoryEdgeAttributes:
            if MA not in edges:
                raise KeyError("Error! the pipe network is missing the " + MA + " edge attribute.")
        # mustHaveOneAttribute = ["pipeRoughness","frictionFactor"]
        if 'frictionFactor' not in edges and 'pipeRoughness' not in edges:
            raise KeyError("Error! the pipe edges must specify a friction factor, or a pipe roughness or both.")

        mandatoryNodeAttributes = ["name", "pressure", "flow", "pboundary","fboundary"]
        for MA in mandatoryNodeAttributes:
            if MA not in nodes:
                raise KeyError("Error! the pipe network is missing the " + MA + " node attribute.")

        #input rationalisation
        if 'pumpHeadGain' not in edges:
            edges['pumpHeadGain'] = 0.0

        #build graphs
        self.digraph = nx.from_pandas_edgelist(edges, 'from', 'to',
                                               edges.columns.values.tolist(),
                                               create_using=nx.DiGraph())
        nx.set_node_attributes(self.digraph, nodes.set_index('name').to_dict('index'))
        self.ugraph = self.digraph.to_undirected(as_view=True)
        self.loops = list(nx.cycle_basis(self.ugraph))
        self.edges = self.digraph.edges
        self.nodes = self.digraph.nodes
        self.uedges = self.ugraph.edges
        self.edgesOnLoop = set(e for loop in self.loops for e in loop)
        self.lastResidual = 9999999
        self.lastDeltaFlow = 9999999

        self.history = []
        self.updateHistory(999999, 9999999)


    def cleanEmpty(self):
        #FIXME: Track down and fix the empty node bug, where some inputs create a graph with a mysterios
        # empty node with no properties.
        """

        Description: A placeholder function that clears the graph of any empty nodes and edges. This is largly to deal with a bug at this time.
        """
        delnodes = []
        for n in self.nodes:
            if self.nodes[n] == {}:
                delnodes.append(n)
        for n in delnodes:
            self.digraph.remove_node(n)

        deledges = []
        for e in self.edges:
            if self.edges[e] == {}:
                deledges.append(e)
        for e in deledges:
            self.digraph.remove_edge(e)


    def visualise(self, filepath, vislevel=1):
        self.cleanEmpty()
        plt.rcParams["figure.figsize"] = (16, 9)
        pos = nx.kamada_kawai_layout(self.digraph, weight='length')
        edgeWidths = [self.digraph[u][v]['flow'] / (0.25 * nx.to_pandas_edgelist(self.digraph)['flow'].max()) for u, v
                      in self.edges()]
        pressures = list(nx.get_node_attributes(self.digraph, 'pressure').values())
        nodeColours = [n / pd.DataFrame.from_dict(self.nodes, orient='index')['pressure'].max() for n in pressures]

        # Plot it, providing a continuous color scale with cmap:
        nx.draw(self.digraph, pos, with_labels=False, node_color=nodeColours, width=edgeWidths, cmap=plt.cm.Blues,
                arrows=True)

        edge_labels = {}
        node_labels = {}

        if vislevel == 1:
            edge_labels = dict(
                [((n1, n2), d['name'])
                 for n1, n2, d in self.digraph.edges(data=True)])

            node_labels = dict([(n, n)
                                for n, d in self.digraph.nodes(data=True)])

        if vislevel == 2:
            edge_labels = dict([((n1, n2), d['name'] + "\nQ:" + "{:.2e}".format(d['flow']) + "\nhl:" + "{:.2e}".format(
                d['headloss']))
                                for n1, n2, d in self.digraph.edges(data=True)])

            node_labels = dict([(n, n + "\nP:" + "{:.2e}".format(d['pressure']))
                                for n, d in self.digraph.nodes(data=True)])

        if vislevel > 0:
            nx.draw_networkx_edge_labels(self.digraph, pos, edge_labels=edge_labels)
            nx.draw_networkx_labels(self.digraph, pos, labels=node_labels)
        plt.savefig(filepath)

    def updateHistory(self, deltaFlow, residual):
        self.history.append([deltaFlow, residual, pd.DataFrame.from_dict(self.nodes, orient='index'),
                             nx.to_pandas_edgelist(self.digraph)])

    def nodeCSV(self):
        nodeDF = pd.DataFrame.from_dict(self.nodes, orient='index')
        nodeDF.reset_index(inplace=True)
        nodeDF.rename(columns={'index': 'name'}, inplace=True)

        return nodeDF.to_csv()

    def edgeCSV(self):
        edgeDF = nx.to_pandas_edgelist(self.digraph)
        edgeDF.drop(columns=['source', 'target'], inplace=True)
        return edgeDF.to_csv()

    def historyCSV(self):
        dfcolumns = ["Iteration", "convcrit", "residual"]
        for node in self.nodes:
            if 'pboundary' in self.nodes[node]:
                dfcolumns.append(str(node) + '_P')
        for edge in self.edges:
            dfcolumns.append(str(self.edges[edge]["name"]) + '_F')

        unpackedData = []

        for i,iteration in enumerate(self.history):
            unpackedIteration = [i, iteration[0], iteration[1]]
            unpackedIteration += iteration[2]["pressure"].tolist() + iteration[3]["flow"].tolist()
            # print(unpackedIteration)
            unpackedData.append(unpackedIteration)
        formatteddata = unpackedData
        historyDataframe = pd.DataFrame(formatteddata, columns=dfcolumns)
        return historyDataframe.to_csv()

