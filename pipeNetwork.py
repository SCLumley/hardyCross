import fluidmech
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

class PipeNetwork:
    def __init__(self, nodes: pd.DataFrame, edges: pd.DataFrame):

        #input checking
        mandatoryEdgeAttributes = ["to","from","name","flow","length","diam"]
        for MA in mandatoryEdgeAttributes:
            if MA not in edges:
                raise KeyError("Error! the pipe network is missing the " + MA + " edge attribute.")
        # mustHaveOneAttribute = ["pipeRoughness","frictionFactor"]
        if 'frictionFactor' not in edges and 'pipeRoughness' not in edges:
            raise KeyError("Error! the pipe edges must specify a friction factor, or a pipe roughness or both.")
        if 'pumpHeadGain' not in edges:
            edges['pumpHeadGain'] = 0.0

        mandatoryNodeAttributes = ["name", "pressure", "flow", "pboundary","fboundary"]
        for MA in mandatoryNodeAttributes:
            if MA not in nodes:
                raise KeyError("Error! the pipe network is missing the " + MA + " node attribute.")



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
        i = 0
        for iteration in self.history:
            unpackedIteration = [i, iteration[0], iteration[1]]
            unpackedIteration += iteration[2]["pressure"].tolist() + iteration[3]["flow"].tolist()
            # print(unpackedIteration)
            unpackedData.append(unpackedIteration)
            i += 1
        formatteddata = unpackedData
        historyDataframe = pd.DataFrame(formatteddata, columns=dfcolumns)
        return historyDataframe.to_csv()

