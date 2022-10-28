import hcsolver as hc
import argparse
import os
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hardy cross solver.')
    parser.add_argument(
        'input',
        help='Filepath to a  settings input file written in python. This must always be entered.',
        type=str
    )

    parser.add_argument(
        '-n',
        '--nodecsv',
        dest="nodeInput",
        help='Filepath to a csv containing a description of the pipe network nodes. Overrides anything specified in the input file',
        type=str
    )

    parser.add_argument(
        '-e',
        '--edgecsv',
        dest="edgeInput",
        help='Filepath to a csv containing a description of the pipe network edges. Overrides anything specified in the input file',
        type=str
    )



    args = parser.parse_args()

    defaultSettings = hc.defaultSettings
    settings={}
    nodes = None
    edges = None


    basedir = os.path.dirname(os.path.realpath(__file__))
    input = open(args.input, "r")
    jobdir = os.path.dirname(os.path.realpath(args.input))
    vars = input.read()
    exec(vars)

    inputSettings = settings.copy()
    settings=defaultSettings
    settings.update(inputSettings)



    if args.nodeInput is not None:
        nodes = pd.read_csv(args.nodeInput)
    if args.edgeInput is not None:
        edges = pd.read_csv(args.edgeInput)

    if settings is {}:
        raise ValueError('Settings have not been specified.')
    if nodes is None:
        raise ValueError('Nodes have not been specified in the supplied inputs.')
    if edges is None:
        raise ValueError('Edges have not been specified in the supplied inputs.')

    pipenetwork = hc.PipeNetwork(nodes, edges)

    hc.solve(pipenetwork,settings)

#    pipenetwork.cleanEmpty()

    historyOutput="out_"+settings["outstring"]+"_hist.csv"
    histcsv=open(os.path.join(jobdir,historyOutput),'w')
    histcsv.write(pipenetwork.historyCSV())
    histcsv.close()

    nodeOutput = "out_" + settings["outstring"] + "_nodes.csv"
    nodecsv=open(os.path.join(jobdir,nodeOutput),'w')
    nodecsv.write(pipenetwork.nodeCSV())
    nodecsv.close()

    edgeOutput = "out_" + settings["outstring"] + "_edges.csv"
    edgecsv=open(os.path.join(jobdir,edgeOutput),'w')
    edgecsv.write(pipenetwork.edgeCSV())
    edgecsv.close()

    if settings["visualise"]:
        pngOutput = "out_" + settings["outstring"] + "_vis.png"
        pipenetwork.visualise(os.path.join(jobdir,pngOutput),settings)