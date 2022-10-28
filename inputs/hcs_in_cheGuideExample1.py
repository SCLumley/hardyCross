#################### Input Data

# Based on example 1 taken from https://cheguide.com/pipe_network.html

nodes = pd.DataFrame({
    'name': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
    'pressure': [0, 100, 0, 0, 0, 0, 0, 0],
    'flow': [+0.3, 0, -0.05, -0.1, 0, 0, 0, -0.15],
    'fboundary': [True, False, True, True, False, False, False, True],
    'pboundary': [False, True, False, False, False, False, False, False]
})


edges = pd.DataFrame({
    'name': ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'],
    'from': ['A', 'A', 'B', ' F', 'B', 'D', 'C', 'G', 'G', 'E'],
    'to':   ['B', 'F', 'D', 'G', 'C', 'E', 'E', 'D', 'H', 'H'],
    'flow': [0.2, 0.1, 0.08, 0.12, 0.07, 0.03, 0.1, 0.05, 0.05, 0.10],
    'length': [300, 250, 350, 125, 350, 125, 300, 125, 350, 125],
    'diam': [0.3, 0.25, 0.2, 0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.15],
    'pipeRoughness': [0.00087 * 0.30, 0.00104 * 0.25, 0.00130 * 0.20, 0.00130 * 0.20, 0.00130 * 0.20, 0.00130 * 0.20,
                      0.00130 * 0.20, 0.00173 * 0.15, 0.00130 * 0.20, 0.00173 * 0.15]
})

# Settings
settings = {
    "deltaConv": 1e-6,
    "residualConv": 1e-6,
    "Printout": True,
    "maxruns": 100000,
    "dynamicFF": True,
    "monotonic": False,
    "outstring": "cheGuide_example_1",
    "visualise": True,
    "accelerator": 0.1
}


# nodes = pd.DataFrame({
#     'name': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
#     'pressure': [0, 100, 0, 0, 0, 0, 0, 0],
#     'flow': [+0.3, 0, -0.05, -0.1, 0, 0, 0, -0.15],
#     'fboundary': [True, False, True, True, False, False, False, True],
#     'pboundary': [False, True, False, False, False, False, False, False]
# })


# edges = pd.DataFrame({
#     'name': ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'],
#     'from': ['A', 'A', 'B', ' F', 'B', 'D', 'C', 'G', 'G', 'E'],
#     'to': ['B', 'F', 'D', 'G', 'C', 'E', 'E', 'D', 'H', 'H'],
#     'flow': [0.2, 0.1, 0.08, 0.12, 0.07, 0.03, 0.1, 0.05, 0.05, 0.10],
#     'length': [300, 250, 350, 125, 350, 125, 300, 125, 350, 125],
#     'diam': [0.3, 0.25, 0.2, 0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.15],
#     'pipeRoughness': [0.00087 * 0.30, 0.00104 * 0.25, 0.00130 * 0.20, 0.00130 * 0.20, 0.00130 * 0.20, 0.00130 * 0.20,
#                       0.00130 * 0.20, 0.00173 * 0.15, 0.00130 * 0.20, 0.00173 * 0.15]
# })