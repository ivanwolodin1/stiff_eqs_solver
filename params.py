def get_simulation_params():
    return {
        't_range': 500,
        'dt': 1e-4,
        'Ca': 0.00072,
        'G0': 0.0000098,
        'phi': 1e-6,
        'Pr': 6.99,
        'Bi': 1,
        'Ma': 0.167,
        'omega': 1.5,
        'b': 2,
        'bc': {'x': 'periodic'},
    }
